"""
Polytomy-aware Newick tree tanglegram with crossing minimization.

Build order:
  1. Newick parser         -> Node, parse_newick
  2. Tree traversal utils  -> subtree_size, iter_leaves, leaf_order
  3. Crossing target collection -> collect_targets
  4. Bitmask DP            -> optimal_child_order   (solves ONE node exactly, with K_MAX cap)
  5. 1-side optimization   -> untangle_cross_1side   (applies step 4 to every node, bottom-up)
  6. 2-side optimization   -> untangle_cross_2side   (alternates both trees)
  7. Crossing count        -> crossings              (reused as-is from the original binary code)
  8. Rendering             -> layout_x, plot_tanglegram
"""

import re, os
import numpy as np
import matplotlib.pyplot as plt

# 1. Newick parser

class Node:
    """A tree node. children is a plain list, so there is no limit on the
    number of children (polytomies are supported natively). A scipy linkage
    matrix, by contrast, has exactly two child slots per row and cannot
    represent this.

    branch_length is the length of the edge connecting this node to its
    PARENT (meaningless for the root). If the Newick string gives an
    explicit ":<number>" for this node, that value is used; otherwise it
    defaults to 1.0.
    """
    __slots__ = ("name", "children", "branch_length")

    def __init__(self, name=None, children=None, branch_length=1.0):
        self.name = name
        self.children = children if children is not None else []
        self.branch_length = branch_length

    @property
    def is_leaf(self):
        return len(self.children) == 0

    def __repr__(self):
        if self.is_leaf:
            return f"Leaf({self.name})"
        return f"Node({len(self.children)} children)"

def load_newick(source):
    if os.path.isfile(source):
        with open(source) as f:
            newick_text = f.read()
    else:
        newick_text = source
    return parse_newick(newick_text)

def parse_newick(s):
    """Parse a Newick string like "(A,B,(C,D,E));" into a Node tree.
 
    When '(' is seen, a child is parsed recursively, then a while loop keeps
    reading more children as long as ',' appears. This while loop is what
    makes the parser handle any number of children (2, 3, 22, ...) instead
    of assuming exactly two.
    """
    s = s.strip().rstrip(";")
    pos = 0
 
    def parse_node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1  # consume '('
            children = [parse_node()]
            while s[pos] == ",":
                pos += 1
                children.append(parse_node())
            assert s[pos] == ")"
            pos += 1
            node = Node(children=children)
        else:
            node = Node()
 
        # read label - stop at '[' too, so a following NHX-style annotation
        # like "[&sequenced=true,has_seq_descendant=true]" is not accidentally
        # absorbed into the label text.
        m = re.match(r"[^,():;\[\]]*", s[pos:])
        label = m.group(0)
        pos += len(label)
        if label:
            node.name = label
 
        # NHX-style annotation "[...]" - may itself contain commas, which
        # must NOT be mistaken for a sibling separator by the outer parser.
        # Consumed here as one atomic block and discarded (not needed for
        # crossing minimization or rendering); nested brackets are not
        # expected in standard NHX and are not handled.
        if pos < len(s) and s[pos] == "[":
            end = s.index("]", pos)
            pos = end + 1
 
        # branch length ":0.123" - use it if present, otherwise the Node
        # default (1.0) from __init__ stays in place.
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":[^,()\[\]]*", s[pos:])
            raw = m.group(0)[1:]  # strip leading ':'
            pos += len(m.group(0))
            try:
                node.branch_length = float(raw)
            except ValueError:
                pass  # malformed length -> keep default of 1.0
 
        # a branch length can itself be followed by another "[...]" block
        # in some NHX variants - handle that too.
        if pos < len(s) and s[pos] == "[":
            end = s.index("]", pos)
            pos = end + 1
 
        return node
 
    return parse_node()

# 2. Tree traversal utils

def subtree_size(node):
    """Number of leaves under this node. Works regardless of child count."""
    if node.is_leaf:
        return 1
    return sum(subtree_size(c) for c in node.children)


def iter_leaves(node):
    """Depth-first traversal (DFS): fully explore one child before moving to
    the next. `for c in node.children` loops over however many children
    exist, so this function needs no change between binary and polytomy
    trees. (The original binary code effectively did the same traversal
    hidden inside scipy's dendrogram() call.)
    """
    if node.is_leaf:
        yield node
    else:
        for c in node.children:
            yield from iter_leaves(c)


def leaf_order(root, as_dict=False):
    """Leaf order (i.e. slot positions) obtained from a DFS walk in the
    tree's current child order."""
    order = [leaf.name for leaf in iter_leaves(root)]
    if as_dict:
        return {lab: i for i, lab in enumerate(order)}
    return order

# 3. Crossing target collection

def _edge_map(edges, key_index):
    """label -> list of partner labels this label is connected to.
    key_index=0 means "tree1 label -> [tree2 labels]"."""
    m = {}
    for e in edges:
        m.setdefault(e[key_index], []).append(e[1 - key_index])
    return m


def collect_targets(node, my_edge_map, other_lindex):
    """Gather, for every leaf under this subtree, the slot positions it
    connects to in the other tree, and return them sorted.

    Sorting matters: step 4 uses np.searchsorted (binary search) on these
    arrays, which requires sorted input.
    """
    targets = []
    for leaf in iter_leaves(node):
        for partner_label in my_edge_map.get(leaf.name, []):
            if partner_label in other_lindex:
                targets.append(other_lindex[partner_label])
    return np.sort(np.asarray(targets, dtype=float))

# 4. Bitmask DP - exactly optimize the child order of ONE node

K_MAX = 50  # Cap on the number of children of a SINGLE node (i.e. the length
            # of children_targets passed into one call of optimal_child_order).
            # This is not a cap on the total number of nodes in a layer or in
            # the whole tree - each node is a separate, independent call.


class TooManyChildrenError(ValueError):
    pass


def _pair_crossing_cost(ta, tb):
    """Crossing count if ta is placed before tb.

    ta, tb are sorted arrays of partner-tree slot positions. Placing ta to
    the left of tb causes a crossing for every pair where a value in ta is
    greater than a value in tb (that connection has to cross back over).

    np.searchsorted(ta, tb, side="right") finds, for each value in tb, how
    many values in ta are <= it. Subtracting that from ta.size gives how
    many values in ta are strictly greater than that tb value - i.e. how
    many crossings that tb value causes. Summing over all of tb gives the
    total crossing count between the two groups.
    """
    if ta.size == 0 or tb.size == 0:
        return 0
    hi = np.searchsorted(ta, tb, side="right")
    return int(np.sum(ta.size - hi))


def optimal_child_order(children_targets, k_max=K_MAX):
    """Bitmask DP (Held-Karp style) that finds the exact crossing-minimizing
    order of a node's children.

    children_targets: list of np.array. children_targets[i] = sorted partner
                       slot positions that the i-th child is connected to.

    k_max: cap on len(children_targets) for THIS call only - i.e. this
           node's child count. Other nodes elsewhere in the tree, or other
           nodes in the same layer, are separate calls and are not affected
           by this cap.

    Returns: (best_cost, order) where order is a permutation of indices
             into children_targets.
    """
    k = len(children_targets)
    if k > k_max:
        raise TooManyChildrenError(
            f"This node has {k} children, exceeding the cap of {k_max}. "
            f"Exactly enumerating 2^{k} subsets is not practical."
        )
    if k <= 1:
        return 0, list(range(k))

    # Precompute the pairwise cost w[i][j] = cost of placing child i
    # immediately before child j, for every ordered pair. This is only
    # O(k^2) and is far cheaper than the DP body below.
    #
    # Larger groupings are not computed separately: if children are placed
    # in order i1, i2, ..., im, the total crossing cost equals the sum of
    # all pairwise costs w[ia][ib] for ia before ib. So triples/quadruples
    # etc. are handled implicitly by summing pairwise costs in the DP loop
    # below, not by building a separate 3-way or 4-way cost table.
    w = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                w[i][j] = _pair_crossing_cost(children_targets[i], children_targets[j])

    FULL = 1 << k  # number of subsets = 2^k
    INF = float("inf")
    dp = [INF] * FULL             # dp[mask] = min cost of having placed the subset `mask`
    parent_choice = [-1] * FULL   # which element was placed last to reach `mask`
    dp[0] = 0

    for mask in range(FULL):
        if dp[mask] == INF:
            continue
        cost_so_far = dp[mask]
        rest = (FULL - 1) ^ mask  # elements not yet placed
        m = rest
        while m:
            j = (m & -m).bit_length() - 1  # next candidate element to place
            m &= m - 1

            # Cost of placing j right after everything already in `mask`:
            # sum of w[i][j] for every i already in mask. This is exactly
            # where pairwise costs are combined to score a larger subset.
            add = 0
            mm = mask
            while mm:
                i = (mm & -mm).bit_length() - 1
                mm &= mm - 1
                add += w[i][j]

            new_mask = mask | (1 << j)
            new_cost = cost_so_far + add
            if new_cost < dp[new_mask]:
                dp[new_mask] = new_cost
                parent_choice[new_mask] = j

    # Backtrack: follow "who was placed last" to reconstruct the order.
    order_rev = []
    mask = FULL - 1
    while mask:
        j = parent_choice[mask]
        order_rev.append(j)
        mask ^= (1 << j)
    order = order_rev[::-1]

    return dp[FULL - 1], order

# 5. 1-side optimization: rotate tree1's nodes to match tree2 (bottom-up)

def untangle_cross_1side(tree1, tree2, edges):
    """tree2 stays fixed. Every polytomy node in tree1 is optimized exactly,
    starting from the leaves and working up.

    `for child in node.children: optimize(child)` is called BEFORE this
    node's own optimization - a post-order (bottom-up) traversal. Children
    must already be settled before the parent can compute correct targets,
    which is the standard bottom-up DP principle.

    This function calls optimal_child_order (step 4) once per polytomy node
    encountered during the traversal - "applying step 4 to every node" means
    exactly this: the traversal (this step) repeatedly invokes the DP
    function built in step 4, once for each node that has 2+ children.
    """
    lindex2 = leaf_order(tree2, as_dict=True)
    my_edge_map = _edge_map(edges, key_index=0)

    def optimize(node):
        for child in node.children:
            optimize(child)  # children first

        if len(node.children) >= 2:
            targets = [collect_targets(c, my_edge_map, lindex2) for c in node.children]
            _, order = optimal_child_order(targets)  # K_MAX cap is enforced here
            node.children = [node.children[i] for i in order]
        # if len(node.children) < 2 there is no ordering decision to make

    optimize(tree1)
    return tree1, tree2

# 6. 2-side optimization: alternate both trees

def untangle_cross_2side(tree1, tree2, edges, max_n_iterations=10):
    """Fix one side and optimize the other, then swap, repeatedly. Moving
    one side can create new optimization opportunities for the other side.
    Stops early if no further improvement is found (current == prev).
    """
    edges_inv = [(e[1], e[0]) for e in edges]

    prev = None
    for _ in range(max_n_iterations):
        tree1, tree2 = untangle_cross_1side(tree1, tree2, edges)
        tree2, tree1 = untangle_cross_1side(tree2, tree1, edges_inv)

        current = crossings(tree1, tree2, edges)
        if current == prev:
            break
        prev = current
        if current == 0:
            break

    return tree1, tree2

# 7. Crossing count (reused as-is from the original binary code)

def _count_inversions(values):
    """Merge-sort based inversion counting, O(n log n). This is pure array
    arithmetic with no dependence on tree structure, so it is reused
    unchanged from the original binary implementation.
    """
    def rec(a):
        if a.size < 2:
            return a, 0
        mid = a.size // 2
        left, cl = rec(a[:mid])
        right, cr = rec(a[mid:])
        cross = int(np.sum(left.size - np.searchsorted(left, right, side="right")))
        return np.sort(np.concatenate([left, right])), cl + cr + cross
    return rec(np.asarray(values, dtype=float))[1]


def crossings(tree1, tree2, edges):
    """Count actual crossings given the trees' CURRENT (already decided)
    leaf order. This is a measurement function, distinct from untangle
    (which decides what order to use), so it stays identical regardless of
    whether the trees are binary or polytomous.
    """
    lindex1 = leaf_order(tree1, as_dict=True)
    lindex2 = leaf_order(tree2, as_dict=True)
    p1 = np.array([lindex1[e[0]] for e in edges], dtype=float)
    p2 = np.array([lindex2[e[1]] for e in edges], dtype=float)
    return _count_inversions(p2[np.lexsort((p2, p1))])

# 8. Rendering (replacement for scipy's dendrogram(), which cannot draw
#    polytomies since it is built around a strictly binary linkage matrix)

def layout_x(root):
    """Assign each leaf an integer slot (y-coordinate in the final plot) via
    DFS. `for c in node.children` loops over any number of children, so this
    needs no special-casing for polytomies.
    """
    leaf_x = {}
    counter = [0]

    def assign(node):
        if node.is_leaf:
            x = counter[0]
            counter[0] += 1
            leaf_x[node.name] = x
            return x, x
        xs, xe = [], []
        for c in node.children:
            a, b = assign(c)
            xs.append(a)
            xe.append(b)
        return min(xs), max(xe)

    assign(root)
    return leaf_x


def _max_cumulative_depth(node, depth=0.0):
    """Distance from the root down to the farthest tip, following actual
    branch_length values (each edge contributes its own branch_length,
    defaulting to 1.0 wherever the Newick input didn't specify one).
    Root itself is at depth 0; a child sits at parent_depth + child's own
    branch_length.
    """
    if node.is_leaf:
        return depth
    return max(_max_cumulative_depth(c, depth + c.branch_length) for c in node.children)


def draw_dendrogram(node, ax, leaf_slot, max_depth, mirror=False, lw=1.0):
    """Recursively draw the actual branch structure of the tree, using each
    edge's branch_length (from the Newick input if given, else the default
    of 1.0) rather than a fixed one-level-per-generation spacing.

    For every internal node:
      - this node's own x-position is its cumulative distance from the root
      - one horizontal segment per child, from this node's x out to each
        child's x (child's x = this node's x + child.branch_length)
      - one vertical segment at this node's x spanning the y-range of
        its children

    A node with 3+ children (a polytomy) is drawn exactly the same way as
    a node with 2 - the vertical segment simply spans more than two
    horizontal segments, and the fork sits at one shared x. That is what
    makes ((A,B,C)) render as a single fork joining three tips instead of
    collapsing into nested binary splits.

    `mirror` flips which side of the axis depth=0 (the root) lands on, so
    that both panels can face outward from a shared root while keeping
    tips oriented consistently - see plot_tanglegram for how it is used.
    """

    def to_draw_x(depth):
        return max_depth - depth if mirror else depth

    def _draw(node, depth):
        if node.is_leaf:
            return depth, leaf_slot[node.name]

        child_coords = [_draw(c, depth + c.branch_length) for c in node.children]
        ys = [c[1] for c in child_coords]
        this_y = (min(ys) + max(ys)) / 2

        draw_this_x = to_draw_x(depth)
        for cx, cy in child_coords:
            draw_cx = to_draw_x(cx)
            ax.plot([draw_this_x, draw_cx], [cy, cy], c="k", lw=lw)
        ax.plot([draw_this_x, draw_this_x], [min(ys), max(ys)], c="k", lw=lw)

        return depth, this_y

    _draw(node, 0.0)


def leaf_draw_positions(node, max_depth, mirror):
    """Each leaf's actual draw_x (after mirroring), so labels can be placed
    right at the tip of its branch instead of a fixed axis edge - branch
    lengths mean tips no longer all line up at the same x.
    """
    positions = {}

    def rec(n, depth):
        if n.is_leaf:
            positions[n.name] = max_depth - depth if mirror else depth
        else:
            for c in n.children:
                rec(c, depth + c.branch_length)

    rec(node, 0.0)
    return positions


def plot_tanglegram(tree1, tree2, edges, figsize=(8, 6)):
    """Draw both trees' actual branch structure (including polytomies) on
    the outer two axes, and connecting lines between matched leaves in the
    middle axis.
    """
    labels1 = leaf_order(tree1)
    labels2 = leaf_order(tree2)
    x1 = layout_x(tree1)
    x2 = layout_x(tree2)
    depth1 = _max_cumulative_depth(tree1)
    depth2 = _max_cumulative_depth(tree2)
 
    fig, (ax1, ax3, ax2) = plt.subplots(1, 3, figsize=figsize)
 
    # Left tree: not mirrored. depth is measured root=0, tips=max_depth
    # (cumulative branch length outward from the root), so draw_x=depth
    # already puts tips (large depth) on the side of ax1 adjacent to the
    # middle panel, and the root (depth=0) on the outer left edge.
    draw_dendrogram(tree1, ax1, x1, depth1, mirror=False)
    pos1 = leaf_draw_positions(tree1, depth1, mirror=False)
    # Tips are aligned at depth1 (the side nearest the middle panel). A tip
    # whose real branch stops short (pos1[name] < depth1) gets a dashed
    # extension out to depth1, so every label lines up at the same x.
    for name, y in x1.items():
        if pos1[name] < depth1 - 1e-9:
            ax1.plot([pos1[name], depth1], [y, y], c="0.6", lw=0.7, ls=(0, (2, 2)))
        ax1.text(depth1 + 0.1, y, name, fontsize=8, ha="left", va="center")
    ax1.set_yticks([])
    ax1.set_xticks([])
    ax1.set_xlim(-0.5, depth1 + 0.7)
    ax1.set_ylim(-1, len(labels1))
    for spine in ax1.spines.values():
        spine.set_visible(False)
 
    # Right tree: mirrored, so tips (depth=max_depth) land at draw_x=0 -
    # the side of ax2 adjacent to the middle panel - and the root
    # (depth=0) lands at draw_x=depth2, the outer right edge. Tips are
    # aligned at 0 the same way.
    draw_dendrogram(tree2, ax2, x2, depth2, mirror=True)
    pos2 = leaf_draw_positions(tree2, depth2, mirror=True)
    for name, y in x2.items():
        if pos2[name] > 1e-9:
            ax2.plot([0, pos2[name]], [y, y], c="0.6", lw=0.7, ls=(0, (2, 2)))
        ax2.text(-0.1, y, name, fontsize=8, ha="right", va="center")
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.set_xlim(-0.7, depth2 + 0.5)
    ax2.set_ylim(-1, len(labels2))
    for spine in ax2.spines.values():
        spine.set_visible(False)
 
    ax3.axis("off")
    ax3.set_xlim(0, 1)
    ax3.set_ylim(-1, max(len(labels1), len(labels2)))
 
    for e in edges:
        if e[0] in x1 and e[1] in x2:
            ax3.plot([0, 1], [x1[e[0]], x2[e[1]]], c="k", lw=1)
 
    fig.tight_layout()
    return fig

def run_polytomy_tanglegram(nwk1, nwk2):
    t1 = load_newick(nwk1) 
    t2 = load_newick(nwk2) 

    edges = [(l, l) for l in leaf_order(t1) if l in leaf_order(t2)]
    
    print("crossings before:", crossings(t1, t2, edges))
    t1, t2 = untangle_cross_2side(t1, t2, edges)
    print("crossings after: ", crossings(t1, t2, edges))
   
    fig = plot_tanglegram(t1, t2, edges)
    plt.show()

    return fig


# Example usage

if __name__ == "__main__":
    # t1 mixes explicit branch lengths (A,B,C,F) with edges that have none
    # (D,E and both internal nodes default to 1.0)
    t1 = parse_newick("((A:2,B:1,C:3),(D,E),F:4);")
    t2 = parse_newick("((C,A,B),(E,D),F);")

    edges = [(l, l) for l in leaf_order(t1) if l in leaf_order(t2)]

    print("crossings before:", crossings(t1, t2, edges))
    t1, t2 = untangle_cross_2side(t1, t2, edges)
    print("crossings after: ", crossings(t1, t2, edges))

    fig = plot_tanglegram(t1, t2, edges)
    plt.show()