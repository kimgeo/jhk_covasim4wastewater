# plotting/plot_mutation_tree.py

import matplotlib.pyplot as plt
import networkx as nx


def tidy_tree_layout(G):
    """
    Simple tidy-style tree layout:
    - Root at top center
    - Parent centered over its children
    - Depth increases downward
    Assumes G is a tree (event-based).
    """

    roots = [n for n in G.nodes if G.in_degree(n) == 0]
    if not roots:
        raise ValueError("No root found.")
    root = roots[0]

    pos = {}
    next_x = {}

    def assign(node, depth):
        children = list(G.successors(node))
        if depth not in next_x:
            next_x[depth] = 0.0

        if not children:
            x = next_x[depth]
            next_x[depth] += 1.0
            pos[node] = (x, -depth)
            return

        for c in children:
            assign(c, depth + 1)

        child_xs = [pos[c][0] for c in children]
        x = sum(child_xs) / len(child_xs)
        pos[node] = (x, -depth)

    assign(root, 0)

    # Center root at 0
    root_x, _ = pos[root]
    for n in pos:
        x, y = pos[n]
        pos[n] = (x - root_x, y)

    return pos


def plot_mutation_tree(
    G,
    savepath=None,
    show=True,
    detection_day=None,
    sequenced_only=False,
    collapse=False,
):
    """
    Plot event-based mutation tree in tidy-style layout.
    Nodes are not drawn except sequenced nodes.
    Branch length is proportional to mutation count.
    """

    from analysis.mutation_tree import extract_sequenced_subtree, collapse_clades

    if sequenced_only:
        G = extract_sequenced_subtree(G)

    if collapse:
        G = collapse_clades(G)

    pos = tidy_tree_layout(G)

    fig, ax = plt.subplots(figsize=(14, 12))

    for src, tgt, data in G.edges(data=True):
        if src not in pos or tgt not in pos:
            continue
        x1, y1 = pos[src]
        x2, y2 = pos[tgt]

        n_mut = data.get("n_mut", 1)
        lw = 0.5 + 0.1 * n_mut

        ax.plot([x1, x2], [y1, y2], color="black", linewidth=lw, alpha=0.8)

    seq_x = []
    seq_y = []
    for n in G.nodes:
        if G.nodes[n].get("sequenced", False):
            if n not in pos:
                continue
            x, y = pos[n]
            seq_x.append(x)
            seq_y.append(y)

    ax.scatter(seq_x, seq_y, color="red", s=20, alpha=0.9)

    ax.set_title("Event-based Tidy Mutation Tree (sequenced = red)", fontsize=14)

    if detection_day is not None:
        ax.text(
            0.5,
            1.02,
            f"Detection day: {detection_day}",
            transform=ax.transAxes,
            ha="center",
            fontsize=11,
            color="blue",
        )

    ax.axis("off")
    plt.tight_layout()

    if savepath:
        plt.savefig(savepath, dpi=200)

    if show:
        plt.show()
    else:
        plt.close(fig)
