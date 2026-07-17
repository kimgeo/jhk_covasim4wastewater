# plotting/plot_mutation_tree.py

import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

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

    pos = graphviz_layout(G, prog="dot")

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

    ax.set_title("Event-based Infection Tree", fontsize=14)

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

def spread_levels(G, pos):
    depths = {}
    root = "None_0"
    for node in G.nodes:
        try:
            depths[node] = nx.shortest_path_length(G, root, node)
        except:
            depths[node] = 0
    level_nodes = {}
    for node, d in depths.items():
        level_nodes.setdefault(d, []).append(node)

    new_pos = {}
    for d, nodes in level_nodes.items():
        for i, node in enumerate(nodes):
            x, y = pos[node]
            new_pos[node] = (i, -d)

    return new_pos


def plot_erase_mutation_tree(
    G,
    savepath=None,
    show=True,
    detection_day=None,
    erase=True,          # True = edge -> delete, False = edge -> transparent
    sequenced_only=False,
    collapse=False,
):
    """
    Plot mutation tree but erase (or fade) edges whose descendants
    do NOT include any sequenced nodes.
    """

    from analysis.mutation_tree import extract_sequenced_subtree, collapse_clades

    # Optionally restrict to sequenced subtree
    if sequenced_only:
        G = extract_sequenced_subtree(G)

    if collapse:
        G = collapse_clades(G)

    pos = graphviz_layout(G, prog="dot")

    fig, ax = plt.subplots(figsize=(14, 12))

    subtree_has_seq = {}
    for n in G.nodes:
        desc = nx.descendants(G, n)
        subtree_has_seq[n] = (
            G.nodes[n].get("sequenced", False)
            or any(G.nodes[d].get("sequenced", False) for d in desc)
        )

    # Draw edges
    for src, tgt, data in G.edges(data=True):
        if src not in pos or tgt not in pos:
            continue

        # If src has no sequenced descendants → erase or fade
        if not subtree_has_seq[tgt]:
            if erase:
                continue
            alpha = 0.4
            color = "black"
        else:
            alpha = 1.0
            color = "red"

        x1, y1 = pos[src]
        x2, y2 = pos[tgt]

        n_mut = data.get("n_mut", 1)
        lw = 0.5 + 0.1 * n_mut

        ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, alpha=alpha)

    # Draw sequenced nodes
    seq_x = []
    seq_y = []
    mut_x = []
    mut_y = []

    for n in G.nodes:
        if G.nodes[n].get("sequenced", False):
            if G.nodes[n].get("mutations", False):
                if n not in pos:
                    continue
                x, y = pos[n]
                mut_x.append(x)
                mut_y.append(y)
            else:
                if n not in pos:
                    continue
                x, y = pos[n]
                seq_x.append(x)
                seq_y.append(y)

    ax.scatter(seq_x, seq_y, color="red", s=20, alpha=0.9)
    ax.scatter(mut_x, mut_y, color="blue", s=20, alpha=0.9)

    ax.set_title("Event-based Infection Tree (sequenced = red)", fontsize=14)

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
