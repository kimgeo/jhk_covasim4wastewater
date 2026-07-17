# analysis/mutation_tree.py

import networkx as nx
from collections import defaultdict


def build_mutation_tree(sim, daily_sequenced_agents):
    """
    Build mutation tree where each node is an infection event, not an agent.
    Node naming: "{agent_id}_{infected_day}".
    This guarantees a tree (each event has exactly one parent).
    """

    infection_log = getattr(sim, "infection_log", None)
    if infection_log is None:
        infection_log = getattr(sim.people, "infection_log", None)
    if infection_log is None:
        raise ValueError("Simulation has no infection_log.")

    sequenced_event_ids = set()
    for day_list in daily_sequenced_agents:
        sequenced_event_ids.update(day_list)

    G = nx.DiGraph()
    last_event_date = {}

    for entry in infection_log:
        src = entry["source"]
        tgt = entry["target"]
        date = entry["date"]

        # child event
        tgt_event = f"{tgt}_{date}"
        if tgt_event not in G:
            G.add_node(tgt_event)
            G.nodes[tgt_event]["sequenced"] = tgt_event in sequenced_event_ids
            G.nodes[tgt_event]["has_seq_descendant"] = G.nodes[tgt_event]["sequenced"]

        if src is None:
            src_event = f"None_0"
            if src_event not in G:
                G.add_node(src_event)
                G.nodes[src_event]["sequenced"] = False
                G.nodes[src_event]["has_seq_descendant"] = False
        else:
            src_date = last_event_date[src]
            src_event = f"{src}_{src_date}"

        G.add_edge(src_event, tgt_event)
        last_event_date[tgt] = date

        if G.nodes[tgt_event]["sequenced"]:
            for ancestor in nx.ancestors(G, tgt_event):
                G.nodes[ancestor]["has_seq_descendant"] = True
    

    return G


def extract_sequenced_subtree(G):
    """
    Extract subtree containing sequenced event nodes and all their ancestors.
    """

    seq_nodes = {n for n in G.nodes if G.nodes[n].get("sequenced", False)}
    if not seq_nodes:
        return G.copy()

    keep = set(seq_nodes)
    for n in seq_nodes:
        keep.update(nx.ancestors(G, n))

    SG = G.subgraph(keep).copy()

    roots = [n for n in SG.nodes if SG.in_degree(n) == 0]
    if len(roots) > 1:
        super_root = "super_root"
        SG.add_node(super_root, sequenced=False)
        for r in roots:
            SG.add_edge(super_root, r)

    return SG


def collapse_clades(G):
    """
    Collapse nodes that share identical mutation paths.
    Mutation path is defined as cumulative mutations from root.
    Works on event-based tree.
    """

    roots = [n for n in G.nodes if G.in_degree(n) == 0]
    if not roots:
        return G.copy()

    root = roots[0]

    def mutation_path(node):
        path = nx.shortest_path(G, root, node)
        muts = []
        for i in range(len(path) - 1):
            src, tgt = path[i], path[i + 1]
            muts.extend(G.edges[src, tgt]["mutations"])
        return tuple(muts)

    groups = {}
    for n in G.nodes:
        mp = mutation_path(n)
        groups.setdefault(mp, []).append(n)

    CG = nx.DiGraph()

    for mp, nodes in groups.items():
        rep = nodes[0]
        CG.add_node(rep)
        CG.nodes[rep]["mutation_path"] = mp
        CG.nodes[rep]["sequenced"] = any(G.nodes[n].get("sequenced", False) for n in nodes)

    for mp, nodes in groups.items():
        rep = nodes[0]
        for n in nodes:
            for child in G.successors(n):
                child_mp = mutation_path(child)
                child_rep = groups[child_mp][0]
                if rep != child_rep:
                    CG.add_edge(rep, child_rep)

    super_root = "super_root"
    CG.add_node(super_root, sequenced=False)
    for n in CG.nodes:
        if n == super_root:
            continue
        if CG.in_degree(n) == 0:
            CG.add_edge(super_root, n)

    return CG
