# analysis/mutation_tree.py

import networkx as nx
from collections import defaultdict


def build_mutation_tree(sim, sequenced_agents=None):
    """
    Build mutation tree where each node is an infection event, not an agent.
    Node naming: "{agent_id}_{event_index}".
    This guarantees a tree (each event has exactly one parent).
    """

    infection_log = getattr(sim, "infection_log", None)
    if infection_log is None:
        infection_log = getattr(sim.people, "infection_log", None)
    if infection_log is None:
        raise ValueError("Simulation has no infection_log.")

    sequenced_agents = set(sequenced_agents or [])

    G = nx.DiGraph()

    # Track latest event node per agent
    last_event_node = {}
    # Track event index per agent
    event_counts = defaultdict(int)

    for entry in infection_log:
        src = entry["source"]
        tgt = entry["target"]

        # Create event node for target
        event_counts[tgt] += 1
        tgt_event = f"{tgt}_{event_counts[tgt]}"

        if tgt_event not in G:
            G.add_node(tgt_event)

        # Mark sequenced if agent is sequenced
        G.nodes[tgt_event]["sequenced"] = tgt in sequenced_agents

        # Variant info if available
        if hasattr(sim.people, "variant"):
            G.nodes[tgt_event]["variant"] = sim.people.variant[tgt]

        # Determine parent event node
        if src is None:
            # Synthetic root event for this target
            src_event = f"root_{tgt_event}"
            if src_event not in G:
                G.add_node(src_event)
                G.nodes[src_event]["sequenced"] = False
                G.nodes[src_event]["variant"] = "root"
        else:
            # Use latest event of source; if none, create first event
            if src not in event_counts:
                event_counts[src] += 1
                src_event = f"{src}_{event_counts[src]}"
                if src_event not in G:
                    G.add_node(src_event)
                    G.nodes[src_event]["sequenced"] = src in sequenced_agents
                    if hasattr(sim.people, "variant"):
                        G.nodes[src_event]["variant"] = sim.people.variant[src]
            else:
                src_event = last_event_node.get(src)
                if src_event is None:
                    event_counts[src] += 1
                    src_event = f"{src}_{event_counts[src]}"
                    if src_event not in G:
                        G.add_node(src_event)
                        G.nodes[src_event]["sequenced"] = src in sequenced_agents
                        if hasattr(sim.people, "variant"):
                            G.nodes[src_event]["variant"] = sim.people.variant[src]

        # Update latest event for source and target
        last_event_node[src] = src_event if src is not None else src_event
        last_event_node[tgt] = tgt_event

        raw_muts = entry.get("branch_mutations", [])
        n_mut = entry.get("n_mutations", len(raw_muts))

        muts = []
        for m in raw_muts:
            if isinstance(m, tuple):
                muts.append("_".join(map(str, m)))
            else:
                muts.append(str(m))

        G.add_edge(src_event, tgt_event, mutations=muts, n_mut=n_mut)

    # Ensure single root via super_root
    roots = [n for n in G.nodes if G.in_degree(n) == 0]
    if len(roots) > 1:
        super_root = "super_root"
        G.add_node(super_root, sequenced=False, variant="root")
        for r in roots:
            G.add_edge(super_root, r, mutations=[], n_mut=0)

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
