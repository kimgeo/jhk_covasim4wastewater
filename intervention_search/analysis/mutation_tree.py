# analysis/mutation_tree.py

import networkx as nx
import covasim as cv
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

def build_infection_seq_tree(sim, daily_sequenced_agents):
    """
    Build a mutation-collapsed infection sequence tree.
    Each node represents a unique mutation state (haplotype).
    Each edge represents a mutation event.
    Also tracks:
        - sequenced: whether any infection event in this mutation state was sequenced
        - has_seq_descendant: whether any descendant mutation state has a sequenced event
    """

    # infection log
    infection_log = getattr(sim, "infection_log", None)
    if infection_log is None:
        infection_log = getattr(sim.people, "infection_log", None)
    if infection_log is None:
        raise ValueError("Simulation has no infection_log.")

    # sequenced event IDs
    sequenced_event_ids = set()
    for day_list in daily_sequenced_agents:
        sequenced_event_ids.update(day_list)

    # mutation tree structures
    mutation_nodes = {}   # mutation_state -> node_id
    node_events = {}      # node_id -> list of event_ids
    G = nx.DiGraph()

    # root mutation state
    root_state = frozenset()
    mutation_nodes[root_state] = "Node_0"
    node_events["Node_0"] = []
    G.add_node("Node_0", mutations=root_state,
               sequenced=False,
               has_seq_descendant=False)

    # agent -> current mutation state
    current_state = {}

    for entry in infection_log:
        if entry["source"] is None:
            tgt = entry["target"]
            date = entry["date"]
            event_id = f"{tgt}_{date}"

            # initial pop_infected is root mutation state
            current_state[tgt] = root_state
            node_events["Node_0"].append(event_id)

            if event_id in sequenced_event_ids:
                G.nodes["Node_0"]["sequenced"] = True
                G.nodes["Node_0"]["has_seq_descendant"] = True

    for entry in infection_log:
        src = entry["source"]
        tgt = entry["target"]
        date = entry["date"]
        event_id = f"{tgt}_{date}"

        branch_mut = tuple(entry.get("branch_mutations", []))

        # parent mutation state
        if src is None:
            parent_state = root_state
        else:
            parent_state = current_state[src]

        if len(branch_mut) == 0:
            current_state[tgt] = parent_state
            node_id = mutation_nodes[parent_state]
            node_events[node_id].append(event_id)

            if event_id in sequenced_event_ids:
                G.nodes[node_id]["sequenced"] = True
                G.nodes[node_id]["has_seq_descendant"] = True
                for ancestor in nx.ancestors(G, node_id):
                    G.nodes[ancestor]["has_seq_descendant"] = True

            continue

        # if mutated -> new mutation state
        new_state = frozenset(parent_state | set(branch_mut))

        # if new mutation state -> new node
        if new_state not in mutation_nodes:
            new_node_id = f"Node_{len(mutation_nodes)}"
            mutation_nodes[new_state] = new_node_id
            node_events[new_node_id] = []

            G.add_node(new_node_id,
                       mutations=new_state,
                       sequenced=False,
                       has_seq_descendant=False)

            # parent -> new_state edge
            parent_node_id = mutation_nodes[parent_state]
            G.add_edge(parent_node_id, new_node_id)

        # agent state update
        current_state[tgt] = new_state

        # event
        node_id = mutation_nodes[new_state]
        node_events[node_id].append(event_id)

        # sequenced flag
        if event_id in sequenced_event_ids:
            G.nodes[node_id]["sequenced"] = True
            G.nodes[node_id]["has_seq_descendant"] = True
            for ancestor in nx.ancestors(G, node_id):
                G.nodes[ancestor]["has_seq_descendant"] = True

    return G, mutation_nodes, node_events

def mutation_tree_to_newick(G, root="Node_0"):

    def dfs(node):
        children = list(G.successors(node))
        label = node

        # annotation
        ann = []
        if G.nodes[node].get("sequenced", False):
            ann.append("sequenced=true")
        if G.nodes[node].get("has_seq_descendant", False):
            ann.append("has_seq_descendant=true")

        if ann:
            label = f"{label}[&{','.join(ann)}]"

        if not children:
            return label

        child_str = ",".join(dfs(c) for c in children)
        return f"({child_str}){label}"

    return dfs(root) + ";"

def extract_sequenced_subtree(G):
    """
    Return a pruned mutation tree containing only nodes
    whose subtree has at least one sequenced event.
    """

    H = G.copy()

    to_remove = [n for n in H.nodes
                 if not H.nodes[n].get("has_seq_descendant", False)]

    H.remove_nodes_from(to_remove)

    return H


def reconstruct_haplotype_from_mutation_state(sim, mutation_state):
    """
    mutation_state: frozenset of mutation IDs (integers)
    returns: nucleotide string
    """
    # reference genome
    ref = sim.sequence_tracker.reference.copy()

    # apply SNPs
    for site, ref_nt, alt_nt in mutation_state:
        ref[site] = alt_nt

    return cv.decode_sequence(ref)

def export_fasta_from_G(G_prune, sim, filepath):
    """
    Write FASTA file where each record corresponds to a mutation-state node in G_prune.
    """

    with open(filepath, "w") as f:
        for node in G_prune.nodes:
            mut_state = G_prune.nodes[node]["mutations"]
            hap = reconstruct_haplotype_from_mutation_state(sim, mut_state)

            f.write(f">{node}\n")
            f.write(hap + "\n")

def infer_tree_jc69(fasta_path):
    """
    Run IQ-TREE with JC69 model.
    Output: fasta_path.treefile
    """
    import subprocess

    cmd = [
        "iqtree",
        "-s", fasta_path,
        "-m", "JC69",
    ]
    subprocess.run(cmd, check=True)

    return fasta_path + ".treefile"

