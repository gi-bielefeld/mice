import copy
from pathlib import Path

def add_node_gfa(u_raw, graph):
    u, u_ori = u_raw[:-1], u_raw[-1]
    if u_raw != "T*" and u+"h" not in graph:
        graph[u+"t"] = set()
        graph[u+"h"] = set()
    if u_raw == "T*" and "T*" not in graph:
        graph["T*"] = set()

def add_edge_gfa(u_raw, v_raw, graph):
    add_node_gfa(u_raw, graph)
    add_node_gfa(v_raw, graph)
    u, u_ori = u_raw[:-1], u_raw[-1]
    v, v_ori = v_raw[:-1], v_raw[-1]
    if u_raw == "T*":
        v_ext = v+"t" if v_ori == "+" else v+"h"
        graph["T*"].add(v_ext)
        graph[v_ext].add("T*")
    elif v_raw == "T*":
        u_ext = u+"h" if u_ori == "+" else u+"t"
        graph["T*"].add(u_ext)
        graph[u_ext].add("T*")
    else:
        u_ext = u+"h" if u_ori == "+" else u+"t"
        v_ext = v+"t" if v_ori == "+" else v+"h"
        graph[u_ext].add(v_ext)
        graph[v_ext].add(u_ext)
    
def genomes_to_graph(genomes, telomer=True):
    graph = dict()
    for genome in genomes:
        for (u_raw, v_raw) in zip(genome[:-1],genome[1:]):
            add_edge_gfa(u_raw, v_raw, graph)
        if telomer:
            add_edge_gfa("T*", genome[0], graph)
            add_edge_gfa(genome[-1], "T*", graph)
    return graph

def get_repeats(genomes):
    repeats = set()

    for genome in genomes:
        genome_repeats = set()
        for el in genome:
            el = el[:-1]
            if el in genome_repeats:
                repeats.add(el)
            genome_repeats.add(el)
    return repeats

def parse_gfa_edges(gfa_path):
    graph = dict()
    with open(gfa_path, 'r') as fh:
        for line in fh:
            if not line or line[0] != 'L':
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            u = parts[1]
            u_ori = parts[2]
            v = parts[3]
            v_ori = parts[4]
            if u+"h" not in graph:
                graph[u+"t"] = set()
                graph[u+"h"] = set()
            if v+"h" not in graph:
                graph[v+"t"] = set()
                graph[v+"h"] = set()
            u_ext = u+"h" if u_ori == "+" else u+"t"
            v_ext = v+"h" if v_ori == "-" else v+"t"
            graph[u_ext].add(v_ext)
            graph[v_ext].add(u_ext)
    return graph

def parse_gfa_paths(gfa_path):
    """
    Read only P-lines from a GFA (v1) file.
    Returns: dict[name] = list like ['1+','2-','3+']
    """
    paths = {}
    with open(gfa_path, 'r') as fh:
        for line in fh:
            if not line or line[0] != 'P':
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[1]
            segs = parts[2].split(',')
            cleaned = [s for s in segs if s and (s.endswith('+') or s.endswith('-'))]
            if cleaned:
                paths[name] = cleaned
    return paths
    
def get_other_ext(node):
    if node[:-1] == "T":
        return node
    return node[:-1]+"h" if node.endswith("t") else node[:-1]+"t"
    
def remove_selfloops(graph):
    remove_node = False
    for node, adj in graph.items():
        for v in adj:
            if node[:-1] == v[:-1]:
                remove_node=True
        if remove_node:
            if node[:-1]+"t" in adj:
                adj.remove(node[:-1]+"t")
            if node[:-1]+"h" in adj:
                adj.remove(node[:-1]+"h")
    
def find_connected_components(graph):
    visited = set()
    components = []
    node_to_component_index = {}

    def dfs(node, component_index):
        visited.add(node)
        node_to_component_index[node] = component_index
        components[component_index].add(node)
        for neighbor in graph.get(node, []):
            if neighbor not in visited:
                dfs(neighbor, component_index)

    for node in graph:
        if node not in visited:
            components.append(set())
            dfs(node, len(components) - 1)

    return components, node_to_component_index

def graph_1_degree_contraction(graph, repeats=set()):
    graph_new = copy.deepcopy(graph)
    partition_tree = dict()
    for k in graph_new.keys():
        partition_tree[k[:-1]] = set()
    
    Q = []
    for k, v in graph_new.items():
        if len(v) == 1: 
            Q.append(k)

    #random.seed(11)
    #random.shuffle(Q)
    #    v--u==u2--w_1
    #            \-w_2
    #    v--w_1
    #     \-w_2
    i = 0
    while i < len(Q):
        u = Q[i]
        i += 1

        if len(graph_new[u]) != 1:
            continue 

        v = list(graph_new[u])[0]
        
        if u == "T*" or v == "T*":
            # telomer should not be merged
            continue 
            
        if v[:-1] in repeats or u[:-1] in repeats:
            #duplicates for now are being left alone
            continue
        
        partition_tree[v[:-1]].add(u[:-1])
        partition_tree[u[:-1]].add(v[:-1])
        
        u2 = get_other_ext(u)
        for w in graph_new[u2]:
            graph_new[w].remove(u2)
            graph_new[w].add(v)
            
            if len(graph_new[w]) == 1:
                Q.append(w)

        graph_new[v].remove(u)
        graph_new[v].update(graph_new[u2])
        
        if len(graph_new[v]) == 1:
            Q.append(v)
        
        #Remove u from graph_new
        graph_new[u] = set()
        graph_new[u2] = set()
        
    (partition, node_to_part) = find_connected_components(partition_tree)
    return (partition, node_to_part, graph_new)

def parse_gfa_paths(gfa_path: Path):
    """Extract only the genome paths (P-lines)."""
    genomes = []
    with gfa_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("P\t"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    segs = fields[2]
                    if segs == "*" or segs == "":
                        genomes.append([])
                    else:
                        genomes.append(segs.split(","))
    return genomes

def main():
    gfa_dir = Path("gfa")
    out_dir = Path("size")
    out_dir.mkdir(parents=True, exist_ok=True)

    for gfa_file in sorted(gfa_dir.glob("*.gfa")):
        genomes = parse_gfa_paths(gfa_file)
        graph = genomes_to_graph(genomes)
        repeats = get_repeats(genomes)
        (partition, node_to_part, new_graph) = graph_1_degree_contraction(graph, repeats)

        out_file = out_dir / (gfa_file.stem + ".txt")
        with out_file.open("w") as out:
            # The -1 is to remove the telomer
            out.write(str(len(partition)-1) + "\n")
    print(f"Expected outputs written to: {out_dir.resolve()}")

if __name__ == "__main__":
    main()

