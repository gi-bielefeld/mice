use crate::io::*;
use std::fs;
use anyhow::Result;
use crate::collections::{HashSet, HashMap};
use std::path;
use crate::cli::Cli;

pub fn run_mice(
    args: &Cli,
) -> Result<()> {
    // Args
    let out_dir = path::Path::new(&args.out_dir);
    if !out_dir.exists() {
        fs::create_dir_all(out_dir)?;
    }

    let graph_input = &args.graph_input;
    let remove_duplicates = args.remove_duplicates;
    let group_by = args.group_by;
    let min_size = args.min_size;
    let dirty = args.dirty;
    let force_ext = None;

    // Lets go
    let (graph_bundle, path_bundle, partition_bundle) = load_graph(graph_input, force_ext, remove_duplicates, group_by, dirty)?;

    let GraphBundle { mut graph, num_nodes, duplicates, } = graph_bundle;
    let GenomeBundle { genomes, num_paths, node_indexer, } = path_bundle;

    let mut partition_bundle = compress_graph(&mut graph, num_nodes, partition_bundle, &duplicates);

    if min_size != 0 {
        filter_min_size(&mut partition_bundle.node_to_part, num_nodes, &genomes, min_size);
        graph = update_graph(&genomes, num_nodes, &partition_bundle.node_to_part);
        partition_bundle = compress_graph(&mut graph, num_nodes, partition_bundle, &duplicates);
    }

    eprintln!("num genomes:\t{}", genomes.len());
    eprintln!("num paths:\t{}", num_paths);
    eprintln!("num nodes:\t{}", num_nodes);
    eprintln!("num partitions:\t{}", partition_bundle.num_parts);
    eprintln!("ratio:\t\t{:.2}", partition_bundle.num_parts as f64 / num_nodes as f64);

    write_paths(out_dir, &genomes, &partition_bundle.node_to_part)?;
    write_partition(out_dir, num_nodes, &partition_bundle.node_to_part, node_indexer)?;

    write_output(graph_input, force_ext, out_dir, &genomes, &partition_bundle.node_to_part)?;

    Ok(())
}

fn filter_min_size(
    node_to_part: &mut [usize], 
    num_nodes: usize, 
    genomes: &HashMap<String, PathBundle>,
    min_size: usize, ) {
    let mut part_counts = vec![0usize; num_nodes];
    for i in 0..num_nodes {
        if node_to_part[i] != FILTERED {
            part_counts[node_to_part[i]] += 1;
        }
    }
    for genome in genomes.values() {
        for (i, path) in genome.paths.iter().enumerate() {
            let path_starts = &genome.path_starts[i];
            let path_ends = &genome.path_ends[i];
            for (j, el) in path.iter().enumerate() {
                //TODO this works for GFF but GFA has no path_ends or path_starts
                if part_counts[el.id] == 1 && (path_ends[j] - path_starts[j] + 1) < min_size {
                    node_to_part[el.id] = FILTERED;
                }
            }
        }
    }
}

fn connected_components(mut node_to_part: Vec<usize>, num_nodes: usize) -> PartitionBundle {
    let mut num_parts = 0usize;

    let mut i = 0usize;
    while i < num_nodes {
        if node_to_part[i] != FILTERED {
            if node_to_part[i] == i {
                num_parts += 1;
            } else {
                let mut j = i;
                let mut stack = Vec::new();
                while node_to_part[node_to_part[j]] != node_to_part[j] {
                    stack.push(j);
                    j = node_to_part[j];
                }

                while let Some(z) = stack.pop() {
                    node_to_part[z] = node_to_part[j];
                }
            }
        }
        i += 1;
    }

    PartitionBundle {
        node_to_part,
        num_parts,
    }
}

// unifying (u_ext, v_ext)
// before:
// v_ext--u_ext==u2_ext--w_1_ext
//                     \-w_2_ext
// after:
// v_ext--w_1_ext
//      \-w_2_ext
pub fn compress_graph(
    graph: &mut [Vec<usize>],
    num_nodes: usize,
    partition_bundle: PartitionBundle,
    duplicates: &HashSet<usize>,
) -> PartitionBundle {
    let mut node_to_part = partition_bundle.node_to_part;
    let mut next_adj = vec![0; 2 * num_nodes + 1];
    let mut queue = Vec::new();
    let telomer = graph.len() - 1;

    for (u_ext, adj) in graph.iter().enumerate() {
        let u = u_ext >> 1;
        if adj.len() == 1 && !duplicates.contains(&u) && u != telomer && u != (adj[0]>>1) {
            queue.push(u_ext);
        }
    }

    while let Some(u_ext) = queue.pop() {
        let u = u_ext >> 1;
        if node_to_part[u] != u || next_adj[u_ext] != graph[u_ext].len() - 1 {
            continue;
        } // Already merged

        let v_ext = graph[u_ext][next_adj[u_ext]];
        let u2_ext = u_ext ^ 1usize;

        if u_ext == telomer || v_ext == telomer {
            continue;
        } // No telomer

        let v = v_ext >> 1;
        if duplicates.contains(&u) || duplicates.contains(&v) {
            continue;
        } // No duplicates

        if node_to_part[v] == u {
            continue;
        } // Already merged to u

        node_to_part[u] = v;

        let adj_u2 = std::mem::take(&mut graph[u2_ext]);
        let adj_u2 = &adj_u2[next_adj[u2_ext]..adj_u2.len()];

        for &w_ext in adj_u2.iter() {
            let w = w_ext >> 1;
            if node_to_part[w] == w {
                graph[w_ext].push(v_ext);
                graph[v_ext].push(w_ext);
                if is_adj_size_one(&mut graph[w_ext], &mut next_adj[w_ext], &node_to_part) && w != (graph[w_ext][next_adj[w_ext]] >> 1) {
                    queue.push(w_ext);
                }
            } else if w_ext == u2_ext { //add self loop
                graph[v_ext].push(v_ext);
            } 
        }

        if is_adj_size_one(&mut graph[v_ext], &mut next_adj[v_ext], &node_to_part) && v != (graph[v_ext][next_adj[v_ext]] >> 1) {
            queue.push(v_ext);
        }

        // Remove u from graph
        graph[u_ext].clear();
    }

    connected_components(node_to_part, num_nodes)
}

fn is_adj_size_one(adj: &mut [usize], next: &mut usize, node_to_part: &[usize]) -> bool {
    let mut i = *next;
    //Find first not merged
    while i < adj.len() && (node_to_part[adj[i] >> 1]) != (adj[i] >> 1) {
        i += 1;
    }
    if i < adj.len() {
        //If found
        let mut j = i + 1;
        while j < adj.len() && (((node_to_part[adj[j] >> 1]) != (adj[j] >> 1)) || adj[j] == adj[i]) {
            j += 1;
        }

        *next = j - 1;
        adj[*next] = adj[i];

        if *next == adj.len() - 1 {
            return true; //Degree 1
        }
    }

    false // Degree != 1
}

#[allow(dead_code)]
pub fn run_mice_test(input: &str, force_ext: Option<&str>) -> Result<usize> {
    let remove_duplicates = 0usize;
    let group_by = false;
    let dirty = false;
    let (graph_bundle, genome_bundle, partition_bundle) = load_graph(input, force_ext, remove_duplicates, group_by, dirty)?;
    let GraphBundle { mut graph, num_nodes, duplicates, } = graph_bundle;

    let partition_bundle = compress_graph(&mut graph, num_nodes, partition_bundle, &duplicates);

    eprintln!("num genomes:\t{}", genome_bundle.genomes.len());
    eprintln!("num paths:\t{}", genome_bundle.num_paths);
    eprintln!("num nodes: {}", num_nodes);
    eprintln!("num partitions: {}", partition_bundle.num_parts);
    eprintln!(
        "ratio: {:.2}",
        partition_bundle.num_parts as f64 / num_nodes as f64
    );

    Ok(partition_bundle.num_parts)
}

