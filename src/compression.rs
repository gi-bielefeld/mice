use crate::io::*;
use std::{fs};
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
    let quorum = args.quorum;
    let group_by = args.group_by;
    let min_size = args.min_size;
    let merge_with_dup = args.merge_with_dup;
    let force_ext = None;

    // Lets go
    let (graph_bundle, path_bundle, partition_bundle) =
        load_graph(graph_input, force_ext, remove_duplicates, quorum, group_by, merge_with_dup)?;

    let GraphBundle { mut graph, num_nodes, duplicates, counts} = graph_bundle;
    let GenomeBundle { genomes, num_paths, node_indexer, } = path_bundle;

    let mut partition_bundle = compress_graph_linear(&mut graph, num_nodes, partition_bundle, &duplicates, &counts);

    if min_size != 0 {
        filter_min_size(&mut partition_bundle.node_to_part, num_nodes, &genomes, min_size);
        graph = update_graph(&genomes, num_nodes, &partition_bundle.node_to_part);
        partition_bundle = compress_graph_linear(&mut graph, num_nodes, partition_bundle, &duplicates, &counts);
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

fn extend_side(
    extremity: usize,
    graph: &mut [Vec<usize>],
    node_to_part: &mut [usize],
    influence: &mut [usize],
    deg: &mut [usize],
    duplicates: &HashSet<usize>,
) {
    let element = extremity >> 1;
    let mut queue: Vec<usize> = [extremity].to_vec();
    let telomere = graph.len()-1;

    while let Some(ux) = queue.pop() {
        let u = ux >> 1;
        for vy in &graph[ux] {
            let v = vy >> 1;
            if vy == &telomere || v == u || node_to_part[v] != v || duplicates.contains(&v) {
                continue;
            }
            if influence[*vy] == UNINITIALIZED || influence[*vy] == extremity {
                influence[*vy] = extremity;
                deg[*vy] -= 1;
                if deg[*vy] == 0 {
                    node_to_part[v] = element;
                    let vx = vy ^ 1usize;
                    queue.push(vx);
                }
            }
        }
    }
}

pub fn compress_graph_linear(
    graph: &mut [Vec<usize>],
    num_nodes: usize,
    partition_bundle: PartitionBundle,
    duplicates: &HashSet<usize>,
    counts: &Vec<usize>,
) -> PartitionBundle {
    let mut node_to_part = partition_bundle.node_to_part;

    let mut counts_dict: HashMap<usize,usize> = HashMap::default();

    // TODO i can prob do this cleaner, in fewer steps, with some clever use of iter, map, collect
    for i in 0..counts.len() {
        counts_dict.insert(i, counts[i]);
    }
    let mut sorted_vec: Vec<(&usize, &usize)> = counts_dict.iter().collect();
    sorted_vec.sort_by(|a,b| b.1.cmp(a.1));

    let mut influence: Vec<usize> = vec![0; 2*num_nodes+1];
    let mut deg: Vec<usize> = vec![0; 2*num_nodes+1];

    for node in 0..graph.len() {
        influence[node] = UNINITIALIZED;
        deg[node] = graph[node].len();
    }

    for (element, _) in sorted_vec {
        if node_to_part[*element] != *element || duplicates.contains(element) {
            continue;
        }
        let elem_tail = 2*element;
        let elem_head = 2*element+1;

        extend_side(elem_tail, graph, &mut node_to_part, &mut influence, & mut deg, duplicates);
        extend_side(elem_head, graph, &mut node_to_part, &mut influence, & mut deg, duplicates);
    }

    let mut num_parts = 0usize;
    let mut i = 0usize;
    // TODO does this need to be a while loop?
    while i < num_nodes {
        if node_to_part[i] != FILTERED {
            if node_to_part[i] == i {
                num_parts += 1;
            }
        }
        i += 1;
    }

    PartitionBundle {
        node_to_part,
        num_parts,
    }
}

#[allow(dead_code)]
pub fn run_mice_hardcoded(input: &str, out_dir: &str, force_ext: Option<&str>) -> Result<usize> {
    let remove_duplicates = 0usize;
    let quorum = 0usize;
    let group_by = true;
    let merge_with_dup = false;
    let (graph_bundle, genome_bundle, partition_bundle) = load_graph(input, force_ext, remove_duplicates, quorum, group_by, merge_with_dup)?;
    let GraphBundle { mut graph, num_nodes, duplicates, counts, } = graph_bundle;
    let GenomeBundle { genomes, num_paths, node_indexer, } = genome_bundle;

    let partition_bundle = compress_graph_linear(&mut graph, num_nodes, partition_bundle, &duplicates, &counts);

    eprintln!("num genomes:\t{}", genomes.len());
    eprintln!("num paths:\t{}", num_paths);
    eprintln!("num nodes: {}", num_nodes);
    eprintln!("num partitions: {}", partition_bundle.num_parts);
    eprintln!(
        "ratio: {:.2}",
        partition_bundle.num_parts as f64 / num_nodes as f64
    );

    let out_dir = path::Path::new(out_dir);
    if !out_dir.exists() {
        fs::create_dir_all(out_dir)?;
    }

    write_paths(out_dir, &genomes, &partition_bundle.node_to_part)?;
    eprintln!("paths written");
    write_partition(out_dir, num_nodes, &partition_bundle.node_to_part, node_indexer)?;
    eprintln!("partitions written");
    write_output(input, force_ext, out_dir, &genomes, &partition_bundle.node_to_part)?;
    eprintln!("output written");
    Ok(partition_bundle.num_parts)

}


#[allow(dead_code)]
pub fn run_mice_test(input: &str, force_ext: Option<&str>) -> Result<usize> {
    let remove_duplicates = 0usize;
    let quorum = 0usize;
    let group_by = false;
    let merge_with_dup = false;
    let (graph_bundle, genome_bundle, partition_bundle) = load_graph(input, force_ext, remove_duplicates, quorum, group_by, merge_with_dup)?;
    let GraphBundle { mut graph, num_nodes, duplicates, counts, } = graph_bundle;

    // let partition_bundle = compress_graph(&mut graph, num_nodes, partition_bundle, &duplicates);
    let partition_bundle = compress_graph_linear(&mut graph, num_nodes, partition_bundle, &duplicates, &counts);

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

