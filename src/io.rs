use anyhow::{bail, Result};
use crate::collections::{HashMap, HashSet};
use std::path;
use std::str;
use std::fs::File;
use std::io::{Read, BufReader,BufWriter, Write};
use flate2::read::MultiGzDecoder;

mod gfa;
mod gff;

//numbers convention for node_to_part array
pub const UNINITIALIZED: usize = usize::MAX;
pub const FILTERED: usize = usize::MAX - 1;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct SignedId {
    pub id: usize,
    pub plus: bool,
}

pub type Path = Vec<SignedId>;

pub struct NodeIndexer {
    pub map: HashMap<Vec<u8>, usize>,
    pub next: usize,
}
impl NodeIndexer {
    fn new() -> Self {
        Self {
            map: HashMap::default(),
            next: 0,
        }
    }
    #[inline]
    fn id_for(&mut self, s: &[u8]) -> usize {
        if let Some(&id) = self.map.get(s) {
            return id;
        }
        let id = self.next;
        self.next += 1;
        self.map.insert(s.to_vec(), id);
        id
    }
}

pub struct GenomeBundle {
    pub genomes: HashMap<String, PathBundle>,
    pub num_paths: usize,
    pub node_indexer: Option<NodeIndexer>,
}

pub struct PathBundle {
    pub paths: Vec<Path>,
    pub path_names: Vec<String>,
    pub path_starts: Vec<Vec<usize>>,
    pub path_ends: Vec<Vec<usize>>,
    pub path_sizes: Vec<usize>,
}

pub struct GraphBundle {
    pub graph: Vec<Vec<usize>>,
    pub num_nodes: usize,
    pub duplicates: HashSet<usize>,
}

pub struct PartitionBundle {
    pub node_to_part: Vec<usize>,
    pub num_parts: usize,
}

// ---------- Public API ----------
pub fn load_graph(
    input: &str,
    force_ext: Option<&str>,
    remove_duplicates: usize,
    group_by: bool,
    dirty: bool,
) -> Result<(GraphBundle, GenomeBundle, PartitionBundle)> {
    find_graph_type(input, force_ext)?.read_graph(input, remove_duplicates, group_by, dirty)
}

pub fn update_graph(
    genomes: &HashMap<String, PathBundle>,
    num_nodes: usize,
    node_to_part: &[usize],
) -> Vec<Vec<usize>> {
    gff::Gff.genomes_to_graph(genomes, num_nodes, node_to_part)
}

pub fn write_output(
    input: &str,
    force_ext: Option<&str>,
    out_dir: &path::Path,
    genomes: &HashMap<String, PathBundle>,
    node_to_part: &[usize],
) -> Result<()> {
    find_graph_type(input, force_ext)?.write_graph(out_dir, genomes, node_to_part)
}

pub fn write_paths(
    out_dir: &path::Path,
    genomes: &HashMap<String, PathBundle>,
    node_to_part: &[usize],
) -> Result<()> {
    let output = out_dir.join("paths.txt");
    std::fs::remove_file(&output).ok();
    let file = File::create(output)?;
    let mut writer = BufWriter::new(file);

    for (genome_name, genome) in genomes.iter() {
        for (i, path) in genome.paths.iter().enumerate() {
            let path_name = &genome.path_names[i];
            writeln!(writer, ">{}#{}", genome_name, path_name)?;
            let mut print_comma = false;
            for el in path {
                if node_to_part[el.id] == el.id {
                    if print_comma {
                        write!(writer, ",")?;
                    }
                    print_comma = true;
                    let part = node_to_part[el.id]+1;
                    let sign = if el.plus { '+' } else { '-' };
                    write!(writer, "{part}{sign}")?;
                }
            }
            writeln!(writer)?;
        }
    }

    Ok(())
}

pub fn write_partition(
    out_dir: &path::Path,
    num_nodes: usize,
    node_to_part: &[usize],
    node_indexer: Option<NodeIndexer>,
) -> Result<()> {
    let output = out_dir.join("partitions.txt");
    std::fs::remove_file(&output).ok();
    let file = File::create(output)?;
    let mut writer = BufWriter::new(file);

    let mut id_to_node_str = vec![Vec::new(); num_nodes];
    if let Some(node_indexer) = node_indexer {
        for (node_str, id) in node_indexer.map {
            id_to_node_str[id] = node_str;
        }
    } else {
        // default 1-based index
        id_to_node_str = (1..=num_nodes).map(|x| x.to_string().into_bytes()).collect();
    }

    let mut partition: Vec<Vec<usize>> = (0..num_nodes).map(|_| Vec::new()).collect();
    for (id, &part) in node_to_part.iter().enumerate() {
        if part != FILTERED && part < num_nodes {
            partition[part].push(id);
        }
    }

    for (core_id, part) in partition.iter().enumerate() {
        if !part.is_empty() {
            write!(writer, "{}:", str::from_utf8(&id_to_node_str[core_id]).unwrap())?;
            for &id in part {
                write!(writer, " {}", str::from_utf8(&id_to_node_str[id]).unwrap())?;
            }
            writeln!(writer)?;
        }
    }

    Ok(())
}

pub fn bufreader_from_compressed_file(file: &str) -> BufReader<Box<dyn Read>> {
    eprintln!("loading graph from {}", &file);
    let f = std::fs::File::open(file).unwrap_or_else(|err| panic!("Error opening file {} (err: {})", &file, err));
    let reader: Box<dyn Read> = if file.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    BufReader::new(reader)
}

//#[allow(dead_code)]
//pub fn load_paths(input: &str, group_by: bool, force_ext: Option<&str>) -> Result<(GenomeBundle, usize)> {
//    find_graph_type(input, force_ext)?.read_paths(input, group_by)
//}

// ---------- Private ----------
fn find_graph_type(input: &str, force_ext: Option<&str>) -> Result<Box<dyn GraphReader>> {
    let ext = force_ext
    .map(|s| s.to_ascii_lowercase())
    .or_else(|| {
        let p = std::path::Path::new(input);
        let last = p.extension()?.to_str()?.to_ascii_lowercase();
        if last == "gz" {
            let stem = p.file_stem()?.to_str()?.to_ascii_lowercase();
            stem.rsplit('.').next().map(|s| s.to_string())
        } else {
            Some(last)
        }
    });

    match ext.as_deref() {
        Some("gfa") => Ok(Box::new(gfa::Gfa)),
        Some("gff") => Ok(Box::new(gff::Gff)),
        Some(other) => bail!("Unsupported input extension: {}", other),
        None => bail!("Cannot infer input type (no extension). Pass force_ext."),
    }
}

// ---------- GraphReader ----------
#[inline]
fn add_edge_to_graph(
    u: usize,
    v: usize,
    graph: &mut [Vec<usize>],
    edge_set: &mut [HashSet<usize>],
) {
    if edge_set[u].insert(v) {
        edge_set[v].insert(u);
        graph[u].push(v);
        graph[v].push(u);
    }
}

pub trait GraphReader {
    fn read_paths(&self, input: &str, group_by: bool) -> Result<(GenomeBundle, usize)>;

    fn write_graph(
        &self,
        out_dir: &path::Path,
        genomes: &HashMap<String, PathBundle>,
        node_to_part: &[usize],
    ) -> Result<()>;

    fn genomes_to_graph(
        &self,
        genomes: &HashMap<String, PathBundle>,
        num_nodes: usize,
        node_to_part: &[usize],
    ) -> Vec<Vec<usize>> {
        let mut edge_set: Vec<HashSet<usize>> =
            (0..num_nodes * 2 + 1).map(|_| HashSet::default()).collect();
        let mut graph: Vec<Vec<usize>> = (0..num_nodes * 2 + 1).map(|_| Vec::new()).collect();
        let telomer = graph.len() - 1;

        for (_, genome) in genomes.iter() {
            for path in genome.paths.iter() {
                let mut i = 0;
                while i < path.len() && node_to_part[path[i].id] != path[i].id {
                    i += 1;
                }
                if i < path.len() {
                    let mut v_extremity = path[i].id * 2 + (!path[i].plus as usize);
                    add_edge_to_graph(telomer, v_extremity, &mut graph, &mut edge_set);
                    let mut u_extremity = path[i].id * 2 + (path[i].plus as usize);
                    while i + 1 < path.len() {
                        i += 1;
                        if node_to_part[path[i].id] != path[i].id {
                            continue;
                        }
                        v_extremity = path[i].id * 2 + (!path[i].plus as usize);
                        add_edge_to_graph(u_extremity, v_extremity, &mut graph, &mut edge_set);
                        u_extremity = path[i].id * 2 + (path[i].plus as usize);
                    }
                    add_edge_to_graph(u_extremity, telomer, &mut graph, &mut edge_set);
                }
            }
        }

        graph
    }

    fn get_genome_duplicates(
        &self,
        genomes: &HashMap<String, PathBundle>,
        remove_duplicates: usize,
    ) -> (HashSet<usize>, HashSet<usize>) {
        if remove_duplicates == 0 || remove_duplicates == 2 {
            let mut duplicates: HashSet<usize> = HashSet::default();
            for genome in genomes.values() {
                let mut duplicates_in_genome: HashSet<usize> = HashSet::default();
                for path in genome.paths.iter() {
                    for &sid in path.iter() {
                        if !duplicates_in_genome.insert(sid.id) {
                            duplicates.insert(sid.id);
                        }
                    }
                }
            }

            eprintln!("num duplicates:\t{}", duplicates.len());

            if remove_duplicates == 0 {
                (duplicates, HashSet::default())
            } else {
                eprintln!("num of filtered duplicates:\t{}", duplicates.len());
                (HashSet::default(), duplicates)
            }
        } else if remove_duplicates >= 3 {
            let mut duplicates: HashSet<usize> = HashSet::default();
            let mut duplicates_to_filter: HashSet<usize> = HashSet::default();
            for genome in genomes.values() {
                let mut count_in_genome: HashMap<usize, usize> = HashMap::default();
                for path in genome.paths.iter() {
                    for &sid in path.iter() {
                        *count_in_genome.entry(sid.id).or_insert(0) += 1;
                    }
                }
                for (&id, &count) in count_in_genome.iter() {
                    if count > 1 {
                        if count >= remove_duplicates {
                            duplicates_to_filter.insert(id);
                        } else {
                            duplicates.insert(id);
                        }
                    }
                }
            }
            eprintln!(
                "num duplicates:\t{}",
                duplicates.len() + duplicates_to_filter.len()
            );
            eprintln!("num of kept duplicates:\t{}", duplicates.len());
            eprintln!(
                "num of filtered duplicates:\t{}",
                duplicates_to_filter.len()
            );
            (duplicates, duplicates_to_filter)
        } else {
            unreachable!("Error, --remove-dup cannot be 1");
        }
    }

    fn read_graph(
        &self,
        input: &str,
        remove_duplicates: usize,
        group_by: bool,
        dirty: bool,
    ) -> Result<(GraphBundle, GenomeBundle, PartitionBundle)> {
        let (genome_bundle, num_nodes) = self.read_paths(input, group_by)?;

        let (mut duplicates, duplicates_to_filter) =
            self.get_genome_duplicates(&genome_bundle.genomes, remove_duplicates);

        if dirty {
            duplicates = HashSet::default();
        }

        let mut node_to_part: Vec<usize> = (0..num_nodes + 1).collect();
        let num_parts = 0usize;

        if remove_duplicates > 0 {
            for el in duplicates_to_filter {
                node_to_part[el] = FILTERED;
            }
        }

        let graph = self.genomes_to_graph(&genome_bundle.genomes, num_nodes, &node_to_part);

        Ok((
            GraphBundle {
                graph,
                num_nodes,
                duplicates,
            },
            genome_bundle,
            PartitionBundle {
            node_to_part,
            num_parts,
            }
        ))
    }
}
