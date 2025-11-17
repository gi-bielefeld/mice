use crate::io::*;
use anyhow::Result;
use std::collections::hash_map;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str;

pub struct Gfa;

impl GraphReader for Gfa {
    fn read_paths(&self, gfa_input: &str, group_by: bool) -> Result<(GenomeBundle, usize)> {
        let (genome_bundle, num_nodes) = Self::parse_gfa_paths(gfa_input, group_by)
            .map_err(|e| anyhow::anyhow!("Error reading GFA: {}", e))?;

        Ok((genome_bundle, num_nodes))
    }

    #[allow(unused_variables)]
    //TODO
    fn write_graph(
        &self,
        out_dir: &path::Path,
        genomes: &HashMap<String, PathBundle>,
        node_to_part: &[usize],
    ) -> Result<()> {
        Ok(())
    }
}

impl Gfa {
    #[inline]
    fn path_name_to_genome_string(path_name: String) -> String {
        match path_name.find("#") {
            Some(idx) => path_name[..idx].to_string(),
            None => path_name,
        }
    }

    #[inline]
    fn parse_path_node(node: &[u8], node_indexer: &mut NodeIndexer) -> SignedId {
        let orient = *node.last().unwrap();
        let name = &node[..node.len() - 1];
        let id = node_indexer.id_for(name);
        SignedId {
            id,
            plus: orient == b'+',
        }
    }

    fn parse_path_seq_to_signed_id_vec(
        data: &[u8],
        node_indexer: &mut NodeIndexer,
    ) -> Result<(String, Path)> {
        let mut it = data.iter();
        let mut start = it.position(|&x| x == b'\t').unwrap() + 1; //skip P
        let end_name = it.position(|&x| x == b'\t').unwrap() + 1;

        let path_name = str::from_utf8(&data[start..end_name + 1])?.to_owned();
        start += end_name; //skip path_name

        let data = &data[start..];
        let end = it
            //.position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
            .rposition(|x| x == &b'\t')
            .unwrap();

        let path: Vec<SignedId> = data[..end]
            .split(|&x| x == b',')
            .map(|node| Self::parse_path_node(node, node_indexer))
            .collect();

        Ok((path_name, path))
    }

    pub fn parse_gfa_paths(filename: &str, group_by: bool) -> Result<(GenomeBundle, usize)> {
        let file = File::open(filename)?;
        let mut reader = BufReader::new(file);

        let mut genomes: HashMap<String, PathBundle> = HashMap::default();
        let mut num_paths = 0usize;
        let mut node_indexer = NodeIndexer::new();

        let mut buf = vec![];
        while reader.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
            if buf[0] == b'P' {
                let (path_name, path) =
                    Self::parse_path_seq_to_signed_id_vec(&buf, &mut node_indexer)?;
                let mut genome_name = path_name.clone();
                if group_by {
                    genome_name = Self::path_name_to_genome_string(genome_name);
                }

                let entry = genomes.entry(genome_name);
                match entry {
                    hash_map::Entry::Occupied(mut p) => {
                        p.get_mut().paths.push(path);
                        p.get_mut().path_names.push(path_name);
                    }
                    hash_map::Entry::Vacant(p) => {
                        p.insert(PathBundle {
                            paths: vec![path],
                            path_names: vec![path_name],
                            path_starts: Vec::new(),
                            path_ends: Vec::new(),
                            path_sizes: Vec::new(),
                        });
                    }
                }
                num_paths += 1;
            }
            buf.clear();
        }

        let num_nodes = node_indexer.next;
        let node_indexer = Some(node_indexer);
        let genome_bundle = GenomeBundle {
            genomes,
            num_paths,
            node_indexer,
        };
        Ok((genome_bundle, num_nodes))
    }
}
