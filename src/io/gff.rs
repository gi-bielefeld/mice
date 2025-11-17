use crate::io::*;
use anyhow::Result;
use std::collections::hash_map;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str;

//The Gff assumed here are from SibeliaZ and Cactus with a ID as attributes and a 1-based index
//example:
//Genome1.Chr2	SibeliaZ	SO:0000856	524938	525114	.	+	.	ID=1
//Genome1.Chr1	SibeliaZ	SO:0000856	521833	522009	.	+	.	ID=2
//Genome2.Chr4	SibeliaZ	SO:0000856	536437	536613	.	+	.	ID=2
pub struct Gff;

impl GraphReader for Gff {
    fn read_paths(&self, gfa_input: &str, group_by: bool) -> Result<(GenomeBundle, usize)> {
        let (genome_bundle, num_nodes) = Self::parse_gff_paths(gfa_input, group_by)
            .map_err(|e| anyhow::anyhow!("Error reading GFF: {}", e))?;

        Ok((genome_bundle, num_nodes))
    }

    fn write_graph(
        &self,
        out_dir: &path::Path,
        genomes: &HashMap<String, PathBundle>,
        node_to_part: &[usize],
    ) -> Result<()> {
        let output = out_dir.join("output.gff");
        std::fs::remove_file(&output).ok();
        let file = File::create(output)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "##gff-version 3")?;

        //Header
        for genome in genomes.values() {
            for (path_size, path_name) in genome.path_sizes.iter().zip(genome.path_names.iter()) {
                writeln!(writer, "##sequence-region {path_name} 1 {path_size}")?;
            }
        }
        //Gff
        for (genome_name, genome) in genomes.iter() {
            for (z, path) in genome.paths.iter().enumerate() {
                let path_name = &genome.path_names[z];
                let path_starts = &genome.path_starts[z];
                let path_ends = &genome.path_ends[z];
                let mut i = 0usize;
                while i < path.len() {
                    while i < path.len() && node_to_part[path[i].id] == FILTERED {
                        i += 1;
                    }
                    if i == path.len() {
                        break;
                    }

                    let part_i = node_to_part[path[i].id]; //First non-FILTERED
                    let mut j = i;
                    while j + 1 < path.len()
                        && (part_i == node_to_part[path[j + 1].id]
                            || node_to_part[path[j + 1].id] == FILTERED)
                    {
                        j += 1;
                    }

                    let end = j;
                    while node_to_part[path[j].id] == FILTERED {
                        j -= 1;
                    }
                    // This is true for the clean mode, but not for the dirty
                    // This is the situation:
                    //               ...i            j      ...
                    // node_to_part: ...pi pi F F pi pi F px...  (pi = part_i, px!=pi)

                    //There are two cases:
                    // node_to_part: ...pi pi F F pi pi...  (pi = part_i)
                    // core:         ...0  1  0 0 0  0 ...
                    // or
                    // core:         ...1  1  0 0 1  1 ...

                    //Count number of core elements (either 1 or many)
                    let mut count_core = 0usize;
                    let mut last_el_core = SignedId {
                        id: UNINITIALIZED,
                        plus: false,
                    };
                    for l in i..=j {
                        if node_to_part[path[l].id] == path[l].id {
                            count_core += 1;
                            last_el_core = path[l];
                        }
                    }

                    match count_core {
                        1 => {
                            let part = node_to_part[last_el_core.id];
                            let strand = if last_el_core.plus { '+' } else { '-' };
                            let first_start = path_starts[i];
                            let last_end = path_ends[j];
                            writeln!(writer, "{path_name}\tmice\tSO:0000856\t{first_start}\t{last_end}\t.\t{strand}\t.\tID={};genome={genome_name}", part+1)?;

                            // Obscured filtered ranges
                            let mut l = i + 1;
                            while l <= j {
                                if node_to_part[path[l].id] == FILTERED {
                                    let start_filter_pos = l;
                                    let mut m = l + 1;
                                    while m <= j && node_to_part[path[m].id] == FILTERED {
                                        m += 1;
                                    }
                                    let end_filter_pos = m - 1;
                                    let left_start_filtered = path_starts[start_filter_pos];
                                    let right_end_filtered = path_ends[end_filter_pos];
                                    let left_end_part = path_ends[start_filter_pos - 1];
                                    let right_start_part = path_starts[end_filter_pos + 1];
                                    if right_start_part > left_end_part + 1 {
                                        let start =
                                            usize::max(left_start_filtered, left_end_part + 1);
                                        let end =
                                            usize::min(right_end_filtered, right_start_part - 1);
                                        writeln!(writer, "{path_name}\tmice\tfiltered\t{start}\t{end}\t.\t{strand}\t.\tgenome={genome_name}")?;
                                    }
                                    l = m;
                                }
                                l += 1;
                            }

                            // Obscured Ns ranges
                            for l in i..j {
                                let left_end = path_ends[l];
                                let right_start = path_starts[l + 1];
                                if right_start > left_end + 1 {
                                    //No overlap
                                    let left_end = left_end + 1;
                                    let right_start = right_start - 1;
                                    writeln!(writer, "{path_name}\tmice\tNs\t{left_end}\t{right_start}\t.\t{strand}\t.\tgenome={genome_name}")?;
                                }
                            }
                        }
                        n if n > 1 => {
                            let mut l = i;
                            while l <= j {
                                let start = path_starts[l];
                                while l <= j && node_to_part[path[l].id] != path[l].id {
                                    l += 1;
                                }
                                if l <= j {
                                    let part = node_to_part[path[l].id];
                                    let strand = if path[l].plus { '+' } else { '-' };
                                    while l < j && node_to_part[path[l+1].id] != path[l+1].id {
                                        l += 1;
                                    }
                                    let end = path_ends[l];
                                    writeln!(writer, "{path_name}\tmice\tSO:0000856\t{start}\t{end}\t.\t{strand}\t.\tID={};genome={genome_name}", part+1)?;
                                    l += 1;
                                }
                            }
                        }
                        _ => {
                            eprintln!("Error: no core found in a substring in {genome_name}:");
                            let part = node_to_part[path[i].id];
                            for x in path[i..=j].iter() {
                                if part != node_to_part[x.id] {
                                    eprintln!("the parts are not correct");
                                }
                            }
                            eprintln!("part: {}", part+1);
                            for x in path[i..=j].iter() {
                                eprint!("{}{} ", x.id+1, if x.plus {'+'} else {'-'});
                            }
                            eprintln!();
                            return Ok(());
                        }
                    }
                    i = end + 1;
                }
            }
        }

        Ok(())
    }
}

struct GffRow {
    start: usize,
    end: usize,
    strand: bool,
    id: usize,
}

struct BoundedPath {
    pub path: Path,
    pub path_starts: Vec<usize>,
    pub path_ends: Vec<usize>,
    pub genome_name: Option<String>,
}

impl Gff {
    fn extract_gff_info_from_row(line: &[u8]) -> Option<(String, String, GffRow)> {
        let mut f = line.splitn(9, |&b| b == b'\t');

        let seqname = str::from_utf8(f.next()?).ok()?.to_owned();
        let _source = f.next()?;
        let _feature = f.next()?;
        let start = str::from_utf8(f.next()?).ok()?.parse::<usize>().ok()?;
        let end = str::from_utf8(f.next()?).ok()?.parse::<usize>().ok()?;
        let _score = f.next()?;
        let strand = f.next()? == b"+";
        let _frame = f.next()?;
        let attributes = str::from_utf8(f.next()?).ok()?;

        let mut id = None;
        let mut genome_name = None;

        for key_val in attributes.split(';') {
            let key_val = key_val.trim();
            if let Some((key, val)) = key_val.split_once('=') {
                if key.eq_ignore_ascii_case("id") {
                    if let Ok(parsed) = val.trim().parse::<usize>() {
                        id = Some(parsed);
                    }
                } else if key.eq_ignore_ascii_case("genome") {
                    genome_name = Some(String::from(val));
                }
            }
        }

        let id = id? - 1; // 1-based index (based on SibeliaZ, Cactus and gfa2gff)
        let genome_name = genome_name?;

        Some((
            seqname,
            genome_name,
            GffRow {
                start,
                end,
                strand,
                id,
            },
        ))
    }

    fn parse_gff_header(line: &[u8]) -> Option<(String, usize, usize)> {
        let mut f = line.splitn(4, |&b| b == b' ');
        let _ = str::from_utf8(f.next()?).ok()?.to_owned();
        let seqname = str::from_utf8(f.next()?).ok()?.to_owned();
        let start = str::from_utf8(f.next()?).ok()?.parse::<usize>().ok()?;
        let end = str::from_utf8(f.next()?)
            .ok()?
            .trim()
            .parse::<usize>()
            .ok()?;

        Some((seqname, start, end))
    }

    pub fn parse_gff_paths(filename: &str, group_by: bool) -> Result<(GenomeBundle, usize)> {
        let mut bounded_paths: HashMap<String, BoundedPath> = HashMap::default();
        let mut header: HashMap<String, usize> = HashMap::default();
        let mut num_nodes = 0usize;

        let file = File::open(filename)?;
        let mut reader = BufReader::new(file);

        // Parsing gff
        let mut buf = vec![];
        while reader.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
            if buf[0] != b'#' {
                if let Some((path_name, genome_name, row)) = Self::extract_gff_info_from_row(&buf) {
                    num_nodes = usize::max(row.id, num_nodes);

                    let genome_name = if group_by { Some(genome_name) } else { None };

                    let el = SignedId {
                        id: row.id,
                        plus: row.strand,
                    };
                    let entry = bounded_paths.entry(path_name);
                    match entry {
                        hash_map::Entry::Occupied(mut p) => {
                            p.get_mut().path.push(el);
                            p.get_mut().path_starts.push(row.start);
                            p.get_mut().path_ends.push(row.end);
                        }
                        hash_map::Entry::Vacant(p) => {
                            p.insert(BoundedPath {
                                path: vec![el],
                                path_starts: vec![row.start],
                                path_ends: vec![row.end],
                                genome_name,
                            });
                        }
                    }
                } else {
                    eprintln!(
                        "Error: unable to read line: {}",
                        str::from_utf8(&buf).unwrap()
                    );
                }
            } else if !buf.starts_with(b"##gff-version") {
                if let Some((path_name, _start, end)) = Self::parse_gff_header(&buf) {
                    header.insert(path_name, end);
                } else {
                    eprintln!(
                        "Error: unable to read header: {}",
                        str::from_utf8(&buf).unwrap()
                    );
                }
            }
            buf.clear();
        }
        num_nodes += 1;

        let num_paths = bounded_paths.len();

        //Populate genomes and sort each bounded_path based on bounded_path.path_starts
        let mut genomes: HashMap<String, PathBundle> = HashMap::default();
        for (path_name, bounded_path) in bounded_paths {
            let path_len = bounded_path.path.len();

            let mut idx: Vec<usize> = (0..path_len).collect();
            idx.sort_by(|&i, &j| bounded_path.path_starts[i].cmp(&bounded_path.path_starts[j]));

            let mut path_new = Vec::with_capacity(path_len);
            let mut path_starts_new = Vec::with_capacity(path_len);
            let mut path_ends_new = Vec::with_capacity(path_len);
            for &i in idx.iter() {
                path_new.push(bounded_path.path[i]);
                path_starts_new.push(bounded_path.path_starts[i]);
                path_ends_new.push(bounded_path.path_ends[i]);
            }

            let entry = genomes.entry(bounded_path.genome_name.unwrap_or(path_name.clone()));
            match entry {
                hash_map::Entry::Occupied(mut p) => {
                    if !header.is_empty() {
                        p.get_mut().path_sizes.push(header[&path_name]);
                    }
                    p.get_mut().paths.push(path_new);
                    p.get_mut().path_names.push(path_name);
                    p.get_mut().path_starts.push(path_starts_new);
                    p.get_mut().path_ends.push(path_ends_new);
                }
                hash_map::Entry::Vacant(p) => {
                    let mut path_sizes = Vec::new();
                    if !header.is_empty() {
                        path_sizes.push(header[&path_name]);
                    }
                    p.insert(PathBundle {
                        paths: vec![path_new],
                        path_names: vec![path_name],
                        path_starts: vec![path_starts_new],
                        path_ends: vec![path_ends_new],
                        path_sizes,
                    });
                }
            }
        }

        let node_indexer = None;
        let genome_bundle = GenomeBundle {
            genomes,
            num_paths,
            node_indexer,
        };
        Ok((genome_bundle, num_nodes))
    }
}
