import random
from pathlib import Path
from itertools import product

def generate_genome_set(genome_len, genome_count, max_id_mult=6, duplicates=False, rng=None):
    rng = rng or random
    max_gene_id = genome_len * max_id_mult
    gene_pool = [f"{i}" for i in range(max_gene_id)]

    genomes = []
    if not duplicates:
        for _ in range(genome_count):
            genome = rng.sample(gene_pool, genome_len)
            rng.shuffle(genome)
            signed_genome = [gene + rng.choice(['+', '-']) for gene in genome]
            genomes.append(signed_genome)
    else:
        for _ in range(genome_count):
            genome = rng.choices(gene_pool, k=genome_len)
            signed_genome = [gene + rng.choice(['+', '-']) for gene in genome]
            genomes.append(signed_genome)
    return genomes

def _parse_signed(g):
    return g[:-1], g[-1]

def write_gfa(genomes, outfile, seq="A", overlap="0M", meta=None):
    meta = meta or {}
    # Segments
    segments = sorted({_parse_signed(g)[0] for path in genomes for g in path})
    # Links
    links = set()
    for path in genomes:
        for a, b in zip(path, path[1:]):
            a_id, a_or = _parse_signed(a)
            b_id, b_or = _parse_signed(b)
            links.add((a_id, a_or, b_id, b_or))

    p = Path(outfile)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as fh:
        tags = ["VN:Z:1.0"]
        tags.append(f"m:i:{meta['m']}")
        tags.append(f"n:i:{meta['n']}")
        tags.append(f"mad_id:i:{meta['mad_id']}")
        tags.append(f"dup:Z:{'T' if meta['dup'] else 'F'}")
        tags.append(f"rep:i:{meta['rep']}")
        fh.write("H\t" + "\t".join(tags) + "\n")

        for sid in segments:
            fh.write(f"S\t{sid}\t{seq}\n")

        for a_id, a_or, b_id, b_or in sorted(links):
            fh.write(f"L\t{a_id}\t{a_or}\t{b_id}\t{b_or}\t{overlap}\n")

        for i, path in enumerate(genomes, start=1):
            segs = ",".join([f"{_parse_signed(g)[0]}{_parse_signed(g)[1]}" for g in path]) if path else "*"
            over = ",".join([overlap] * (len(path) - 1)) if len(path) > 1 else "*"
            fh.write(f"P\tg{i}\t{segs}\t{over}\n")

if __name__ == "__main__":
    out_dir = Path("gfa")
    out_dir.mkdir(parents=True, exist_ok=True)

    m_list = [2, 3, 5, 7, 100, 1000]
    n_list = [2, 3, 5, 15]
    mad_list = [1, 2, 3, 4]
    dup_list = [False, True]

    for m, n, mad_id, dup in product(m_list, n_list, mad_list, dup_list):
        for rep in range(1, 11):  # 10 replicates
            seed = hash((m, n, mad_id, dup, rep)) & 0xFFFFFFFF
            rng = random.Random(seed)

            genomes = generate_genome_set(
                genome_len=m,
                genome_count=n,
                max_id_mult=mad_id,
                duplicates=dup,
                rng=rng
            )

            fname = f"m{m}_n{n}_mad{mad_id}_dup{'T' if dup else 'F'}_rep{rep}.gfa"
            path = out_dir / fname

            write_gfa(
                genomes,
                path,
                seq="A",
                overlap="0M",
                meta={"m": m, "n": n, "mad_id": mad_id, "dup": dup, "rep": rep}
            )

    print(f"Wrote all GFAs to: {out_dir.resolve()}")

