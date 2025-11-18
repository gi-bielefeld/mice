from Bio import AlignIO
import sys
from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument("maf")
args = parser.parse_args()

alignments = AlignIO.parse(args.maf, "maf")
bid=0
for alignment in alignments:
    bid+=1
    for record in alignment:
        sid = ".".join(record.id.split(".")[1:])
        if sid.startswith("Anc"):
           continue
        source = "cactus"
        feature = "SO:0000856"
        score="."
        frame="."
        attribute=f"ID={bid}"
        an = record.annotations
        start = an["start"]
        end = an["start"]+an["size"]-1
        strand = "+" if  an["strand"] == 1 else "-"
        print(f"{sid}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attribute}")
