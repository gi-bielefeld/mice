#!/bin/bash

workdir=/path/to/dir/
filelist=$workdir/filelist.txt
reference="TODO"


rm -r js; /usr/bin/time -f "Time\n%e" -o $workdir/benchmark.txt cactus-pangenome ./js $filelist --outDir $workdir/cactus --reference $reference --outName $workdir # --filter 0 --clip 0

rm -r js; cactus-hal2maf ./js $workdir/cactus/$workdir.full.hal $workdir/maf.maf --refGenome $reference

python3 ./maf2gff.py $workdir/maf.maf > $workdir/blocks.gff
