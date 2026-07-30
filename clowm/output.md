# Output

MICE outputs three files into your chosen `--out_dir` (named `mice_output` by default).

## `output.gff`

This GFF file contains the annotations of the resulting synteny blocks. It follows the standard GFF3 format.

A section of an example `output.gff` file looks like this:

```
##gff-version 3
##sequence-region NC_000913.3 1 4641652
...
NC_000913.3     mice    SO:0000856      1       5593    .       +       .       ID=4051;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5564    5594    .       +       .       ID=3289;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5565    5597    .       +       .       ID=3223;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5568    5599    .       +       .       ID=3167;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5570    5600    .       +       .       ID=3263;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5571    5601    .       -       .       ID=3213;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5572    5602    .       -       .       ID=3296;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5573    5658    .       -       .       ID=253315;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5629    5671    .       -       .       ID=277551;genome=GCF_000005845.2_ASM584v2_genomic.fna
NC_000913.3     mice    SO:0000856      5642    9967    .       -       .       ID=1727;genome=GCF_000005845.2_ASM584v2_genomic.fna
...
```

The first section of the file, commented with `##`, contains the metadata header. Following this, each entry has the following columns:

1. Sequence region (chromosome, contig etc.) in which the synteny block is located
2. `mice` feature source identifier
3. "Conserved" ontology tag (feature type) 
4. Start position of block
5. End position of block
7. _Score (ignored)_
8. Strand
9. _Frame (ignored)_
10. Attributes `ID`, which gives the unique numerical ID of the synteny block (corresponding to its canonical anchor); and `genome`, the block's genome of origin.


## `paths.txt`

This text file contains the input genomes re-spelled with block instances, rather than nucleotides.

A section of an example `paths.txt` file looks like this:

```
>GCF_000005845.2_ASM584v2_genomic.fna#NC_000913.3
4051+,3289+,3223+,3167+,3263+,3213-,3296-,253315-,277551-,1727-,7202+,136610-,132800-,945+,22460-,22453+,285692+,23731-,105587-,107243-,36801-,36806-,36817+,234996+,210660-,39203-,212279+,175468+,175799+,153393-,17748+,232769+,230480+,133073+,37466-,193224-,23750-,257200+,138875+,83614+,83624-,226872-,297677+,290125-,287142-,285542-,281148+,25726-,283953-,18910+,18878-,147687-,112540-,44299+,44330-,44334+,183770-,241689+,242415+,150534-,17870+,101067-,101828-,45916-,241879+,48739+,946-,49654+,29476+,150670+,17980-,18012+,11550+,524-,106095-,20055+,20108+,20172-,20110+,20073-,293561-,3559-,60082+,45268+,45349+,45406+,45220+,2463-,10577+,4717+,13293-,47565+,6023+,2762-,485-,34645-,1488-,20057+,19994+,293593+,13196-,13183-,60557-,3268+,3121-,3256-,3263-,3167-,3223-,3289-,59851-,6080-,290476-,290500-,20057+,19994+,293593+,60498-,3278-,3213+,3263-,3167-,3223-,3289-,59851-,6080-,290476-,290500-,20057+,19994+,293593+,13196-,13183-,242527+,45751-,71775-,6069-,1406+,277412+,3247-,3159-,3318+,...
```

The file will be `2n` lines long, where `n` is the amount of individual input genomes. Odd-numbered lines will have a FASTA-style header with the name of the genome about to be re-spelled; the line underneath each header will be the entire genome, respelled as a comma-separated sequence of block instances in the format `<ID><SIGN>`, where `<ID>` is the block's unique `ID` attribute (as given in `output.gff`) and `<SIGN>` is either `+` or `-`, denoting whether this instance of the block is in its "forward" or "backward" orientation, respectively.

## `partitions.txt`

This text file contains the elements which where merged into each block's canonical anchor.

A section of an example `partitions.txt` file looks like this:

```
2: 2
4: 4 422 1565 2417 2453 2629 2632 3889 ...
```

The file will have as many lines as there are synteny blocks. Each line will start with `<ID>: `, the ID of the canonical anchor which represents that synteny block, followed by a space-separated list of element IDs which were merged into this synteny block. The first element in this list always matches the canonical anchor ID (as, by definition, every block contains its canonical anchor); if this is the only element present in the list, the synteny block is a singleton.