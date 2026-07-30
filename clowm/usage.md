# Parameters used by MICE

## Input / Output Options

| Name | Datatype | Default | Description |
| :-: | :-: | :-: | :- | 
| `input_fasta` | Filepath | - | A filepath to a single GFF-format file, where each entry has a 1-based `ID=` attribute. <br></br> The file must be plain-text (not "rich text", eg. `.rtf` or `.docx`), but can have any name and any suffix (ie. it doesn't necessarily need to end in `.gff` or similar). <br></br> Do not submit raw text or any other non-filepath input into this parameter - only a string representing a filepath, pointing to a file that already exists in your chosen S3 bucket.|
| `out_dir` | Filepath | - | A filepath to a directory where you want MICE to save its three output files. This can be just the root of an S3 bucket, or a directory inside it (if it doesn't exist, it will be created). See the "Output" tab for more information on the output format.|

## Execution Options

| Name | Datatype | Default | Description |
| :-: | :-: | :-: | :- | 
| `remove_dup` | Integer | 0 | Remove an element if it occurs at least `remove_dup` times in any single genome. Set to 0 to disable this functionality. Minimum value 0. A value of 1 will result in an empty input; if you want this option enabled, use 2 or above. |
| `quorum` | Integer | 0 | Remove an element if it occurs in fewer that `quorum` of the input genomes. Set to 0 to disable this functionality. Minimum value 0. If you set this option to a value higher than the amount of input genomes, it will result in an empty input. |
| `min_size` | Integer | 0 | After the first compaction step, drop yet-unmerged elements shorter than `min_size` base pairs, then compact again. Set to 0 to disable this functionality. Minimum value 0. |
| `no_group_by` | Boolean | False | If true, every unique sequence region (chromosome, contig etc.) is treated as its own genome. Used eg. for the purpose of determining whether an element is duplicated (it could be duplicated within a whole genome, but in two different chromosomes). |
| `merge_with_dup` | Boolean | False | When `True`, avoids marking duplicated elements in each genome and merges them as if they were unique. <br></br> Setting this option to `True` effectively activates the MICE (dup) mode described in the paper. Keeping it `False` instead keeps MICE (Bp bij) mode active. |
