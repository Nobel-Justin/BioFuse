# psl_seq
`psl_seq` loads PSL file (tpsl/gpsl) and reference genome, and generates the FASTA file of transcripts (tpsl) or genes (gpsl).

## Usage

```
Usage:   perl BioFuse.pl psl_seq <[Options]>

Options:

   # Input and Output #
    -p   [s]  PSL file. <required>
    -f   [s]  whole genome ref fasta file. <required>
              NOTE: please make sure that the ref_seg in the PSL file are all included in the
                    whole genome ref fasta file.
    -o   [s]  output fa file. <required>

   # Options #
    -t   [s]  type of PSL file, 'gene' or 'trans'. <required>
    -ex  [i]  extend length of gene/trans body. [0]
    -ab  [s]  the refseg you want to avoid, allow mutilple input. [NULL]
    -nm  [s]  the gene/trans name you want to keep only, allow mutilple input. [NULL]
    -bt  [s]  the biotype you want to keep only, allow mutilple input. [NULL]

    -h        show this help

Version:
   0.85 at 2019-05-04

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
