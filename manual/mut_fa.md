# mut_fa
`mut_fa` loads mutation list and reference genome, and generates FASTA file (local region) of given mutations.

## Usage

```
Usage:   perl BioFuse.pl mut_fa <[Options]>

Options:
    -f  [s]  reference genome fasta file. <required>
    -m  [s]  mutation list. <required>
    -o  [s]  output mutated fasta file. <required>
    -s  [i]  flanking size around mutation, minimum:200. [1000]
    -d  [i]  minimum distance allowed between neighbor mutations. [0, disabled]
             check snv/ins/del
    -da [i]  action when found distance less than '-d'. [x]
             'x': alert and eXit; 's': Skip latter mutation and continue.
    -b  [s]  region BED file to filter mutations.
    -one     given region is one-based coordinate, not BED format. [disabled]
    -h       show this help

Version:
   0.05 at 2019-10-22

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
