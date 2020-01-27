# pslToBed
`pslToBed` loads PSL file (tpsl/gpsl), and generates BED of genetic regions.

## Usage

```
Usage:   perl BioFuse.pl pslToBed <[Options]>

Options:

   # Input and Output #
    -tp [s]  transcript PSL file. <required>
    -o  [s]  output bed file. <required>

   # Options #
    -pz [i]  promoter size extended to upstream from the first exon. [1000]
    -tz [i]  terminator size extended to downstream from the last exon. [1000]
    -ab [s]  the refseg you want to avoid, allow mutilple input. [NULL]
    -gn [s]  the gene name you want to keep only, allow mutilple input. [NULL]
    -bt [s]  the biotype you want to keep only, allow mutilple input. [NULL]
    -rt [s]  the region type you want to keep only, allow mutilple input. [NULL]
              NOTE: region types: CDS/UTR5/UTR3/EXON/INTRON/PROMOTER/TERMINATOR
    -pm [s]  mode to select 'protein_coding' transcript for CDS/UTR5/UTR3. [trNO]
              NOTE: available mode: trNO, longest
    -em [s]  mode to select transcript for EXON/INTRON/PROMOTER/TERMINATOR. [protein_trNO]
              NOTE: available mode: protein_trNO, trNO, longest, merge
    -one     output region is one-based coordinate, not BED format. [disabled]
    -abp     include abnormal protein_coding transcript. [disabled]

    -h       show this help

Version:
   0.03 at 2019-10-07

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
