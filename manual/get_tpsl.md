# get_tpsl
`get_tpsl` loads GTF gene annotation file and reference genome, and generates the transcripts PSL file.

## Usage

```
Usage:   perl BioFuse.pl get_tpsl <[Options]>

Options:
     -g   [s]  gtf file. <required>
     -f   [s]  whole genome ref fasta file. <required>
     -o   [s]  trans PSL file to output. <required>
     -rft [s]  TABLE-separated list of refseg symbols relationship. [optional]
                NOTE: Generaly, the ref_seg symbols in GTP file (normally, the first column) is different from
                      that in reference file (the '-f' parameter). So, you should specify the corresponding relationship
                      of refseg symbols between GTF file and reference file.
                      Format:  refseg_symbol_in_gtf  refseg_symbol_in_reference
                      e.g.,  10  chr10
     -sor [s]  select the source of gtf, can be used in multiple times. Default to accept all sources.
                NOTE: From v75, ensembl provides the source_info stored by 'gene_source' tag.
                      It will accept the transcript whose gene_source tag contains your input.
                      For instance, once you input 'havana', it will accept both 'havana' and 'ensembl_havana'.
                                    if you input 'ensembl', it means both 'ensembl' and 'ensembl_havana' will be accepted.
     -stc [s]  the standard start codon sequence, forced to upper cases, can be used in multiple times. [ATG]
                NOTE: Although, to our knowledge, the start condon sequence ATG is extensive in all species,
                      we still provide this para for user to state other sequences. Actually, some genes have
                      different sequence, such as, MYC-001 with CTG, TEAD4-001 with TTG, and some genes in
                      mitochondria and chloroplast. So, if you want to save these genes as 'protein-coding',
                      but not 'protein-coding-with-abnormal-start_codon', please state their special start_codon
                      sequences to take into account, like ' -stc CTG -stc TTG '.
     -cbd [s]  the cytoBand file. [optional]
                e.g., Download human file from UCSC:
                   For hg19:  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
                   For hg38:  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
     -h        show this help

Version:
   0.83 at 2019-05-04

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
