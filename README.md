# BioFuse

Util PERL module of general functions and objects applied in bioinformatics software development.

- Author: Wenlong Jia
- Email:  wenlongkxm@gmail.com

## Version
0.16

`2022-04-03`

## Installation

BioFuse is written in PERL. All of its functions are packaged into a standalone PERL module. Besides, it also requires other additional modules.

To install BioFuse, you need to download BioFuse and add the current directory to the `PERL5LIB` path.
```bash
git clone https://github.com/Nobel-Justin/BioFuse.git
PERL5LIB=$PERL5LIB:$PWD; export PERL5LIB
```
List of additional PERL modules required:
- [Carp](https://metacpan.org/pod/Carp)
- [JSON](https://metacpan.org/pod/JSON)
- [List::Util](https://metacpan.org/pod/List::Util)
- [Math::Trig](https://metacpan.org/pod/Math::Trig)
- [POSIX](https://metacpan.org/pod/distribution/perl/ext/POSIX/lib/POSIX.pod)
- [LWP::UserAgent](https://metacpan.org/pod/LWP::UserAgent)

If you encounter problems, please open an issue at the [project on Github](https://github.com/Nobel-Justin/BioFuse/issues).

## Commands

Besides of functions as perl modules, BioFuse also provides several commands.

- [get_tpsl](./manual/get_tpsl.md)

  `get_tpsl` generates the transcripts PSL file from GTF file.

- [get_gpsl](./manual/get_gpsl.md)

  `get_gpsl` generates the gene PSL file from GTF file.

- [psl_seq](./manual/psl_seq.md)

  `psl_seq` generates the FASTA file of transcripts (tpsl) or genes (gpsl).

- [pslToBed](./manual/pslToBed.md)

  `pslToBed` generates BED of genetic regions from PSL file.

- [mut_fa](./manual/mut_fa.md)

  `mut_fa` generates FASTA file based on given mutations.

- [prot_dom](./manual/prot_dom.md)

  `prot_dom` fetches domain information from online database.
