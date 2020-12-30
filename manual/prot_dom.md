# prot_dom
`prot_dom` loads gene/trans/protein id, and fetches domain information from online database.

## Usage

```
 Usage:   perl BioFuse.pl prot_dom <[Options]>

 Options:

    # Input and Output #
     -i  [s]  list of query id. <required>
              for ncbi: ENSG/P/T or [NX]M/P id.
              for ensm: ONLY ENST id.
     -o  [s]  output domain list. <required>

    # Options #
     -d  [s]  select database. [ncbi/ensm]
              ncbi: NCBI-GENE; ensm: Ensembl-Domains
     -t  [i]  maximum try times to connect database website for one query. [10]

     -h       show this help

 Version:
    0.03 at 2020-01-15

 Author:
    Wenlong Jia (wenlongkxm@gmail.com)
```
