#!/usr/bin/perl -w
use strict;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::BioInfo::GeneAnno::GTFtoGenePSL;

my ($VERSION, $DATE, $AUTHOR, $EMAIL);

#----- version --------
$VERSION = "1.0";
$DATE = '2019-08-15';
#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

my $module = 'BioFuse::BioInfo::GeneAnno::GTFtoGenePSL';
my $usage = "
     Usage:   perl GetGpsl.pl <[Options]>

     Options:
         -g   [s]  gtf file. <required>
         -o   [s]  gene PSL file to output. <required>
         -rft [s]  list of refseg symbols relationship. [optional]
                    NOTE: Generaly, the ref_seg symbols in GTP file (normally, the first column) is different from 
                          that in reference file (the '-f' parameter). So, you should specify the corresponding relationship
                          of refseg symbols between GTF file and reference file.
                          Format:  refseg_symbol_in_gtf  \\t  refseg_symbol_in_reference
                          e.g.,  10 \\t chr10
         -sor [s]  select the source of gtf, can be used in multiple times. Default to accept all sources.
                    NOTE: From v75, ensembl provides the source_info stored by 'gene_source' tag.
                          It will accept the transcript whose gene_source tag contains your input.
                          For instance, once you input 'havana', it will accept both 'havana' and 'ensembl_havana'.
                                        if you input 'ensembl', it means both 'ensembl' and 'ensembl_havana' will be accepted.
         -cbd [s]  the cytoBand file. [optional]
                    e.g., Download human file from UCSC:
                       For hg19:  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
                       For hg38:  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
         -h        show this help

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
\n";

# load module variants
"$module"->Load_moduleVar_to_pubVarPool;
# get options
"$module"->Get_Cmd_Options;

if("$module"->para_alert){
	warn_and_exit $usage;
}

"$module"->GTFtoGenePSL;
