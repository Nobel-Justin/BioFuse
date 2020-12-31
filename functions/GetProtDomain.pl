#!/usr/bin/perl -w
use strict;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::BioInfo::GeneAnno::GetProtDomain;

my ($VERSION, $DATE, $AUTHOR, $EMAIL);

#----- version --------
$VERSION = "1.0";
$DATE = '2020-11-22';
#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

my $module = 'BioFuse::BioInfo::GeneAnno::GetProtDomain';
my $usage = "
     Usage:   perl GetProtDomain.pl <[Options]>

     Options:

        # Input and Output #
         -i  [s]  list of query id. <required>
                  for ncbi: ENSG/P/T or [NX]M/P id.
                  for ensm: ONLY ENST id.
         -o  [s]  output domain list. <required>

        # Options #
         -d  [s]  select database. <required>
                  ncbi: NCBI-GENE; ensm: Ensembl-Domains
         -t  [i]  maximum try times to connect database website for one query. [10]

         -h       show this help

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

"$module"->GetProtDomain;
