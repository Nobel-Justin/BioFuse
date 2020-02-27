#!/usr/bin/perl -w
use strict;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::BioInfo::FASTA::GetMutFA;

my ($VERSION, $DATE, $AUTHOR, $EMAIL);

#----- version --------
$VERSION = "1.0";
$DATE = '2019-10-22';
#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

my $module = 'BioFuse::BioInfo::FASTA::GetMutFA';
my $usage = "
     Usage:   perl GetMutFA.pl <[Options]>

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

"$module"->GetMutFasta;
