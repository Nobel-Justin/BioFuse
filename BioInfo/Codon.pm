package BioFuse::BioInfo::Codon;

use strict;
use warnings;
use BioFuse::Util::Log qw/ stout_and_sterr /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              Load_Codon
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Codon';
#----- version --------
$VERSION = "0.31";
$DATE = '2018-10-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @function_list = qw/
                        Load_Codon
                     /;

#----------- load codon table ---------
sub Load_Codon{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my ($Codon_to_Amino_Href) = @_;

    $Codon_to_Amino_Href = {
        'TTT','F',    'CTT','L',    'ATT','I',    'GTT','V',
        'TTC','F',    'CTC','L',    'ATC','I',    'GTC','V',
        'TTA','L',    'CTA','L',    'ATA','I',    'GTA','V',
        'TTG','L',    'CTG','L',    'ATG','M',    'GTG','V',
        'TCT','S',    'CCT','P',    'ACT','T',    'GCT','A',
        'TCC','S',    'CCC','P',    'ACC','T',    'GCC','A',
        'TCA','S',    'CCA','P',    'ACA','T',    'GCA','A',
        'TCG','S',    'CCG','P',    'ACG','T',    'GCG','A',
        'TAT','Y',    'CAT','H',    'AAT','N',    'GAT','D',
        'TAC','Y',    'CAC','H',    'AAC','N',    'GAC','D',
        'TAA','*',    'CAA','Q',    'AAA','K',    'GAA','E',
        'TAG','*',    'CAG','Q',    'AAG','K',    'GAG','E',
        'TGT','C',    'CGT','R',    'AGT','S',    'GGT','G',
        'TGC','C',    'CGC','R',    'AGC','S',    'GGC','G',
        'TGA','*',    'CGA','R',    'AGA','R',    'GGA','G',
        'TGG','W',    'CGG','R',    'AGG','R',    'GGG','G'
    };
    # inform
    stout_and_sterr "[INFO]\tLoad codon ok!\n";
}

1; ## tell the perl script the successful access of this module.
