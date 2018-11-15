package BioFuse::BioInfo::GeneAnno::PSLtoFASTA;

use strict;
use warnings;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::GeneAnno::PSL qw/ read_unit_region_from_PSL extract_exon_seq_from_genome_and_output /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              PSLtoFASTA
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::PSLtoFASTA';
#----- version --------
$VERSION = "0.82";
$DATE = '2018-11-15';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        PSLtoFASTA
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} get_tpsl <[Options]>

     Options:
         -p  [s]  PSL database file. <required>
         -f  [s]  whole genome ref fasta file. <required>
                  NOTE: please make sure that the ref_seg in the PSL file are all included in the 
                        whole genome ref fasta file.
         -o  [s]  output fa file. <required>
         -ab [s]  the refseg you want to avoid. [NULL]
                  NOTE: could input in mutilple times
         -h       show this help

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
    $V_Href->{ $_->[0] } = $_->[1] for
        map {
            if( !exists $V_Href->{$_->[0]} ){
                ( $_ );
            }
            else{
                warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            [ psl => undef ],
            [ out => undef ],
            [ whole_genome => undef ],
            # option
            [ abandon_refseg => [] ],

            # intermediate variants
            [ Chr_Exon_Info => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['out'],
                                  ['psl'],
                                  ['whole_genome']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-p:s"  => \$V_Href->{psl},
        "-o:s"  => \$V_Href->{out},
        "-f:s"  => \$V_Href->{whole_genome},
        # options
        "-ab:s" => \$V_Href->{abandon_refseg},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{psl})
             || !file_exist(filePath=>$V_Href->{whole_genome})
             || !defined $V_Href->{out}
            );
}

#--- get fa files based on PSL file ---
sub PSLtoFASTA{

    # load psl file
    my %Chr_Exon_Info;
    read_unit_region_from_PSL( Refseg_Exon_Info_Href => $V_Href->{Chr_Exon_Info},
                               psl_file => $V_Href->{psl},
                               avoid_refseg_Aref => $V_Href->{abandon_refseg}
                            );
    # output fasta
    extract_exon_seq_from_genome_and_output( Refseg_Exon_Info_Href => $V_Href->{Chr_Exon_Info},
                                             out => $V_Href->{out},
                                             whole_genome => $V_Href->{whole_genome}
                                            );
}

#--- 
1; ## tell the perl script the successful access of this module.
