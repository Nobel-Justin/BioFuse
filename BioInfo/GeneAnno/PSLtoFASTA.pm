package BioFuse::BioInfo::GeneAnno::PSLtoFASTA;

use strict;
use warnings;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::Objects::GeneAnno::PSL_OB;

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
$VERSION = "0.85";
$DATE = '2019-05-04';

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
     Usage:   perl $V_Href->{MainName} psl_seq <[Options]>

     Options:

        # Input and Output #
         -p   [s]  PSL file. <required>
         -f   [s]  whole genome ref fasta file. <required>
                   NOTE: please make sure that the ref_seg in the PSL file are all included in the 
                         whole genome ref fasta file.
         -o   [s]  output fa file. <required>

        # Options #
         -t   [s]  type of PSL file, 'gene' or 'trans'. <required>
         -ex  [i]  extend length of gene/trans body. [0]
         -ab  [s]  the refseg you want to avoid, allow mutilple input. [NULL]
         -nm  [s]  the gene/trans name you want to keep only, allow mutilple input. [NULL]
         -bt  [s]  the biotype you want to keep only, allow mutilple input. [NULL]

         -h        show this help

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
            ## use 'psl' in BioFuse::LoadOn
            ## use 'whole_genome' in BioFuse::LoadOn
            [ seqFA => undef ],
            # option
            [ psl_type => undef ],
            [ extendLength => 0 ],
            [ abandon_refseg => [] ],
            [ require_objname => [] ],
            [ require_biotype => [] ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['seqFA'],
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
        "-o:s"  => \$V_Href->{seqFA},
        "-f:s"  => \$V_Href->{whole_genome},
        # options
        "-t:s"  => \$V_Href->{psl_type},
        "-ex:i" => \$V_Href->{extendLength},
        "-ab:s" => \@{$V_Href->{abandon_refseg}},
        "-nm:s" => \@{$V_Href->{require_objname}},
        "-bt:s" => \@{$V_Href->{require_biotype}},
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
             || $V_Href->{extendLength} < 0
             || !defined $V_Href->{seqFA}
             || !defined $V_Href->{psl_type}
            );
}

#--- get fa files based on PSL file ---
sub PSLtoFASTA{
    # psl object
    my $psl = BioFuse::BioInfo::Objects::GeneAnno::PSL_OB->new(filePath => $V_Href->{psl}, psl_type => $V_Href->{psl_type});
    # load psl file
    $psl->load_GeneOrTrans_from_PSL( skipRefSegHref => { map {($_,1)} @{$V_Href->{abandon_refseg}} },
                                     requireObjHref => { map {($_,1)} @{$V_Href->{require_objname}} },
                                     requireBtyHref => { map {($_,1)} @{$V_Href->{require_biotype}} }
                                   );
    # output fasta
    $psl->extract_GeneOrTrans_seq( outputFasta => $V_Href->{seqFA},
                                   whole_genome => $V_Href->{whole_genome},
                                   extendLength => $V_Href->{extendLength}
                                 );
}

#--- 
1; ## tell the perl script the successful access of this module.
