package BioFuse::BioInfo::GeneAnno::GTFtoGenePSL;

use strict;
use warnings;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::GeneAnno::GTF qw/ read_GTF add_refseg_cytoband create_gene_PSL /;
use BioFuse::BioInfo::CytoBand qw/ load_cytoband /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              GTFtoGenePSL
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::GTFtoGenePSL';
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
                        GTFtoGenePSL
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} get_gpsl <[Options]>

     Options:
         -g   [s]  gtf database file. <required>
         -o   [s]  output gene PSL file. <required>
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
         -cbd [s]  the cytoBand database file. [optional]
                    e.g., Download human file from UCSC:
                       For hg18:  http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz
                       For hg19:  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
                       For hg38:  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
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
            [ out => undef ],
            [ gtf => undef ],
            [ refseg_transform_list => undef ],
            [ cytoBand_file => undef ],
            # option
            [ gtf_source => [] ],

            # intermediate variants
            [ GTF_Info => {} ],
            [ cytoBand => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['out'],
                                  ['gtf']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-o:s"  => \$V_Href->{out},
        "-g:s"  => \$V_Href->{gtf},
        "-rft:s"=> \$V_Href->{refseg_transform_list},
        "-cbd:s"=> \$V_Href->{cytoBand_file},
        # options
        "-sor:s"=> \@{$V_Href->{gtf_source}},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{gtf})
             || !defined $V_Href->{out}
            );
}

#--- get gene PSL file from GTF file---
sub GTFtoGenePSL{

    # read GTF
    read_GTF( GTF_Info_Href => $V_Href->{GTF_Info},
              gtf_file => $V_Href->{gtf},
              refseg_transform_list => $V_Href->{refseg_transform_list},
              gtf_gene_source_Aref => $V_Href->{gtf_source}
            );
    # read cytoBand
    if( defined $V_Href->{cytoBand_file} ){
        load_cytoband( cytoBand_Href => $V_Href->{cytoBand},
                       cytoBand_file => $V_Href->{cytoBand_file}
                     );
        add_refseg_cytoband( GTF_Info_Href => $V_Href->{GTF_Info},
                             cytoBand_Href => $V_Href->{cytoBand}
                            );
    }
    # output gene psl
    create_gene_PSL( GTF_Info_Href => $V_Href->{GTF_Info},
                     PSL_file => $V_Href->{out}
                    );
    
}

#--- 
1; ## tell the perl script the successful access of this module.
