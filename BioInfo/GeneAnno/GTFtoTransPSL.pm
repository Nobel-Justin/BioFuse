package BioFuse::BioInfo::GeneAnno::GTFtoTransPSL;

use strict;
use warnings;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::GeneAnno::GTF qw/ read_GTF create_trans_PSL mark_abnormal_Start_codon /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              GTFtoTransPSL
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::GTFtoTransPSL';
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
                        GTFtoTransPSL
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} get_tpsl <[Options]>

     Options:
          -g   [s]  gtf file. <required>
          -f   [s]  whole genome ref fasta file. <required>
          -o   [s]  trans PSL file to output. <required>
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
            ## use 'gtf' in BioFuse::LoadOn
            ## use 'whole_genome' in BioFuse::LoadOn
            ## use 'cytoBand_file' in BioFuse::LoadOn
            [ refseg_transform_list => undef ],
            # option
            [ gtf_source => [] ],

            # intermediate variants
            [ GTF_gene => {} ],
            [ cytoBand => {} ],
            [ _refseg_transform => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['psl'],
                                  ['whole_genome'],
                                  ['gtf']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-o:s"  => \$V_Href->{psl},
        "-g:s"  => \$V_Href->{gtf},
        "-rft:s"=> \$V_Href->{refseg_transform_list},
        "-cbd:s"=> \$V_Href->{cytoBand_file},
        "-f:s"  => \$V_Href->{whole_genome},
        # options
        "-sor:s"=> \@{$V_Href->{gtf_source}},
        "-stc:s"=> \@{$V_Href->{Start_codon}},
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
             || !file_exist(filePath=>$V_Href->{whole_genome})
             || !defined $V_Href->{psl}
            );
}

#--- get trans PSL file from GTF file---
sub GTFtoTransPSL{
    # read GTF
    read_GTF;
    # check start codon
    mark_abnormal_Start_codon;
    # output trans psl
    create_trans_PSL;
}

#--- 
1; ## tell the perl script the successful access of this module.
