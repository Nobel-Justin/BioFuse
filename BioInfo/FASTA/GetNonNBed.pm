package BioFuse::BioInfo::FASTA::GetNonNBed;

use strict;
use warnings;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Write /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;
use BioFuse::Util::String qw/ getRegexRegion /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              GetNonNBedfromFasta
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::FASTA::GetNonNBed';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-11-29';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        GetNonNBedfromFasta
                        get_nonN_region_from_segseq
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} get_nonN <[Options]>

     Options:
         -f  [s]  fasta file. <required>
         -o  [s]  output BED file. <required>
         -one     output region is started from one, not BED format. [disabled]
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
            ## use 'whole_genome' in BioFuse::LoadOn
            [ nonN_bed => undef ],
            # option
            [ oneBase => 0 ],
            # list to abs-path
            [ ToAbsPath_Aref => [ ['nonN_bed'],
                                  ['whole_genome']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-o:s"  => \$V_Href->{nonN_bed},
        "-f:s"  => \$V_Href->{whole_genome},
        # option
        "-one"  => \$V_Href->{oneBase},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{whole_genome})
             || !defined $V_Href->{nonN_bed}
            );
}

#--- get non-N region bed file from reference fasta ---
sub GetNonNBedfromFasta{
    # read fasta file
    open (my $bedfh, Try_GZ_Write($V_Href->{nonN_bed})) || die "fail write output nonN_bed file: $!\n";
    read_fasta_file( FaFile => $V_Href->{whole_genome},
                     subrtRef => \&get_nonN_region_from_segseq,
                     subrtParmAref => [bedfh => $bedfh, oneBase => $V_Href->{oneBase}]
                   );
    close $bedfh;
}

#--- get non-N region from given chr-seq ---
## get zero-start(BED) interval.
sub get_nonN_region_from_segseq{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $segName = $parm{segName};
    my $segSeq_Sref = $parm{segSeq_Sref};
    my $bedfh = $parm{bedfh};
    my $oneBase = $parm{oneBase} || 0; # not bed, but start from one

    my $nonN_region_Aref = getRegexRegion( StrTest => $$segSeq_Sref,
                                           regex => 'N',
                                           dismatch => 1,
                                           ignore_case => 1,
                                           zerobase => !$oneBase
                                        );
    print {$bedfh} join("\t", $segName, @$_)."\n" for @$nonN_region_Aref;
    # inform
    stout_and_sterr "[INFO]\tget non-N region of refseg $segName OK.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
