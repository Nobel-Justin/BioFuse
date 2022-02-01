package BioFuse::BioInfo::Alignment::SAMtools;

use strict;
use warnings;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              SAMtool_calmd
              SAMtools_mpileup
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Alignment::SAMtools';
#----- version --------
$VERSION = "0.01";
$DATE = '2021-12-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        SAMtool_calmd
                        SAMtools_mpileup
                     /;

#--- samtools calmd pipeline ---
sub SAMtool_calmd{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ref = $parm{ref};
    my $sBamPath = $parm{sBamPath};
    my $oBamPath = $parm{oBamPath};
    my $samtools = $parm{samtools};

    # samtools calmd
    my $cmd = "($samtools calmd -b -r $sBamPath $ref >$oBamPath 2>/dev/null) && "
             ."($samtools index $oBamPath)";
    trible_run_for_success($cmd, 'calmd', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tSAMtools calmd bam finished.\n"
                         ."\tI=$sBamPath\n"
                         ."\tO=$oBamPath\n";
}

#--- samtools mpileup to vcf pipeline ---
sub SAMtools_mpileup{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ref = $parm{ref};
    my $sBamPath = $parm{sBamPath};
    my $vcfPath = $parm{vcfPath};
    my $samtools = $parm{samtools};

    # samtools calmd
    my $cmd = "($samtools mpileup -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR' -A -B -d 1000000 -v -O -s -m 3 -L 1000000 -f $ref -o $vcfPath $sBamPath)";
    trible_run_for_success($cmd, 'mpileup', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tSAMtools mpileup to vcf finished.\n"
                         ."\tI=$sBamPath\n"
                         ."\tO=$vcfPath\n";
}

1; ## tell the perl script the successful access of this module.
