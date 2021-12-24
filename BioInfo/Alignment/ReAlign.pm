package BioFuse::BioInfo::Alignment::ReAlign;

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
              GATK_ReAlign
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Alignment::ReAlign';
#----- version --------
$VERSION = "0.01";
$DATE = '2021-12-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        GATK_ReAlign
                     /;

#--- GATK-realn pipeline ---
sub GATK_ReAlign{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ref = $parm{ref};
    my $sBamPath = $parm{sBamPath};
    my $oBamPath = $parm{oBamPath};
    my $gatk = $parm{gatk};
    my $java = $parm{java} || 'java';
    my $jmem = $parm{jmem} || '1g';

    # GATK realn
    my $intervals = "$oBamPath.intervals";
    my $cmd = "($java -Xmx$jmem -jar $gatk -T RealignerTargetCreator -R $ref -I $sBamPath -o $intervals 2>/dev/null) && "
             ."($java -Xmx$jmem -jar $gatk -T IndelRealigner -R $ref -I $sBamPath -o $oBamPath -targetIntervals $intervals -maxInMemory 300000 -maxReads 100000 -l INFO 2>/dev/null)";
    trible_run_for_success($cmd, 'GATK_realn', {esdo_Nvb=>0});
    # sweep
    `rm -rf $intervals`;
    # inform
    stout_and_sterr "[INFO]\tGATK realn bam finished.\n"
                         ."\tI=$sBamPath\n"
                         ."\tO=$oBamPath\n";
}

1; ## tell the perl script the successful access of this module.
