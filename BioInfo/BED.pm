package BioFuse::BioInfo::BED;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ stout_and_sterr warn_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Array qw/ merge /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              read_bed_file
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::BED';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-12-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        read_bed_file
                     /;

#--- read bed file ---
## load and convert to one-start interval
sub read_bed_file{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $bedFile = $parm{bedFile};
    my $nonName = $parm{nonName} || 0;
    my $oneBase = $parm{oneBase} || 0; # not bed, but start from one
    my $loadAsBED = $parm{loadAsBED} || 0; # record interval in zero-start format
    my $skipMerge = $parm{skipMerge} || 0;

    # check bed file existence
    file_exist(filePath => $bedFile, alert => 1);

    my $ItvHref = {};
    # offset to record interval
    my $offset = 0;
    $offset = -1 if  $oneBase &&  $loadAsBED;
    $offset =  1 if !$oneBase && !$loadAsBED;
    # read bed file and record intervals
    open (BED, Try_GZ_Read($bedFile)) || die "fail to read bed file: $!\n";
    while(<BED>){
        next if /^#/;
        my @ele = split;
        my ($refseg, $stpos, $edpos) = @ele[0,1,2];
        # record as one-start interval
        if($nonName){
            push @{$ItvHref->{$refseg}}, [$stpos+$offset, $edpos+0];
        }
        else{
            my $name = $ele[3] || '_UNDEF_';
            push @{$ItvHref->{$name}->{$refseg}}, [$stpos+$offset, $edpos+0];
        }
    }
    close BED;

    # merge interval
    unless($skipMerge){
        # prepare
        my @ArefToMerge =   $nonName
                          ? values %$ItvHref
                          : map {values %$_} values %$ItvHref;
        # sort and merge
        @{$_} = @{merge(itvAfListAf => [$_], mergeAdjacent => !$loadAsBED)} for @ArefToMerge;
    }

    return $ItvHref;
}

1; ## tell the perl script the successful access of this module.
