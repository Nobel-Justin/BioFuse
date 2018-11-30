package BioFuse::BioInfo::BED;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ stout_and_sterr warn_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Array qw/ mergeOverlap /;

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
$DATE = '2018-11-29';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        read_bed_file
                     /;

#--- read bed file ---
sub read_bed_file{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $BedFile = $parm{BedFile};
    my $ItvHref = $parm{ItvHref};
    my $nonName = $parm{nonName} || 0;
    my $oneBase = $parm{oneBase} || 0; # not bed, but start from one
    my $skipMerge = $parm{skipMerge} || 0;

    # check bed file existence
    file_exist(filePath => $BedFile, alert => 1);

    # read bed file and record region
    my $offset = $oneBase ? 0 : 1;
    open (BED, Try_GZ_Read($BedFile)) || die "fail to read bed file: $!\n";
    while(<BED>){
        next if /^#/;
        my @ele = split;
        my ($refseg, $stpos, $edpos) = @ele[0,1,2];
        # record
        if($nonName){
            push @{$ItvHref->{$refseg}}, [$stpos+$offset, $edpos+0];
        }
        else{
            my $name = $ele[3] || '_UNDEF_';
            push @{$ItvHref->{$name}->{$refseg}}, [$stpos+$offset, $edpos+0];
        }
    }
    close BED;

    return if $skipMerge;
    # prepare
    my @ArefToMerge =   $nonName
                      ? values %$ItvHref
                      : map {values %$_} values %$ItvHref;
    # sort and merge
    ## bed format cannot merge adjacent interval, such as (0,4) and (5,6).
    ## (0,4): pos 0 to 3; while (5,6): pos 5 to 5
    mergeOverlap(regionAref => $_, mergeAdjacent => $oneBase) for @ArefToMerge;
}

1; ## tell the perl script the successful access of this module.
