package BioFuse::BioInfo::Position;

use strict;
use warnings;
use List::Util qw/ min /;
use BioFuse::BioInfo::BED qw/ read_bed_file /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              load_region_for_randPos
              get_random_pos
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Position';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-12-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        load_region_for_randPos
                        get_random_pos
                     /;

#--- load region for random position selection ---
sub load_region_for_randPos{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ItvFile = $parm{ItvFile};
    my $ItvAref = $parm{ItvAref};
    my $oneBase = $parm{oneBase} || 0;
    my $winSize = $parm{winSize} || 1000;

    # as bed file with one-base
    my $ItvHref = read_bed_file( bedFile => $ItvFile,
                                 nonName => 1,
                                 oneBase => $oneBase,
                                 loadAsBED => 0
                               );
    # convert to index region array for next random selection equally
    for my $refseg (keys %$ItvHref){
        for my $itvAf (@{$ItvHref->{$refseg}}){
            my ($stpos, $edpos) = @$itvAf;
            if($edpos-$stpos+1 <= $winSize){
                push @$ItvAref, [$refseg, $stpos, $edpos];
            }
            else{
                for (; $stpos < $edpos; $stpos += $winSize){
                    my $this_ed = min($stpos+$winSize-1, $edpos);
                    push @$ItvAref, [$refseg, $stpos, $this_ed];
                }
            }
        }
    }
}

#--- randomly pick certain positions from given interval ---
## such as from nonN region of reference
sub get_random_pos{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ItvAref = $parm{ItvAref};
    my $randCnt = $parm{randCnt};

    my @randPos;
    my $IdxCnt = scalar @$ItvAref;
    for (1 .. $randCnt){
        my $randItvAf = $ItvAref->[int(rand($IdxCnt))];
        my $randPos = $randItvAf->[1] + int(rand($randItvAf->[2]-$randItvAf->[1]+1));
        push @randPos, [$randItvAf->[0], $randPos];
    }
    return \@randPos;
}

1; ## tell the perl script the successful access of this module.
