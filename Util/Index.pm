package BioFuse::Util::Index;

use strict;
use warnings;
use Data::Dumper;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              Pos2Idx
              IndexRegion
              FindOverlapIdxRegion
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Index';
#----- version --------
$VERSION = "0.05";
$DATE = '2018-11-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        Pos2Idx
                        IndexRegion
                        FindOverlapIdxRegion
                     /;

#--- get posWinIdx for given pos ---
sub Pos2Idx{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $winSize = $parm{winSize} || 1000;
    return int( $parm{pos} / $winSize );
}

#--- return idx-region Href of given region Aref ---
sub IndexRegion{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ItvAref = $parm{ItvAref};
    my $winSize = $parm{winSize} || 1000;

    my $IdxItvHref = {};
    for my $Itv (@$ItvAref){
        my @PosIdx = map { ( &Pos2Idx(pos => $_, winSize => $winSize) ) } @$Itv;
        shift @PosIdx if($PosIdx[0] == $PosIdx[1]);
        push @{$IdxItvHref->{$_}}, $Itv for ($PosIdx[0] .. $PosIdx[-1]);
    }
    return $IdxItvHref;
}

#--- query idx-region overlap with given region ---
## note that 'queryMode' could be one of pos/itv/ob_pos/ob_itv
sub FindOverlapIdxRegion{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $IdxItvHref = $parm{IdxItvHref};
    my $regionAref = $parm{regionAref};
    my $winSize = $parm{winSize}; # should give by user
    my $queryMode = $parm{queryMode};

    my @overlapObj = ();
    my @PosIdx = map { ( &Pos2Idx(pos => $_, winSize => $winSize) ) } @$regionAref;
    shift @PosIdx if($PosIdx[0] == $PosIdx[1]);
    # inner PosIdx
    push @overlapObj, @{$IdxItvHref->{$_}} for grep exists $IdxItvHref->{$_}, ($PosIdx[0]+1 .. $PosIdx[-1]-1);
    # boundary PosIdx
    my $posMode = ($queryMode =~ /pos$/); # pos or itv
    my $objMode = ($queryMode =~ /^ob_/); # obj-func ?
    for my $PosIdx ( grep exists $IdxItvHref->{$_}, @PosIdx ){
        # coordinates-sorted by pos or itv_st in ascending order
        for my $Obj ( @{$IdxItvHref->{$PosIdx}} ){
            if($posMode){ # single-pos
                my $pos = $objMode ? $Obj->get_pos : $Obj; # scalar
                next if $pos < $regionAref->[0];
                last if $pos > $regionAref->[1];
            }
            else{ # interval
                my $itv = $objMode ? $Obj->get_itv : $Obj; # array ref
                next if $itv->[1] < $regionAref->[0];
                last if $itv->[0] > $regionAref->[1];
            }
            # record
            push @overlapObj, $Obj;
        }
    }
    # return
    return \@overlapObj;
}

#--- 
1; ## tell the perl script the successful access of this module.
