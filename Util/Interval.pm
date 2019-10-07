package BioFuse::Util::Interval;

use strict;
use warnings;
use List::Util qw/ min max /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use Data::Dumper;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              Get_Two_Seg_Olen
              merge
              intersect
              exclude
              arrange_region_with_clip
              deal_circular_extended_part
              clip_region_edge
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Interval';
#----- version --------
$VERSION = "0.34";
$DATE = '2019-10-07';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        Get_Two_Seg_Olen
                        merge
                        intersect
                        exclude
                        arrange_region_with_clip
                        deal_circular_extended_part
                        clip_region_edge
                     /;

#--- Get the overlapped length of two intervals ---
# Note: only works for intervals with integer boundaries
# if want to use it for demical boundaries, inter the fifth option
sub Get_Two_Seg_Olen{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my ($s1,$e1,$s2,$e2,$offset) = @_;
    $offset = 1 unless defined $offset;
    my $overlap_len = min($e1,$e2) - max($s1,$s2) + $offset;
    return $overlap_len < 0 ? 0 : $overlap_len;
}

#--- merge intervals ---
## works on integer interval !!!
# default max_dist is 0, which means will not merge adjacent intervals, such as [10,19],[20,30]
# when set mergeAdjacent, max_dist is set as at least one, then adjacent intervals will be merged
sub merge{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $itvAfListAf = $parm{itvAfListAf};
    my $max_dist = $parm{max_dist} || 0;
    my $mergeAdjacent = $parm{mergeAdjacent} || 0;

    # set max_dist as at least one
    $max_dist = max(1, $max_dist) if $mergeAdjacent;
    # gather given intervals
    my @mergeItv;
    for my $itvAf (@$itvAfListAf){
        push @mergeItv, [ @$_ ] for @$itvAf; # copy value
    }
    # sort by end-pos
    @mergeItv = sort {$a->[1] <=> $b->[1]} @mergeItv;
    # merge overlap
    for (my $i = $#mergeItv; $i > 0; $i--){
        # check credible interval
        if($mergeItv[$i]->[0] > $mergeItv[$i]->[1]){
            warn_and_exit "<ERROR>\tindex $i has reversed interval: [@{$mergeItv[$i]}]\n";
        }
        # check overlap
        if($mergeItv[$i]->[0] <= $mergeItv[$i-1]->[1] + $max_dist){
            my $new_st = min($mergeItv[$i]->[0],$mergeItv[$i-1]->[0]);
            splice(@mergeItv, $i-1, 2, [$new_st, $mergeItv[$i]->[1]]);
        }
    }
    # sort by st-pos
    @mergeItv = sort {$a->[0] <=> $b->[0]} @mergeItv;
    # return ref of sorted merged intervals
    return \@mergeItv;
}

#--- get intersction of given intervals ---
## works on integer interval !!!
sub intersect{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $itvAfListAf = $parm{itvAfListAf};
    my $bedFormat = $parm{bedFormat} || 0;
    my $skipMerge = $parm{skipMerge} || 0;

    my $lastItscAf = [];
    # get intersection one by one
    for my $i (0 .. scalar(@$itvAfListAf)-1){
        # get sorted merged interval
        my $mergeItvAf =  $skipMerge
                        ? [$itvAfListAf->[$i]] # directly use refer
                        : &merge(itvAfListAf => [$itvAfListAf->[$i]], mergeAdjacent => !$bedFormat); # copy value in `merge` func
        # zeroBase (BED) to oneBase, if set
        if($bedFormat){
            $_->[0]++ for @$mergeItvAf;
        }
        # update intersection
        if($i == 0){
            $lastItscAf = [map{[@$_]} @{$mergeItvAf}]; # copy value
        }
        else{
            my @intersection;
            my ($i_a, $iCeil_a) = (0, scalar @$lastItscAf);
            my ($i_b, $iCeil_b) = (0, scalar @$mergeItvAf);
            while($i_a < $iCeil_a && $i_b < $iCeil_b){
                my ($st_a, $ed_a) = @{$lastItscAf->[$i_a]}[0,1];
                my ($st_b, $ed_b) = @{$mergeItvAf->[$i_b]}[0,1];
                # have overlap?
                push @intersection, [max($st_a, $st_b), min($ed_a, $ed_b)] if &Get_Two_Seg_Olen($st_a, $ed_a, $st_b, $ed_b);
                # move index
                $ed_a > $ed_b ? $i_b++ : $i_a++;
            }
            # update
            $lastItscAf = \@intersection;
        }
        # oneBase to zeroBase (BED), if set
        if($bedFormat && $skipMerge){
            $_->[0]-- for @$mergeItvAf; # as use refer
        }
        # no intersection found
        last unless scalar @$lastItscAf;
    }
    # oneBase to zeroBase (BED), if set
    if($bedFormat){
        $_->[0]-- for @$lastItscAf;
    }
    # return array ref of intersection interval
    return $lastItscAf;
}

#--- get intervals after exclusion ---
## works on integer interval !!!
## exclude [itvAfListAf_e] from [itvAfListAf_s]
sub exclude{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $itvAfListAf_s = $parm{itvAfListAf_s}; # source
    my $itvAfListAf_e = $parm{itvAfListAf_e}; # exclude
    my $bedFormat = $parm{bedFormat} || 0;

    # get sorted merged interval
    my $mergeItvAf_s = &merge(itvAfListAf => $itvAfListAf_s, mergeAdjacent => !$bedFormat);
    my $mergeItvAf_e = &merge(itvAfListAf => $itvAfListAf_e, mergeAdjacent => !$bedFormat);
    # zeroBase (BED) to oneBase, if set
    if($bedFormat){
        $_->[0]++ for @$mergeItvAf_s;
        $_->[0]++ for @$mergeItvAf_e;
    }
    # get intersection
    my $itsctItvAf_e = &intersect(itvAfListAf => [$mergeItvAf_s, $mergeItvAf_e]);
    # do exclusion
    my @remainItv;
    my ($i_s, $iCeil_s) = (0, scalar @$mergeItvAf_s);
    my ($i_e, $iCeil_e) = (0, scalar @$itsctItvAf_e);
    while($i_s < $iCeil_s && $i_e < $iCeil_e){
        my ($st_s, $ed_s) = @{$mergeItvAf_s->[$i_s]}[0,1];
        my ($st_e, $ed_e) = @{$itsctItvAf_e->[$i_e]}[0,1];
        # have overlap?
        my $overlap = &Get_Two_Seg_Olen($st_s, $ed_s, $st_e, $ed_e);
        if($st_s <= $st_e){
            if($overlap){
                push @remainItv, [$st_s, $st_e - 1] if $st_s < $st_e;
                $mergeItvAf_s->[$i_s]->[0] = $ed_e + 1 if $ed_e < $ed_s;
            }
            else{
                push @remainItv, [$st_s, $ed_s];
            }
        }
        # move index
        $ed_s > $ed_e ? $i_e++ : $i_s++;
    }
    # if source has left intervals
    push @remainItv, [@{$mergeItvAf_s->[$_]}] for ($i_s .. $iCeil_s-1);
    # oneBase to zeroBase (BED), if set
    if($bedFormat){
        $_->[0]-- for @remainItv;
    }
    # return array ref of remained interval
    return \@remainItv;
}

#--- arrange one region in expected region ---
sub arrange_region_with_clip{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $lftPos = $parm{lftPos}; # left-edge of the original region
    my $rgtpos = $parm{rgtpos}; # right-edge of the original region
    my $min_Pos = $parm{min_Pos}; # minimum for clip
    my $max_Pos = $parm{max_Pos}; # maximum for clip
    my $cirl_oLen = $parm{cirl_oLen}; # circular position

    # deal with circular extended part, if possible
    my @candReg_Aref = ( defined $cirl_oLen )
                      ? &deal_circular_extended_part( lftPos => $lftPos, rgtpos => $rgtpos, origLen => $cirl_oLen )
                      : ( [ $lftPos, $rgtpos ] );
    # clip region
    my @Region;
    for my $Aref ( @candReg_Aref ){
        $Aref->[0] = &clip_region_edge( LR_Edge => 'L', edgePos => $Aref->[0], min_Pos => $min_Pos, max_Pos => $max_Pos );
        $Aref->[1] = &clip_region_edge( LR_Edge => 'R', edgePos => $Aref->[1], min_Pos => $min_Pos, max_Pos => $max_Pos );
        # record
        push @Region, $Aref if( $Aref->[0] < $Aref->[1] );
    }

    return @Region;
}

#--- deal with the extended region ---
## e.g., for the circular virus
sub deal_circular_extended_part{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $lftPos = $parm{lftPos};
    my $rgtpos = $parm{rgtpos};
    my $origLen = $parm{origLen};

    if( $lftPos > $origLen ){
        return [ $lftPos - $origLen, $rgtpos - $origLen ];
    }
    elsif(   $lftPos < $origLen
          && $rgtpos > $origLen
    ){
        return ( [ $lftPos,          $origLen ],
                 [      1, $rgtpos - $origLen ]  );
    }
    else{
        return [ $lftPos, $rgtpos ];
    }
}

#--- clip region edge based on allowed min/max position ---
sub clip_region_edge{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $edgePos = $parm{edgePos}; # the original edge position
    my $LR_Edge = $parm{LR_Edge}; # which edge
    my $min_Pos = $parm{min_Pos}; # minimum for clip
    my $max_Pos = $parm{max_Pos}; # maximum for clip

    if( $LR_Edge =~ /^L/i ){
        return min( max( $edgePos, $min_Pos ), $max_Pos );
    }
    else{
        return max( min( $edgePos, $max_Pos ), $min_Pos );
    }
}

1; ## tell the perl script the successful access of this module.
