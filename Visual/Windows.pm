package BioFuse::Visual::Windows;

use strict;
use warnings;
use POSIX qw/ceil/;
use List::Util qw/ max min sum /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              load_regResolForWindowlize
              get_pos_windowNO
              windowNO_to_posInterval
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::Windows';
#----- version --------
$VERSION = "0.09";
$DATE = '2018-12-03';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        load_regResolForWindowlize
                        get_pos_windowNO
                        windowNO_to_posInterval
                     /;

# assume regResol_Href has such structure:
# regResol_Href = {}
# regResol_Href -> $refseg -> $st_pos -> 'st_pos' = $st_pos
# regResol_Href -> $refseg -> $st_pos -> 'ed_pos' = $ed_pos
# regResol_Href -> $refseg -> $st_pos -> 'resol'  = $resol

#--- load region resolution for further windowlize ---
sub load_regResolForWindowlize{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $reg_resol_file = $parm{reg_resol_file};
    my $regResol_Href = $parm{regResol_Href}; # to store resolution of region
    # destReg_Href -> {}
    # destReg_Href -> refseg = $refseg
    # destReg_Href -> st_pos = $st_pos
    # destReg_Href -> ed_pos = $ed_pos
    my $destReg_Href = $parm{destReg_Href}; # for check completeness
    my $destRegAutoFigureSize = $parm{destRegAutoFigureSize}; # || (defined $destReg_Href ? 500 : 0);
    my $autoAdjustResol = $parm{autoAdjustResol} || $destRegAutoFigureSize;

    my %drawReg_edgeShowBool = (st_pos=>0, ed_pos=>0);

    open (RES, Try_GZ_Read($reg_resol_file)) || die "fail read resolution file: $!\n";
    while (<RES>){
        next if(/^#/);
        my ($reg_seg, $reg_st, $reg_ed, $resolution) = (split)[0..3];
        # required segment of region
        if(    defined $destReg_Href
            && $reg_seg ne $destReg_Href->{refseg}
        ){
            next;
        }
        # simple check
        if(    $reg_st > $reg_ed
            || $resolution <= 0
        ){
            warn_and_exit "<ERRO>\tWrong data in resolution file.\n"
                                 ."\t$_"
                                 ."\treg_resol_file\n";
        }
        # clip region with destReg
        $reg_st = max( $reg_st, $destReg_Href->{st_pos} );
        $reg_ed = min( $reg_ed, $destReg_Href->{ed_pos} );
        next if($reg_st >= $reg_ed);
        # adjust resolution
        my $orig_resolution = $resolution;
        if( $autoAdjustResol ){
            my $adjust_resol = $resolution;
            my $regSize = $reg_ed - $reg_st + 1;
            my $doTime = 1;
            my $ratio = 3;
            # find proper resol
            if( $resolution >= 1 ){
                $adjust_resol = $regSize;
                $adjust_resol = int( sqrt( $adjust_resol ) * $ratio ) for ( 1 .. $doTime );
            }
            # minimum regShowSize
            my $MinRegShowSize = 15;
            $adjust_resol = min( $adjust_resol, (sprintf "%.2f", $regSize / $MinRegShowSize) );
            # minimum resol is at least 1.0
            $adjust_resol = max( $adjust_resol, 1 );
            # update
            if( $resolution != $adjust_resol ){
                $resolution = $adjust_resol;
            }
        }
        # record
        if( defined $destReg_Href ){
            # only record for destReg
            if( Get_Two_Seg_Olen( $reg_st, $reg_ed, $destReg_Href->{st_pos}, $destReg_Href->{ed_pos} ) ){
                $regResol_Href->{$reg_seg}->{$reg_st} = { st_pos => $reg_st, ed_pos=>$reg_ed, resol=>$resolution, orig_resol=>$orig_resolution };
            }
            else{
                warn "<WARN>\tIgnore $reg_seg:$reg_st-$reg_ed from resolution file.\n";
            }
            # mark destReg edge found or not
            for my $edge_key (qw/ st_pos ed_pos /){
                if(    $reg_st  <= $destReg_Href->{$edge_key}
                    && $reg_ed  >= $destReg_Href->{$edge_key}
                ){
                    $drawReg_edgeShowBool{$edge_key} = 1;
                }
            }
        }
        else{
            $regResol_Href->{$reg_seg}->{$reg_st} = { st_pos => $reg_st, ed_pos=>$reg_ed, resol=>$resolution, orig_resol=>$orig_resolution };
        }
    }
    close RES;

    # destReg work
    if( defined $destReg_Href ){
        # check destReg is whole covered
        if( $drawReg_edgeShowBool{st_pos} * $drawReg_edgeShowBool{ed_pos} == 0 ){
            warn_and_exit "<ERROR>\tCannot find draw region edge-positions in resolution file.\n";
        }
        # resolution auto adjustment for whole destReg
        if( $autoAdjustResol ){
            # displayed pixel interval size
            my $destRegEdPos_winNO = get_pos_windowNO( 
                                                        pos => $destReg_Href->{ed_pos},
                                                        refseg => $destReg_Href->{refseg},
                                                        regResol_Href => $regResol_Href,
                                                        stpos_1stwin => $destReg_Href->{st_pos}
                                                    );
            # adjust resol in proportion
            my $ratio = $destRegEdPos_winNO / $destRegAutoFigureSize;
            for my $reg_seg (sort keys %$regResol_Href){
                $regResol_Href->{$reg_seg}->{$_}->{resol} = max( 1, sprintf "%.2f", $regResol_Href->{$reg_seg}->{$_}->{resol} * $ratio ) for keys %{$regResol_Href->{$reg_seg}};
            }
        }
    }

    # check region linkage, region should connected with each other precisely
    for my $reg_seg (sort keys %$regResol_Href){
        my @reg_st = sort {$a<=>$b} keys %{$regResol_Href->{$reg_seg}};
        for my $i (0 .. $#reg_st){
            # check connection with previous region
            my ($st_now, $ed_now) = ( $reg_st[$i],   $regResol_Href->{$reg_seg}->{$reg_st[$i]}->{ed_pos} );
            if( $i != 0 ){
                my ($st_pre, $ed_pre) = ( $reg_st[$i-1], $regResol_Href->{$reg_seg}->{$reg_st[$i-1]}->{ed_pos} );
                if( $st_now != $ed_pre + 1 ){
                    warn_and_exit "<ERROR>\tfind non-linked regions on $reg_seg: [$st_pre, $ed_pre] and [$st_now, $ed_now]\n"
                                         ."\t'pre_end_Pos($ed_pre)' must equal to 'now_start_Pos($st_now) - 1'.\n";
                }
            }
            # report resolution if changes.
            my $used_resol = $regResol_Href->{$reg_seg}->{$reg_st[$i]}->{resol};
            my $orig_resol = $regResol_Href->{$reg_seg}->{$reg_st[$i]}->{orig_resol};
            if( $used_resol != $orig_resol ){
                stout_and_sterr "[INFO]\tAdjust resolution of region ($reg_seg:$st_now-$ed_now) from $orig_resol to $used_resol\n";
            }
        }
    }
}

#--- get window NO of given position ---
sub get_pos_windowNO{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $refseg = $parm{refseg};
    my $pos = $parm{pos};
    my $regResol_Href = $parm{regResol_Href};
    my $stpos_1stwin = $parm{stpos_1stwin};
    my $max_winNO = $parm{maxWinNOtoDraw};
    my $give_winLastPos = $parm{give_winLastPos}; # will give the last pos of the window
    my $winResol_Sref = $parm{winResol_Sref};

    if( $pos < $stpos_1stwin ){
        # warn "The pos to determine location is smaller than start position.\n$pos < $stpos_1stwin\n";
        return 1;
    }

    my $window_NO = 0;
    my $sameWinLastPos = 0;
    for my $reg_st (sort {$a<=>$b} keys %{$regResol_Href->{$refseg}} ){
        my $reg_ed = $regResol_Href->{$refseg}->{$reg_st}->{ed_pos};
        my $resol  = $regResol_Href->{$refseg}->{$reg_st}->{resol};
        next if($reg_ed < $stpos_1stwin);
        my $cal_st = max($reg_st, $stpos_1stwin);
        my $cal_ed = min($reg_ed, $pos);
        my $window_count = POSIX::ceil( ($cal_ed - $cal_st + 1) / $resol );
        $window_NO += $window_count;
        # found the point
        if($cal_ed == $pos){
            $sameWinLastPos = min( $reg_ed, int($cal_st+$window_count*$resol-1) );
            $$winResol_Sref = $resol if(defined $winResol_Sref);
            last;
        }
    }

    # reset to the first window
    $window_NO = 1 if($window_NO <= 0);
    # max window NO control
    $window_NO = $max_winNO if(defined($max_winNO) && $window_NO > $max_winNO);

    return $give_winLastPos
           ? ($sameWinLastPos, $window_NO)
           : ($window_NO);
}

#--- get position range represented by given window NO ---
sub windowNO_to_posInterval{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $refseg = $parm{refseg};
    my $windowNO = $parm{windowNO};
    my $regResol_Href = $parm{regResol_Href};
    my $stpos_1stwin = $parm{stpos_1stwin};

    if( $windowNO <= 0 ){
        warn_and_exit "Negative window_NO is not allowed in 'windowNO_to_posInterval' func.\n";
    }

    my $window_sum = 0;
    for my $reg_st (sort {$a<=>$b} keys %{$regResol_Href->{$refseg}} ){
        my $reg_ed = $regResol_Href->{$refseg}->{$reg_st}->{ed_pos};
        my $resol  = $regResol_Href->{$refseg}->{$reg_st}->{resol};
        next if($reg_ed < $stpos_1stwin);
        my $cal_st = max($reg_st, $stpos_1stwin);
        my $cal_ed = $reg_ed;
        my $window_count = POSIX::ceil( ($cal_ed - $cal_st + 1) / $resol );
        if( $window_sum + $window_count >= $windowNO ){
            my $itval_stPos = max( int( $reg_st + ($windowNO - $window_sum - 1) * $resol ), $stpos_1stwin );
            my $itval_edPos = max( int( $itval_stPos + $resol - 1 ), $itval_stPos );
            return ( $itval_stPos, $itval_edPos );
        }
        else{
            $window_sum += $window_count;
        }
    }

    warn_and_exit "<ERROR>\tCannot find pos-interval of given window_NO ($windowNO) via region-resol table in 'windowNO_to_posInterval' func.\n"
                        ."\treached sum window: $window_sum\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
