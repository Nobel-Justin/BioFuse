package BioFuse::Visual::ArrangeObj;

use strict;
use warnings;
use POSIX qw/ceil/;
use List::Util qw/ max min sum /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Visual::SVG_Util::Color qw/ %COLOR_DB /;
use BioFuse::Visual::SVG_Util::Font qw/ show_text_in_line /;
use BioFuse::Visual::Windows qw/ get_pos_windowNO  /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              allocate_object_bilater
              allocate_object_vertical
              test_overlap_at_one_layer
              get_ObjSpanPixItval
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::ArrangeObj';
#----- version --------
$VERSION = "0.10";
$DATE = '2020-11-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        allocate_object_bilater
                        allocate_object_vertical
                        test_overlap_at_one_layer
                        get_ObjSpanPixItval
                     /;

#--- arrange position of object bilaterally to avoid overlap ---
sub allocate_object_bilater{

    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $obj_Href = $parm{obj_Href};
    my $objIdx_key = $parm{objIdx_key};
    my $newLoc_key = $parm{newLoc_key};
    my $objWid_key = $parm{objWid_key};
    my $obj_gap  = $parm{obj_gap} || 5;
    my $adj_ratio = $parm{adj_ratio} || 0.25;
    my $timeToRpt = $parm{timeToRpt} || 1000;
    my $maxTryTime = $parm{maxTryTime} || 3E5;
    my $disable_L = $parm{disable_L} || 0;
    my $disable_R = ($disable_L) ? 0 : ($parm{disable_R} || 0);

    my $adjust_step = $obj_gap * $adj_ratio;
    my @sorted_objIdx_by_loc = sort {$obj_Href->{$a}->{$objIdx_key}<=>$obj_Href->{$b}->{$objIdx_key}} keys %$obj_Href;
    my $start_idx = 0;
    my $arrange_time = 0;
    ARRANGE: {
        if(++$arrange_time % $timeToRpt == 1){
            stout_and_sterr "$arrange_time time arrangement operation pass. Start from index $start_idx\n";
        }
        if($arrange_time > $maxTryTime){
            stout_and_sterr "$arrange_time time exceeds maximum ($maxTryTime), stop allocation operation now.\n";
            return;
        }
        for my $idx ($start_idx  .. $#sorted_objIdx_by_loc){
            my $objIdx = $sorted_objIdx_by_loc[$idx];
            my $newLoc = $obj_Href->{$objIdx}->{$newLoc_key};
            my $objWid = $obj_Href->{$objIdx}->{$objWid_key};
            my ($L_edge, $R_edge) = ($newLoc - $objWid / 2, $newLoc + $objWid / 2);
            # check the overlap with left object
            if(    !$disable_L
                && $idx - 1 >= 0
            ){
                my $left_objIdx = $sorted_objIdx_by_loc[$idx-1];
                my $left_newLoc = $obj_Href->{$left_objIdx}->{$newLoc_key};
                my $left_locWid = $obj_Href->{$left_objIdx}->{$objWid_key};
                my ($left_L_edge, $left_R_edge) = ($left_newLoc - $left_locWid / 2, $left_newLoc + $left_locWid / 2);
                if($L_edge - $left_R_edge < $obj_gap){
                    $obj_Href->{$objIdx}->{$newLoc_key} += $adjust_step;
                    $obj_Href->{$left_objIdx}->{$newLoc_key} -= $adjust_step;
                    $start_idx = $idx - 1;
                    redo ARRANGE;
                }
            }
            # check the overlap with right object
            if(    !$disable_R
                && $idx + 1 <= $#sorted_objIdx_by_loc
            ){
                my $right_objIdx = $sorted_objIdx_by_loc[$idx+1];
                my $right_newLoc = $obj_Href->{$right_objIdx}->{$newLoc_key};
                my $right_locWid = $obj_Href->{$right_objIdx}->{$objWid_key};
                my ($right_L_edge, $right_R_edge) = ($right_newLoc - $right_locWid / 2, $right_newLoc + $right_locWid / 2);
                if($right_L_edge - $R_edge < $obj_gap){
                    $obj_Href->{$objIdx}->{$newLoc_key} -= $adjust_step;
                    $obj_Href->{$right_objIdx}->{$newLoc_key} += $adjust_step;
                    $start_idx = $idx + 1;
                    redo ARRANGE;
                }
            }
        }
    }
}

#--- arrange layer of object vertically to avoid overlap ---
sub allocate_object_vertical{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $layer_href = $parm{layer_href};
       $layer_href = {1=>[]} if( !defined $layer_href );
    my $obj_Href = $parm{obj_Href};
    my $obj_itval_key = $parm{obj_itval_key};
    my $obj_layer_key = $parm{obj_layer_key};
    my $obj_gap = $parm{obj_gap} || 5;
    my $max_layer_Sref = $parm{max_layer_Sref};

    # sort by left position
    my @sorted_obj = sort { 
                            $obj_Href->{$a}->{$obj_itval_key}->[0]->[0]
                            <=>
                            $obj_Href->{$b}->{$obj_itval_key}->[0]->[0]
                          } keys %$obj_Href;
    # no objects
    return if( scalar(@sorted_obj) == 0 );
    # allocation
    for my $obj_name ( @sorted_obj ){
        my $obj_info_Href = $obj_Href->{$obj_name};

        #--- find the proper layer
        my $ok_layer_NO = 0;
        for my $layer_NO (sort {$a<=>$b} keys %$layer_href){
            my $ok_bool = ! &test_overlap_at_one_layer(
                                                        Sbjct_Reg_Aref => $layer_href->{$layer_NO},
                                                        Query_Reg_Aref => $obj_info_Href->{$obj_itval_key},
                                                        obj_gap => $obj_gap
                                                      );
            if( $ok_bool ){
                $ok_layer_NO = $layer_NO;
                last;
            }
            else{
                $ok_layer_NO = $layer_NO + 1;
            }
        }
        # recode the layer of this obj
        $obj_info_Href->{$obj_layer_key} = $ok_layer_NO;
        # recode usage at this layer
        push @{ $layer_href->{$ok_layer_NO} }, @{ $obj_info_Href->{$obj_itval_key} };
    }
    # update
    $$max_layer_Sref = max( keys %$layer_href );
}

#--- test object overlap at one layer ---
## if find overlap, return TURE
sub test_overlap_at_one_layer{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $Sbjct_Reg_Aref = $parm{Sbjct_Reg_Aref};
    my $Query_Reg_Aref = $parm{Query_Reg_Aref};
    my $obj_gap = $parm{obj_gap} || 5;

    my $no_ovp_bool = 1;
    for my $occu_reg_Aref ( @$Sbjct_Reg_Aref ){
        my $occu_reg_lftEdge = $occu_reg_Aref->[0] - $obj_gap;
        my $occu_reg_rgtEdge = $occu_reg_Aref->[1] + $obj_gap;
        for my $obj_reg_Aref ( @$Query_Reg_Aref ){
            my $ovp_len = Get_Two_Seg_Olen( $occu_reg_lftEdge, $occu_reg_rgtEdge, $obj_reg_Aref->[0], $obj_reg_Aref->[1] );
            if( $ovp_len != 0 ){ # find overlap at this layer
                $no_ovp_bool = 0;
                last;
            }
        }
        # stop
        last if( $no_ovp_bool == 0 );
    }

    return !$no_ovp_bool;
}

#--- calculate the span interval of each object ---
sub get_ObjSpanPixItval{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ObjInfoPool_Href = $parm{ObjInfoPool_Href};
    my $pos2xAxisPos_para_Aref = $parm{pos2xAxisPos_para_Aref};
    my $axisZX = $parm{axisZX};
       $axisZX = 100 if( !defined $axisZX );
    my $labelFtsz = $parm{labelFtsz} || 12;
    my $labelFtfm = $parm{labelFtfm} || 'Arial';
    my $labelHeight = $parm{labelHeight} || 15;
    my $bodyLabelGap = $parm{bodyLabelGap} || 5;
    my $bodyRegAref_key = $parm{bodyRegAref_key};
    my $bodyPixItval_key = $parm{bodyPixItval_key};
    my $bodyWithLabelPixItval_key = $parm{bodyWithLabelPixItval_key};
    my $labelText_key = $parm{labelText_key};
    my $bodyCol_key = $parm{bodyCol_key};
    my $bodyColRef_key = $parm{bodyColRef_key}; # for assign color if set, or else auto

    my %col_usage;
    for my $obj_name ( sort keys %$ObjInfoPool_Href ){
        my $objInfo_Href = $ObjInfoPool_Href->{$obj_name};
        # get object X-axis display interval
        for my $bodyRegAref (@{$objInfo_Href->{$bodyRegAref_key}}){
            # draw interval of object body
            my $body_lftPix = $axisZX + get_pos_windowNO( pos => $bodyRegAref->[0], @$pos2xAxisPos_para_Aref );
            my $body_rgtPix = $axisZX + get_pos_windowNO( pos => $bodyRegAref->[1], @$pos2xAxisPos_para_Aref );
            push @{$objInfo_Href->{$bodyPixItval_key}}, [ $body_lftPix, $body_rgtPix ];
            # label interval
            if(    defined $bodyWithLabelPixItval_key
                && defined $labelText_key
            ){
                my $temp_SVG_obj = 'temp';
                my $text_size_Aref = show_text_in_line(
                                                        svg_obj => \$temp_SVG_obj,
                                                        text_x => $body_lftPix - $bodyLabelGap,
                                                        text_y => 100, # temp
                                                        text => $objInfo_Href->{$labelText_key},
                                                        font_family => $labelFtfm,
                                                        font_size => $labelFtsz,
                                                        text_col => 'black',
                                                        text_anchor => 'end',
                                                        height_adjust => 1,
                                                        height_limit => $labelHeight,
                                                        width_limit => 0,
                                                        rotate_degree => 0,
                                                        draw_bool => 0, # do not draw
                                                    );
                # full width of this object info
                my $bodyLabel_lftPix = $body_lftPix - $bodyLabelGap - $text_size_Aref->[1];
                my $bodyLabel_rgtPix = $body_rgtPix;
                push @{$objInfo_Href->{$bodyWithLabelPixItval_key}}, [ $bodyLabel_lftPix, $bodyLabel_rgtPix ];
            }
        }
        # set object color
        if( defined $bodyCol_key ){
            my $objBody_col;
            if ( defined $bodyColRef_key ){
                my $bodyColRef = $objInfo_Href->{$bodyColRef_key};
                if( exists $col_usage{$bodyColRef} ){
                    $objBody_col = $col_usage{$bodyColRef};
                }
                else{
                    my %bodyColUsed  = map { ($_,1) } values %col_usage;
                    my @COLOR_DB_idx = map { $COLOR_DB{$_} } sort {$a<=>$b} grep /^\d+$/ && !exists($bodyColUsed{$COLOR_DB{$_}}), keys %COLOR_DB;
                    if( @COLOR_DB_idx == 0 ){
                        @COLOR_DB_idx = map { $COLOR_DB{$_} } sort {$a<=>$b} grep /^\d+$/, keys %COLOR_DB;
                        $objBody_col = $COLOR_DB_idx[int(rand(scalar(@COLOR_DB_idx)))];
                    }
                    else{
                        $objBody_col = $COLOR_DB_idx[0];
                    }
                    # record
                    $col_usage{$bodyColRef} = $objBody_col;
                }
            }
            else{
                my @COLOR_DB_idx = grep /^\d+$/, sort keys %COLOR_DB;
                my $count = sum( values %col_usage ) || 0;
                $objBody_col = $COLOR_DB{ $COLOR_DB_idx[ $count % (scalar(@COLOR_DB_idx) || 1) ] };
                # record
                $col_usage{$objBody_col} ++;
            }
            # assign color
            $objInfo_Href->{$bodyCol_key} = $objBody_col;
        }
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
