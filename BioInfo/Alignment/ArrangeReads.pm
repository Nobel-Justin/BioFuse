package BioFuse::BioInfo::Alignment::ArrangeReads;

use strict;
use warnings;
use List::Util qw/ max first /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              ArrangeReadsAmongRefseg
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ArrangeReads';
#----- version --------
$VERSION = "0.02";
$DATE = '2021-12-13';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        ArrangeReadsAmongRefseg
                     /;

#--- arrange reads among multiple refseg ---
# !!! the source bam should be Nsort (after merged by multiple refseg alignments)
sub ArrangeReadsAmongRefseg{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $nSortBam = $parm{nSortBam}; # bam object
    my $outBam = $parm{outBam}; # bam object
    my $onlyPropMap = $parm{onlyPropMap} || 0; # only properly mapped
    my $skipSuppMap = $parm{skipSuppMap} || 0; # skip supplement/second alignment
    my $skipDupRead = $parm{skipDupRead} || 0; # skip PCR duplicated reads
    my $skipPEunmap = $parm{skipPEunmap} || 0; # skip PE-reads that both unmap
    # refseg need to keep SQ of bam header
    ## for mapped reads, only keep alignment on this refseg
    my $keepRefsegID = $parm{keepRefsegID} || undef;
    # to allow the keepRefsegID is not the only refseg with max AS
    my $keepRefNsolo = $parm{keepRefNsolo} || 0;
    # swap when PE ends map to diff-refseg respectively
    my $swapAlignEnds = $parm{swapAlignEnds} || 0;
    # !! effective when keepRefsegID is not provided
       $swapAlignEnds = 0 if defined $keepRefsegID;

    # start writing
    $outBam->start_write;
    # header
    if(defined $keepRefsegID){
        $outBam->write(content=>$_) for grep !/^\@SQ/ || /\sSN:$keepRefsegID\s/, @{$nSortBam->header_Af};
    }
    else{
        $outBam->write(content=>$_) for @{$nSortBam->header_Af};
    }
    # arrange reads
    my $viewOpt = '';
       $viewOpt .= "-f 0x2" if $onlyPropMap;
       $viewOpt .= "-F 0x900" if $skipSuppMap; # -F 0x800 + 0x100
       $viewOpt .= "-F 0x400" if $skipDupRead;
    $nSortBam->smartBam_PEread(viewOpt => $viewOpt, deal_peOB_pool => 1, subrtRef => \&SelectRefsegForReads, subrtParmAref => \@_, quiet => 1);
    # stop writing
    $outBam->stop_write;
}

#--- arrange reads among multiple refseg ---
sub SelectRefsegForReads{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $outBam = $parm{outBam}; # bam object
    my $skipPEunmap = $parm{skipPEunmap} || 0; # skip PE-reads that both unmap
    # for mapped reads, only keep alignment on this refseg
    my $keepRefsegID = $parm{keepRefsegID} || undef;
    # to allow the keepRefsegID is not the only refseg with max AS
    my $keepRefNsolo = $parm{keepRefNsolo} || 0;
    # swap when PE ends map to diff-refseg respectively
    my $swapAlignEnds = $parm{swapAlignEnds} || 0;
    # !! effective when keepRefsegID is not provided
       $swapAlignEnds = 0 if defined $keepRefsegID;

    for my $pe_OB (@$pe_OB_poolAf){
        my %rOB_Af = map{ ($_, $pe_OB->rOB_Af(reads_end=>$_)) } (1,2);
        # mapped rOB
        my %map_rOB = map { ($_, [grep !$_->is_unmap, @{$rOB_Af{$_}}]) } (1,2);
        my %endsMap = map { ($_, scalar(@{$map_rOB{$_}})) } (1,2);
        # different scenarios
        if( !($endsMap{1} or $endsMap{2}) ){ # _1 and _2 are all unmap
            next if $skipPEunmap;
            for my $end (1,2){ # unmap, just output one time (use first)
                $outBam->write(content=>$_->printSAM."\n") for first {$_->is_unmap} @{$rOB_Af{$end}};
            }
        }
        else{
            # alignment score of each refseg
            my %refseg2AS;
            $refseg2AS{$_->mseg} += $_->alignScore for map {@{$rOB_Af{$_}}} (1,2);
            # avail refseg with max alignment score
            my $max_AS = max(values %refseg2AS);
            my @max_AS_refseg = grep $refseg2AS{$_} == $max_AS, keys %refseg2AS;
            # if keepRefsegID is set
            if(defined $keepRefsegID){
                if($keepRefNsolo){ # allow the keepRefsegID is not the only refseg with max AS
                    @max_AS_refseg = grep $_ eq $keepRefsegID, @max_AS_refseg;
                    next if @max_AS_refseg == 0;
                }
                else{ # the keepRefsegID must be the only refseg with max AS
                    next if @max_AS_refseg != 1 || $max_AS_refseg[0] ne $keepRefsegID;
                }
            }
            # go on: different scenarios
            if( $endsMap{1} xor $endsMap{2} ){ # _1/_2, one has map (may have unmap simultaneously), another are all unmap
                $outBam->write(content=>$_->printSAM."\n")
                    for grep {$_->mseg eq $max_AS_refseg[0] || $_->p_mseg eq $max_AS_refseg[0]}
                        map {@{$rOB_Af{$_}}} (1,2);
            }
            else{ # _1 and _2 both have map (may have unmap simultaneously)
                # firstly, find refseg that both _1 and _2 mapped
                my %end2mseg;
                $end2mseg{1}{$_->mseg} = 1 for @{$map_rOB{1}};
                $end2mseg{2}{$_->mseg} = 1 for @{$map_rOB{2}};
                my @msegBiEnd = sort {$refseg2AS{$b}<=>$refseg2AS{$a}} grep exists $end2mseg{2}{$_}, keys %{$end2mseg{1}};
                if(@msegBiEnd != 0){
                    $outBam->write(content=>$_->printSAM."\n")
                        for grep {$_->mseg eq $msegBiEnd[0] || $_->p_mseg eq $msegBiEnd[0]}
                            map {@{$rOB_Af{$_}}} (1,2);
                }
                else{ # each refseg has just one end map, swap ends?
                    if(defined $keepRefsegID || !$swapAlignEnds){ # not allowed when keepRefsegID is set, or disabled
                        $outBam->write(content=>$_->printSAM."\n")
                            for grep {$_->mseg eq $max_AS_refseg[0] || $_->p_mseg eq $max_AS_refseg[0]}
                                map {@{$rOB_Af{$_}}} (1,2);
                    }
                    else{ # do swap!
                        # find each end with max AS
                        my %maxAS_rOB = (1=>undef,2=>undef);
                        for my $end (1,2){
                            for my $rOB (@{$map_rOB{$end}}){
                                $maxAS_rOB{$end} = $rOB if !defined $maxAS_rOB{$end} || $maxAS_rOB{$end}->alignScore < $rOB->alignScore;
                            }
                        }
                        # swap: update rOB attributes, and output
                        for my $tEnd (1,2){ # [t]his end
                            my $mEnd = $tEnd % 2 + 1; # [m]ate end
                            my $t_rOB = $maxAS_rOB{$tEnd};
                            my $m_rOB = $maxAS_rOB{$mEnd};
                            # update
                            ## mated end mapped refseg and position
                            $t_rOB->update_attr(attr_id=>'p_mseg', value=>$m_rOB->mseg);
                            $t_rOB->update_attr(attr_id=>'p_mpos', value=>$m_rOB->mpos);
                            ## set Tlen as zero
                            $t_rOB->update_attr(attr_id=>'tlen', value=>0);
                            ## flags
                            $t_rOB->  set_flag(flag=>0x1 );
                            $t_rOB->  set_flag(flag=>0x20) if $m_rOB->is_rv_map;
                            $t_rOB->unset_flag(flag=>0x20) if $m_rOB->is_fw_map;
                            $t_rOB->unset_flag(flag_Af=>[0x2,0x4,0x8,0x200,0x400]);
                            ## mate-end MQ
                            $t_rOB->update_optfd(tag=>'MQ:i:', value=>$m_rOB->mapQ);
                            # output
                            $outBam->write(content=>$t_rOB->printSAM."\n");
                        }
                    }
                }
            }
        }
    }
}

1; ## tell the perl script the successful access of this module.
