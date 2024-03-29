package BioFuse::BioInfo::Objects::SeqData::Reads_OB;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ max /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::BioInfo::Objects::Allele::AlleleOnReads_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::SeqData::Reads_OB';
#----- version --------
$VERSION = "0.23";
$DATE = '2022-03-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        pid
                        endNO
                        mseg
                        p_mseg
                        mpos
                        p_mpos
                        mapQ
                        rg_id
                        rg_OB
                        rlen
                        mReadLen
                        mRefLen
                        cigar
                        rseq
                        lenFromCigar
                        barcode
                        add_BarcToPid
                        optfd_str
                        add_str_to_optfd
                        optfd_has_regex
                        update_optfd
                        mmCount
                        editDist
                        alignScore
                        judgeAlign
                        is_proper_map
                        is_fw_map
                        is_p_fw_map
                        is_rv_map
                        is_p_rv_map
                        is_unmap
                        is_p_unmap
                        is_2ndmap
                        free_2ndmap
                        is_suppmap
                        free_suppmap
                        is_dup
                        is_mltmap
                        is_good_cigar
                        is_softclip
                        is_hardclip
                        is_clip
                        has_MDtag
                        is_closeAlign
                        foreClipLen
                        hindClipLen
                        biClipLen
                        clipMatch
                        maxClipLen
                        tlen_FixTlen_with_S
                        digestMDtag
                        mDetailAf
                        fuseCigarMD
                        getNearAltDist
                        get_pos_allele
                        printSAM
                        printFQ
                        printFA
                        update_attr
                        set_flag
                        unset_flag
                        find_rgOB
                        add_RGrIDprefToPid
                        update_rgOB_maxRlen
                        update_rid_RGprefix
                        get_pos_allele_v1_BaseOnCigar
                     /;

#--- structure of object
# reads_OB -> SAMtext = $SAMtext # for simple load
# reads_OB -> pid = $pid
# reads_OB -> flag = $flag
# reads_OB -> mseg = $mseg
# reads_OB -> mpos = $mpos
# reads_OB -> mapQ = $mapQ
# reads_OB -> cigar = $cigar
# reads_OB -> MDtag = $MDtag, might lack
# reads_OB -> p_mseg = $p_mseg
# reads_OB -> p_mpos = $p_mpos
# reads_OB -> tlen = $tlen
# reads_OB -> rlen = $rlen
# reads_OB -> mReadLen = $mReadLen
# reads_OB -> mRefLen = $mRefLen
# reads_OB -> rseq = $rseq
# reads_OB -> baseQ = $baseQ
# reads_OB -> optfd = $optfd
# reads_OB -> rg_id = $rg_id
# reads_OB -> rg_OB = $rg_OB
# reads_OB -> endNO = $endNO
# reads_OB -> AS = $AS, alignment score
# reads_OB -> mDetailAf = [ {type => 'M/R/V/D/I/S/H/N', len => $len, str=>$str}, ... ] 
#    note: 'str' only for 'V' and 'D'(undef if no MD_tag)

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;
    my $SAMtext = $parm{ReadsLineText};
    my $_rc_optfd  = $parm{_rc_optfd} || 0;
    my $_rc_rdstr  = $parm{_rc_rdstr} || 0;
    my $simpleLoad = $parm{simpleLoad} || 0; # only record basic-info, will disable '_rc_optfd' and '_rc_rdstr'
    my $reads_end  = $parm{reads_end}; # for R1 and R2 aligned respectively, e.g., Hi-C PE-reads

    my $reads_OB = {};
    # basic infos
    chomp($SAMtext);
    my @term = split /\t+/, $SAMtext;

    # basic alignment info
    $reads_OB->{pid}   = $term[0];
    $reads_OB->{flag}  = $term[1];
    $reads_OB->{mseg}  = $term[2];
    $reads_OB->{mpos}  = $term[3];
    ## reads end NO.
    if( defined $reads_end ){
        $reads_OB->{endNO} = $reads_end;
        $reads_OB->{flag} |= ( $reads_end == 1 ? 0x40 : 0x80 ); # update flag
    }
    else{
        $reads_OB->{endNO} = ($reads_OB->{flag} & 0x40) ? 1 : 2; # 0x80 is '2'
    }

    # simple load
    if($simpleLoad){
        $reads_OB->{SAMtext}  = $SAMtext;
        bless($reads_OB);
        return $reads_OB;
    }

    # normal load
    $reads_OB->{mapQ}  = $term[4];
    $reads_OB->{cigar} = $term[5];
    ($reads_OB->{MDtag}) = ($SAMtext =~ /\sMD:Z:(\S+)/); # sometimes, SAM lacks MD-tag
    ## paired-end info
    $reads_OB->{p_mseg} = $term[6];
    $reads_OB->{p_mpos} = $term[7];
    $reads_OB->{tlen}   = $term[8];
    ## read string
    $reads_OB->{rlen}   = length($term[9]); # note that Hard-clip is short!
    if( $_rc_rdstr ){
        $reads_OB->{rseq}  = $term[9];
        $reads_OB->{baseQ} = $term[10];
    }
    ## optional fields
    if( $_rc_optfd ){
        $reads_OB->{optfd} = join("\t", @term[11 .. $#term]);
    }
    ### reads group
    $reads_OB->{rg_id} = ($SAMtext =~ /RG\:Z\:(\S+)/) ? $1 : undef;
    ### make up cigar 'H' for GATK realn
    if(    $reads_OB->{cigar} !~ /H/
        && $SAMtext =~ /OC\:Z\:(\S+H\S*)/
    ){
        my $origCigar = $1;
        if( $origCigar =~ /^(\d+H)/ ){ $reads_OB->{cigar} = $1.$reads_OB->{cigar} }
        if( $origCigar =~ /(\d+H)$/ ){ $reads_OB->{cigar} = $reads_OB->{cigar}.$1 }
    }

    bless($reads_OB);
    return $reads_OB;
}

#--- return pid ---
sub pid{
    my $reads_OB = shift;
    return $reads_OB->{pid};
}

#--- return endNO ---
sub endNO{
    my $reads_OB = shift;
    return $reads_OB->{endNO};
}

#--- return mapped segment name ---
sub mseg{
    my $reads_OB = shift;
    return $reads_OB->{mseg};
}

#--- return mapped segment name of the paired end ---
sub p_mseg{
    my $reads_OB = shift;
    return $reads_OB->{p_mseg};
}

#--- return mapped position ---
sub mpos{
    my $reads_OB = shift;
    return $reads_OB->{mpos};
}

#--- return mapped position of the paired end ---
sub p_mpos{
    my $reads_OB = shift;
    return $reads_OB->{p_mpos};
}

#--- return mapping quality ---
sub mapQ{
    my $reads_OB = shift;
    return $reads_OB->{mapQ};
}

#--- return reads group ID ---
sub rg_id{
    my $reads_OB = shift;
    if(!defined $reads_OB->{rg_id}){
        cluck_and_exit "<ERROR>\tCan not find rg_id of reads $reads_OB->{pid}/$reads_OB->{endNO}\n";
    }
    else{
        return $reads_OB->{rg_id};
    }
}

#--- return reads group OB ---
sub rg_OB{
    my $reads_OB = shift;
    if(!defined $reads_OB->{rg_OB}){
        cluck_and_exit "<ERROR>\tCan not find rg_OB of reads $reads_OB->{pid}/$reads_OB->{endNO}\n".
                              "\tPlease use 'find_rgOB' func of reads OB.\n";
    }
    else{
        return $reads_OB->{rg_OB};
    }
}

#--- return read length ---
## based on col.10 in SAM
sub rlen{
    my $reads_OB = shift;
    return $reads_OB->{rlen};
}

#--- return length of read's mapped part ---
## mReadLen: just remove length of bilateral soft-clip part(s); hard-clip is rlen
sub mReadLen{
    my $reads_OB = shift;
    unless( exists $reads_OB->{mReadLen} ){
        $reads_OB->{mReadLen} = $reads_OB->rlen - $reads_OB->biClipLen;
    }
    return $reads_OB->{mReadLen};
}

#--- return the length of mapped part on reference genome ---
## mRefLen: consider InDel and skip
sub mRefLen{
    my $reads_OB = shift;
    unless( exists $reads_OB->{mRefLen} ){
        $reads_OB->{mRefLen} =   $reads_OB->mReadLen
                               + $reads_OB->lenFromCigar(type=>'D')
                               + $reads_OB->lenFromCigar(type=>'N')
                               - $reads_OB->lenFromCigar(type=>'I');
    }
    return $reads_OB->{mRefLen};
}

#--- return cigar ---
sub cigar{
    my $reads_OB = shift;
    return $reads_OB->{cigar};
}

#--- return read seq ---
sub rseq{
    my $reads_OB = shift;
    return $reads_OB->{rseq};
}

#--- get length of given type from Cigar ---
## M,D,I,H,S,N
sub lenFromCigar{
    my $reads_OB = shift;
    my %parm = @_;
    my $type = $parm{type};

    my $len = 0;
    my $cigar = $reads_OB->{cigar};
    while( $cigar =~ /\d+$type/ ){ # still has
        if( $cigar =~ /^(\d+)(\D)/ ){ # trim from left, might not the required type
            $len += $1 if( $2 eq $type );
        }
        $cigar =~ s/^\d+\D//;
    }

    return $len;
}

#--- return barcode of this reads ---
sub barcode{
    my $reads_OB = shift;
    my %parm = @_;
    my $barcSAMtag = $parm{barcSAMtag} || 'BX:Z:'; # default is 10x-barcode
    my $barcPidPrf = $parm{barcPidPrf} || ':__BARC__:';

    if(    defined $reads_OB->{optfd}
        && $reads_OB->{optfd} =~ /$barcSAMtag(\S+)/
    ){
        return $1;
    }
    elsif( $reads_OB->{pid} =~ /$barcPidPrf([^:]+)/ ){
        return $1;
    }
    else{
        cluck_and_exit "<ERROR>\tcannot confirm 10x-barcode of reads\n".Dumper($reads_OB);
    }
}

#--- update pid with barcode postfix
sub add_BarcToPid{
    my $reads_OB = shift;
    my %parm = @_;
    my $barcSAMtag = $parm{barcSAMtag} || 'BX:Z:'; # default is 10x-barcode
    my $barcPidPrf = $parm{barcPidPrf} || ':__BARC__:';

    if($reads_OB->{pid} !~ /$barcPidPrf/){
        $reads_OB->{pid} = $reads_OB->{pid} . $barcPidPrf . $reads_OB->barcode(barcSAMtag=>$barcSAMtag, barcPidPrf=>$barcPidPrf);
    }
    else{
        cluck_and_exit "<ERROR>\tthis reads likely already has barcode in pid\n".Dumper($reads_OB);
    }
}

#--- return optfd_str ---
sub optfd_str{
    my $reads_OB = shift;

    if( defined $reads_OB->{optfd} ){
        return $reads_OB->{optfd};
    }
    else{
        cluck_and_exit "<ERROR>\tCannot find optfd_string of reads $reads_OB->{pid} end-$reads_OB->{endNO}.\n";
    }
}

#--- add string to optfd_str ---
sub add_str_to_optfd{
    my $reads_OB = shift;
    my %parm = @_;
    $reads_OB->{optfd} .= $parm{str};
}

#--- test whether optfd_str match given regex ---
## generally, it is string
sub optfd_has_regex{
    my $reads_OB = shift;
    my %parm = @_;
    return ($reads_OB->optfd_str =~ /$parm{regex}/);
}

#--- update optfd ---
## !! note: Z, Printable string, including space ([\S ]+)
sub update_optfd{
    my $reads_OB = shift;
    my %parm = @_;
    my $tag = $parm{tag};
    my $value = $parm{value};
    if($reads_OB->{optfd} =~ /\b($tag[^\t]*)/){ # /\b($tag[\S ]+)\b/
        my $aim = quotemeta $1;
        $reads_OB->{optfd} =~ s/\b$aim/$tag$value/;
    }
    else{
        $reads_OB->{optfd} .= "\t$tag$value";
    }
}

#--- return mismatch count ---
sub mmCount{
    my $reads_OB = shift;
    if($reads_OB->optfd_str =~ /XM:i:(\d+)/){
        return $1;
    }
    elsif($reads_OB->optfd_str =~ /MD:Z:(\w+)/){
        (my $md = $1) =~ s/\d//g;
        return length($md);
    }
    else{
        cluck_and_exit "<ERROR>\tCannot find mismatch count info of reads $reads_OB->{pid} end-$reads_OB->{endNO}.\n";
    }
}

#--- return edit distance ---
sub editDist{
    my $reads_OB = shift;
    if($reads_OB->optfd_str =~ /NM:i:(\d+)/){
        return $1;
    }
    else{
        cluck_and_exit "<ERROR>\tCannot find edit distance info of reads $reads_OB->{pid} end-$reads_OB->{endNO}.\n";
    }
}

#--- return alignment score ---
sub alignScore{
    my $reads_OB = shift;
    if(defined $reads_OB->{AS}){
        return $reads_OB->{AS};
    }
    if($reads_OB->optfd_str =~ /AS:i:(\d+)/){
        $reads_OB->{AS} = $1;
        return $1;
    }
    else{
        return 0;
    }
}

#--- do judgement on alignment ---
## SP: Supplementary alignment
## SD: Secondary alignment
## MM: Multiple-Mapped
## LM: Low Mapping quality
sub judgeAlign{
    my $reads_OB = shift;
    my %parm = @_;
    my $type = $parm{type};
    my $min_mapQ = $parm{min_mapQ};

    if($type eq 'SP'){
        return $reads_OB->is_suppmap;
    }
    elsif($type eq 'SD'){
        return $reads_OB->is_2ndmap;
    }
    elsif($type eq 'MM'){
        return $reads_OB->is_mltmap;
    }
    elsif($type eq 'LM'){
        return ($reads_OB->mapQ < $min_mapQ);
    }
}

#--- test whether it is proper-mapped ---
sub is_proper_map{
    my $reads_OB = shift;
    return $reads_OB->{flag} & 0x2;
}

#--- test whether it is forward-mapped ---
sub is_fw_map{
    my $reads_OB = shift;
    return !($reads_OB->{flag} & 0x4) && !($reads_OB->{flag} & 0x10);
}

#--- test whether its paired end is forward-mapped ---
sub is_p_fw_map{
    my $reads_OB = shift;
    return !($reads_OB->{flag} & 0x8) && !($reads_OB->{flag} & 0x20);
}

#--- test whether it is reversed-mapped ---
sub is_rv_map{
    my $reads_OB = shift;
    return !($reads_OB->{flag} & 0x4) && ($reads_OB->{flag} & 0x10);
}

#--- test whether its paired end is reversed-mapped ---
sub is_p_rv_map{
    my $reads_OB = shift;
    return !($reads_OB->{flag} & 0x8) && ($reads_OB->{flag} & 0x20);
}

#--- test whether it is un-mapped ---
sub is_unmap{
    my $reads_OB = shift;
    return $reads_OB->{flag} & 0x4;
}

#--- test whether its paired end is un-mapped ---
sub is_p_unmap{
    my $reads_OB = shift;
    return $reads_OB->{flag} & 0x8;
}

#--- test whether it is 0x100, secondary alignment ---
sub is_2ndmap{
    my $reads_OB = shift;
    return $reads_OB->{flag} & 0x100;
}

#--- free it from 0x100 if it is secondary alignment ---
sub free_2ndmap{
    my $reads_OB = shift;
    $reads_OB->{flag} -= 0x100 if $reads_OB->is_2ndmap;
}

#--- test whether it is 0x800, supplementary alignemnt ---
sub is_suppmap{
    my $reads_OB = shift;
    return $reads_OB->{flag} & 0x800;
}

#--- free it from 0x800 if it is supplementary alignemnt ---
sub free_suppmap{
    my $reads_OB = shift;
    $reads_OB->{flag} -= 0x800 if $reads_OB->is_suppmap;
}

#--- test whether it is 0x400, duplicated alignemnt ---
sub is_dup{
    my $reads_OB = shift;
    return $reads_OB->{flag} & 0x400;
}

#--- test whether it is multiple-mapped (mapQ is 0) ---
sub is_mltmap{
    my $reads_OB = shift;
    return ( $reads_OB->{mapQ} == 0 );
}

#--- test whether it has good cigar ---
# good: lacks 'N', 'H' tag
# 'N' may affect RNAseq-alignment
sub is_good_cigar{
    my $reads_OB = shift;
    return $reads_OB->{cigar} !~ /[NH]/;
}

#--- test whether it is soft-clip ---
sub is_softclip{
    my $reads_OB = shift;
    return $reads_OB->{cigar} =~ /S/;
}

#--- test whether it is hard-clip ---
sub is_hardclip{
    my $reads_OB = shift;
    return $reads_OB->{cigar} =~ /H/;
}

#--- test whether it is soft/hard-clip ---
sub is_clip{
    my $reads_OB = shift;
    return ( $reads_OB->is_softclip || $reads_OB->is_hardclip );
}

#--- test whether have MD-tag ---
sub has_MDtag{
    my $reads_OB = shift;
    return defined $reads_OB->{MDtag};
}

#--- test whether two read_OB are close alingment ---
sub is_closeAlign{
    my $reads_OB = shift;
    my %parm = @_;
    my $test_rOB = $parm{test_rOB};
    my $distance = $parm{distance} || 1000;

    if(    $reads_OB->is_unmap
        || $test_rOB->is_unmap
    ){
        return 0;
    }

    if(        $reads_OB->mseg eq $test_rOB->mseg
        && abs($reads_OB->mpos -  $test_rOB->mpos) < $distance
    ){
        return 1;
    }
    else{
        return 0;
    }
}

#--- return fore-clip length of given type ---
## default 'S'
## memory debug: 1) avoid variable in regex is OK
## 2018-11-13    2) avoid variable in '[]' in regex is BAD
sub foreClipLen{
    my $reads_OB = shift;
    my %parm = @_;
    my $clipType = $parm{clipType} || 'S';
    # my $clipRegx = "[$clipType]";
    # if($reads_OB->{cigar} =~ /^(\d+)$clipRegx/i){
    #     return $1;
    # }
    # else{
    #     return 0;
    # }
    if(   $clipType =~ /^s$/i
       && $reads_OB->{cigar} =~ /^(\d+)S/i
    ){
        return $1;
    }
    if(   $clipType =~ /^h$/i
       && $reads_OB->{cigar} =~ /^(\d+)H/i
    ){
        return $1;
    }
    if(   $clipType =~ /s/i
       && $clipType =~ /h/i
       && $reads_OB->{cigar} =~ /^(\d+)[SH]/i
    ){
        return $1;
    }
    # lack the required clip type
    return 0;
}

#--- return hind-clip length of given type ---
## default 'S'
sub hindClipLen{
    my $reads_OB = shift;
    my %parm = @_;
    my $clipType = $parm{clipType} || 'S';
    # my $clipRegx = "[$clipType]";
    # if($reads_OB->{cigar} =~ /(\d+)$clipRegx$/i){
    #     return $1;
    # }
    # else{
    #     return 0;
    # }
    if(   $clipType =~ /^s$/i
       && $reads_OB->{cigar} =~ /(\d+)S$/i
    ){
        return $1;
    }
    if(   $clipType =~ /^h$/i
       && $reads_OB->{cigar} =~ /(\d+)H$/i
    ){
        return $1;
    }
    if(   $clipType =~ /s/i 
       && $clipType =~ /h/i
       && $reads_OB->{cigar} =~ /(\d+)[SH]$/i
    ){
        return $1;
    }
    # lack the required clip type
    return 0;
}

#--- get bi-clip length of given type ---
## default 'S'
sub biClipLen{
    my $reads_OB = shift;
    my %parm = @_;
    my $clipType = $parm{clipType} || 'S';
    return (  $reads_OB->foreClipLen( clipType => $clipType )
            + $reads_OB->hindClipLen( clipType => $clipType ) );
}

#--- check whether other_rOB has clipped part matched to this reads ---
## theoretically， this should be the same reads aligned to different positions/refsegs
##                 e.g., split-reads of rearrangement
## compare rlen with sum of selected unilateral clipped parts of two rOB
## consider inner-insertion (1+r) and micro-homology (1-r)
sub clipMatch{
    my $reads_OB = shift;
    my %parm = @_;
    my $other_rOB = $parm{other_rOB};
    my $CLsumRLthreR = $parm{CLsumRLthreR} || 0.2; # [1-r,1+r]

    # both must be clipped aligned
    return 0 unless $reads_OB->is_clip && $other_rOB->is_clip;
    # this rOB ClipLen & MapLen*
    my $t_foreCL = $reads_OB->foreClipLen;
    my $t_hindCL = $reads_OB->hindClipLen;
    my $t_mapLen = $reads_OB->mReadLen;
    my $t_foreML = $t_mapLen - $t_foreCL;
    my $t_hindML = $t_mapLen - $t_hindCL;
    # other rOB ClipLen & MapLen*
    my $o_foreCL = $other_rOB->foreClipLen;
    my $o_hindCL = $other_rOB->hindClipLen;
    my $o_mapLen = $other_rOB->mReadLen;
    my $o_foreML = $o_mapLen - $o_foreCL;
    my $o_hindML = $o_mapLen - $o_hindCL;
    # whether clipped parts match with each other
    my $CLsumMaxDiff = $reads_OB->rlen * $CLsumRLthreR * 2;
    if(   ($reads_OB->is_fw_map && $other_rOB->is_fw_map)
       || ($reads_OB->is_rv_map && $other_rOB->is_rv_map)
    ){
        # return 1 if   abs($reads_OB->foreClipLen + $other_rOB->hindClipLen - $reads_OB->rlen) <= $CLsumMaxDiff
        #            || abs($reads_OB->hindClipLen + $other_rOB->foreClipLen - $reads_OB->rlen) <= $CLsumMaxDiff;
        return 1 if   abs($t_foreML-$o_foreCL) + abs($t_hindCL-$o_hindML) <= $CLsumMaxDiff
                   || abs($t_foreCL-$o_foreML) + abs($t_hindML-$o_hindCL) <= $CLsumMaxDiff;
    }
    else{
        # return 1 if   abs($reads_OB->foreClipLen + $other_rOB->foreClipLen - $reads_OB->rlen) <= $CLsumMaxDiff
        #            || abs($reads_OB->hindClipLen + $other_rOB->hindClipLen - $reads_OB->rlen) <= $CLsumMaxDiff;
        return 1 if   abs($t_foreML-$o_hindCL) + abs($t_hindCL-$o_foreML) <= $CLsumMaxDiff
                   || abs($t_foreCL-$o_hindML) + abs($t_hindML-$o_foreCL) <= $CLsumMaxDiff;
    }
    # last, not match
    return 0;
}

#--- return max clip length of given type ---
sub maxClipLen{
    my $reads_OB = shift;
    my %parm = @_;
    my $clipType = $parm{clipType} || 'S';
    return max( $reads_OB->foreClipLen( clipType => $clipType ),
                $reads_OB->hindClipLen( clipType => $clipType ) );
}

#--- get the length to fix Tlen based on soft-clip Edge part ---
sub tlen_FixTlen_with_S{
    my $reads_OB = shift;

    if( $reads_OB->is_fw_map ){
        return $reads_OB->foreClipLen;
    }
    elsif( $reads_OB->is_rv_map ){
        return $reads_OB->hindClipLen;
    }
}

#--- return digested MD tag ---
sub digestMDtag{

    my $reads_OB = shift;

    return [] unless( $reads_OB->has_MDtag );

    my $MDtag = $reads_OB->{MDtag}; # might be empty
    my @MDtag;
    while( length($MDtag) != 0 ){
        if( $MDtag =~ /^(\d+)/ ){ # match
            push @MDtag, ['ref',$1];
        }
        elsif( $MDtag =~ /^([^\^\d])/ ){ # snv, take single letter each time
            push @MDtag, ['snv',$1];
        }
        elsif( $MDtag =~ /^\^(\D+)/ ){ # del
            push @MDtag, ['del',$1];
        }
        # update
        $MDtag =~ s/^\^?$MDtag[-1]->[1]//;
        # discard solo '0'
        pop @MDtag if $MDtag[-1]->[1] eq '0';
    }

    return \@MDtag;
}

#--- return the mDetailAf (complete map info) obtained from fuseCigarMD ---
sub mDetailAf{
    my $reads_OB = shift;
    $reads_OB->fuseCigarMD unless exists $reads_OB->{mDetailAf};
    return $reads_OB->{mDetailAf};
}

#--- fuse Cigar and MD tag to get complete map info ---
# R: Ref, matched
# V: snV, mismatched
# I: Ins, insertion
# D: Del, deletion
# S: Soft-clipped
# H: Hard-clipped
# N: skipped ref
sub fuseCigarMD{
    my $reads_OB = shift;
    my $CIGAR = $reads_OB->{cigar};
    my $hasMD = $reads_OB->has_MDtag;
    my $MDAf  = $reads_OB->digestMDtag;

    my @mDetail;
    while( length($CIGAR) != 0 ){
        my ($len, $Ctag) = ($CIGAR =~ /^(\d+)(\D)/);
        # check
        unless( defined $len && defined $Ctag ){
            cluck_and_exit "<ERROR>\tcannot recognize cigar of this reads.\n".Dumper($reads_OB);
        }
        # shorten cigar
        $CIGAR =~ s/^\d+\D//;
        # different scenarios
        if( $Ctag eq 'M' ){ # Mapped with possible mismatch ONLY
            # sometimes, MD-tag is empty
            if( !$hasMD ){
                push @mDetail, { type => 'M', len => $len }; # 'refORsnv'
            }
            else{
                while($len > 0){
                    my $MDtype = $MDAf->[0]->[0];
                    # MDtype should be 'ref' or 'snv'
                    if( $MDtype eq 'ref' ){ # 'R'
                        my $MDlen = $MDAf->[0]->[1];
                        if( $MDlen > $len ){
                            # record
                            push @mDetail, { type => 'R', len => $len };
                            # update
                            $MDAf->[0]->[1] -= $len;
                            $len = 0;
                        }
                        else{
                            # record
                            push @mDetail, { type => 'R', len => $MDlen };
                            # update
                            shift @$MDAf;
                            $len -= $MDlen;
                        }
                    }
                    elsif( $MDtype eq 'snv' ){ # 'V', length must be 1
                        # record
                        push @mDetail, { type => 'V', len => 1, str => $MDAf->[0]->[1] };
                        # update
                        shift @$MDAf;
                        $len--;
                    }
                    else{
                        cluck_and_exit "<ERROR>\tbad 'M' part: Cigar($len$Ctag), MD tag (@{$MDAf->[0]}).\n".Dumper($reads_OB);
                    }
                }
            }
        }
        elsif( $Ctag eq 'D' ){ # Deletion
            # record
            push @mDetail, { type => 'D', len => $len, str => undef };
            # has MD-tag
            if( $hasMD ){
                # MD tag should be 'del'
                if(           $MDAf->[0]->[0]  ne 'del'
                    || length($MDAf->[0]->[1]) != $len
                ){
                    cluck_and_exit "<ERROR>\tbad 'D' part: Cigar($len$Ctag), MD tag (@{$MDAf->[0]}).\n".Dumper($reads_OB);
                }
                # update
                $mDetail[-1]->{str} = $MDAf->[0]->[1];
                shift @$MDAf;
            }
        }
        elsif(    $Ctag eq 'I' # Insertion
               || $Ctag eq 'S' # Soft-Clip
               || $Ctag eq 'H' # Hard-Clip
               || $Ctag eq 'N' # Skip-Ref, esp. for RNA-seq reads
        ){ # no action on MDtag
            push @mDetail, { type => $Ctag, len => $len };
        }
    }

    # exchange a[DN]bI to bIa[DN]
    for my $i ( 0 .. $#mDetail-1 ){
        my $j = $i + 1;
        if(    $mDetail[$i]->{type} =~ /[DN]/
            && $mDetail[$j]->{type} eq 'I'
        ){
            @mDetail[$i,$j] = @mDetail[$j,$i];
        }
    }

    # connect aIbI or aDbD
    for my $i ( reverse 1 .. $#mDetail ){
        my $j = $i - 1;
        if(    $mDetail[$i]->{type} =~ /[DI]/
            && $mDetail[$i]->{type} eq $mDetail[$j]->{type}
        ){
            # merge to element with lower index
            $mDetail[$j]->{len} += $mDetail[$i]->{len};
            $mDetail[$j]->{str} .= $mDetail[$i]->{str} if( $mDetail[$i]->{type} eq 'D' );
            # discard element with higher index
            splice( @mDetail, $i, 1 );
        }
    }

    # record
    $reads_OB->{mDetailAf} = \@mDetail;
}

#--- get distance to the nearest alteration on the reads ---
# from the given index and prime
# if no alteration encountered, return reads-edge distance
sub getNearAltDist{
    my $reads_OB = shift;
    my %parm = @_;
    my $Bidx = $parm{Bidx};
    my $prime = $parm{prime};

    my $mDetailAf = $reads_OB->{mDetailAf};
    my $mLastBidx = scalar(@$mDetailAf)-1;

    my @idx = ( $prime == 3 ? ($Bidx+1 .. $mLastBidx) : (reverse 0 .. $Bidx-1) );
    my $dist = $parm{distBase} || 0;
    for my $i (@idx){
        if( $mDetailAf->[$i]->{type} !~ /[MR]/ ){ # 'M' cannot be accurate
            return $dist;
        }
        else{
            $dist += $mDetailAf->[$i]->{len};
        }
    }
    return $dist;
}

#--- get allele at given position based on Cigar and MDtag ---
# return allele_OB
# snv: single base matched ref/alt base
# ins: seq inserted behind 'pos'
# del: offset means moving N bases from 1st-del-pos is the 'pos', so seq deleted behind 'pos' has offset -1.
# priority: 'preferINDEL' > snv > ins/del
sub get_pos_allele{
    my $reads_OB = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $refBase = $parm{refBase} || undef; # only use for 'refORsnv'
    my $preferINDEL = $parm{preferINDEL} || 0;
    my $no_warn = $parm{no_warn};

    # initialize allele_OB
    my $allele_OB = BioFuse::BioInfo::Objects::Allele::AlleleOnReads_OB->new( chr => $chr, pos => $pos, rOBmLen => $reads_OB->mReadLen );
       $allele_OB->loadInfo( refBase => $refBase ) if( defined $refBase );

    # unmap check
    if( $reads_OB->is_unmap ){
        unless( $no_warn ){
            warn "<WARN>\treads $reads_OB->{pid} end-$reads_OB->{endNO} is unmapped.\n".Dumper($reads_OB);
        }
        $allele_OB->setMiss;
        return $allele_OB;
    }
    # chr check
    if( $reads_OB->{mseg} ne $chr ){
        unless( $no_warn ){
            warn "<WARN>\treads $reads_OB->{pid} end-$reads_OB->{endNO} doesn't match query chr ($chr).\n".Dumper($reads_OB);
        }
        $allele_OB->setMiss;
        return $allele_OB;
    }

    # attributes
    my $mpos = $reads_OB->{mpos};
    my $rSeq = uc($reads_OB->{rseq});
    my $mlen = $reads_OB->mReadLen;
    my $qSeq = $reads_OB->{baseQ};
    my $CIGAR = $reads_OB->{cigar};
    # the complete map info
    ## fusing Cigar and MD-tag info, the later might lack
    my $mDetailAf = $reads_OB->mDetailAf;
    my $mLastBidx = scalar(@$mDetailAf)-1;

    # find the pos
    my $dist_shift = 0;
    my $rSpn_refLen = 0;
    my $rCov_lastPos = -1; # M/D/N last covered pos
    my $rCov_lastBase; # M/D/N last covered base
    my $rCov_lastBaseQ; # M/D/N last covered base quality
    for my $Bidx ( 0 .. $mLastBidx ){
        my $type = $mDetailAf->[$Bidx]->{type};
        my $len  = $mDetailAf->[$Bidx]->{len};
        # different scenarios
        if(    $type eq 'R' # ref matched
            || $type eq 'V' # snv
            || $type eq 'M' # ref OR snv
        ){
            $rCov_lastPos = $mpos + $rSpn_refLen + $len - 1;
            # M-part cover the want-pos
            if(    $rCov_lastPos > $pos
                || (   $rCov_lastPos == $pos
                       # ignore INDEL next to, or next is not INDEL
                    && $Bidx != $mLastBidx
                    && ( !$preferINDEL || $mDetailAf->[$Bidx+1]->{type} !~ /[ID]/ )
                   )
            ){
                my $sub_i = $len + $pos - $rCov_lastPos - 1;
                my $allele = substr($rSeq, $sub_i, 1);
                my $rEdgeDist5 = $rSpn_refLen + $sub_i + 1 + $dist_shift; # this pos included
                my $rEdgeDist3 = $mlen - $rEdgeDist5 + 1; # this pos included
                my $alleleType = ( $type eq 'R' ? 'ref' : ($type eq 'V' ? 'snv' : 'refORsnv') );
                   $alleleType = ( uc($allele) eq uc($refBase) ? 'ref' : 'snv' ) if( $type eq 'M' && defined $refBase );
                $allele_OB->loadInfo(
                    type => $alleleType,
                    allele => $allele,
                    baseQ => substr($qSeq, $sub_i, 1),
                    rEdgeDist5 => $rEdgeDist5,
                    rEdgeDist3 => $rEdgeDist3,
                    nAltDist5 => $reads_OB->getNearAltDist( Bidx => $Bidx, prime => 5, distBase => $sub_i ),
                    nAltDist3 => $reads_OB->getNearAltDist( Bidx => $Bidx, prime => 3, distBase => $len-$sub_i-1 )
                );
                return $allele_OB;
            }
            # update for next block
            $rCov_lastBase  = substr($rSeq, $len-1, 1);
            $rCov_lastBaseQ = substr($qSeq, $len-1, 1);
            $rSpn_refLen += $len;
            $rSeq = substr($rSeq, $len);
            $qSeq = substr($qSeq, $len);
        }
        elsif($type eq 'I'){
            # inserts just behind the want-pos
            if( $rCov_lastPos == $pos ){
                my $rEdgeDist5 = $rSpn_refLen + $dist_shift; # forebase
                my $rEdgeDist3 = $mlen - $rEdgeDist5 - $len; # nextbase
                $allele_OB->loadInfo(
                    type => 'ins',
                    allele => substr($rSeq, 0, $len),
                    baseQ => substr($qSeq, 0, $len),
                    rEdgeDist5 => $rEdgeDist5,
                    rEdgeDist3 => $rEdgeDist3,
                    nAltDist5 => $reads_OB->getNearAltDist( Bidx => $Bidx, prime => 5 ),
                    nAltDist3 => $reads_OB->getNearAltDist( Bidx => $Bidx, prime => 3 ),
                    forebase  => $rCov_lastBase,
                    forebaseQ => $rCov_lastBaseQ,
                    nextbase  => substr($rSeq, $len, 1),
                    nextbaseQ => substr($qSeq, $len, 1)
                );
                return $allele_OB;
            }
            # update for next block
            $rSeq = substr($rSeq, $len);
            $qSeq = substr($qSeq, $len);
            $dist_shift += $len;
        }
        elsif(   $type eq 'D'
              || $type eq 'N' # for RNA-seq reads
        ){
            my $del1stPos = $mpos + $rSpn_refLen;
            # deletion covers the want-pos
            if(    $type eq 'D'
                && $pos <= $del1stPos + $len - 1
            ){
                my $rEdgeDist5 = $rSpn_refLen + $dist_shift; # forebase
                my $rEdgeDist3 = $mlen - $rEdgeDist5; # nextbase
                $allele_OB->loadInfo(
                    type => 'del',
                    allele => $mDetailAf->[$Bidx]->{str}, # might be undef if lacks MD-tag
                    baseQ => undef,
                    delSize => $len,
                    offset => $pos - $del1stPos, # -1 is 'forebase'
                    rEdgeDist5 => $rEdgeDist5,
                    rEdgeDist3 => $rEdgeDist3,
                    nAltDist5 => $reads_OB->getNearAltDist( Bidx => $Bidx, prime => 5 ),
                    nAltDist3 => $reads_OB->getNearAltDist( Bidx => $Bidx, prime => 3 ),
                    forebase  => $rCov_lastBase,
                    forebaseQ => $rCov_lastBaseQ,
                    nextbase  => substr($rSeq, 0, 1),
                    nextbaseQ => substr($qSeq, 0, 1)
                );
                return $allele_OB;
            }
            # update for next block
            $rCov_lastPos = $mpos + $rSpn_refLen + $len - 1;
            $rCov_lastBase = undef; # fine, as 'D' must be followed by 'R/V', then get updated
            $rCov_lastBaseQ = undef;
            $rSpn_refLen += $len;
            $dist_shift -= $len;
        }
        elsif($type eq 'S'){
            # just skip S part, update for next block
            $rSeq = substr($rSeq, $len);
            $qSeq = substr($qSeq, $len);
        }
        elsif($type eq 'H'){
            1; # nothing to do
        }
    }

    # reach here, means this read doesn't cover the pos.
    unless( $no_warn ){
        warn "<WARN>\treads $reads_OB->{pid} end-$reads_OB->{endNO} doesn't cover query pos ($chr:$pos).\n".Dumper($reads_OB);
    }
    $allele_OB->setMiss;
    return $allele_OB;
}

#--- print reads in SAM format ---
sub printSAM{
    my $reads_OB = shift;
    # simple load
    return $reads_OB->{SAMtext} if defined $reads_OB->{SAMtext};
    # normal load
    my @SAM_fields = qw/ pid flag mseg mpos mapQ cigar p_mseg p_mpos tlen rseq baseQ optfd /;
    return join("\t", map {($reads_OB->{$_})} @SAM_fields );
}

#--- print reads in FastQ format ---
sub printFQ{
    my $reads_OB = shift;
    my $rSeq = $reads_OB->{rseq}; # initialize
    my $bQul = $reads_OB->{baseQ};
    if($reads_OB->is_rv_map){
        ($rSeq = reverse uc($rSeq)) =~ tr/ACGT/TGCA/;
        $bQul  = reverse $bQul;
    }
    return join("\n", "\@$reads_OB->{pid}/$reads_OB->{endNO}", $rSeq, '+', $bQul);
}

#--- print reads in FastA format ---
sub printFA{
    my $reads_OB = shift;
    my $rSeq = $reads_OB->{rseq}; # initialize
    if($reads_OB->is_rv_map){
        ($rSeq = reverse uc($rSeq)) =~ tr/ACGT/TGCA/;
    }
    return join("\n", "\>$reads_OB->{pid}/$reads_OB->{endNO}", $rSeq);
}

#---- update attribute ---
sub update_attr{
    my $reads_OB = shift;
    my %parm = @_;
    $reads_OB->{ $parm{attr_id} } = $parm{value};
}

#---- set flag ---
sub set_flag{
    my $reads_OB = shift;
    my %parm = @_;
    $parm{flag_Af} = [$parm{flag}] if !defined $parm{flag_Af};
    $reads_OB->{flag} |= $_ for @{$parm{flag_Af}};
}

#---- unset flag ---
sub unset_flag{
    my $reads_OB = shift;
    my %parm = @_;
    $parm{flag_Af} = [$parm{flag}] if !defined $parm{flag_Af};
    for (@{$parm{flag_Af}}){
        $reads_OB->{flag} -= $_ if $reads_OB->{flag} & $_;
    }
}

#--- find rg_OB based on rg_id ---
sub find_rgOB{
    my $reads_OB = shift;
    my %parm = @_;
    my $rgid2rgOB_Href = $parm{rgid2rgOB_Href};

    if(    !defined $reads_OB->{rg_id}
        || !exists $rgid2rgOB_Href->{$reads_OB->{rg_id}}
      ){
        cluck_and_exit "<ERROR>\tCan not get available RG_ID from reads $reads_OB->{pid}/$reads_OB->{endNO}\n";
    }
    else{
        $reads_OB->{rg_OB} = $rgid2rgOB_Href->{$reads_OB->{rg_id}};
        return $rgid2rgOB_Href->{$reads_OB->{rg_id}};
    }
}

#--- add rg rid prefix to pid ---
sub add_RGrIDprefToPid{
    my $reads_OB = shift;
    my $rg_OB = $reads_OB->{rg_OB};
    $rg_OB->{pid} = $rg_OB->rID_prefix . $rg_OB->{pid};
}

#--- update the maximum read length of its rg_OB ---
sub update_rgOB_maxRlen{
    my $reads_OB = shift;
    $reads_OB->{rg_OB}->update_maxRlen(endNO => $reads_OB->{endNO}, rLen => $reads_OB->{rlen});
}

#--- extract and record new read-id FS prefix ---
sub update_rid_RGprefix{
    my $reads_OB = shift;
    my %parm = @_;
    my $Tool_Tag = $parm{Tool_Tag};

    my $rg_OB = $reads_OB->{rg_OB};
    if( $reads_OB->{pid} =~ /^(.*_${Tool_Tag}rgNO(\d+)${Tool_Tag}_)(.+)$/ ){
        $rg_OB->{rID_prefix_prev} = $1;
        $rg_OB->{RG_NO} = $2;
        # add tissue tag
        my $curr_tissue_tag = ucfirst( $rg_OB->{tissue} );
        my $prev_tissue_tag = '';
        my $rg_info_prefix = '';
        if( $rg_OB->{rID_prefix_prev} =~ /_${Tool_Tag}ts(.+)${Tool_Tag}(.+)$/ ){
            $prev_tissue_tag = $1;
            $rg_info_prefix = $2;
        }
        else{ # no tissue tag added before
            $rg_info_prefix = $rg_OB->{rID_prefix_prev};
        }
        # construct current read-id prefix
        $rg_OB->{rID_prefix} = "_${Tool_Tag}ts$curr_tissue_tag$prev_tissue_tag${Tool_Tag}$rg_info_prefix";
    }
    else{
        cluck_and_exit "<ERROR>\tUnrecognized read-id prefix ($reads_OB->{pid}) from $rg_OB->{tissue} bam file:\n"
                            ."\t$rg_OB->{bam_OB}->{filepath}\n";
    }

    # inform
    stout_and_sterr "[INFO]\tExtract and update reads-id prefix of reads group $rg_OB->{RG_ID}:\n"
                         ."\tfrom $rg_OB->{rID_prefix_prev} to $rg_OB->{rID_prefix}\n";
}

#!!!!!!!!!!!!!#
# deprecated! #
#!!!!!!!!!!!!!#

#--- get allele at given position ---
# snv: single base matched ref/alt base
# ins: seq inserted behind 'pos'
# del: offset means moving N bases from 1st-del-pos is the 'pos', so seq deleted behind 'pos' has offset -1.
# priority: 'preferINDEL' > snv > ins/del
#* this version is done during visiting-study in Umich *#
#* deprecated on 2018/10/17 *#
sub get_pos_allele_v1_BaseOnCigar{

    my $reads_OB = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $preferINDEL = $parm{preferINDEL} || 0;
    my $no_warn = $parm{no_warn};

    # attributes
    my $mseg = $reads_OB->{mseg};
    my $mpos = $reads_OB->{mpos};
    my $rseq = uc($reads_OB->{rseq});
    my $mlen = $reads_OB->mReadLen;
    my $baseq = $reads_OB->{baseQ};
    my $CIGAR = $reads_OB->{cigar};

    # cigar unmap
    if( $reads_OB->is_unmap ){
        warn "<WARN>\treads $reads_OB->{pid} end-$reads_OB->{endNO} is unmapped.\n" unless $no_warn;
        return { type => 'miss', allele => undef };
    }
    # chr check
    if( $mseg ne $chr ){
        warn "<WARN>\treads $reads_OB->{pid} end-$reads_OB->{endNO} ($mseg) doesn't match query chr ($chr).\n" unless $no_warn;
        return { type => 'miss', allele => undef };
    }

    # exchange a[DN]bI to bIa[DN] in cigar
    if( $CIGAR =~ /\d+[DN]\d+I/ ){
        while( $CIGAR =~ /(\d+[DN])(\d+I)/    ){ my ($DNp,$Ip)=($1,$2);      $CIGAR =~ s/$DNp$Ip/$Ip$DNp/                };
        while( $CIGAR =~ /(\d+)([ID])(\d+)\2/ ){ my ($sum,$type)=($1+$3,$2); $CIGAR =~ s/${1}$type${3}$type/${sum}$type/ }; # connect aIbI or aDbD
    }

    # find the pos
    my $i = 0;
    my $dist_shift = 0;
    my $rCov_lastPos = -1; # M/D/N last covered pos
    my $rCov_lastBase; # M/D/N last covered base
    my $rCov_lastBaseQ; # M/D/N last covered base quality
    while( length($CIGAR) != 0 ){
        my ($len, $type) = ($CIGAR =~ /^(\d+)(\D)/);
        # check
        unless( defined $len && defined $type ){
            cluck_and_exit "<ERROR>\tcannot recognize cigar of this reads.\n"
                                ."\t$reads_OB->{pid},end-$reads_OB->{endNO},$reads_OB->{mseg},$reads_OB->{mpos},$reads_OB->{cigar}\n";
        }
        # shorten cigar
        $CIGAR =~ s/^\d+\D//;
        # different scenarios
        if($type eq 'M'){ # mapped with possible mismatch ONLY
            $rCov_lastPos = $mpos + $i + $len - 1;
            # M-part cover the want-pos
            if(    $rCov_lastPos > $pos
                || (   $rCov_lastPos == $pos
                    && (   !$preferINDEL        # ignore INDEL next to
                        || $CIGAR !~ /^\d+[ID]/ # or, next is not INDEL
                       )
                   )
            ){
                my $sub_i = $len + $pos - $rCov_lastPos - 1;
                my $allele = substr($rseq, $sub_i, 1);
                my $baseQ = substr($baseq, $sub_i, 1);
                my $dist5 = $i + $sub_i + 1 + $dist_shift;
                my $dist3 = $mlen - $dist5 + 1;
                return { type => 'refORsnv', allele => $allele, baseQ => $baseQ, dist5 => $dist5, dist3 => $dist3 };
            }
            # update for next block
            $rCov_lastBase = substr($rseq, $len-1, 1);
            $rCov_lastBaseQ = substr($baseq, $len-1, 1);
            $i += $len;
            $rseq = substr($rseq, $len);
            $baseq = substr($baseq, $len);
        }
        elsif($type eq 'I'){
            # inserts just behind the want-pos
            if( $rCov_lastPos == $pos ){
                my $allele = substr($rseq, 0, $len);
                my $baseQ = substr($baseq, 0, $len);
                my $dist5 = $i + $dist_shift; # forebase
                my $dist3 = $mlen - $dist5 - $len; # nextbase
                my $forebase = $rCov_lastBase;
                my $forebaseQ = $rCov_lastBaseQ;
                my $nextbase = substr($rseq, $len, 1);
                my $nextbaseQ = substr($baseq, $len, 1);
                return { type => 'ins', allele => $allele, baseQ => $baseQ, dist5 => $dist5, dist3 => $dist3,
                         forebase => $forebase, forebaseQ => $forebaseQ, nextbase => $nextbase, nextbaseQ => $nextbaseQ };
            }
            # update for next block
            $rseq  = substr($rseq, $len);
            $baseq = substr($baseq, $len);
            $dist_shift += $len;
        }
        elsif(   $type eq 'D'
              || $type eq 'N' # for RNA-seq reads
        ){
            my $del1stPos = $mpos + $i;
            # deletion covers the want-pos
            if(    $type eq 'D'
                && $pos <= $del1stPos + $len - 1
            ){
                my $delSize = $len;
                my $dist5 = $i + $dist_shift; # forebase
                my $dist3 = $mlen - $dist5; # nextbase
                my $forebase = $rCov_lastBase;
                my $forebaseQ = $rCov_lastBaseQ;
                my $nextbase = substr($rseq, 0, 1);
                my $nextbaseQ = substr($baseq, 0, 1);
                my $offset = $pos - $del1stPos; # -1 is 'forebase'
                return { type => 'del', allele => undef, baseQ => undef, dist5 => $dist5, dist3 => $dist3,
                         delSize => $delSize, offset => $offset,
                         forebase => $forebase, forebaseQ => $forebaseQ, nextbase => $nextbase, nextbaseQ => $nextbaseQ };
            }
            # update for next block
            $rCov_lastPos = $mpos + $i + $len - 1;
            $rCov_lastBase = undef; # fine, as 'D' must be followed by 'M', then get updated
            $rCov_lastBaseQ = undef;
            $i += $len;
            $dist_shift -= $len;
        }
        elsif($type eq 'S'){
            # just skip S part, update for next block
            $rseq  = substr($rseq, $len);
            $baseq = substr($baseq, $len);
        }
        elsif($type eq 'H'){
            1; # nothing to do
        }
    }

    # reach here, means this read doesn't cover the pos.
    warn "<WARN>\treads $reads_OB->{pid} end-$reads_OB->{endNO} ($mseg:$mpos:$reads_OB->{cigar}) doesn't cover query pos ($chr:$pos).\n" unless $no_warn;
    return { type => 'miss', allele => undef };
}

1; ## tell the perl script the successful access of this module.
