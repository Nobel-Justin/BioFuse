package BioFuse::BioInfo::Objects::Allele::AlleleOnReads_OB;

use strict;
use warnings;
use List::Util qw/ sum /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::BioInfo::Quality qw/ baseQ_char2score /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Allele::AlleleOnReads_OB';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-11-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        setMiss
                        loadInfo
                        get_chr
                        get_pos
                        get_refBase
                        get_rOBmLen
                        get_type
                        get_qual
                        get_rEdgeDist5
                        get_rEdgeDist3
                        get_nAltDist5
                        get_nAltDist3
                        get_alleleSeq
                        get_deloffset
                        get_score
                     /;

#--- structure of object
# allele_OB -> chr = $chr
# allele_OB -> pos = $pos
# allele_OB -> refBase = $refBase, might have
# allele_OB -> rOBmLen = $rOBmLen
# allele_OB -> type = $type (ref/snv/ins/del)
# allele_OB -> allele = $allele
# allele_OB -> baseQ = $baseQ
# allele_OB -> rEdgeDist5 = $rEdgeDist5
# allele_OB -> rEdgeDist3 = $rEdgeDist3
# allele_OB -> nAltDist5 = $nAltDist5, distance to nearest alteration 5-prime
# allele_OB -> nAltDist3 = $nAltDist3, distance to nearest alteration 3-prime
# allele_OB -> forebase = $forebase, for indel
# allele_OB -> forebaseQ = $forebaseQ, for indel
# allele_OB -> nextbase = $nextbase, for indel
# allele_OB -> nextbaseQ = $nextbaseQ, for indel
# allele_OB -> offset = $offest, for del
# allele_OB -> delSize = $delSize, for del

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $allele_OB = {};
    $allele_OB->{chr} = $parm{chr};
    $allele_OB->{pos} = $parm{pos};
    $allele_OB->{refBase} = $parm{refBase} || undef;
    $allele_OB->{rOBmLen} = $parm{rOBmLen};
    $allele_OB->{type} = undef;
    $allele_OB->{allele} = undef;

    bless($allele_OB);
    return $allele_OB;
}

#--- set this allele is missed ---
sub setMiss{
    my $allele_OB = shift;
    $allele_OB->{type} = 'miss';
    $allele_OB->{allele} = undef;
}

#--- record given info ---
sub loadInfo{
    my $allele_OB = shift;
    my %parm = @_;
    $allele_OB->{$_} = $parm{$_} for keys %parm;
}

#--- return allele chr ---
sub get_chr{
    my $allele_OB = shift;
    return $allele_OB->{chr};
}

#--- return allele pos ---
sub get_pos{
    my $allele_OB = shift;
    return $allele_OB->{pos};
}

#--- return allele ref-Base ---
sub get_refBase{
    my $allele_OB = shift;
    return $allele_OB->{refBase};
}

#--- return rOB mapped length ---
sub get_rOBmLen{
    my $allele_OB = shift;
    return $allele_OB->{rOBmLen};
}

#--- return allele type ---
sub get_type{
    my $allele_OB = shift;
    return $allele_OB->{type};
}

#--- return score of this allele ---
## score = NORM_dist_score + NORM_qual_score
## NORM_dist_score = sqrt( rEdgeDist_score * aNearDist_score ) / r_mReadLen
## rEdgeDist_score = sqrt( rEdgeDist5 * rEdgeDist3 );
## aNearDist_score = sqrt( nAltDist5 * nAltDist3 );
## NORM_qual_score = get_qual / 40
## PS. base-quality range approximates to 40
sub get_score{
    my $allele_OB = shift;
    my %parm = @_;
    my $Q_offset = $parm{Q_offset} || 33; # Sanger and Illumina 1.8+.

    return sqrt(   sqrt( $allele_OB->get_rEdgeDist5 * $allele_OB->get_rEdgeDist3 )
                 * sqrt( $allele_OB->get_nAltDist5  * $allele_OB->get_nAltDist3  )
               ) / $allele_OB->get_rOBmLen
           + $allele_OB->get_qual(Q_offset => $Q_offset) / 40;
}

#--- return allele quality ---
sub get_qual{
    my $allele_OB = shift;
    my %parm = @_;
    my $Q_offset = $parm{Q_offset} || 33; # Sanger and Illumina 1.8+.

    if(    $allele_OB->{type} eq 'ref'
        || $allele_OB->{type} eq 'snv'
    ){
        return baseQ_char2score(Q_char => $allele_OB->{baseQ}, Q_offset => $Q_offset);
    }
    elsif(   $allele_OB->{type} eq 'ins'
          || $allele_OB->{type} eq 'del'
    ){
        my @baseQ = ( $allele_OB->{forebaseQ}, $allele_OB->{nextbaseQ} );
        push @baseQ, split //, $allele_OB->{baseQ} if( defined $allele_OB->{baseQ} ); # baseQ of inserted-seq
        return sum( map { ( baseQ_char2score(Q_char=>$_, Q_offset=>$Q_offset) ) } @baseQ ) / scalar(@baseQ);
    }
    # not available
    return -1;
}

#--- return distance to 5-prime reads-edge ---
sub get_rEdgeDist5{
    my $allele_OB = shift;
    return $allele_OB->{rEdgeDist5};
}

#--- return distance to 3-prime reads-edge ---
sub get_rEdgeDist3{
    my $allele_OB = shift;
    return $allele_OB->{rEdgeDist3};
}

#--- return distance to 5-prime nearest alteration ---
sub get_nAltDist5{
    my $allele_OB = shift;
    return $allele_OB->{nAltDist5};
}

#--- return distance to 3-prime nearest alteration ---
sub get_nAltDist3{
    my $allele_OB = shift;
    return $allele_OB->{nAltDist3};
}

#--- return allele sequence ---
sub get_alleleSeq{
    my $allele_OB = shift;
    return $allele_OB->{allele};
}

#--- return allele del-Size ---
sub get_delSize{
    my $allele_OB = shift;
    return $allele_OB->{delSize};
}

#--- return allele del-offset ---
sub get_deloffset{
    my $allele_OB = shift;
    return $allele_OB->{offset};
}

1; ## tell the perl script the successful access of this module.
