package BioFuse::BioInfo::Objects::Trans_OB;

use strict;
use warnings;
use Data::Dumper;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Trans_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        get_ENSid
                        get_ref_seg
                        get_strand
                        get_use_name
                        get_ori_name
                        get_biotype
                        get_exon_region
                        get_find_seq_mark
                        mark_find_seq
                     /;

#--- structure of object
## for gene PSL format, refer to module BioFuse::BioInfo::GeneAnno::GTF_transOB
# trans -> ENSid = $ENSid
# trans -> ref_seg = $ref_seg
# trans -> strand = $strand
# trans -> cytoband = $cytoband
# trans -> ori_name = $trans_ori_name
# trans -> use_name = $trans_use_name
# trans -> biotype = $trans_biotype
# trans -> gene_use_name = $gene_use_name
# trans -> exon = [ [st1,ed1], [st2,ed2], ... ]
# trans -> find_seq = [0]/1
##--- for protein coding trans ---##
# trans -> CDS  = [ [st1,ed1], [st2,ed2], ... ]
# trans -> protein_ENSid = $protein_ENSid
# trans -> start_codon = [ $pos1, $pos2, $pos3 ]
# trans -> stop_codon  = [ $pos1, $pos2, $pos3 ]
##--- specific for FuseSV virus meta-info ---##
# trans -> locus_tag = $locus_tag ()
# trans -> note = $note (specific for FuseSV virus meta-info)
# trans -> product = $product (specific for FuseSV virus meta-info)
## deprecated
# trans -> source = $gtf_source
# trans -> version = $trans_version
# trans -> exon_len = $exon_len (string)
# trans -> exon_smallEdge = $exon_smallEdge (string)
# trans -> exon_count = $exon_count

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;
    my $pslLine = $parm{pslLine};
    my $useExts = $parm{useExts} || 0; # load info specific for FuseSV virus meta-info

    my $trans = {};
    ## basic attributes
    my @ele = split /\s+/, $pslLine;
    $trans->{ref_seg}  = $ele[13];
    $trans->{strand}   = $ele[8];
    $trans->{cytoband} = $ele[14];
    $trans->{ENSid}    = $ele[10];
    $trans->{use_name} = $ele[9];
    $trans->{ori_name} = $ele[11];
    $trans->{small_edge} = $ele[15] + 1;
    $trans->{large_edge} = $ele[16];
    $trans->{biotype}  = $ele[12];
    $trans->{gene_use_name} = $ele[21];
    ## exon
    # $trans->{exon_len} = $ele[18]; # string
    # $trans->{exon_smallEdge} = $ele[20]; # string
    # $trans->{exon_count} = $ele[17];
    my @exon_len = split /,/, $ele[18];
    my @exon_smp = split /,/, $ele[20];
    $trans->{exon} = [ map { [$exon_smp[$_]+1, $exon_smp[$_]+$exon_len[$_]] } (0 .. $#exon_len) ];
    ## CDS and protein
    $trans->{CDS}  = [ map { /(\d+)\((\d+)\)/; [$1+1, $1+$2] } split /,/, $ele[19] ];
    $trans->{protein_ENSid} = $ele[6] if $ele[6] ne 'NA';
    $trans->{start_codon} = [ sort {$a<=>$b} split /,/, $ele[4] ] if $ele[4] ne ',';
    $trans->{stop_codon}  = [ sort {$a<=>$b} split /,/, $ele[5] ] if $ele[5] ne ',';
    ## other info
    # $trans->{source}   = $ele[0];
    # $trans->{version}  = $ele[1];
    $trans->{find_seq} = 0;
    ## specific for FuseSV virus meta-info
    if($useExts){
        $trans->{locus_tag} = $ele[2];
        $trans->{note}      = $ele[3];
        $trans->{product}   = $ele[7];
    }

    bless($trans);
    return $trans;
}

#--- return the ENSid ---
sub get_ENSid{
    my $trans = shift;
    return $trans->{ENSid};
}

#--- return the ref_seg ---
sub get_ref_seg{
    my $trans = shift;
    return $trans->{ref_seg};
}

#--- return the strand ---
sub get_strand{
    my $trans = shift;
    return $trans->{strand};
}

#--- return the trans use name ---
sub get_use_name{
    my $trans = shift;
    return $trans->{use_name};
}

#--- return the trans original name ---
sub get_ori_name{
    my $trans = shift;
    return $trans->{ori_name};
}

#--- return trans biotype ---
sub get_biotype{
    my $trans = shift;
    return $trans->{biotype};
}

#--- return Aref of exon region ---
sub get_exon_region{
    my $trans = shift;
    return $trans->{exon};
}

#--- return mark of find seq ---
sub get_find_seq_mark{
    my $trans = shift;
    return $trans->{find_seq};
}

#--- find seq, then mark as 1 ---
sub mark_find_seq{
    my $trans = shift;
    $trans->{find_seq} = 1;
}

1; ## tell the perl script the successful access of this module.
