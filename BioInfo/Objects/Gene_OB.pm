package BioFuse::BioInfo::Objects::Gene_OB;

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

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Gene_OB';
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
                        get_exon_region
                        get_find_seq_mark
                        mark_find_seq
                     /;

#--- structure of object
## for gene PSL format, refer to module BioFuse::BioInfo::GeneAnno::GTF_geneOB
# gene -> ENSid = $ENSid
# gene -> ref_seg = $ref_seg
# gene -> strand = $strand
# gene -> cytoband = $cytoband
# gene -> ori_name = $gene_ori_name
# gene -> use_name = $gene_use_name
# gene -> biotype = $gene_biotype
# gene -> small_edge = $small_edge
# gene -> large_edge = $large_edge
# gene -> merged_exon = [ [st1,ed1], [st2,ed2], ... ]
# gene -> find_seq = [0]/1
## deprecated
# gene -> source = $gtf_source
# gene -> version = $gene_version
# gene -> exon_len = $exon_len (string)
# gene -> exon_smallEdge = $exon_smallEdge (string)
# gene -> exon_count = $exon_count

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;
    my $pslLine = $parm{pslLine};

    my $gene = {};
    ## basic attributes
    my @ele = split /\s+/, $pslLine;
    $gene->{ref_seg}  = $ele[13];
    $gene->{strand}   = $ele[8];
    $gene->{cytoband} = $ele[14];
    $gene->{ENSid}    = $ele[10];
    $gene->{use_name} = $ele[9];
    $gene->{ori_name} = $ele[11];
    $gene->{small_edge} = $ele[15] + 1;
    $gene->{large_edge} = $ele[16];
    $gene->{biotype}  = $ele[12];
    ## all uniq exon
    # $gene->{exon_count} = $ele[17];
    # $gene->{exon_len} = $ele[18]; # string
    # $gene->{exon_smallEdge} = $ele[20]; # string
    ## merged exon
    $gene->{merged_exon} = [ map { /(\d+)\((\d+)\)/; [$1+1, $1+$2] } split /,/, $ele[19] ];
    ## other info
    # $gene->{source}   = $ele[0];
    # $gene->{version}  = $ele[1];
    $gene->{find_seq} = 0;

    bless($gene);
    return $gene;
}

#--- return the ENSid ---
sub get_ENSid{
    my $gene = shift;
    return $gene->{ENSid};
}

#--- return the ref_seg ---
sub get_ref_seg{
    my $gene = shift;
    return $gene->{ref_seg};
}

#--- return the strand ---
sub get_strand{
    my $gene = shift;
    return $gene->{strand};
}

#--- return the gene use name ---
sub get_use_name{
    my $gene = shift;
    return $gene->{use_name};
}

#--- return the gene original name ---
sub get_ori_name{
    my $gene = shift;
    return $gene->{ori_name};
}

#--- return Aref of exon region ---
sub get_exon_region{
    my $gene = shift;
    return $gene->{merged_exon};
}

#--- return mark of find seq ---
sub get_find_seq_mark{
    my $gene = shift;
    return $gene->{find_seq};
}

#--- find seq, then mark as 1 ---
sub mark_find_seq{
    my $gene = shift;
    $gene->{find_seq} = 1;
}

1; ## tell the perl script the successful access of this module.
