package BioFuse::BioInfo::GeneAnno::GTF_geneOB;

use strict;
use warnings;
use List::Util qw/ min max /;
use BioFuse::BioInfo::GeneAnno::GTF_transOB;
use BioFuse::BioInfo::CytoBand qw/ get_cytoband /;
use BioFuse::Util::Array qw/ mergeOverlap /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::GTF_geneOB';
#----- version --------
$VERSION = "0.07";
$DATE = '2018-11-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        load_gtf_info
                        add_cytoband_info
                        get_ENSid
                        get_ref_seg
                        get_strand
                        get_gene_use_name
                        get_gene_ori_name
                        get_trans_use_name
                        get_trans_ori_name
                        get_transOB_Aref
                        refine_gtf_info
                        refined_use_name
                        get_gene_psl_line
                     /;

#--- structure of object
# gene -> ENSid = $ENSid
# gene -> ref_seg = $ref_seg
# gene -> strand = $strand
# gene -> cytoband = $cytoband
# gene -> ori_name = $gene_ori_name
# gene -> use_name = $gene_use_name
# gene -> source = $gtf_source
# gene -> version = $gene_version
# gene -> biotype = $gene_biotype
# gene -> isoforms = {trans_ENSid_1=>$GTF_transOB_1, trans_ENSid_2=>$GTF_transOB_2, .. }

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;
    my $gtfLineOB = $parm{gtfLineOB};

    my $gene = {};
    # load basic attributes of gene
    $gene->{ENSid}    = $gtfLineOB->get_ENSid(type => 'gene');
    $gene->{source}   = $gtfLineOB->get_source(type => 'gene');
    $gene->{version}  = $gtfLineOB->get_version(type => 'gene');
    $gene->{biotype}  = $gtfLineOB->get_biotype(type => 'gene');
    $gene->{ori_name} = $gtfLineOB->get_name(type => 'gene');
    $gene->{use_name} = $gene->{ori_name};
    $gene->{ref_seg}  = $gtfLineOB->get_refSeg;
    $gene->{strand}   = $gtfLineOB->get_strand;
    $gene->{isoforms} = {};

    bless($gene);
    return $gene;
}

#--- load the gtf info to the GTF_geneOB object ---
# e.g., transcript, exon, CDS, codon, ...
# these all should belong to one transcript
sub load_gtf_info{
    my $gene = shift;
    my %parm = @_;
    my $gtfLineOB = $parm{gtfLineOB};

    # initialize one GTF_transOB object
    my $trans_ENSid = $gtfLineOB->get_ENSid(type => 'trans');
    unless(exists $gene->{isoforms}->{$trans_ENSid}){
        $gene->{isoforms}->{$trans_ENSid} = BioFuse::BioInfo::GeneAnno::GTF_transOB->new(gtfLineOB => $gtfLineOB);
    }
    # load other information, exon, CDS, codon, ...
    my $trans = $gene->{isoforms}->{$trans_ENSid};
    my $region_type = $gtfLineOB->get_regionType;
    my @region_parm = (st => $gtfLineOB->get_stpos, ed => $gtfLineOB->get_edpos);
    # exon info of transcript
    $trans->load_exon_region(@region_parm) if $region_type =~ /exon/i;
    # CDS info of transcript, if available
    $trans->load_CDS_region(@region_parm, gtfLineOB => $gtfLineOB) if $region_type =~ /CDS/i;
    # Codon
    $trans->load_codon_info(@region_parm, codon_type => lc($region_type)) if $region_type =~ /codon/i;
}

#--- add cytoband info ---
sub add_cytoband_info{
    my $gene = shift;
    my %parm = @_;
    $gene->{cytoband} = get_cytoband(cytoBand_Href => $parm{cytoBand_Href}, refseg => $gene->{ref_seg},
                                     regionStP => $gene->{small_edge}, regionEdP => $gene->{large_edge});
    # trans
    $_->add_cytoband_info(cytoBand_Href => $parm{cytoBand_Href}) for values %{$gene->{isoforms}};
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
sub get_gene_use_name{
    my $gene = shift;
    return $gene->{use_name};
}

#--- return the gene original name ---
sub get_gene_ori_name{
    my $gene = shift;
    return $gene->{ori_name};
}

#--- return Aref of the trans use names ---
sub get_trans_use_name{
    my $gene = shift;
    return [map {$_->get_trans_use_name} values %{$gene->{isoforms}}];
}

#--- return Aref of the trans original names ---
sub get_trans_ori_name{
    my $gene = shift;
    return [map {$_->get_trans_ori_name} values %{$gene->{isoforms}}];
}

#--- return Aref of GTF trans OB of this gene OB ---
sub get_transOB_Aref{
    my $gene = shift;
    return [values %{$gene->{isoforms}}];
}

#--- note trans_biotype with abnormal number of codon ---
sub refine_gtf_info{
    my $gene = shift;
    my @gene_edge;
    for my $trans (values %{$gene->{isoforms}}){
        # refine protein's info
        $trans->note_abnor_codon_num;
        # exon info
        $trans->refine_exon;
        push @gene_edge, @{$trans->get_edge};
    }
    # gene region
    $gene->{small_edge} = min(@gene_edge);
    $gene->{large_edge} = max(@gene_edge);
}

#--- add ENSid to use name to refine duplication ---
sub refined_use_name{
    my $gene = shift;
    $gene->{use_name} = $gene->{ori_name} . '_' . $gene->{ENSid};
}

#--- return the gene PSL line string ---
sub get_gene_psl_line{
    my $gene = shift;
    # exon region info
    my %exon;
    for my $trans (values %{$gene->{isoforms}}){
        $exon{$_->[0]}{($_->[1]-$_->[0]+1)} = 1 for @{$trans->get_exon_Aref};
    }
    # de-dup
    my @exon_small_edge;
    my @exon_len;
    my @all_exon;
    for my $exon_small_edge (sort {$a<=>$b} keys %exon){
        for my $exon_len (sort {$a<=>$b} keys %{$exon{$exon_small_edge}}){
            push @exon_small_edge, $exon_small_edge - 1;
            push @exon_len, $exon_len;
            push @all_exon, [$exon_small_edge, $exon_small_edge+$exon_len-1];
        }
    }
    # exon count
    my $exon_number = scalar(@exon_small_edge);
    # merged exon region
    mergeOverlap(regionAref => \@all_exon, mergeAdjacent => 1);
    my @merged_exon = map {($_->[0]-1).'('.($_->[1]-$_->[0]+1).')'} @all_exon;

    #------- Format defination of gene PSL file ------------#
    #------- customized by SOAPfuse ------------------------#
    # Note: all the smallest edge position of region has been subtracted by 1, other info, such as largest edge
    #       or point position remains no changes.
    #  column_NO.    Description
    #           1    gene source
    #           2    gene version
    #           3           Reserved
    #           4           Reserved
    #           5           Reserved
    #           6           Reserved
    #           7           Reserved
    #           8           Reserved
    #           9    sense strand
    #          10    gene name for bioinformatics analysis
    #          11    gene ENS_id
    #          12    gene name (original) from the GTF file
    #          13    gene biotype
    #          14    refseg
    #          15    cytoband on refseg
    #          16    the smallest position on the plus strand of the refseg, and subtract 1
    #          17    the largest  position on the plus strand of the refseg, no modification
    #          18    the number of exon
    #          19    length of each exon, corresponding to the NO.21 column
    #          20    merged exon region info, format: smallest_position(len), sorted by small->large
    #          21    the smallest position on the plus strand of the refseg of each exon, sorted by small->large
    #-------------------------------------------------------#
    my @psl_line;
    push @psl_line, $gene->{source};
    push @psl_line, $gene->{version};
    push @psl_line, '0' for (3 .. 8);
    push @psl_line, $gene->{strand};
    push @psl_line, $gene->{use_name};
    push @psl_line, $gene->{ENSid};
    push @psl_line, $gene->{ori_name};
    push @psl_line, $gene->{biotype};
    push @psl_line, $gene->{ref_seg};
    push @psl_line, $gene->{cytoband} || 'NA';
    push @psl_line, $gene->{small_edge} - 1;
    push @psl_line, $gene->{large_edge};
    push @psl_line, $exon_number;
    push @psl_line, join(',',@exon_len).',';
    push @psl_line, join(',',@merged_exon).',';
    push @psl_line, join(',',@exon_small_edge).',';

    return join("\t", @psl_line);
}

#--- 
1; ## tell the perl script the successful access of this module.
