package BioFuse::BioInfo::Objects::GeneAnno::GTF_transOB;

use strict;
use warnings;
use List::Util qw/ min max /;
use BioFuse::BioInfo::CytoBand qw/ get_cytoband /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::GeneAnno::GTF_transOB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        load_exon_region
                        load_CDS_region
                        load_codon_info
                        add_cytoband_info
                        get_ENSid
                        get_ref_seg
                        get_strand
                        get_biotype
                        get_trans_use_name
                        get_trans_ori_name
                        get_edge
                        get_startCodonPos
                        get_stopCodonPos
                        get_exon_Aref
                        get_CDS_Aref
                        update_biotype
                        note_abnor_codon_num
                        refine_exon
                        refined_use_name
                        get_trans_psl_line
                     /;

#--- structure of object
# trans -> ENSid = $ENSid
# trans -> ref_seg = $ref_seg
# trans -> strand = $strand
# trans -> cytoband = $cytoband
# trans -> ori_name = $trans_ori_name
# trans -> use_name = $trans_use_name
# trans -> source = $gtf_source
# trans -> version = $trans_version
# trans -> biotype = $trans_biotype
# trans -> exon = { $st1=>$ed1, $st2=>$ed2, .. }
##--- for protein coding trans ---##
# trans -> CDS  = { $st1=>$ed1, $st2=>$ed2, .. }
# trans -> protein_ENSid = $protein_ENSid
# trans -> start_codon = { $pos1=>1, $pos2=>1, $pos3=>1 }
# trans -> stop_codon  = { $pos1=>1, $pos2=>1, $pos3=>1 }
##--- specific for FuseSV virus meta-info ---##
# trans -> locus_tag = $locus_tag ()
# trans -> note = $note (specific for FuseSV virus meta-info)
# trans -> product = $product (specific for FuseSV virus meta-info)

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;
    my $gtfLineOB = $parm{gtfLineOB};

    my $trans = {};
    # load basic attributes
    $trans->{ENSid}     = $gtfLineOB->get_ENSid(type => 'trans');
    $trans->{source}    = $gtfLineOB->get_source(type => 'trans');
    $trans->{version}   = $gtfLineOB->get_version(type => 'trans');
    $trans->{biotype}   = $gtfLineOB->get_biotype(type => 'trans');
    $trans->{ori_name}  = $gtfLineOB->get_name(type => 'trans');
    $trans->{use_name}  = $trans->{ori_name};
    $trans->{ref_seg}   = $gtfLineOB->get_refSeg;
    $trans->{strand}    = $gtfLineOB->get_strand;
    ## specific for FuseSV virus meta-info
    $trans->{locus_tag} = $gtfLineOB->get_transLocusTag;
    $trans->{note}      = $gtfLineOB->get_transNote;
    $trans->{product}   = $gtfLineOB->get_transProduct;

    bless($trans);
    return $trans;
}

#--- load exon region ---
sub load_exon_region{
    my $trans = shift;
    my %parm = @_;
    $trans->{exon}->{$parm{st}} = $parm{ed};
}

#--- load CDS region ---
sub load_CDS_region{
    my $trans = shift;
    my %parm = @_;
    $trans->{CDS}->{$parm{st}} = $parm{ed};
    $trans->{protein_ENSid} = $parm{gtfLineOB}->get_proteinENSid unless defined $trans->{protein_ENSid};
}

#--- load codon info ---
sub load_codon_info{
    my $trans = shift;
    my %parm = @_;
    $trans->{$parm{codon_type}}->{$_} = 1 for ($parm{st} .. $parm{ed});
}

#--- add cytoband info ---
sub add_cytoband_info{
    my $trans = shift;
    my %parm = @_;
    $trans->{cytoband} = get_cytoband(cytoBand_Href => $parm{cytoBand_Href}, refseg => $trans->{ref_seg},
                                      regionStP => $trans->{small_edge}, regionEdP => $trans->{large_edge});
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

#--- return trans biotype ---
sub get_biotype{
    my $trans = shift;
    return $trans->{biotype};
}

#--- return the trans use name ---
sub get_trans_use_name{
    my $trans = shift;
    return $trans->{use_name};
}

#--- return the trans original name ---
sub get_trans_ori_name{
    my $trans = shift;
    return $trans->{ori_name};
}

#--- return [small_edge, large_edge] ---
sub get_edge{
    my $trans = shift;
    return [$trans->{small_edge}, $trans->{large_edge}];
}

#--- return Aref of sorted start codon position ---
sub get_startCodonPos{
    my $trans = shift;
    return [sort {$a<=>$b} keys %{$trans->{start_codon}}];
}

#--- return Aref of sorted stop codon position ---
sub get_stopCodonPos{
    my $trans = shift;
    return [sort {$a<=>$b} keys %{$trans->{stop_codon}}];
}

#--- return Aref of trans exon region ---
sub get_exon_Aref{
    my $trans = shift;
    return [map {[$_,$trans->{exon}->{$_}]} sort {$a<=>$b} keys %{$trans->{exon}}];
}

#--- return Aref of trans CDS region ---
sub get_CDS_Aref{
    my $trans = shift;
    return [map {[$_,$trans->{CDS}->{$_}]} sort {$a<=>$b} keys %{$trans->{CDS}}];
}

#--- update biotype ---
sub update_biotype{
    my $trans = shift;
    my %parm = @_;
    $trans->{biotype} = $parm{biotype};
}

#--- note the abnormal start/stop codon num ---
sub note_abnor_codon_num{
    my $trans = shift;
    if($trans->get_biotype eq 'protein_coding'){
        my $st_codon = scalar keys %{$trans->{start_codon}};
        my $sp_codon = scalar keys %{$trans->{stop_codon}};
        if($st_codon != 3 && $sp_codon != 3){
            $trans->update_biotype(biotype => "protein_coding_with_abnor_codon_num($st_codon,$sp_codon)");
        }
        elsif($st_codon != 3 && $sp_codon == 3){
            $trans->update_biotype(biotype => "protein_coding_with_abnor_st_codon_num($st_codon)");
        }
        elsif($st_codon == 3 && $sp_codon != 3){
            $trans->update_biotype(biotype => "protein_coding_with_abnor_sp_codon_num($sp_codon)");
        }
    }
}

#--- refine exon info ----
sub refine_exon{
    my $trans = shift;
    # exon number
    $trans->{exon_number} = scalar keys %{$trans->{exon}};
    # get trans bilateral positions
    $trans->{small_edge} = min(keys   %{$trans->{exon}});
    $trans->{large_edge} = max(values %{$trans->{exon}});
}

#--- add ENSid to use name to refine duplication ---
sub refined_use_name{
    my $trans = shift;
    $trans->{use_name} = $trans->{ori_name} . '_' . $trans->{ENSid};
}

#--- return the trans PSL line string ---
sub get_trans_psl_line{
    my $trans = shift;
    my %parm = @_;
    my $show_metainfo = $parm{show_metainfo} || 0;
    my $gene_use_name = $parm{gene_use_name};

    #------- Format defination of trans PSL file -----------#
    #------- customized by SOAPfuse ------------------------#
    # Note: all the smallest edge position of region has been subtracted by 1, other info, such as largest edge
    #       or point position (e.g., start_codon) remains no changes.
    #  column_NO.    Description
    #           1    transcript source
    #           2    transcript version
    #           3           Reserved (maybe 'locus_tag' specific for FuseSV virus meta-info)
    #           4           Reserved (maybe 'note' specific for FuseSV virus meta-info)
    #           5    start_codon info, if available
    #           6    stop_codon info, if available
    #           7    protein ENS_id, if available
    #           8           Reserved (maybe 'product' specific for FuseSV virus meta-info)
    #           9    sense strand
    #          10    transcript name for bioinformatics analysis
    #          11    transcript ENS_id
    #          12    transcript name (original) from the GTF file
    #          13    biotype
    #          14    refseg
    #          15    cytoband on refseg
    #          16    the smallest position on the plus strand of the refseg, and subtract 1
    #          17    the largest  position on the plus strand of the refseg, no modification
    #          18    the number of exon
    #          19    length of each exon, corresponding to the NO.21 column
    #          20    CDS region info, if available format: smallest_position(len), sorted by small->large
    #          21    the smallest position on the plus strand of the refseg of each exon, sorted by small->large
    #          22    gene name for bioinformatics analysis
    #-------------------------------------------------------#
    my @psl_line;
    push @psl_line, $trans->{source};
    push @psl_line, $trans->{version};
    push @psl_line, ($show_metainfo ? $trans->{locus_tag} : '0');
    push @psl_line, ($show_metainfo ? $trans->{note} : '0');
    push @psl_line, join(',', @{$trans->get_startCodonPos}).',';
    push @psl_line, join(',', @{$trans->get_stopCodonPos}) .',';
    push @psl_line, $trans->{protein_ENSid} || 'NA';
    push @psl_line, ($show_metainfo ? $trans->{product} : '0');
    push @psl_line, $trans->{strand};
    push @psl_line, $trans->{use_name};
    push @psl_line, $trans->{ENSid};
    push @psl_line, $trans->{ori_name};
    push @psl_line, $trans->{biotype};
    push @psl_line, $trans->{ref_seg};
    push @psl_line, $trans->{cytoband} || 'NA';
    push @psl_line, $trans->{small_edge} - 1;
    push @psl_line, $trans->{large_edge};
    push @psl_line, $trans->{exon_number};
    push @psl_line, join(',', map {$_->[1]-$_->[0]+1} @{$trans->get_exon_Aref}) . ',';
    push @psl_line, join(',', map {($_->[0]-1).'('.($_->[1]-$_->[0]+1).')'} @{$trans->get_CDS_Aref}) . ',';
    push @psl_line, join(',', map {$_->[0]-1} @{$trans->get_exon_Aref}) . ',';
    push @psl_line, $gene_use_name;

    return join("\t", @psl_line);
}

#--- 
1; ## tell the perl script the successful access of this module.
