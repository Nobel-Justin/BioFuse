package BioFuse::BioInfo::Objects::GeneAnno::PSL_OB;

use strict;
use warnings;
use Data::Dumper;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ stout_and_sterr /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;
use BioFuse::BioInfo::Objects::GeneAnno::Gene_OB;
use BioFuse::BioInfo::Objects::GeneAnno::Trans_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::GeneAnno::PSL_OB';
#----- version --------
$VERSION = "0.05";
$DATE = '2019-05-04';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        load_GeneOrTrans_from_PSL
                        get_objPOOLHf
                        extract_GeneOrTrans_seq
                        output_exon_seq
                     /;

#--- structure of object
# psl -> filePath = $filePath
# psl -> psl_type = $psl_type (gene/trans)
# psl -> objPOOL = {}
# psl -> objPOOLRecKey = $key (refseg/name)

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $psl = {};
    # load basic attributes of psl
    $psl->{filePath} = $parm{filePath};
    $psl->{psl_type} = $parm{psl_type} || 'gene'; # or trans
    $psl->{objPOOL} = {}; # main container

    bless($psl);
    return $psl;
}

#--- read PSL file ---
sub load_GeneOrTrans_from_PSL{
    my $psl = shift;
    my %parm = @_;
    my $useExts = $parm{useExts} || 0; # load info specific for FuseSV virus meta-info
    my $recordKey = $parm{recordKey} || 'refseg'; # genename
    # filter
    my $skipRefSegHref = $parm{skipRefSegHref}; # refseg
    my $requireObjHref = $parm{requireObjHref}; # original name
    my $requireBtyHref = $parm{requireBtyHref}; # biotype
    # reset
    $skipRefSegHref = undef if defined $skipRefSegHref && ! scalar keys %$skipRefSegHref;
    $requireObjHref = undef if defined $requireObjHref && ! scalar keys %$requireObjHref;
    $requireBtyHref = undef if defined $requireBtyHref && ! scalar keys %$requireBtyHref;

    # reset attributes
    $psl->{objPOOL} = {};
    $psl->{objPOOLRecKey} = $recordKey;
    # gene or trans object
    my $objModule = $psl->{psl_type} eq 'gene' ? 'BioFuse::BioInfo::Objects::GeneAnno::Gene_OB' : 'BioFuse::BioInfo::Objects::GeneAnno::Trans_OB';
    stout_and_sterr "[INFO]\tapply object module $objModule for $psl->{psl_type} PSL file.\n";
    # read PSL
    open (PSL, Try_GZ_Read($psl->{filePath})) || die"fail psl_file: $!\n";
    while(<PSL>){
        next if /^#/;
        # create object
        my $object = $objModule->new(pslLine => $_, useExts => $useExts);
        # filter
        next if defined $skipRefSegHref &&  exists $skipRefSegHref->{$object->get_ref_seg};
        next if defined $requireObjHref && !exists $requireObjHref->{$object->get_ori_name};
        next if defined $requireBtyHref && !exists $requireBtyHref->{$object->get_biotype};
        # record
        my $key;
        if($recordKey eq 'refseg'){
            $key = $object->get_ref_seg;
        }
        elsif($recordKey eq 'genename'){
            $key = $psl->{psl_type} eq 'gene' ? $object->get_use_name : $object->get_gene_use_name;
        }
        push @{$psl->{objPOOL}->{$key}}, $object;
    }
    close PSL;
    # inform
    stout_and_sterr "[INFO]\tload information from $psl->{psl_type} PSL OK.\n"
                         ."\t$psl->{filePath}\n";
}

#--- return objPOOL Href ---
sub get_objPOOLHf{
    my $psl = shift;
    return $psl->{objPOOL};
}

#--- output base ---
sub extract_GeneOrTrans_seq{
    my $psl = shift;
    my %parm = @_;
    my $outputFasta = $parm{outputFasta};
    my $whole_genome = $parm{whole_genome};
    my $extendLength = $parm{extendLength} || 0;

    # map refseg to objects
    my $refSegToObjHf = {};
    if($psl->{objPOOLRecKey} eq 'refseg'){
        $refSegToObjHf = $psl->{objPOOL};
    }
    else{
        for my $objAf (values %{$psl->{objPOOL}}){
            push @{$refSegToObjHf->{$_->get_ref_seg}}, $_ for @$objAf;
        }
    }
    # read fasta file
    open (my $outfh, Try_GZ_Write($outputFasta)) || die "fail write output Fasta: $!\n";
    read_fasta_file( FaFile => $whole_genome, needSeg_Href => $refSegToObjHf,
                     subrtRef => \&output_exon_seq,
                     subrtParmAref => [fh => $outfh, refSegToObjHf => $refSegToObjHf, extendLength => $extendLength] );
    close $outfh;
    # check the remained unit
    for my $ref_seg (sort keys %$refSegToObjHf){
        stout_and_sterr "<WARN>\tseqeuence of ".$_->get_use_name." from Refseg $ref_seg remains undetected.\n" for grep !$_->get_find_seq_mark, @{$refSegToObjHf->{$ref_seg}};
    }
    # inform
    stout_and_sterr "[INFO]\tcreate Exon_Seq file OK!\n"
                         ."\t$outputFasta\n";
}

#--- extract exon region from given chr-seq ---
sub output_exon_seq{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $segName = $parm{segName};
    my $segSeq_Sref = $parm{segSeq_Sref};
    my $fh = $parm{fh};
    my $refSegToObjHf = $parm{refSegToObjHf};
    my $extendLength = $parm{extendLength} || 0;

    my $line_base = 50;
    for my $object (sort {$a->get_use_name cmp $b->get_use_name} @{$refSegToObjHf->{$segName}}){
        # get exon sequence along plus strand orientation
        my $exon_reg = $object->get_exon_region;
        my $exon_seq = join('', map {substr($$segSeq_Sref, $_->[0]-1, $_->[1]-$_->[0]+1)} @$exon_reg);
        # extend?
        if($extendLength){
            $exon_seq = substr($$segSeq_Sref, $exon_reg->[ 0]->[0]-$extendLength-1, $extendLength)
                       .$exon_seq
                       .substr($$segSeq_Sref, $exon_reg->[-1]->[1], $extendLength);
        }
        # capital letter
        $exon_seq = uc $exon_seq;
        # reversed complementary?
        ($exon_seq = reverse $exon_seq) =~ tr/ACGT/TGCA/ if $object->get_strand eq '-';
        # output
        print {$fh} ">".$object->get_use_name."\n";
        print {$fh} substr($exon_seq,$_*$line_base,$line_base)."\n" for (0 .. int((length($exon_seq)-1)/$line_base));
        # once ok delete the unit
        $object->mark_find_seq;
    }
}

1; ## tell the perl script the successful access of this module.
