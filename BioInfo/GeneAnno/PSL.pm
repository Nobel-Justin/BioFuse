package BioFuse::BioInfo::GeneAnno::PSL;

use strict;
use warnings;
use Data::Dumper;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ stout_and_sterr /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;
use BioFuse::BioInfo::Objects::Gene_OB;
use BioFuse::BioInfo::Objects::Trans_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              load_GeneOrTrans_from_PSL
              extract_GeneOrTrans_seq
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::PSL';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-11-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        load_GeneOrTrans_from_PSL
                        extract_GeneOrTrans_seq
                        output_exon_seq
                     /;

#--- read PSL file ---
sub load_GeneOrTrans_from_PSL{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ObjectPoolHref = $parm{ObjectPoolHref};
    my $psl_file = $parm{psl_file};
    my $psl_type = $parm{psl_type} || 'gene'; # or trans
    my $useExts = $parm{useExts} || 0; # load info specific for FuseSV virus meta-info

    # gene or trans object
    my $objModule = $psl_type eq 'gene' ? 'BioFuse::BioInfo::Objects::Gene_OB' : 'BioFuse::BioInfo::Objects::Trans_OB';
    stout_and_sterr "[INFO]\tapply object module $objModule for $psl_type PSL file.\n";
    # read PSL
    open (PSL, Try_GZ_Read($psl_file)) || die"fail $psl_file: $!\n";
    while(<PSL>){
        next if /^#/;
        # create object
        my $object = $objModule->new(pslLine => $_, useExts => $useExts);
        # filter
        next if defined $parm{skipRefSegHref} && exists $parm{skipRefSegHref}->{$object->get_ref_seg};
        next if defined $parm{requireObjHref} && exists $parm{requireObjHref}->{$object->get_ori_name};
        # record
        push @{$ObjectPoolHref->{$object->get_ref_seg}}, $object;
    }
    close PSL;
    # inform
    stout_and_sterr "[INFO]\tread exon region from $psl_type PSL ok.\n"
                         ."\t$psl_file\n";
}

#--- output base ---
sub extract_GeneOrTrans_seq{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ObjectPoolHref = $parm{ObjectPoolHref};
    my $outputFasta = $parm{outputFasta};
    my $whole_genome = $parm{whole_genome};
    my $extendLength = $parm{extendLength} || 0;

    # read fasta file
    open (my $outfh, Try_GZ_Write($outputFasta)) || die "fail write output Fasta: $!\n";
    read_fasta_file( FaFile => $whole_genome, needSeg_Href => $ObjectPoolHref,
                     subrtRef => \&output_exon_seq,
                     subrtParmAref => [fh => $outfh, ObjectPoolHref => $ObjectPoolHref, extendLength => $extendLength] );
    close $outfh;

    # check the remained unit
    for my $ref_seg (sort keys %$ObjectPoolHref){
        stout_and_sterr "<WARN>\tseqeuence of ".$_->get_use_name." from Refseg $ref_seg remains undetected.\n" for grep !$_->get_find_seq_mark, @{$ObjectPoolHref->{$ref_seg}};
    }
    # inform
    stout_and_sterr "[INFO]\tcreate Exon_Seq file ok!\n"
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
    my $ObjectPoolHref = $parm{ObjectPoolHref};
    my $extendLength = $parm{extendLength} || 0;

    my $line_base = 50;
    for my $object (sort {$a->get_use_name cmp $b->get_use_name} @{$ObjectPoolHref->{$segName}}){
        # get exon sequence along plus strand orientation
        my $exon_reg = $object->get_exon_region;
        my $exon_seq = join('', map {substr($$segSeq_Sref, $_->[0]-1, $_->[1]-$_->[0]+1)} @$exon_reg);
        # extend?
        if($extendLength){
            $exon_seq = substr($$segSeq_Sref, $exon_reg->[ 0]->[0]-1, $exon_reg->[ 0]->[1]-$exon_reg->[ 0]->[0]+1)
                       .$exon_seq
                       .substr($$segSeq_Sref, $exon_reg->[-1]->[0]-1, $exon_reg->[-1]->[1]-$exon_reg->[-1]->[0]+1);
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
