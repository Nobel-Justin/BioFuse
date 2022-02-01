package BioFuse::BioInfo::Objects::GeneAnno::GTF_OB;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::BioInfo::Objects::GeneAnno::GTF_geneOB;
use BioFuse::BioInfo::Objects::GeneAnno::GTF_lineOB;
use BioFuse::BioInfo::CytoBand qw/ load_cytoband /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;

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

$MODULE_NAME = 'BioFuse::BioInfo::Objects::GeneAnno::GTF_OB';
#----- version ---
$VERSION = "0.08";
$DATE = '2019-05-04';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#----- functions in this pm --------#
my @function_list = qw/
                        new
                        load_GTF
                        read_refseg_transform
                        refine_GTF_info
                        refine_names
                        add_refseg_cytoband
                        create_gene_PSL
                        create_trans_PSL
                        mark_abnormal_Start_codon
                        check_Start_codon_Seq
                     /;

#--- structure of object
# gtf -> filePath = $filePath
# gtf -> GTF_gene = {}
# gtf -> _refseg_transform = {}
# gtf -> gtf_source = "$s1,$s2,$s3,..,$sX"

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $gtf = {};
    # load basic attributes of gtf
    $gtf->{filePath} = $parm{filePath};
    $gtf->{GTF_gene} = {}; # main container

    bless($gtf);
    return $gtf;
}

#--- load the initial gtf info to GTF_geneOB type from the gtf file ---
## parm -> refseg_transform_list
## parm -> gtf_source (array-ref)
## parm -> cytoBand_file
sub load_GTF{
    my $gtf = shift;
    my %parm = @_;

    # read refseg transforming list
    if(defined $parm{refseg_transform_list}){
        $gtf->read_refseg_transform(@_);
    }
    else{
        undef $gtf->{_refseg_transform};
    }
    # gtf source selection
    if(defined $parm{gtf_source}){
        $gtf->{gtf_source} = join(',', @{$parm{gtf_source}});
    }
    else{
        undef $gtf->{gtf_source};
    }

    # read gtf file to GTF_geneOB type
    open (GTF, Try_GZ_Read($gtf->{filePath})) || die "cannot read gtf_file: $!\n";
    while(<GTF>){
        next if(/^\#/); # comments, e.g., "#!genome-build GRCh38.p2" in v80
        chomp;
        my $gtfLineOB = BioFuse::BioInfo::Objects::GeneAnno::GTF_lineOB->new(gtf_lineinfo => $_);
        my $refSeg = $gtfLineOB->get_refSeg;
        # transform refSeg
        if(defined $gtf->{_refseg_transform}){
            if(exists $gtf->{_refseg_transform}->{$refSeg}){
                $gtfLineOB->update_refSeg(refSeg => $gtf->{_refseg_transform}->{$refSeg});
            }
            else{ # all unspecified refSeg will be skipped.
                next;
            }
        }
        # source selection
        if(    $gtf->{gtf_source}
            && !$gtfLineOB->match_gtf_source(gtf_source => $gtf->{gtf_source})
        ){
            next;
        }
        # record info
        my ($gene_ENSid) = $gtfLineOB->get_ENSid(type => 'gene');
        if($gtfLineOB->get_regionType eq 'gene'){
            # to initialize GTF_geneOB
            $gtf->{GTF_gene}->{$gene_ENSid} = BioFuse::BioInfo::Objects::GeneAnno::GTF_geneOB->new(gtfLineOB => $gtfLineOB);
        }
        else{
            # the transcript might come without no 'gene' info ahead, then warn. (before v75).
            unless(exists $gtf->{GTF_gene}->{$gene_ENSid}){
                # make up the GTF_geneOB object
                $gtf->{GTF_gene}->{$gene_ENSid} = BioFuse::BioInfo::Objects::GeneAnno::GTF_geneOB->new(gtfLineOB => $gtfLineOB);
                stout_and_sterr "<WARN>\tgtf file lacks the line of 'gene' information of $gene_ENSid. Reconstructed it.\n";
            }
            # load more info, must belong to one transcript
            $gtf->{GTF_gene}->{$gene_ENSid}->load_gtf_info(gtfLineOB => $gtfLineOB);
        }
    }
    close GTF;
    # inform
    stout_and_sterr "[INFO]\tload gtf OK.\n"
                         ."\t$gtf->{filePath}\n";

    # refine the codon info
    $gtf->refine_GTF_info;
    # refine the gene/transcript name
    $gtf->refine_names;
    # add cytoBand info if set
    if(defined $parm{cytoBand_file}){
        load_cytoband(cytoBand_Href => $gtf->{cytoBand}, cytoBand_file => $parm{cytoBand_file});
        $gtf->add_refseg_cytoband;
    }
}

#--- read the file contains the refseg transform information ---
sub read_refseg_transform{
    my $gtf = shift;
    my %parm = @_;
    my $refseg_transform_list = $parm{refseg_transform_list};

    open (RT, Try_GZ_Read($refseg_transform_list)) || die "cannot read refseg_transform_list: $!\n";
    while(<RT>){
        chomp;
        my ($refseg_in_gtf,$refseg_for_use) = split /\t+/;
        $gtf->{_refseg_transform}->{$refseg_in_gtf} = $refseg_for_use;
    }
    close RT;
    # inform
    stout_and_sterr "[INFO]\tread refseg transform list ok.\n"
                         ."\t$refseg_transform_list\n";
}

#--- refine the codon number of transcript ---
sub refine_GTF_info{
    my $gtf = shift;
    $_->refine_gtf_info for values %{$gtf->{GTF_gene}};
    # inform
    stout_and_sterr "[INFO]\trefine gtf data ok.\n";
}

#--- refine duplicated use_names of gene and transcript ---
sub refine_names{
    my $gtf = shift;
    # get name count
    my (%_genename_count,%_transname_count);
    for my $gene (values %{$gtf->{GTF_gene}}){
        $_genename_count{$gene->get_gene_use_name}++;
        $_transname_count{$_}++ for @{$gene->get_trans_use_name};
    }
    # refine duplicated gene/trans use_name
    for my $gene (values %{$gtf->{GTF_gene}}){
        # this gene use name appears more than 1 time
        $gene->refined_use_name if $_genename_count{$gene->get_gene_use_name} > 1;
        for my $trans (@{$gene->get_transOB_Aref}){
            # this trans use name appears more than 1 time
            $trans->refined_use_name if $_transname_count{$trans->get_trans_use_name} > 1;
        }
    }
    # inform
    stout_and_sterr "[INFO]\trefine gene/trans names ok.\n";
}

#--- add the refseg cytoband info ---
sub add_refseg_cytoband{
    my $gtf = shift;
    $_->add_cytoband_info(cytoBand_Href => $gtf->{cytoBand}) for values %{$gtf->{GTF_gene}};
    # inform
    stout_and_sterr "[INFO]\tadd the refseg cytoband ok.\n";
}

#--- creat gene psl file ---
sub create_gene_PSL{
    my $gtf = shift;
    my %parm = @_;
    my $gPSLpath = $parm{gPSLpath};
    # output gene psl
    open (GPSL, Try_GZ_Write($gPSLpath)) || die "fail create PSL_file: $!\n";
    for my $_gene_ENSid (sort keys %{$gtf->{GTF_gene}}){
        print GPSL $gtf->{GTF_gene}->{$_gene_ENSid}->get_gene_psl_line . "\n";
    }
    close GPSL;
    # inform
    stout_and_sterr "[INFO]\twrite gene PSL ok!\n"
                         ."\t$gPSLpath\n";
}

#--- creat trans psl file ---
sub create_trans_PSL{
    my $gtf = shift;
    my %parm = @_;
    my $tPSLpath = $parm{tPSLpath};
    my $show_metainfo = $parm{show_metainfo} || 0;

    # check and mark start codon
    if(defined $parm{Start_codon}){
        $gtf->mark_abnormal_Start_codon(@_);
    }
    # output trans psl
    open (TPSL, Try_GZ_Write($tPSLpath)) || die "fail create PSL_file: $!\n";
    for my $_gene_ENSid (sort keys %{$gtf->{GTF_gene}}){
        my $gene = $gtf->{GTF_gene}->{$_gene_ENSid};
        my $gene_use_name = $gene->get_gene_use_name;
        print TPSL $_->get_trans_psl_line(show_metainfo => $show_metainfo, gene_use_name => $gene_use_name) . "\n" for @{$gene->get_transOB_Aref};
    }
    close TPSL;
    # inform
    stout_and_sterr "[INFO]\twrite trans PSL ok!\n"
                         ."\t$tPSLpath\n";
}

#--- refine the start codon seq according to the genome ---
sub mark_abnormal_Start_codon{
    my $gtf = shift;
    my %parm = @_;

    # standard start codon allowed
    my %Start_codon = map {(uc($_),1)} @{$parm{Start_codon}};
    stout_and_sterr "[INFO]\tStandard start codon sequence: ".join(',',sort keys %Start_codon)."\n";
    # select the 'protein_coding' transcript
    my %transToCheck;
    for my $gene (values %{$gtf->{GTF_gene}}){
        push @{$transToCheck{$gene->get_ref_seg}}, $_ for grep $_->get_biotype eq 'protein_coding', @{$gene->get_transOB_Aref};
    }
    # get start codon sequence and check
    read_fasta_file( FaFile => $parm{whole_genome}, needSeg_Href => \%transToCheck,
                     subrtRef => \&check_Start_codon_Seq,
                     subrtParmAref => [transToCheck_Href => \%transToCheck, Start_codon_Href => \%Start_codon] );
    # inform
    stout_and_sterr "[INFO]\tcheck start codon seq ok!\n";
}

#--- check start codon from given chr-seq ---
sub check_Start_codon_Seq{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $segName = $parm{segName};
    my $segSeq_Sref = $parm{segSeq_Sref};
    my $transToCheck_Href = $parm{transToCheck_Href};
    my $Start_codon_Href = $parm{Start_codon_Href};

    for my $trans (@{$transToCheck_Href->{$segName}}){
        my $StartCodon_Seq = join( '', map {uc(substr($$segSeq_Sref,$_-1,1))} @{$trans->get_startCodonPos} );
        # for minus strand gene
        ($StartCodon_Seq = reverse $StartCodon_Seq) =~ tr/ACGT/TGCA/ if($trans->get_strand eq '-');
        # judge the start codon seq
        unless(exists $Start_codon_Href->{$StartCodon_Seq}){
            $trans->update_biotype(biotype => "protein_coding_with_abnor_st_codon_seq($StartCodon_Seq)");
            stout_and_sterr "[INFO]\t".$trans->get_ENSid.' '.$trans->get_trans_ori_name." has abnormal start codon sequence: $StartCodon_Seq\n";
        }
    }
}

1; ## tell the perl script the successful access of this module.
