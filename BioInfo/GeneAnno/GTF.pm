package BioFuse::BioInfo::GeneAnno::GTF;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::GeneAnno::GTF_geneOB;
use BioFuse::BioInfo::GeneAnno::GTF_lineOB;
use BioFuse::BioInfo::CytoBand qw/ load_cytoband /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              read_GTF
              create_trans_PSL
              create_gene_PSL
              mark_abnormal_Start_codon
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::GTF';
#----- version ---
$VERSION = "0.07";
$DATE = '2018-11-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#----- functions in this pm --------#
my @functoion_list = qw/
                        read_GTF
                        read_refseg_transform
                        load_GTF_gene
                        refine_GTF_info
                        refine_names
                        add_refseg_cytoband
                        create_gene_PSL
                        create_trans_PSL
                        mark_abnormal_Start_codon
                        check_Start_codon_Seq
                     /;

#--- read GTF file ---
sub read_GTF{
    # read refseg transforming list
    &read_refseg_transform(@_);

    # read gtf file to GTF_geneOB type
    &load_GTF_gene(@_);

    # refine the codon info
    &refine_GTF_info(@_);

    # refine the gene/transcript name
    &refine_names(@_);

    # add cytoBand info if set
    &add_refseg_cytoband(@_);
}

#--- read the file contains the refseg transform information ---
sub read_refseg_transform{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    return unless defined $V_Href->{refseg_transform_list};
    open (RT, Try_GZ_Read($V_Href->{refseg_transform_list})) || die "cannot read refseg_transform_list: $!\n";
    while(<RT>){
        chomp;
        my ($refseg_in_gtf,$refseg_for_use) = split /\t+/;
        $V_Href->{_refseg_transform}->{$refseg_in_gtf} = $refseg_for_use;
    }
    close RT;
    # inform
    stout_and_sterr "[INFO]\tread refseg transform list ok.\n"
                         ."\t$V_Href->{refseg_transform_list}\n";
}

#--- load the initial gtf info to GTF_geneOB type from the gtf file ---
sub load_GTF_gene{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    my $gtf_source = join(',', @{$V_Href->{gtf_source}});
    open (GTF, Try_GZ_Read($V_Href->{gtf})) || die "cannot read gtf_file: $!\n";
    while(<GTF>){
        next if(/^\#/); # comments, e.g., "#!genome-build GRCh38.p2" in v80
        chomp;
        my $gtfLineOB = BioFuse::BioInfo::GeneAnno::GTF_lineOB->new(gtf_lineinfo => $_);
        my $refSeg = $gtfLineOB->get_refSeg;
        # transform refSeg
        if(defined $V_Href->{refseg_transform_list}){
            if(exists $V_Href->{_refseg_transform}->{$refSeg}){
                $gtfLineOB->update_refSeg(refSeg => $V_Href->{_refseg_transform}->{$refSeg});
            }
            else{ # all unspecified refSeg will be skipped.
                next;
            }
        }
        # source selection
        if(    $gtf_source
            && !$gtfLineOB->match_gtf_source(gtf_source => $gtf_source)
        ){
            next;
        }
        # record info
        my ($gene_ENSid) = $gtfLineOB->get_ENSid(type => 'gene');
        if($gtfLineOB->get_regionType eq 'gene'){
            # to initialize GTF_geneOB
            $V_Href->{GTF_gene}->{$gene_ENSid} = BioFuse::BioInfo::GeneAnno::GTF_geneOB->new(gtfLineOB => $gtfLineOB);
        }
        else{
            # the transcript might come without no 'gene' info ahead, then warn. (before v75).
            unless(exists $V_Href->{GTF_gene}->{$gene_ENSid}){
                # make up the GTF_geneOB object
                $V_Href->{GTF_gene}->{$gene_ENSid} = BioFuse::BioInfo::GeneAnno::GTF_geneOB->new(gtfLineOB => $gtfLineOB);
                stout_and_sterr "<WARN>\tgtf file lacks the line of 'gene' information of $gene_ENSid. Reconstructed it.\n";
            }
            # load more info, must belong to one transcript
            $V_Href->{GTF_gene}->{$gene_ENSid}->load_gtf_info(gtfLineOB => $gtfLineOB);
        }
    }
    close GTF;
    # inform
    stout_and_sterr "[INFO]\tload initial gtf info OK.\n"
                         ."\t$V_Href->{gtf}\n";
}

#--- refine the codon number of transcript ---
sub refine_GTF_info{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    $_->refine_gtf_info for values %{$V_Href->{GTF_gene}};
    # inform
    stout_and_sterr "[INFO]\trefine the gtf info ok.\n";
}

#--- refine duplicated use_names of gene and transcript ---
sub refine_names{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    # get name count
    my (%_genename_count,%_transname_count);
    for my $gene (values %{$V_Href->{GTF_gene}}){
        $_genename_count{$gene->get_gene_use_name}++;
        $_transname_count{$_}++ for @{$gene->get_trans_use_name};
    }
    # refine duplicated gene/trans use_name
    for my $gene (values %{$V_Href->{GTF_gene}}){
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
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    return unless defined $V_Href->{cytoBand_file};
    load_cytoband(cytoBand_Href => $V_Href->{cytoBand}, cytoBand_file => $V_Href->{cytoBand_file});
    $_->add_cytoband_info(cytoBand_Href => $V_Href->{cytoBand}) for values %{$V_Href->{GTF_gene}};
    # inform
    stout_and_sterr "[INFO]\tadd the refseg cytoband ok.\n";
}

#--- creat gene psl file ---
sub create_gene_PSL{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    open (GPSL, Try_GZ_Write($V_Href->{psl})) || die "fail create PSL_file: $!\n";
    for my $_gene_ENSid (sort keys %{$V_Href->{GTF_gene}}){
        print GPSL $V_Href->{GTF_gene}->{$_gene_ENSid}->get_gene_psl_line . "\n";
    }
    close GPSL;
    # inform
    stout_and_sterr "[INFO]\twrite gene PSL ok!\n"
                         ."\t$V_Href->{psl}\n";
}

#--- creat trans psl file ---
sub create_trans_PSL{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $show_metainfo = $parm{show_metainfo} || 0;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    open (TPSL, Try_GZ_Write($V_Href->{psl})) || die "fail create PSL_file: $!\n";
    for my $_gene_ENSid (sort keys %{$V_Href->{GTF_gene}}){
        my $gene = $V_Href->{GTF_gene}->{$_gene_ENSid};
        my $gene_use_name = $gene->get_gene_use_name;
        print TPSL $_->get_trans_psl_line(show_metainfo => $show_metainfo, gene_use_name => $gene_use_name) . "\n" for @{$gene->get_transOB_Aref};
    }
    close TPSL;
    # inform
    stout_and_sterr "[INFO]\twrite trans PSL ok!\n"
                         ."\t$V_Href->{psl}\n";
}

#--- refine the start codon seq according to the genome ---
sub mark_abnormal_Start_codon{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $V_Href = exists $parm{V_Href} ? $parm{V_Href} : BioFuse::LoadOn->load_variants_dict;

    # standard start codon
    my %Start_codon = map {(uc($_),1)} @{$V_Href->{Start_codon}};
    stout_and_sterr "[INFO]\tStandard start codon sequence: ".join(',',sort keys %Start_codon)."\n";
    # select the 'protein_coding' transcript
    my %transToCheck;
    for my $gene (values %{$V_Href->{GTF_gene}}){
        push @{$transToCheck{$gene->get_ref_seg}}, $_ for grep $_->get_biotype eq 'protein_coding', @{$gene->get_transOB_Aref};
    }
    # get start codon sequence and check
    read_fasta_file( FaFile => $V_Href->{whole_genome}, needSeg_Href => \%transToCheck,
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
