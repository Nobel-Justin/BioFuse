package BioFuse::BioInfo::GeneAnno::GTF;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::BioInfo::GeneAnno::GTF_geneOB;
use BioFuse::BioInfo::CytoBand qw/ load_cytoband /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              read_GTF
              add_refseg_cytoband
              create_trans_PSL
              create_gene_PSL
              read_refseg_transform
              load_GTF_Info_type_from_gtf
              refine_GTF_info
              refine_names
              mark_abnormal_Start_codon
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::GTF';
#----- version ---
$VERSION = "0.06";
$DATE = '2018-01-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#----- functions in this pm --------#
my @functoion_list = qw/
                        read_GTF
                        add_refseg_cytoband
                        create_trans_PSL
                        create_gene_PSL
                        read_refseg_transform
                        load_GTF_Info_type_from_gtf
                        refine_GTF_info
                        refine_names
                        mark_abnormal_Start_codon
                        check_Start_codon_Seq
                     /;

#--- read GTF file ---
sub read_GTF{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href,$gtf_file,$refseg_transform_list,$gtf_gene_source_Aref) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};
    my $gtf_file = $parm{gtf_file};
    my $refseg_transform_list = $parm{refseg_transform_list};
    my $gtf_gene_source_Aref = $parm{gtf_gene_source_Aref};

    # read refseg transforming list
    my %_refseg_transform;
    if($refseg_transform_list){
        &read_refseg_transform(refseg_transform_list=>$refseg_transform_list, _refseg_transform_Href=>\%_refseg_transform);
    }

    # gtf source selection
    my $_gtf_gene_source_string = join(',', @$gtf_gene_source_Aref);
    my $_gtf_gene_source_bool_or_Sref = (scalar(@$gtf_gene_source_Aref)==0) ? 0 : \$_gtf_gene_source_string;

    # read gtf file to GTF_geneOB type
    &load_GTF_Info_type_from_gtf(GTF_Info_Href => $GTF_Info_Href,
                                 gtf_file => $gtf_file,
                                 _refseg_transform_Href => \%_refseg_transform,
                                 user_set_transform_list_bool => $refseg_transform_list,
                                 gtf_gene_source_Sref => $_gtf_gene_source_bool_or_Sref );

    # refine the codon info
    &refine_GTF_info(GTF_Info_Href => $GTF_Info_Href);

    # refine the gene/transcript name
    &refine_names(GTF_Info_Href => $GTF_Info_Href);
}

#--- add the refseg cytoband info ---
# optional
sub add_refseg_cytoband{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href,$cytoBand_Href) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};
    my $cytoBand_Href = $parm{cytoBand_Href};

    for my $_gene_ENSid (keys %$GTF_Info_Href){
        $GTF_Info_Href->{$_gene_ENSid}->add_cytoband_info(cytoBand_Href => $cytoBand_Href);
    }

    stout_and_sterr "[INFO]\tadd the refseg cytoband ok.\n";
}

#--- creat trans psl file ---
sub create_trans_PSL{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href,$PSL_file,$show_metainfo) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};
    my $PSL_file = $parm{PSL_file};
    my $show_metainfo = $parm{show_metainfo} || 0;

    open (my $PSL_fh, Try_GZ_Write($PSL_file)) || die "fail create $PSL_file: $!\n";
    for my $_gene_ENSid (sort keys %$GTF_Info_Href){
        $GTF_Info_Href->{$_gene_ENSid}->write_trans_psl_info(psl_fh => $PSL_fh, show_metainfo => $show_metainfo);
    }
    close $PSL_fh;

    stout_and_sterr "[INFO]\twrite trans_PSL ok!\n";
}

#--- creat gene psl file ---
sub create_gene_PSL{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href,$PSL_file) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};
    my $PSL_file = $parm{PSL_file};

    open (my $PSL_fh, Try_GZ_Write($PSL_file)) || die "fail create $PSL_file: $!\n";
    for my $_gene_ENSid (sort keys %$GTF_Info_Href){
        $GTF_Info_Href->{$_gene_ENSid}->write_gene_psl_info(psl_fh => $PSL_fh);
    }
    close $PSL_fh;

    stout_and_sterr "[INFO]\twrite gene_PSL ok!\n";
}

#--- read the file contains the refseg transform information ---
sub read_refseg_transform{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($refseg_transform_list,$_refseg_transform_Href) = @_;
    my %parm = @_;
    my $refseg_transform_list = $parm{refseg_transform_list};
    my $_refseg_transform_Href = $parm{_refseg_transform_Href};

    open (RT, Try_GZ_Read($refseg_transform_list)) || die "cannot read $refseg_transform_list: $!\n";
    while(<RT>){
        chomp;
        my ($refseg_in_gtf,$refseg_for_use) = split /\t+/;
        $$_refseg_transform_Href{$refseg_in_gtf} = $refseg_for_use;
    }
    close RT;

    stout_and_sterr "[INFO]\tread refseg transform list ok.\n";
}

#--- load the initial gtf info to GTF_geneOB type from the gtf file ---
sub load_GTF_Info_type_from_gtf{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href, $gtf_file, $_refseg_transform_Href, $user_set_transform_list_bool, $gtf_gene_source_Sref) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};
    my $gtf_file = $parm{gtf_file};
    my $_refseg_transform_Href = $parm{_refseg_transform_Href};
    my $user_set_transform_list_bool = $parm{user_set_transform_list_bool};
    my $gtf_gene_source_Sref = $parm{gtf_gene_source_Sref};

    open (my $GTF_fh, Try_GZ_Read($gtf_file)) || die "cannot read $gtf_file: $!\n";
    while(<$GTF_fh>){
        next if(/^\#/); # commend, such as "#!genome-build GRCh38.p2" in Ensemble Released v80

        my ($ref_seg, $region_type) = (split)[0,2];

        # when the refseg trans_list is set, all unspecified ref_seg will be skipped.
        next if($user_set_transform_list_bool && !exists($_refseg_transform_Href->{$ref_seg}));
        $ref_seg = $_refseg_transform_Href->{$ref_seg} || $ref_seg;

        # source selection
        if($gtf_gene_source_Sref){ # $gtf_gene_source_Sref may be 0 when users do not set any seletion.
            my ($gene_source) = (/gene_source\s\"([^\"]+)\"/);
            next if( $gene_source && $$gtf_gene_source_Sref !~ /$gene_source/ ); # gene_source may not exist in ensembl before v75
        }

        # create a GTF_geneOB type
        my ($gene_id) = (/gene_id\s\"([^\"]+)\"/);
        if($region_type eq 'gene'){ # firstly constrcut the gene_id key as the GTF_geneOB object
            $GTF_Info_Href->{$gene_id} = {}; # initialize
            $GTF_Info_Href->{$gene_id} = BioFuse::BioInfo::GeneAnno::GTF_geneOB->new(ref_seg => $ref_seg, gtf_lineinfo => $_);
        }
        else{
            if( exists($GTF_Info_Href->{$gene_id}) ){
                $GTF_Info_Href->{$gene_id}->load_gtf_info(gtf_lineinfo => $_);
            }
            else{ # if the transcript comes, but no its gene_id GTF_geneOB object found, it will be warned. (version before v75).
                warn ("The gtf source lacks the line of 'gene' information of $gene_id. Reconstructed it.\n");
                # reconstructed GTF_geneOB object
                $GTF_Info_Href->{$gene_id} = {}; # initialize
                $GTF_Info_Href->{$gene_id} = BioFuse::BioInfo::GeneAnno::GTF_geneOB->new(ref_seg => $ref_seg, gtf_lineinfo => $_);
                # load more info
                $GTF_Info_Href->{$gene_id}->load_gtf_info(gtf_lineinfo => $_);
            }
        }
    }
    close $GTF_fh;

    stout_and_sterr "[INFO]\tload initial gtf info ok.\n";
}

#--- refine the codon number of transcript ---
sub refine_GTF_info{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};

    for my $_gene_ENSid (keys %$GTF_Info_Href){
        $GTF_Info_Href->{$_gene_ENSid}->refine_gtf_info;
    }

    stout_and_sterr "[INFO]\trefine the gtf info ok.\n";
}

#--- refine the name of gene and transcript ---
sub refine_names{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};

    my (%_genename_count,%_transname_count);
    for my $_gene_ENSid (keys %$GTF_Info_Href){
        my $gene = $GTF_Info_Href->{$_gene_ENSid};
        my $_gene_use_name = $gene->get_gene_use_name;
        $_genename_count{$_gene_use_name}++;
        my @_trans_use_name = $gene->get_trans_use_name;
        $_transname_count{$_}++ for @_trans_use_name;
    }

    for my $_gene_ENSid (keys %$GTF_Info_Href){
        my $gene = $GTF_Info_Href->{$_gene_ENSid};
        my $_gene_use_name = $gene->get_gene_use_name;
        if($_genename_count{$_gene_use_name} > 1){ # this gene use name appears more than 1 time
            $gene->{gene_use_name} = $gene->{gene_ori_name} . '_' . $_gene_ENSid;
        }
        for my $_trans_ENSid (keys %{$gene->{trans_info}}){
            my $trans = $gene->{trans_info}->{$_trans_ENSid};
            my $_trans_use_name = $trans->{trans_use_name};
            if($_transname_count{$_trans_use_name} > 1){ # this trans use name appears more than 1 time
                $trans->{trans_use_name} = $trans->{trans_ori_name} . '_' . $_trans_ENSid;
            }
        }
    }

    stout_and_sterr "[INFO]\trefine gene/trans names ok.\n";
}

#--- refine the start codon seq according to the genome ---
sub mark_abnormal_Start_codon{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href,$Start_codon_Aref,$whole_genome) = @_;
    my %parm = @_;
    my $GTF_Info_Href = $parm{GTF_Info_Href};
    my $Start_codon_Aref = $parm{Start_codon_Aref};
    my $whole_genome = $parm{whole_genome};

    stout_and_sterr "[INFO]\tcheck strat codon sequence now...\n";

    my $Start_codon_Href = {};
    %$Start_codon_Href= map {(uc($_),1)} @$Start_codon_Aref;
    stout_and_sterr "[INFO]\tStandard strat codon sequence: ".join(',',sort keys %$Start_codon_Href)."\n";

    # select the protein coding transcript
    my $Refseg_to_StartCodon_Href = {};
    for my $_gene_ENSid (keys %$GTF_Info_Href){
        my $gene_Href = $GTF_Info_Href->{$_gene_ENSid};
        my $ref_seg = $gene_Href->{ref_seg};
        for my $_trans_ENSid (keys %{$gene_Href->{trans_info}}){
            # only deal protein_coding
            my $trans_biotype = $gene_Href->{trans_info}->{$_trans_ENSid}->{trans_biotype};
            if($trans_biotype eq 'protein_coding'){
                $Refseg_to_StartCodon_Href->{$ref_seg}->{$_gene_ENSid} = 1; # store the gene Href to save the memory
                last; # find one protein_coding, no need to next loop
            }
        }
    }

    # get start codon sequence and check
    read_fasta_file( FaFile => $whole_genome, needSeg_Href => $Refseg_to_StartCodon_Href,
                     subrtRef => \&check_Start_codon_Seq,
                     subrtParmAref => [GTF_Info_Href => $GTF_Info_Href, Refseg_to_StartCodon_Href => $Refseg_to_StartCodon_Href, Start_codon_Href => $Start_codon_Href] );

    stout_and_sterr "[INFO]\tcheck start codon seq ok!\n";
}

#--- check start codon from given chr-seq ---
sub check_Start_codon_Seq{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($GTF_Info_Href,$Start_codon_Aref,$whole_genome) = @_;
    my %parm = @_;
    my $segName = $parm{segName};
    my $segSeq_Sref = $parm{segSeq_Sref};
    my $GTF_Info_Href = $parm{GTF_Info_Href};
    my $Refseg_to_StartCodon_Href = $parm{Refseg_to_StartCodon_Href};
    my $Start_codon_Href = $parm{Start_codon_Href};

    foreach my $_gene_ENSid (keys %{$Refseg_to_StartCodon_Href->{$segName}}) {
        my $gene_Href = $GTF_Info_Href->{$_gene_ENSid};
        my $strand = $gene_Href->{strand};
        for my $_trans_ENSid (keys %{$gene_Href->{trans_info}}){
            my $trans_Href = $gene_Href->{trans_info}->{$_trans_ENSid};
            # only deal protein_coding
            my $trans_biotype = $trans_Href->{trans_biotype};
            if($trans_biotype eq 'protein_coding'){
                my $StartCodon_Seq = '';
                for my $StartCodon_Pos (sort {$a<=>$b} keys %{$trans_Href->{start_codon}}){
                    $StartCodon_Seq .= uc(substr($$segSeq_Sref,$StartCodon_Pos-1,1));
                }
                # for minus strand gene
                ($StartCodon_Seq = reverse $StartCodon_Seq) =~ tr/ACGT/TGCA/ if($strand eq '-');
                # judge the startcodon seq
                unless(exists($Start_codon_Href->{$StartCodon_Seq})){
                    $trans_Href->{trans_biotype} = "protein_coding_with_abnor_st_codon_seq($StartCodon_Seq)";
                    stout_and_sterr "[INFO]\t$trans_Href->{ENSid} $trans_Href->{trans_ori_name} has abnormal start codon sequence: $StartCodon_Seq\n";
                }
            }
        }
    }
}

1; ## tell the perl script the successful access of this module.
