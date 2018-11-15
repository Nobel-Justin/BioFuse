package BioFuse::BioInfo::GeneAnno::PSL;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Array qw/ mergeOverlap /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              read_unit_region_from_PSL
              modify_pos_info
              extract_exon_seq_from_genome_and_output
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::PSL';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-11-15';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        read_unit_region_from_PSL
                        modify_pos_info
                        extract_exon_seq_from_genome_and_output
                        output_exon_seq
                     /;

#--- read GTF file ---
sub read_unit_region_from_PSL{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($Refseg_Exon_Info_Href,$psl_file,$avoid_refseg_Aref) = @_;
    my %parm = @_;
    my $Refseg_Exon_Info_Href = $parm{Refseg_Exon_Info_Href};
    my $psl_file = $parm{psl_file};
    my $avoid_refseg_Aref = $parm{avoid_refseg_Aref};

    # refseg user want to avoid
    my %avoid_refseg;
    @avoid_refseg{@$avoid_refseg_Aref} = () if($avoid_refseg_Aref);

    open (PSL, Try_GZ_Read($psl_file)) || die"fail $psl_file: $!\n";
    while(<PSL>){
        my ($strand,$unit,$refseg,$length,$positive_st) = (split)[8,9,13,18,20];
        next if(exists($avoid_refseg{$refseg})); # refseg ignored by user
        #--- change the pos to non-overlapped regions ---
        my @length = split /,/,$length;
        my @positive_st = split /,/,$positive_st;
        &modify_pos_info(exLen_Aref => \@length, stPos_Aref => \@positive_st);
        #-------- record the trans pos info --------
        @{$$Refseg_Exon_Info_Href{$refseg}{$unit}} = ($strand,\@length,\@positive_st);
    }
    close PSL;

    stout_and_sterr "[INFO]\tRead exon region from PSL ok.\n";
}

#--- modify the psl pos info ---
sub modify_pos_info{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($length_array,$positive_st_array) = @_;
    my %parm = @_;
    my $exLen_Aref = $parm{exLen_Aref};
    my $stPos_Aref = $parm{stPos_Aref};

    my @exonItv = map { [ $stPos_Aref->[$_] + 1,
                          $stPos_Aref->[$_] + $exLen_Aref->[$_] ]
                      } (0 .. scalar(@$exLen_Aref)-1);
    # merge overlap
    mergeOverlap(regionAref => \@exonItv, mergeAdjacent => 1);
    # update
    @$stPos_Aref = map { $_->[0] - 1 } @exonItv;
    @$exLen_Aref = map { $_->[1] - $_->[0] + 1 } @exonItv;
}

#--- output base ---
sub extract_exon_seq_from_genome_and_output{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    # my ($Refseg_Exon_Info_Href,$out,$whole_genome) = @_;
    my %parm = @_;
    my $Refseg_Exon_Info_Href = $parm{Refseg_Exon_Info_Href};
    my $out = $parm{out};
    my $whole_genome = $parm{whole_genome};

    open (my $outfh, Try_GZ_Write($out)) || die "fail $out: $!\n";
    read_fasta_file( FaFile => $whole_genome, needSeg_Href => $Refseg_Exon_Info_Href,
                     subrtRef => \&output_exon_seq,
                     subrtParmAref => [fh => $outfh, Refseg_Exon_Info_Href => $Refseg_Exon_Info_Href] );

    close $outfh;

    # check the remained unit
    for my $refseg (keys %$Refseg_Exon_Info_Href){
        for my $unit (sort keys %{$$Refseg_Exon_Info_Href{$refseg}}){
            stout_and_sterr "<WARN>\tUnit $unit from Refseg $refseg remains undetected.\n";
        }
    }

    stout_and_sterr "[INFO]\tCreate Exon_Seq file $out ok!\n";
}

#--- extract exon region from given chr-seq ---
sub output_exon_seq{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $segName = $parm{segName};
    my $segSeq_Sref = $parm{segSeq_Sref};
    my $fh = $parm{fh};
    my $Refseg_Exon_Info_Href = $parm{Refseg_Exon_Info_Href};

    my $line_base = 50;
    foreach my $unit (sort keys %{$$Refseg_Exon_Info_Href{$segName}}) {
        my ($strand,$exon_len_Aref,$exon_plus_st_Aref) = @{$$Refseg_Exon_Info_Href{$segName}{$unit}};
        my $unit_seq = '';
        for my $i (1 .. scalar(@$exon_len_Aref)){
            my $exon_len = $$exon_len_Aref[$i-1];
            my $exon_plus_st = $$exon_plus_st_Aref[$i-1] + 1;
            my $exon_seq = uc(substr($$segSeq_Sref,$exon_plus_st-1,$exon_len));
            if($strand eq '-'){ # reversed complementary
                ($exon_seq = reverse $exon_seq) =~ tr/ACGT/TGCA/;
                $unit_seq = $exon_seq.$unit_seq;
            }
            else{
                $unit_seq .= $exon_seq;
            }
        }
        print {$fh} ">$unit\n";
        print {$fh} substr($unit_seq,$_*$line_base,$line_base)."\n" for (0 .. int((length($unit_seq)-1)/$line_base));
        # once ok delete the unit
        delete $$Refseg_Exon_Info_Href{$segName}{$unit};
    }
}

1; ## tell the perl script the successful access of this module.
