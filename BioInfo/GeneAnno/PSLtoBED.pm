package BioFuse::BioInfo::GeneAnno::PSLtoBED;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ max /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Write /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::GeneAnno::PSL qw/ load_GeneOrTrans_from_PSL /;
use BioFuse::Util::Array qw/ merge /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              PSLtoBED
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::PSLtoBED';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        PSLtoBED
                        output_region_of_transOB
                        record_protein_coding_regions
                        record_exon_regions
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} pslToBed <[Options]>

     Options:

        # Input and Output #
         -tp [s]  transcript PSL file. <required>
         -o  [s]  output bed file. <required>

        # Options #
         -pz [i]  promoter size extended to upstream from the first exon. [1000]
         -tz [i]  terminator size extended to downstream from the last exon. [1000]
         -ab [s]  the refseg you want to avoid, allow mutilple input. [NULL]
         -gn [s]  the gene name you want to keep only, allow mutilple input. [NULL]
         -bt [s]  the biotype you want to keep only, allow mutilple input. [NULL]
         -rt [s]  the region type you want to keep only, allow mutilple input. [NULL]
                   NOTE: region types: CDS/UTR5/UTR3/EXON/INTRON/PROMOTER/TERMINATOR
         -pm [s]  mode to select 'protein_coding' transcript for CDS/UTR5/UTR3. [trNO]
                   NOTE: available mode: trNO, longest
         -em [s]  mode to select transcript for EXON/INTRON/PROMOTER/TERMINATOR. [protein_trNO]
                   NOTE: available mode: protein_trNO, trNO, longest, merge
         -one     output region is started from one, not BED format. [disabled]
         -abp     include abnormal protein_coding transcript. [disabled]

         -h       show this help

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
    $V_Href->{ $_->[0] } = $_->[1] for
        map {
            if( !exists $V_Href->{$_->[0]} ){
                ( $_ );
            }
            else{
                warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            ## use 'psl' in BioFuse::LoadOn
            [ tpsl => undef ],
            [ bed_file => undef ],
            # option
            [ promoter_size => 1000 ],
            [ terminator_size => 1000 ],
            [ abandon_refseg => [] ],
            [ require_genename => [] ],
            [ require_biotype => [] ],
            [ require_regtype => [] ],
            [ prcd_mode => 'trNO' ],
            [ exon_mode => 'protein_trNO' ],
            [ include_ab_protein => 0 ],
            [ oneBase => 0 ], # one-start region, not BED

            # intermediate variants
            [ GeneToTransOB => {} ],
            [ Regions_Hf => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['tpsl'],
                                  ['bed_file']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-tp:s" => \$V_Href->{tpsl},
        "-o:s"  => \$V_Href->{bed_file},
        # options
        "-pz:i" => \$V_Href->{promoter_size},
        "-tz:i" => \$V_Href->{terminator_size},
        "-ab:s" => \@{$V_Href->{abandon_refseg}},
        "-gn:s" => \@{$V_Href->{require_genename}},
        "-bt:s" => \@{$V_Href->{require_biotype}},
        "-rt:s" => \@{$V_Href->{require_regtype}},
        "-pm:s" => \$V_Href->{prcd_mode},
        "-abp"  => \$V_Href->{include_ab_protein},
        "-em:s" => \$V_Href->{exon_mode},
        "-one"  => \$V_Href->{oneBase},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{tpsl})
             || $V_Href->{promoter_size} <= 0
             || $V_Href->{terminator_size} <= 0
             || !defined $V_Href->{bed_file}
            );
}

#--- get genetic regions based on PSL file ---
sub PSLtoBED{
    # load trans psl file
    load_GeneOrTrans_from_PSL( ObjectPoolHref => $V_Href->{GeneToTransOB},
                               psl_file => $V_Href->{tpsl},
                               psl_type => 'trans',
                               recordKey => 'genename',
                               skipRefSegHref => { map {($_,1)} @{$V_Href->{abandon_refseg}} },
                               requireBtyHref => { map {($_,1)} @{$V_Href->{require_biotype}} }
                             );
    # output region
    &output_region_of_transOB;
}

#--- output series of regions of transcript ---
sub output_region_of_transOB{

    # genename filter
    if(scalar @{$V_Href->{require_genename}}){
        my %require_genename = map {($_,1)} @{$V_Href->{require_genename}};
        delete $V_Href->{GeneToTransOB}->{$_} for grep ! exists $require_genename{$_}, keys %{$V_Href->{GeneToTransOB}};
    }

    # record series of regions
    for my $transOB_Af (values %{$V_Href->{GeneToTransOB}}){
        # protein_coding? CDS/UTR
        my @transOB_prcd =   $V_Href->{include_ab_protein}
                           ? (grep $_->get_biotype =~ /^protein_coding/, @$transOB_Af)
                           : (grep $_->get_biotype eq  'protein_coding', @$transOB_Af);
        &record_protein_coding_regions(transOB_prcd_Af => \@transOB_prcd) if scalar @transOB_prcd;
        # EXON (protein, trNO, longest, merge)
        &record_exon_regions(transOB_Af => $transOB_Af, transOB_prcd_Af => \@transOB_prcd);
    }

    # region filter
    if(scalar @{$V_Href->{require_regtype}}){
        my %require_regtype = map {(uc($_),1)} @{$V_Href->{require_regtype}};
        delete $V_Href->{Regions_Hf}->{$_} for grep ! exists $require_regtype{$_}, keys %{$V_Href->{Regions_Hf}};
    }

    # output region
    my $offset = $V_Href->{oneBase} ? 0 : 1;
    open (BED, Try_GZ_Write($V_Href->{bed_file})) || die "fail write output region file: $!\n";
    for my $region_type (sort keys %{$V_Href->{Regions_Hf}}){
        my $regType_Hf = $V_Href->{Regions_Hf}->{$region_type};
        for my $refseg (sort keys %$regType_Hf){
            # merge
            $regType_Hf->{$refseg} = merge(itvAfListAf => [$regType_Hf->{$refseg}], mergeAdjacent => 1);
            print BED join("\t", $refseg, $_->[0]-$offset, $_->[1], $region_type)."\n" for @{$regType_Hf->{$refseg}};
            # inform
            stout_and_sterr "[INFO]\toutput $region_type $refseg region OK.\n";
        }
    }
    close BED;
}

#--- record protein coding relevant regions ---
## CDS, UTR3, UTR5
sub record_protein_coding_regions{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $transOB_prcd_Af = $parm{transOB_prcd_Af};

    # select protein_coding transcript for CDS region
    if($V_Href->{prcd_mode} eq 'trNO'){
        @$transOB_prcd_Af = sort {$a->get_ori_name cmp $b->get_ori_name} @$transOB_prcd_Af;
    }
    elsif($V_Href->{prcd_mode} eq 'longest'){
        @$transOB_prcd_Af = sort {$b->get_exon_sumLen <=> $a->get_exon_sumLen} @$transOB_prcd_Af;
    }

    # record
    my $refseg = $transOB_prcd_Af->[0]->get_ref_seg;
    push @{$V_Href->{Regions_Hf}->{CDS}->{$refseg}},  [ @$_ ] for @{$transOB_prcd_Af->[0]->get_CDS_region};
    push @{$V_Href->{Regions_Hf}->{UTR5}->{$refseg}}, [ @$_ ] for @{$transOB_prcd_Af->[0]->get_UTR_region(p=>5)};
    push @{$V_Href->{Regions_Hf}->{UTR3}->{$refseg}}, [ @$_ ] for @{$transOB_prcd_Af->[0]->get_UTR_region(p=>3)};
}

#--- record regions derived from exon ---
## EXON, INTRON, PROMOTER, TERMINATOR
sub record_exon_regions{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $transOB_Af = $parm{transOB_Af};
    my $transOB_prcd_Af = $parm{transOB_prcd_Af};

    # select transcript for EXON region
    my @transOB_forExon;
    if(    $V_Href->{exon_mode} eq 'protein_trNO'
        && scalar @$transOB_prcd_Af
    ){
        @transOB_forExon = ($transOB_prcd_Af->[0]);
    }
    elsif(   $V_Href->{exon_mode} eq 'trNO'
          || (   $V_Href->{exon_mode} eq 'protein_trNO'
              && scalar @transOB_forExon == 0
             )
    ){
        @transOB_forExon = (sort {$a->get_ori_name cmp $b->get_ori_name} @$transOB_Af)[0];
    }
    elsif($V_Href->{exon_mode} eq 'longest'){
        @transOB_forExon = (sort {$b->get_exon_sumLen <=> $a->get_exon_sumLen} @$transOB_Af)[0];
    }
    elsif($V_Href->{exon_mode} eq 'merge'){
        @transOB_forExon = @$transOB_Af;
    }

    # get exon region
    my @EXON;
    for my $transOB (@transOB_forExon){
        push @EXON, [ @$_ ] for @{$transOB->get_exon_region};
    }
    ## merge exon regions (for the 'merge' mode)
    ## after this, regions are sorted in ascending order
    @EXON = @{ merge(itvAfListAf => [\@EXON], mergeAdjacent => 1) };

    my $refseg = $transOB_Af->[0]->get_ref_seg;
    # record EXON
    push @{$V_Href->{Regions_Hf}->{EXON}->{$refseg}},  [ @$_ ] for @EXON;
    # record INTRON
    push @{$V_Href->{Regions_Hf}->{INTRON}->{$refseg}}, [ $EXON[$_]->[1]+1, $EXON[$_+1]->[0]-1 ] for (0 .. $#EXON-1);
    # record PROMOTER and TERMINATOR
    if($transOB_Af->[0]->get_strand eq '+'){
        push @{$V_Href->{Regions_Hf}->{PROMOTER}->{$refseg}},   [ max($EXON[0]->[0]-$V_Href->{promoter_size}, 1), max($EXON[0]->[0]-1, 1) ];
        push @{$V_Href->{Regions_Hf}->{TERMINATOR}->{$refseg}}, [ $EXON[-1]->[1]+1, $EXON[-1]->[1]+$V_Href->{terminator_size} ];
    }
    else{
        push @{$V_Href->{Regions_Hf}->{PROMOTER}->{$refseg}},   [ $EXON[-1]->[1]+1, $EXON[-1]->[1]+$V_Href->{promoter_size} ];
        push @{$V_Href->{Regions_Hf}->{TERMINATOR}->{$refseg}}, [ max($EXON[0]->[0]-$V_Href->{terminator_size}, 1), max($EXON[0]->[0]-1, 1) ];
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
