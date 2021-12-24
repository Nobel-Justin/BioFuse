package BioFuse::BioInfo::Objects::Variant::VCF_OB;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use Data::Dumper;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::VCF qw/ check_I16_alt /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Variant::VCF_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2021-12-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        filepath
                        addTool
                        load2refpos
                        toNormVCF
                     /;

#--- structure of object
# vcf -> filepath = $filepath
# vcf -> tag = $tag
# vcf -> tools = {bcftools => $bcftools}

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $vcf = {};
    $vcf->{filepath} = $parm{filepath} || undef;
    $vcf->{tag} = $parm{tag} || undef;

    bless($vcf);
    return $vcf;
}

#--- return file path ---
sub filepath{
    my $vcf = shift;
    return $vcf->{filepath};
}

#--- add tool ---
sub addTool{
    my $vcf = shift;
    my %parm = @_;
    $vcf->{tools}->{lc($_)} = $parm{$_} for keys %parm;
}

#--- read VCF and load each pos info to given Hash ---
sub load2refpos{
    my $vcf = shift;
    my %parm = @_;
    my $refposHf = $parm{refposHf};
    my $bcftools = $parm{bcftools} || $vcf->{tools}->{bcftools};
    my $useIDV = $parm{useIDV} || 0; # only when alt count == 1
    my $I16minBQ = $parm{I16minBQ} || 20; # vcf I16 min mean baseQ
    my $I16minMQ = $parm{I16minMQ} || 20; # vcf I16 min mean mapQ
    my $I16minTD = $parm{I16minTD} || 5;  # vcf I16 min mean tail dist

    my @theme_tag;
    open (VCF, "$bcftools view $vcf->{filepath} |") || die "Cannot bcftools read $vcf->{filepath}: $!\n";
    while(<VCF>){
        if(/^\#\#/){ # description
            next;
        }
        elsif(/^\#/){ # theme
            (my $theme_line = lc($_)) =~ s/^#//;
            @theme_tag = split /\s+/, $theme_line;
        }
        else{
            my @info = split;
            my %rgOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);

            # information
            my $pos = $rgOb{pos};
            (my $alt = $rgOb{alt}) =~ s/,?\<\*\>$//; # v4.2; has the annoying <*>
            my $info = $rgOb{info};
            my $is_indel = ($info =~ /INDEL;/) ? 1 : 0;

            my $refpos_OB = $refposHf->{$pos}; # BioFuse::BioInfo::Objects::Allele::RefPos_OB
            cluck_and_exit "<ERROR>\tcannot find refpos_OB (pos=$pos).\n".Dumper($vcf) if !defined $refpos_OB;

            # allels
            my (@Dwh, @Dfw, @Drv, @I16);
            if($info =~ /AD=([^;]+);/) { @Dwh = split /,/, $1; } # depth of whole   mapping
            if($info =~ /ADF=([^;]+);/){ @Dfw = split /,/, $1; } # depth of forward mapping
            if($info =~ /ADR=([^;]+);/){ @Drv = split /,/, $1; } # depth of reverse mapping
            if($info =~ /I16=([^;]+);/){ @I16 = split /,/, $1; }

            # whole depth
            ## donot use DP tag, its raw reads depth.
            ## one pos may be in two lines, DP is duplicated.
            ## my ($DP) = ($info =~ /DP=(\d+);/);
            # only consider the non-indel line
            $refpos_OB->addDepth(add=>sum(@Dwh)) if !$is_indel;

            # ref-allel
            my $refWdepth = shift @Dwh;
            my $refFdepth = shift @Dfw;
            my $refRdepth = shift @Drv;
            $refpos_OB->setRefDepth(refdepthAf=>[$refWdepth||0,$refFdepth||0,$refRdepth||0]);

            # mutations
            my @alt = split /,/, $alt;
            if(    scalar(@alt) != 0
                && (  # consider IDV especially for single INDEL
                      # still filter the distance from nearest read-end
                      (    $is_indel
                        && $useIDV
                        && scalar(@alt) == 1
                        && check_I16_alt(I16Af=>\@I16, min_mBQ=>$I16minBQ, min_mMQ=>$I16minMQ, min_mTD=>$I16minTD, only_TD=>1)
                      )
                      # full filtration
                    || check_I16_alt(I16Af=>\@I16, min_mBQ=>$I16minBQ, min_mMQ=>$I16minMQ, min_mTD=>$I16minTD)
                   )
               ){
                if(!$is_indel){ # SNP, i donot want to use unless-else; once it isn't noted as INDEL, it must be one letter.
                    # record
                    $refpos_OB->addMut(mut_id=>'snp,'.uc($alt[$_]), depthAf=>[$Dwh[$_], $Dfw[$_], $Drv[$_], $refWdepth]) for 0 .. $#alt;
                }
                else{ # indels
                    # IDV
                    my $IDV_value = 0;
                    my $use_IDV_bool = ( $useIDV && scalar(@alt) == 1 );
                    if( $use_IDV_bool ){ ($IDV_value) = (/IDV=(\d+);/) }
                    # this is the string around this pos
                    my $ref_allel = $rgOb{ref};
                    for my $i (0 .. $#alt){
                        my $alt_allel = $alt[$i];
                        if(length($ref_allel) < length($alt_allel)){ # insertion
                            my @ref_bases = split //, $ref_allel;
                            my @alt_bases = split //, $alt_allel;
                            # check postfixed identical base
                            my $postfix_ibase = 0;
                            while( scalar(@ref_bases) != 0 && uc($ref_bases[-1]) eq uc($alt_bases[-1]) ){
                                pop @ref_bases;
                                pop @alt_bases;
                                $postfix_ibase ++;
                            }
                            # check prefixed identical base
                            my $prefix_ibase = 0;
                            while( scalar(@ref_bases) != 0 && uc($ref_bases[0]) eq uc($alt_bases[0]) ){
                                shift @ref_bases;
                                shift @alt_bases;
                                $prefix_ibase ++;
                            }
                            # get the inserted seq
                            my $insert_seq = join('', @alt_bases);
                            # check
                            if( scalar(@ref_bases) != 0 || scalar(@alt_bases) == 0 ){
                                cluck_and_exit "<ERROR>:\twrong insertion info from VCF.\n"
                                                      ."\tref_allel:$ref_allel; alt_allel:$alt_allel\n".Dumper($vcf);
                            }
                            # assign the insertion to the position before its first base
                            my $ins_refpos = $pos + $prefix_ibase - 1;
                            my $ins_refpos_OB = $refposHf->{$ins_refpos}; # BioFuse::BioInfo::Objects::Allele::RefPos_OB
                            cluck_and_exit "<ERROR>\tcannot find refpos_OB (pos=$ins_refpos).\n".Dumper($vcf) if !defined $ins_refpos_OB;
                            # record
                            my $depthAf =   $use_IDV_bool
                                          ? [ $IDV_value, int($IDV_value/2), int($IDV_value/2), $refWdepth ] # use IDV to pretend
                                          : [ $Dwh[$i], $Dfw[$i], $Drv[$i], $refWdepth ];
                            $ins_refpos_OB->addMut(mut_id=>'ins,'.uc($insert_seq), depthAf=>$depthAf);
                        }
                        elsif(length($ref_allel) > length($alt_allel)){ # deletion
                            my @ref_bases = split //, $ref_allel;
                            my @alt_bases = split //, $alt_allel;
                            # check postfixed identical base
                            my $postfix_ibase = 0;
                            while( scalar(@alt_bases) != 0 && uc($ref_bases[-1]) eq uc($alt_bases[-1]) ){
                                pop @ref_bases;
                                pop @alt_bases;
                                $postfix_ibase ++;
                            }
                            # check prefixed identical base
                            my $prefix_ibase = 0;
                            while( scalar(@alt_bases) != 0 && uc($ref_bases[0]) eq uc($alt_bases[0]) ){
                                shift @ref_bases;
                                shift @alt_bases;
                                $prefix_ibase ++;
                            }
                            # get the deleted seq
                            my $delete_seq = join('', @ref_bases);
                            # check
                            if( scalar(@alt_bases) != 0 || scalar(@ref_bases) == 0 ){
                                cluck_and_exit "<ERROR>:\twrong deletion info from VCF.\n"
                                                      ."\tref_allel:$ref_allel; alt_allel:$alt_allel\n".Dumper($vcf);
                            }
                            # assign the deletion to the first base of itself
                            my $del_refpos = $pos + $prefix_ibase;
                            my $del_refpos_OB = $refposHf->{$del_refpos}; # BioFuse::BioInfo::Objects::Allele::RefPos_OB
                            cluck_and_exit "<ERROR>\tcannot find refpos_OB (pos=$del_refpos).\n".Dumper($vcf) if !defined $del_refpos_OB;
                            # count the deleted part into depth also
                            $del_refpos_OB->addDepth(add=>$Dwh[$i]);
                            # record, last element is for the aimPos, not this one
                            my $depthAf =   $use_IDV_bool
                                          ? [ $IDV_value, int($IDV_value/2), int($IDV_value/2), 0 ] # use IDV to pretend
                                          : [ $Dwh[$i], $Dfw[$i], $Drv[$i], 0 ];
                            $del_refpos_OB->addMut(mut_id=>'del,'.uc($delete_seq), depthAf=>$depthAf);
                        }
                    }
                }
            }
        }
    }
    close VCF;
}

#--- output norm vcf ---
sub toNormVCF{
    my $vcf = shift;
    my %parm = @_;
    my $ref = $parm{ref};
    my $norm_vcf = $parm{norm_vcf};
    my $bcftools = $parm{bcftools} || $vcf->{tools}->{bcftools};

    # bcftools norm
    my $cmd = "($bcftools norm -f $ref -o $norm_vcf->{filepath} $vcf->{filepath})";
    trible_run_for_success($cmd, 'norm', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tBCFtools norm vcf finished.\n"
                         ."\tI=$vcf->{filepath}\n"
                         ."\tO=$norm_vcf->{filepath}\n";
}

1; ## tell the perl script the successful access of this module.
