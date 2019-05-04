package BioFuse::BioInfo::Objects::GeneAnno::GTF_lineOB;

use strict;
use warnings;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::GeneAnno::GTF_lineOB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        get_refSeg
                        get_strand
                        get_regionType
                        get_stpos
                        get_edpos
                        get_ENSid
                        get_name
                        get_source
                        get_version
                        get_biotype
                        get_transLocusTag
                        get_transNote
                        get_transProduct
                        get_transExonNO
                        get_proteinENSid
                        update_refSeg
                        match_gtf_source
                     /;

#--- structure of object
# gtfLine -> refSeg = $refSeg
# gtfLine -> sourceORbiotype = $sourceORbiotype, the 2nd column of gtf
# gtfLine -> regionType = $regionType
# gtfLine -> regionFstp = $regionFstp, forward strand start pos
# gtfLine -> regionFedp = $regionFedp, forward strand end pos
# gtfLine -> strand = $strand
# gtfLine -> metaInfo = $metaInfo, most information of gene/trans/exon/protein comes from this string

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my @ele = split /\t+/, $parm{gtf_lineinfo};
    my $gtfLine = {};
    $gtfLine->{refSeg} = $ele[0];
    $gtfLine->{sourceORbiotype}  = $ele[1];
    $gtfLine->{regionType} = $ele[2];
    $gtfLine->{regionFstp} = $ele[3];
    $gtfLine->{regionFedp} = $ele[4];
    $gtfLine->{strand} = $ele[6];
    $gtfLine->{metaInfo} = $ele[8];

    bless($gtfLine);
    return $gtfLine;
}

#--- return refSeg ---
sub get_refSeg{
    my $gtfLine = shift;
    return $gtfLine->{refSeg};
}

#--- return gene strand ---
sub get_strand{
    my $gtfLine = shift;
    return $gtfLine->{strand};
}

#--- return region type ---
sub get_regionType{
    my $gtfLine = shift;
    return $gtfLine->{regionType};
}

#--- return region start pos ---
sub get_stpos{
    my $gtfLine = shift;
    return $gtfLine->{regionFstp};
}

#--- return region end pos ---
sub get_edpos{
    my $gtfLine = shift;
    return $gtfLine->{regionFedp};
}

#--- return gene/trans ens_id ---
sub get_ENSid{
    my $gtfLine = shift;
    my %parm = @_;
    my $type = $parm{type} || 'gene';

    my $regexKeyword = $type eq 'gene' ? 'gene_id' : 'transcript_id';
    if($gtfLine->{metaInfo} =~ /$regexKeyword\s\"([^\"]+)\"/){
        return $1;
    }
    else{
        warn_and_exit "<ERROR>\tcannot find $type ENS-id from gtf-lineInfo object.\n".Dumper($gtfLine);
    }
}

#--- return gene/trans name ---
sub get_name{
    my $gtfLine = shift;
    my %parm = @_;
    my $type = $parm{type} || 'gene';

    my $regexKeyword = $type eq 'gene' ? 'gene_name' : 'transcript_name';
    my $name = ($gtfLine->{metaInfo} =~ /$regexKeyword\s\"([^\"]+)\"/ ? $1 : 'NA-'.$gtfLine->get_ENSid(type => $type));
       $name =~ s/[\(\)\/\\\s]+/_/g;
    return $name;
}

#--- return gene/trans source ---
# former versions (before v75) do not contain.
sub get_source{
    my $gtfLine = shift;
    my %parm = @_;
    my $type = $parm{type} || 'gene';

    my $regexKeyword = $type eq 'gene' ? 'gene_source' : 'transcript_source';
    return ($gtfLine->{metaInfo} =~ /$regexKeyword\s\"([^\"]+)\"/ ? $1 : 'NA');
}

#--- return gene/trans version ---
# former versions (before v77) do not contain.
sub get_version{
    my $gtfLine = shift;
    my %parm = @_;
    my $type = $parm{type} || 'gene';

    my $regexKeyword = $type eq 'gene' ? 'gene_version' : 'transcript_version';
    return ($gtfLine->{metaInfo} =~ /$regexKeyword\s\"([^\"]+)\"/ ? "v$1" : 'NA');
}

#--- return gene/trans biotype ---
# former versions (long-long ago..) do not contain, but list it at column NO.2
sub get_biotype{
    my $gtfLine = shift;
    my %parm = @_;
    my $type = $parm{type} || 'gene';

    my $regexKeyword = $type eq 'gene' ? 'gene_biotype' : 'transcript_biotype';
    return ($gtfLine->{metaInfo} =~ /$regexKeyword\s\"([^\"]+)\"/ ? $1 : $gtfLine->{sourceORbiotype});
}

#--- return trans locus_tag ---
# specific for FuseSV virus meta-info
sub get_transLocusTag{
    my $gtfLine = shift;
    return ($gtfLine->{metaInfo} =~ /locus_tag\s\"([^\"]+)\"/ ? $1 : 'N/A');
}

#--- return trans note ---
# specific for FuseSV virus meta-info
sub get_transNote{
    my $gtfLine = shift;
    return ($gtfLine->{metaInfo} =~ /note\s\"([^\"]+)\"/ ? $1 : 'N/A');
}

#--- return trans product ---
# specific for FuseSV virus meta-info
sub get_transProduct{
    my $gtfLine = shift;
    return ($gtfLine->{metaInfo} =~ /product\s\"([^\"]+)\"/ ? $1 : 'N/A');
}

#--- return trans exon NO ---
sub get_transExonNO{
    my $gtfLine = shift;
    return ($gtfLine->{metaInfo} =~ /exon_number\s\"([^\"]+)\"/ ? $1 : 'N/A');
}

#--- return protein ENSid ---
sub get_proteinENSid{
    my $gtfLine = shift;
    return ($gtfLine->{metaInfo} =~ /protein_id\s\"([^\"]+)\"/ ? $1 : 'NA');
}

#--- update refSeg ---
sub update_refSeg{
    my $gtfLine = shift;
    my %parm = @_;
    $gtfLine->{refSeg} = $parm{refSeg};
}

#--- match gtf source ---
sub match_gtf_source{
    my $gtfLine = shift;
    my %parm = @_;
    if(    $gtfLine->{metaInfo} =~ /gene_source\s\"([^\"]+)\"/
        && $parm{gtf_source} !~ /$1/
    ){
        return 0;
    }
    else{ # gene_source may not exist in ensembl before v75
        return 1;
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
