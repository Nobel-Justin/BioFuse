package BioFuse::BioInfo::GeneAnno::GetProtDomain;

use strict;
use warnings;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::LoadOn;
use BioFuse::Util::Web;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              GetProtDomain
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::GeneAnno::GetProtDomain';
#----- version --------
$VERSION = "0.02";
$DATE = '2020-01-09';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        GetProtDomain
                        ncbi_gene_query
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     # fetch domain information from NCBI Gene database #

     Usage:   perl $V_Href->{MainName} prot_dom <[Options]>

     Options:

        # Input and Output #
         -i  [s]  list of id. <required>
                  ENSG/P/T id for ensembl and NM/P id for refgene.
         -o  [s]  output domain list. <required>

        # Options #
         -t  [i]  maximum times to connect ncbi website. [10]
         -v  [i]  minimum interval in second. [3]

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
            [ id_list => undef ],
            [ domain_file => undef ],
            # option
            [ max_try => 10 ],
            [ min_interval => 3 ],

            # setting
            [ url => { prefix => 'https://www.ncbi.nlm.nih.gov/gene/?term=',
                       postfix => '&report=full_report&format=text'  } ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['gene_list'],
                                  ['domain_file']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-i:s"  => \$V_Href->{id_list},
        "-o:s"  => \$V_Href->{domain_file},
        # options
        "-t:i"  => \$V_Href->{max_try},
        "-v:i"  => \$V_Href->{min_interval},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{id_list})
             || $V_Href->{max_try} <= 0
             || $V_Href->{min_interval} <= 0
             || !defined $V_Href->{domain_file}
            );
}

#--- get protein domain information for genetic id ---
sub GetProtDomain{
    # open fh
    open (IDLIST, Try_GZ_Read($V_Href->{id_list})) or die "can not read id list: $!\n";
    open (DOM, Try_GZ_Write($V_Href->{domain_file})) or die "cannot write domain files: $!\n";
    open (DOME,Try_GZ_Write("$V_Href->{domain_file}.error")) or die "cannot write domain error files: $!\n";
    # header
    my @header = qw/ id ncbi_gid ens_t ens_p rfg_t rfg_p dom_st dom_ed db_id dom_name dom_desp /;
    print DOM join("\t", @header)."\n";
    # id
    while(<IDLIST>){
        next if /^#/;
        chomp;
        &ncbi_gene_query(id=>$_);
        sleep int(rand(3)) + $V_Href->{min_interval};
    }
    # close fh
    close IDLIST;
    close DOME;
    close DOM;
}

#--- fetch data from website ---
sub ncbi_gene_query{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $id = $parm{id};

    my $url = $V_Href->{url}->{prefix} . $id . $V_Href->{url}->{postfix};
    my $html;
    my $time = $V_Href->{max_try};
    while($time){
        # warn "$time\t$url\n" if $V_Href->{in_debug};
        $html = urlToHtmlText(url=>$url);
        warn "$html\n" if $V_Href->{in_debug};
        if($html eq 'fail' || $html =~ /An error has occured/i){
            sleep int(rand($V_Href->{max_try}-$time)) + $V_Href->{min_interval};
            $time--;
        }
        else{
            last;
        }
    }
    warn "$id\n$html\n" if $V_Href->{in_debug};
    # error
    unless($time){
        warn "<ERROR>\tfail to fetch ncbi gene info for $id\n";
        print DOME "$id\t$url\n";
        return;
    }
    # analysis
    my @c = split /\n/, $html;
    my $gene_id;
    my ($rfgT, $rfgP) = qw/ - - /;
    my $last_line = '-';
    for (my $i=0; $i<=$#c; $i++){
        if($c[$i] =~ /Gene ID:\s*(\d+)/){
            $gene_id = $1;
        }
        # [NX]M -> [NX]P, always have
        # XM_ (mRNA), XR_ (non-coding RNA), and XP_ (protein): RefSeq pipeline
        # NM_ (mRNA), NR_ (non-coding RNA), and NP_ (protein): RefSeq curated records
        if($c[$i] =~ /([NX]M_[\d\.]+)\s\S+\s([NX]P_[\d\.]+)/i){
            ($rfgT, $rfgP) = ($1, $2);
        }
        # Conserved Domains (3)summary
        if($c[$i] =~ /Conserved\s*Domains\s*\((\d+)\)\s*summary/i){
            my $domC = $1;
            # ensembl id, some transcript may lack
            my ($ensP, $ensT) = qw/ - - /;
            if($last_line =~ /Related\s*Ensembl:\s*(ENSP[\d\.]+),\s*(ENST[\d\.]+)/i){
                ($ensP, $ensT) = ($1, $2);
            }
            warn "$domC\n$last_line\n$ensP, $ensT\n" if $V_Href->{in_debug};
            # each domain
            while($domC--){
                1 while $c[++$i] =~ /^\s+$/; # pfam00870:
                my ($db_id) = ($c[$i] =~ /(\S+):/);
                1 while $c[++$i] =~ /^\s+$/; # Location:95 - 289
                my ($stp, $edp) = ($c[$i] =~ /Location:\s*(\d+)\s*-\s*(\d+)/i);
                1 while $c[++$i] =~ /^\s+$/; # P53; P53 DNA-binding domain
                my ($dom_name, $dom_desp) = ($c[$i] =~ /([^\s;]+);?\s*(.*)/);
                $dom_desp ||= '-';
                # output
                print DOM join("\t", $id, $gene_id, $ensT, $ensP, $rfgT, $rfgP, $stp, $edp, $db_id, $dom_name, $dom_desp)."\n";
            }
        }
        $last_line = $c[$i];
    }
    # inform
    stout_and_sterr "[INFO]\t$_ done\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
