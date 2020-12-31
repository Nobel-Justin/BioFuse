package BioFuse::BioInfo::GeneAnno::GetProtDomain;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
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
$VERSION = "0.03";
$DATE = '2020-11-22';

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
                        query_by_idList
                        get_db_info
                        ncbi_gene_query_protDom
                        ensm_tran_query_protDom
                        write_domain_tsv
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     # fetch domain information from database #

     Usage:   perl $V_Href->{MainName} prot_dom <[Options]>

     Options:

        # Input and Output #
         -i  [s]  list of query id. <required>
                  for ncbi: ENSG/P/T or [NX]M/P id.
                  for ensm: ONLY ENST id.
         -o  [s]  output domain list. <required>

        # Options #
         -d  [s]  select database. <required>
                  ncbi: NCBI-GENE; ensm: Ensembl-Domains
         -t  [i]  maximum try times to connect database website for one query. [10]

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
            [ db => 'null' ], # ncbi, ensm
            [ max_try => 10 ],
            [ min_interval => 3 ],

            # setting
            [ db_list => {ncbi => 1, ensm => 1} ],
            [ url => { ncbi => { prefix => 'https://www.ncbi.nlm.nih.gov/gene/?term=',
                                 postfix => '&report=full_report&format=text'  },
                       ensm => { prefix => 'http://ensembl.org/Homo_sapiens/Component/Transcript/Domains/domains?t=',
                                 postfix => '' } } ],
            [ ensm_header => [ 'source', 'dom_st', 'dom_ed', 'dom_desp', 'dom_acce', 'interpro' ] ],
            # container
            [ ncbi_gene => {} ],
            [ ensm_tran => {} ],
            [ failquery => [] ],
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
        "-d:s"  => \$V_Href->{db},
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
             || !exists $V_Href->{db_list}->{$V_Href->{db}}
             || $V_Href->{max_try} <= 0
             || $V_Href->{min_interval} <= 0
             || !defined $V_Href->{domain_file}
            );
}

#--- get protein domain information for genetic id ---
sub GetProtDomain{
    # query by list
    &query_by_idList;
    # output
    &write_domain_tsv;
}

#--- read id list and query one by one ---
sub query_by_idList{
    open (IDLIST, Try_GZ_Read($V_Href->{id_list})) or die "can not read id list: $!\n";
    while(<IDLIST>){
        next if /^#/;
        chomp;
        if($V_Href->{db} eq 'ncbi'){
            &ncbi_gene_query_protDom(query_id => $_);
        }
        elsif($V_Href->{db} eq 'ensm'){
            &ensm_tran_query_protDom(query_id => $_);
        }
        sleep int(rand(3)) + $V_Href->{min_interval};
    }
    close IDLIST;
}

#--- fetch data via url ---
sub get_db_info{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $url = $parm{url};
    my $query_id = $parm{query_id};

    my $dbInfo;
    my $time = $V_Href->{max_try};
    while($time){
        warn "$time\t$url\n" if $V_Href->{in_debug};
        $dbInfo = urlToHtmlText(url=>$url);
        warn "$dbInfo\n" if $V_Href->{in_debug};
        if(   $dbInfo eq 'fail'
           || $dbInfo =~ /An error has occured/i
           || $dbInfo =~ /Server Exception/i
           || $dbInfo =~ /AJAX error/i
        ){
            sleep int(rand($V_Href->{max_try}-$time)) + $V_Href->{min_interval};
            $time--;
        }
        else{
            last;
        }
    }
    warn "$query_id\n$dbInfo\n" if $V_Href->{in_debug};
    # error
    unless($time){
        warn "<ERROR>\tfail to fetch $V_Href->{db} info for $query_id\n";
        push @{$V_Href->{failquery}}, [$query_id, $url];
        return 0;
    }
    return $dbInfo;
}

#--- fetch data from ncbi gene database ---
sub ncbi_gene_query_protDom{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $query_id = $parm{query_id};

    # fetch ncbi-gene data
    my $url = $V_Href->{url}->{ncbi}->{prefix} . $query_id . $V_Href->{url}->{ncbi}->{postfix};
    my $dbInfo = &get_db_info(query_id => $query_id, url => $url);
    return unless $dbInfo;
    # analysis
    my @c = split /\n/, $dbInfo;
    my $gene_id;
    my ($rfgT, $rfgP) = qw/ - - /;
    my $last_line = '-';
    for (my $i=0; $i<=$#c; $i++){
        if($c[$i] =~ /Gene ID:\s*(\d+)/){
            $gene_id = $1;
            if(exists $V_Href->{ncbi_gene}->{$gene_id}){
                warn "<WARN>\t$gene_id gene has been processed.\n";
                # record the duplicated
                push @{$V_Href->{failquery}}, [$query_id, "$gene_id-has-been-processed(duplicated)"];
                return;
            }
        }
        # [NX]M -> [NX]P, always have
        # XM_ (mRNA), XR_ (non-coding RNA), and XP_ (protein): RefSeq pipeline
        # NM_ (mRNA), NR_ (non-coding RNA), and NP_ (protein): RefSeq curated records
        if($c[$i] =~ /([NX]M_[\d\.]+)\s\S+\s([NX]P_[\d\.]+)/i){
            ($rfgT, $rfgP) = ($1, $2);
            $V_Href->{ncbi_gene}->{$gene_id}->{$rfgT} = { rfgP => $rfgP, dom => [], query_id => $query_id };
        }
        # Conserved Domains (x)summary
        if($c[$i] =~ /Conserved\s*Domains\s*\((\d+)\)\s*summary/i){
            my $domC = $1;
            my $rfgT_hf = $V_Href->{ncbi_gene}->{$gene_id}->{$rfgT};
            # ensembl id, some transcript may lack
            my ($ensP, $ensT) = qw/ - - /;
            if($last_line =~ /Related\s*Ensembl:\s*(ENSP[\d\.]+),\s*(ENST[\d\.]+)/i){
                ($ensP, $ensT) = ($1, $2);
            }
            $rfgT_hf->{ensT} = $ensT;
            $rfgT_hf->{ensP} = $ensP;
            warn "$domC\n$last_line\n$ensP, $ensT\n" if $V_Href->{in_debug};
            # each domain
            while($domC--){
                1 while $c[++$i] =~ /^\s+$/; # pfam00870:
                my ($dom_acce) = ($c[$i] =~ /(\S+):/);
                1 while $c[++$i] =~ /^\s+$/; # Location:95 - 289
                my ($stp, $edp) = ($c[$i] =~ /Location:\s*(\d+)\s*-\s*(\d+)/i);
                1 while $c[++$i] =~ /^\s+$/; # P53; P53 DNA-binding domain
                my ($dom_name, $dom_desp) = ($c[$i] =~ /([^\s;]+);?\s*(.*)/);
                $dom_desp ||= '-';
                # record
                push @{$rfgT_hf->{dom}}, { dom_name => $dom_name, dom_desp => $dom_desp,
                                           stp => $stp, edp => $edp, dom_acce => $dom_acce };
            }
        }
        $last_line = $c[$i];
    }
    # inform
    stout_and_sterr "[INFO]\t$query_id done\n";
}

#--- fetch data from ensembl database ---
sub ensm_tran_query_protDom{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $query_id = uc($parm{query_id});

    # check
    if($query_id !~ /^ENST/i){
        warn_and_exit "<ERROR>\tONLY ENST id allowed for database ensm. $query_id\n";
        exit 1;
    }
    # fetch ensembl data
    my $url = $V_Href->{url}->{ensm}->{prefix} . $query_id . $V_Href->{url}->{ensm}->{postfix};
    my $dbHtml = &get_db_info(query_id => $query_id, url => $url);
    return unless $dbHtml;
    # analysis
    $dbHtml =~ s/(<t[hrd])/\n$1/g;
    $dbHtml =~ s/(<\/?tr>)/\n$1/g;
    my @c = split /\n/, $dbHtml;
    for (my $i=0; $i<=$#c; $i++){
        last if $c[$i] =~ /Other features/i;
        if($c[$i] =~ /<tr>/){
            my %info;
            for my $j (1..6){
                1 while $c[++$i] =~ /^\s+$/;
                my ($info) = ($c[$i] =~ /<td[^>]+>(.+)<\/td>/);
                $info = $1 if $info =~ /<a.+href=.+>(\S+)<\/a>/;
                $info{$V_Href->{ensm_header}->[$j-1]} = $info;
            }
            push @{$V_Href->{ensm_tran}->{$query_id}->{dom}}, \%info;
        }
    }
    # inform
    stout_and_sterr "[INFO]\t$query_id done\n";
}

#--- output domain info ---
sub write_domain_tsv{
    open (DOM, Try_GZ_Write($V_Href->{domain_file})) or die "cannot write domain files: $!\n";
    if($V_Href->{db} eq 'ncbi'){
        my @header = qw/ query_id ncbi_gid ens_t ens_p rfg_t rfg_p dom_st dom_ed dom_acce dom_name dom_desp /;
        print DOM join("\t", @header)."\n";
        for my $gene_id (sort keys %{$V_Href->{ncbi_gene}}){
            for my $rfgT (sort keys %{$V_Href->{ncbi_gene}->{$gene_id}}){
                my $rfgT_Hf = $V_Href->{ncbi_gene}->{$gene_id}->{$rfgT};
                for my $dom_Hf (@{$rfgT_Hf->{dom}}){
                    print DOM join("\t", $rfgT_Hf->{query_id}, $gene_id,
                                         $rfgT_Hf->{ensT}, $rfgT_Hf->{ensP},
                                         $rfgT, $rfgT_Hf->{rfgP},
                                         $dom_Hf->{stp}, $dom_Hf->{edp},
                                         $dom_Hf->{dom_acce}, $dom_Hf->{dom_name}, $dom_Hf->{dom_desp}
                                  )."\n";
                }
            }
        }
    }
    elsif($V_Href->{db} eq 'ensm'){
        print DOM join("\t", 'ens_t', @{$V_Href->{ensm_header}})."\n";
        for my $ensT (sort keys %{$V_Href->{ensm_tran}}){
            for my $dom_hf (@{$V_Href->{ensm_tran}->{$ensT}->{dom}}){
                print DOM join("\t", $ensT, map {$dom_hf->{$_}} @{$V_Href->{ensm_header}})."\n";
            }
        }
    }
    close DOM;
    # query failed
    open (DOME,Try_GZ_Write("$V_Href->{domain_file}.error")) or die "cannot write domain error files: $!\n";
    print DOME join("\t",@$_)."\n" for @{$V_Href->{failquery}};
    close DOME;
}

#--- 
1; ## tell the perl script the successful access of this module.
