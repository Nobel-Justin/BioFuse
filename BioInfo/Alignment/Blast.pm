package BioFuse::BioInfo::Alignment::Blast;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ min /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::String qw/ mapSeqWithMM /;
use BioFuse::BioInfo::FASTA qw/ FAfhFindSeq read_fasta_file /;
use BioFuse::BioInfo::Objects::Allele::RefPos_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              blast_formatdb
              blastn_query
              BlastM0ToRefPos
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Alignment::Blast';
#----- version --------
$VERSION = "0.03";
$DATE = '2022-01-19';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        blast_formatdb
                        blastn_query
                        BlastM0ToRefPos
                        dealOneQueryMapIF
                        digest_m0MapBlock
                        adjust_INSbehMM
                        fillGapForLargeInDel
                        mapIF2RefPos
                     /;

#--- blast formatdb ---
sub blast_formatdb{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ref_fa = $parm{ref_fa};
    my $formatdb = $parm{formatdb};

    # construct blast reference
    my $cmd = "$formatdb -i $ref_fa -p F -l $ref_fa.log";
    trible_run_for_success($cmd, 'formatdb_virusRef', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tblast formatdb done.\n"
                         ."\tR=$ref_fa\n";
}

#--- blastn query database ---
sub blastn_query{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ref_fa = $parm{ref_fa};
    my $query_fa = $parm{query_fa};
    my $output = $parm{output};
    my $blastall = $parm{blastall};
    my $format_m = $parm{format_m} || 0;

    # blastn
    my $cmd = "$blastall -m $format_m -a 8 -p blastn -F F -e 1e-06 -i $query_fa -d $ref_fa -o $output";
    trible_run_for_success($cmd, 'blastn_virusRef', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tblastn query done.\n"
                         ."\tR=$ref_fa\n"
                         ."\tQ=$query_fa\n"
                         ."\tO=$output\n";
}

#--- read query mapping info from blast m0 result ---
sub BlastM0ToRefPos{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $refposHf = $parm{refposHf}; # refpos hash ref
    my $ref_fa   = $parm{ref_fa}; # reference, i.e., subject
    my $query_fa = $parm{query_fa};
    my $blast_m0 = $parm{blast_m0};
    my $qyIDlist = $parm{qyIDlist}; # only deal with query id in this list
    my $skipGapF = $parm{skipGapF} || 0; # not fill gap for large InDel
    my $firstMap = $parm{firstMap} || 0; # only use the first map of each query_id
    my $debug    = $parm{debug};

    # load query id list if given
    my %qyID;
    if(defined $qyIDlist){
        open (QID,Try_GZ_Read($qyIDlist)) || die "fail read $qyIDlist: $!\n";
        while(<QID>){
            my ($query_id) = (split)[0];
            $qyID{$query_id} = 1;
        }
        close QID;
    }

    # to load query sequence one by one when use later, if not skipGapF
    open (my $qyFAfh,Try_GZ_Read($query_fa)) || die "fail read $query_fa: $!\n";
    # load all subject sequence
    my %sbjct_seq;
    read_fasta_file(FaFile=>$ref_fa, Seg2Seq_Href=>\%sbjct_seq);

    # initialize the subject's refpos hash
    for my $sbjct_id (keys %sbjct_seq){
        my @sbjct_seq = split //, $sbjct_seq{$sbjct_id};
        $refposHf->{$sbjct_id}->{$_+1} = BioFuse::BioInfo::Objects::Allele::RefPos_OB->new(pos=>$_+1, refAllele=>$sbjct_seq[$_]) for 0..$#sbjct_seq;
    }

    my %mapIF; # mapinfo hash
    # option array pre-defined
    my @opt = ( mapIF_Hf => \%mapIF,     refposHf => $refposHf,
                sjSeq_Hf => \%sbjct_seq,   qyFAfh => $qyFAfh,
                firstMap => $firstMap,   skipGapF => $skipGapF  );
    # load m0 file
    my $last_query_id = '_fs_init_undef_fs_';
    open (my $m0fh,Try_GZ_Read($blast_m0)) || die "fail read $blast_m0: $!\n";
    while (<$m0fh>){
        if(!/^Query=/){
            next;
        }
        else{
            # result(s) of one query_id starts here
            my ($query_id) = (/Query=\s(\S+)/);
            # required query_id?
            next if defined $qyIDlist && !exists $qyID{$query_id};
            # deal the last query_id
            if(    !exists $mapIF{$query_id}
                &&  exists $mapIF{$last_query_id}
            ){
                # deal map info of last query_id
                &dealOneQueryMapIF(query_id=>$last_query_id, @opt);
                # sweep, release memory
                delete $mapIF{$last_query_id};
            }
            # deal with this maping block
            &digest_m0MapBlock(m0fh=>$m0fh, query_id=>$query_id, mapIF_Hf=>\%mapIF, debug=>$debug);
            # update
            $last_query_id = $query_id;
        }
    }
    close $m0fh;
    # deal map info of last query_id
    &dealOneQueryMapIF(query_id=>$last_query_id, @opt) if exists $mapIF{$last_query_id};
    close $qyFAfh;

    # inform
    stout_and_sterr "[INFO]\tload blast m0 result to refpos hash done.\n"
                         ."\tI=$blast_m0\n"
                         ."\tR=$ref_fa\n"
                         ."\tQ=$query_fa\n";
}

#--- deal map info of given query_id, and info to refpos hash ---
sub dealOneQueryMapIF{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $mapIF_Hf = $parm{mapIF_Hf};
    my $refposHf = $parm{refposHf}; # refpos hash ref
    my $query_id = $parm{query_id};
    my $sjSeq_Hf = $parm{sjSeq_Hf};
    my $qyFAfh = $parm{qyFAfh};
    my $firstMap = $parm{firstMap} || 0; # only use the first map
    my $skipGapF = $parm{skipGapF} || 0; # not fill gap for large InDel

    # Try to concatenate devision by possible large-InDel, if not skipGapF
    &fillGapForLargeInDel(
        mapIF_Hf => $mapIF_Hf,
        query_id => $query_id,
        sjSeq_Hf => $sjSeq_Hf,
          qyFAfh => $qyFAfh
    ) if !$skipGapF;
    # summary map info of last query_id to refpos
    &mapIF2RefPos(
        mapIF_Hf => $mapIF_Hf,
        refposHf => $refposHf,
        query_id => $query_id,
        firstMap => $firstMap
    );
}

#--- one mapping block in m0 format ---
## donot record in mapIF_Hf when 'No hits found'
sub digest_m0MapBlock{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $m0fh = $parm{m0fh};
    my $query_id = $parm{query_id};
    my $mapIF_Hf = $parm{mapIF_Hf}; # ref of mapinfo hash
    my $debug = $parm{debug};

    # my $IUPAC_code = 'acgtryswkmndhvn'; # \., -, use later

    # subject mapped by seq of this query_id
    my ($sbjct_id);
    # SIGNAL for control
    ## 1, initial, 'wait result'
    ## 2, meet Score
    my $signal = 1;
    while($signal){
        my $line_info = <$m0fh>;
        if($line_info !~ /^\s+Score\s+=/){
            # no hits, just skip
            if( $line_info =~ /No\s+hits\s+found/i ){
                last;
            }
            # mapped subject id
            elsif( $line_info =~ /^>(.+)\s$/ ){
                $sbjct_id = $1;
                # reset the signal as 'wait result'.
                $signal = 1;
                next;
            }
            # search next line for certain query result
            elsif( $signal == 1 ){
                next;
            }
            # empty lines
            elsif( $line_info =~ /^\s+$/ ){
                next;
            }
            # other not-allowed info, means finish of this query_id's result
            else{
                last;
            }
        }
        # meet Score:
        $signal = 2;
        # instance: Score =  149 bits (75), Expect = 6e-40
        my ($score, $expect_E) = ($line_info =~ /Score\s+=\s+([\d\.e\+]+)\s+bits.*Expect\s+=\s+(\S+)/);
        #die "$line_info\n$score, $expect_E\n"; ## debug
        # instance: Identities = 84/87 (96%)
        my $Identities_line_info = <$m0fh>;
        my ($mapped_length, $identities) = ($Identities_line_info =~ /Identities\s+=\s+\d+\/(\d+)\s+\(([^\(\)]+)\)/);
        #die "$Identities_line_info\n$mapped_length, $identities\n"; ## debug
        ($identities = '0.' . $identities) =~ s/\%$//;
        # instance: Strand = Plus / Plus
        my $Strand_line_info = <$m0fh>;
        my ($query_strd, $sbjct_strd) = ($Strand_line_info =~ /Strand\s+=\s+(\S+)\s+\/\s+(\S+)/);
        #die "$Strand_line_info\n$query_strd, $sbjct_strd\n"; ## debug
        $query_strd = ($query_strd =~ /Plus/) ? 1 : -1;
        $sbjct_strd = ($sbjct_strd =~ /Plus/) ? 1 : -1;
        # read query mapping details repeatedly
        my $cat_map_info;
        my @cat_query_info;
        my @cat_sbjct_info;
        while (<$m0fh>){
            if(/^Query:/){ # find Query
                chomp(my $query_info = $_);
                chomp(my $map_info = <$m0fh>);
                chomp(my $sbjct_info = <$m0fh>);
                # split info for details
                my ($prefix, $query_st, $query_seq, $query_ed) = ($query_info =~ /(Query:\s+(\d+)\s+)(\S+)\s+(\d+)/i); # /(Query:\s+(\d+)\s+)([acgt\-n]+)\s+(\d+)/i
                my $pf_len = length($prefix);
                $map_info =~ s/^\s{$pf_len,$pf_len}//;
                my ($sbjct_st, $sbjct_seq, $sbjct_ed) = ($sbjct_info =~ /(\d+)\s+(\S+)\s+(\d+)/i); # /(\d+)\s+([$IUPAC_code\.\-]+)\s+(\d+)/i
                # check
                if(length($query_seq) != length($sbjct_seq) || length($query_seq) != length($map_info)){
                    cluck_and_exit  "query_id: $query_id\nquery seq length is not equal to sbjct seq length:\n"
                                   ."query: $query_seq\n"
                                   ."align: $map_info\n"
                                   ."sbjct: $sbjct_seq\n"
                                   ."$query_info\n"
                                   ."$sbjct_info\n";
                }
                # try to concatnate the mapping details
                if(!defined($cat_map_info)){ # first meets
                    # just load on
                    $cat_map_info = $map_info;
                    @cat_query_info = ($query_st, $query_seq, $query_ed);
                    @cat_sbjct_info = ($sbjct_st, $sbjct_seq, $sbjct_ed);
                }
                else{
                    # smart check the continuous attribute
                    if( ($query_st-$cat_query_info[2]) != $query_strd || ($sbjct_st-$cat_sbjct_info[2]) != $sbjct_strd ){
                        cluck_and_exit  "the mapping is not continuous, Query_ID: $query_id\n"
                                       ."query: prev:$cat_query_info[0]-$cat_query_info[2], now:$query_st-$query_ed\n"
                                       ."sbjct: prev:$cat_sbjct_info[0]-$cat_sbjct_info[2], now:$sbjct_st-$sbjct_ed\n";
                    }
                    # concatnate
                    $cat_map_info      .= $map_info;
                    $cat_query_info[1] .= $query_seq;
                    $cat_query_info[2]  = $query_ed; # update the end
                    $cat_sbjct_info[1] .= $sbjct_seq;
                    $cat_sbjct_info[2]  = $sbjct_ed; # update the end
                }
                # check the map_info length, and trigger breaking
                last if length($cat_map_info) == $mapped_length;
            }
        }
        # specifically, adjust the map-info for 'insertion' after 'mismatch' position
        # only change 'sbjct_seq'
        &adjust_INSbehMM(
            querySeq_Sref => \$cat_query_info[1],
             mapInfo_Sref => \$cat_map_info,
            sbjctSeq_Sref => \$cat_sbjct_info[1],
                    debug => $debug
        );
        # store this mapping for this Query
        push @{$mapIF_Hf->{$query_id}}, {        score => $score,
                                              expect_E => $expect_E,
                                                 ident => $identities,
                                            query_strd => $query_strd,
                                            sbjct_strd => $sbjct_strd,
                                            align_strd => $sbjct_strd * $query_strd,
                                              sbjct_id => $sbjct_id,
                                              map_info => $cat_map_info,
                                            query_info => \@cat_query_info,
                                            sbjct_info => \@cat_sbjct_info
                                        };
        #print "$cat_query_info[1]\n$cat_map_info\n$cat_sbjct_info[1]\n"; # debug
    }
}

#--- specifically, adjust the map-info for 'INSertion' behind 'MisMatch' ---
## only change 'sbjct_seq'
sub adjust_INSbehMM{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $querySeq_Sref = $parm{querySeq_Sref};
    my $mapInfo_Sref  = $parm{mapInfo_Sref};
    my $sbjctSeq_Sref = $parm{sbjctSeq_Sref};
    my $debug = $parm{debug};

    # find the insertion pos
    my $update_sign = 0;
    my @sbjctSeq = split //, $$sbjctSeq_Sref;
    my @querySeq = split //, $$querySeq_Sref;
    for ( my $i = 0; $i <= $#sbjctSeq; $i++ ){
        # insertion follows mismatch
        if(    $sbjctSeq[$i] eq '-' # insertion start here
            && $sbjctSeq[$i-1] ne '-' # last pos is not insertion, {never be}.
            && $querySeq[$i-1] ne '-' # last pos is not deletion.
            && lc($sbjctSeq[$i-1]) ne lc($querySeq[$i-1]) # last pos is mismatch
        ){
            # 1st, check the mismatch extending to 5-prime
            my $mism_i = $i-1;
            while(    $mism_i > 0 # cannot reach the 5-prime end
                   && $sbjctSeq[$mism_i-1] ne '-' # last pos is not insertion, {never be}.
                   && $querySeq[$mism_i-1] ne '-' # last pos is not deletion.
                   && lc($sbjctSeq[$mism_i-1]) ne lc($querySeq[$mism_i-1]) # last one is still mismatch
            ){
                $mism_i --;
            }
            # 2nd, check the insertion extending to 3-prime
            my $ins_len = 1; # already is one-bp long
            while(    $i+$ins_len <= $#sbjctSeq
                   && $sbjctSeq[$i+$ins_len] eq '-'
            ){
                $ins_len ++;
            }
            # 3rd, use map_info to validate
            my $substr_idx = $mism_i;
            my $substr_len = $i - $mism_i + $ins_len;
            my $thisReg_map_Info = substr($$mapInfo_Sref,  $substr_idx, $substr_len);
            my $thisReg_querySeq = substr($$querySeq_Sref, $substr_idx, $substr_len);
            my $thisReg_sbjctSeq = substr($$sbjctSeq_Sref, $substr_idx, $substr_len);
            ## map info of this region should be all blanks
            if( $thisReg_map_Info =~ /\S/ ){
                cluck_and_exit `date`."<ERROR>:\tmeet wrong map-info when adjust insertion [idx:$i] after mismatch from Blast result:\n"
                                     ."\tquerySeq: >$thisReg_querySeq<\n"
                                     ."\tmap_info: >$thisReg_map_Info<\n"
                                     ."\tsbjctSeq: >$thisReg_sbjctSeq<\n";
            }
            # 4th, exchange location of insertion and mismatch
            ## only need to change 'sbjctSeq'
            my $thisReg_new_sbjctSeq = substr($thisReg_sbjctSeq, $i-$mism_i) . substr($thisReg_sbjctSeq, 0, $i-$mism_i);
            for my $j ( 1 .. $substr_len ){
                $sbjctSeq[$mism_i+$j-1] = substr($thisReg_new_sbjctSeq, $j-1, 1);
            }
            ## record update sign
            $update_sign = 1;
            ## infrom this exchange in debug mode
            if( $debug ){
                stout_and_sterr  "Adjust insertion [idx:$i] after mismatch from Blast result.\n"
                                ." From:\n"
                                ."\tquerySeq: >$thisReg_querySeq<\n"
                                ."\tmap_info: >$thisReg_map_Info<\n"
                                ."\tsbjctSeq: >$thisReg_sbjctSeq<\n"
                                ." To:\n"
                                ."\tquerySeq: >$thisReg_querySeq<\n"
                                ."\tmap_info: >$thisReg_map_Info<\n"
                                ."\tsbjctSeq: >$thisReg_new_sbjctSeq<\n";
            }
            # update the $i
            $i += $ins_len;
        }
    }
    # update the 'sbjct_seq'
    $$sbjctSeq_Sref = join('', @sbjctSeq) if $update_sign;
}

#--- try to fill the gap induced by large insertion or deletion ---
sub fillGapForLargeInDel{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $mapIF_Hf = $parm{mapIF_Hf};
    my $query_id = $parm{query_id};
    my $sjSeq_Hf = $parm{sjSeq_Hf};
    my $qyFAfh = $parm{qyFAfh};

    # prepare sequence of this query_id
    my $qySeq_Hf = { $query_id => FAfhFindSeq(FAfh=>$qyFAfh, needSegName=>$query_id) };

    my $MapInfo_Aref = $mapIF_Hf->{$query_id};
    my %discard_idx;
    # find the connection iteratively until no more found
    while(1){
        my $map_count = scalar(@$MapInfo_Aref);
        # no more mappings to deal with
        if( $map_count - scalar(keys %discard_idx) == 1 ){
            last;
        }
        # update sign
        my $update_sign = 0;
        # loop to select fore mapping
        for my $fore_map_idx ( 0 .. $map_count-1 ){
            # already discard, old mapping
            next if( exists($discard_idx{$fore_map_idx}) );
            # fore mapping information
            my $fore_mapInfo_Href = $MapInfo_Aref->[$fore_map_idx];
            my $fore_map_sbjct = $fore_mapInfo_Href->{sbjct_id};
            my $fore_map_Qstrd = $fore_mapInfo_Href->{query_strd};
            my $fore_map_Sstrd = $fore_mapInfo_Href->{sbjct_strd};
            my $fore_map_QedP  = $fore_mapInfo_Href->{query_info}->[2];
            my $fore_map_SedP  = $fore_mapInfo_Href->{sbjct_info}->[2];
            my $fore_map_SpanL = length $fore_mapInfo_Href->{map_info};
            # loop to test each other mappings as follower mapping
            for my $cand_map_idx ( 0 .. $map_count-1 ){
                # already discard, old mapping
                next if( exists($discard_idx{$cand_map_idx}) );
                # cand mapping information
                my $cand_mapInfo_Href = $MapInfo_Aref->[$cand_map_idx];
                my $cand_map_sbjct = $cand_mapInfo_Href->{sbjct_id};
                my $cand_map_Qstrd = $cand_mapInfo_Href->{query_strd};
                my $cand_map_Sstrd = $cand_mapInfo_Href->{sbjct_strd};
                my $cand_map_QstP  = $cand_mapInfo_Href->{query_info}->[0];
                my $cand_map_SstP  = $cand_mapInfo_Href->{sbjct_info}->[0];
                my $cand_map_SpanL = length $cand_mapInfo_Href->{map_info};
                # filters
                if(    $cand_map_idx == $fore_map_idx # not itself
                    || $fore_map_sbjct ne $cand_map_sbjct # must map to same subject
                    || $fore_map_Qstrd ne $cand_map_Qstrd # query mapping strand
                    || $fore_map_Sstrd ne $cand_map_Sstrd # sbjct mapping strand
                ){
                    next;
                }
                # the two segment must be continuous at position
                my $map_Qdist = abs($cand_map_QstP - $fore_map_QedP);
                my $map_Sdist = abs($cand_map_SstP - $fore_map_SedP);
                ## large DEL: 'none' OR 'slight' gap on query_seq, and 'large' gap on sbjct_seq
                ## large INS: 'none' OR 'slight' gap on sbjct_seq, and 'large' gap on query_seq
                ## large INS is from 'unkown source', Here
                if(    $fore_map_Qstrd * $fore_map_QedP < $cand_map_Qstrd * $cand_map_QstP
                    && $fore_map_Sstrd * $fore_map_SedP < $cand_map_Sstrd * $cand_map_SstP
                    && (  # query_seq has large deletion
                          (   $map_Qdist <= 30 # empirical distance in practise
                           && $map_Sdist >= $map_Qdist
                           && $map_Sdist <= $map_Qdist + 100 # maximum two 50bp deletions allowed
                          )
                        || # query_seq has large insertion (unkown source)
                          (   $map_Qdist >= $map_Sdist
                           && $map_Sdist <= 30
                           && $map_Qdist <= $map_Sdist + 100 # maximum two 50bp insertions allowed
                          )
                       )
                ){
                    # prepare makeup seq for query_seq gap
                    my $makeup_Qidx = min($fore_map_QedP, $cand_map_QstP);
                    my $makeup_Qlen = $map_Qdist - 1;
                    my $makeup_Qseq = uc( substr($qySeq_Hf->{$query_id}, $makeup_Qidx, $makeup_Qlen) );
                    if( $fore_map_Qstrd == -1 && $makeup_Qlen > 0 ){
                        ($makeup_Qseq = reverse uc($makeup_Qseq)) =~ tr/ACGT/TGCA/;
                    }
                    # prepare makeup seq for sbjct_seq gap
                    my $makeup_Sidx = min($fore_map_SedP, $cand_map_SstP);
                    my $makeup_Slen = $map_Sdist - 1;
                    my $makeup_Sseq = uc( substr($sjSeq_Hf->{$fore_map_sbjct}, $makeup_Sidx, $makeup_Slen) );
                    if( $fore_map_Sstrd == -1 && $makeup_Slen > 0 ){
                        ($makeup_Sseq = reverse uc($makeup_Sseq)) =~ tr/ACGT/TGCA/;
                    }
                    # compare and makeup map-info
                    ## remained part, maybe mismatch(s)
                    my ($remPartMMct, $remPartMidx, $remPartMdet) = mapSeqWithMM(str_a=>$makeup_Qseq, str_b=>$makeup_Sseq);
                    my $headGapLen = $remPartMidx;
                    my $tailGapLen = abs($makeup_Slen - $makeup_Qlen) - $headGapLen;
                    my $makeup_MapI = (' ' x $headGapLen) . $remPartMdet . (' ' x $tailGapLen);
                    # gap part
                    if( $makeup_Slen > $makeup_Qlen ){ # query_seq has large deletion
                        $makeup_Qseq = ('-' x $headGapLen) . $makeup_Qseq . ('-' x $tailGapLen);
                    }
                    else{ # query_seq has large insertion
                        $makeup_Sseq = ('-' x $headGapLen) . $makeup_Sseq . ('-' x $tailGapLen);
                    }
                    # concatenate 'makeup'+'cand' to 'fore'
                    ## still record in 'fore_mapInfo_Href'
                    $fore_mapInfo_Href->{score} += $cand_mapInfo_Href->{score};
                    $fore_mapInfo_Href->{ident} = (   $fore_mapInfo_Href->{ident}*$fore_map_SpanL
                                                    + $cand_mapInfo_Href->{ident}*$cand_map_SpanL
                                                  ) / ($fore_map_SpanL+$cand_map_SpanL);
                    $fore_mapInfo_Href->{map_info} .= ($makeup_MapI . $cand_mapInfo_Href->{map_info});
                    $fore_mapInfo_Href->{query_info}->[1] .= lc($makeup_Qseq . $cand_mapInfo_Href->{query_info}->[1]);
                    $fore_mapInfo_Href->{query_info}->[2] = $cand_mapInfo_Href->{query_info}->[2];
                    $fore_mapInfo_Href->{sbjct_info}->[1] .= lc($makeup_Sseq . $cand_mapInfo_Href->{sbjct_info}->[1]);
                    $fore_mapInfo_Href->{sbjct_info}->[2] = $cand_mapInfo_Href->{sbjct_info}->[2];
                    # inform
                    stout_and_sterr  "[INFO]:\tConnect two blast mapping blocks (large DEL/uks-INS) of query_id($query_id):\n"
                                    ."\tfore_mapping interval, idx=$fore_map_idx:\n"
                                    ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_map_QedP;"
                                    .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_map_SedP\n"
                                    ."\ttail_mapping interval, idx=$cand_map_idx:\n"
                                    ."\t Qst:$cand_map_QstP; Qed:$cand_mapInfo_Href->{query_info}->[2];"
                                    .  " Sst:$cand_map_SstP; Sed:$cand_mapInfo_Href->{sbjct_info}->[2]\n"
                                    ."\tnew_mapping interval, idx=$fore_map_idx:\n"
                                    ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_mapInfo_Href->{query_info}->[2];"
                                    .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_mapInfo_Href->{sbjct_info}->[2]\n"
                                    # debug, show the new connected mapping block
                                    # ."\t '$fore_mapInfo_Href->{query_info}->[1]'\n"
                                    # ."\t '$fore_mapInfo_Href->{map_info}'\n"
                                    # ."\t '$fore_mapInfo_Href->{sbjct_info}->[1]'\n"
                                    ."\tmakeup_mapInfo:\n"
                                    ."\t Q: '$makeup_Qseq'\n"
                                    ."\t M: '$makeup_MapI'\n"
                                    ."\t S: '$makeup_Sseq'\n";
                    # record cand_map_idx as old_discarded
                    $discard_idx{$cand_map_idx} = 1;
                    # reocrd update sign
                    $update_sign = 1;
                    # stop to start the whole new detection loop
                    last;
                }
                ## large INS: 'none' OR 'slight' overlap on query_seq, and 'large' overlap on sbjct_seq
                ## large INS is tandem duplicated, Here
                elsif(    $fore_map_Qstrd * $fore_map_QedP >= $cand_map_Qstrd * $cand_map_QstP - 1
                       && $fore_map_Sstrd * $fore_map_SedP >  $cand_map_Sstrd * $cand_map_SstP
                       &&  # query_seq has large insertion (tandem duplicated)
                          (   $map_Qdist <= 10 # empirical distance in practise
                           && $map_Sdist >= 10 # $map_Qdist
                           && $map_Sdist <= 50 # it is the max ins-size detected from traditional PE seq
                           && $map_Sdist <= $fore_map_SpanL / 5 # fore maps long enough, empirical value
                           && $map_Sdist <= $cand_map_SpanL / 5 # cand maps long enough, empirical value
                          )
                ){
                    # 1st, discard possible extended end of cand query_seq to make it continuous
                    my $QextLen = ( $fore_map_QedP - $cand_map_QstP ) * $fore_map_Qstrd + 1;
                    # modify the candidate info
                    my @cand_Qseq = split //, $cand_mapInfo_Href->{query_info}->[1];
                    my $cand_Qidx = -1;
                    while( $QextLen > 0 ){
                        $cand_Qidx ++;
                        if( $cand_Qseq[$cand_Qidx] ne '-' ){
                            $QextLen --;
                        }
                    }
                    if( $cand_Qidx >= 0 ){
                        my $DiscLen = $cand_Qidx + 1;
                        # map info
                        $cand_mapInfo_Href->{map_info} = substr( $cand_mapInfo_Href->{map_info}, $DiscLen );
                        # query_seq
                        $cand_mapInfo_Href->{query_info}->[1] = substr( $cand_mapInfo_Href->{query_info}->[1], $DiscLen );
                        # sbjct_seq
                        my $Slen = grep !/-/, ( split //, substr( $cand_mapInfo_Href->{sbjct_info}->[1], 0, $DiscLen ) );
                        $cand_mapInfo_Href->{sbjct_info}->[1] = substr( $cand_mapInfo_Href->{sbjct_info}->[1], $DiscLen );
                        $cand_map_SstP += $Slen * $fore_map_Sstrd;
                    }
                    # 2nd, introduce insertion to make subject_seq continuous
                    my $SextLen = ( $fore_map_SedP - $cand_map_SstP ) * $fore_map_Sstrd + 1;
                    # modify the candidate info
                    my @cand_Sseq = split //, $cand_mapInfo_Href->{sbjct_info}->[1];
                    my $cand_Sidx = -1;
                    while( $SextLen > 0 ){
                        $cand_Sidx ++;
                        if( $cand_Sseq[$cand_Sidx] ne '-' ){
                            $SextLen --;
                        }
                    }
                    if( $cand_Sidx >= 0 ){ # must be
                        my $ReplcLen = $cand_Sidx + 1;
                        # map info
                        $cand_mapInfo_Href->{map_info}        = join('', (' ' x abs($ReplcLen))) . substr($cand_mapInfo_Href->{map_info},        $ReplcLen);
                        # sbjct_seq
                        $cand_mapInfo_Href->{sbjct_info}->[1] = join('', ('-' x abs($ReplcLen))) . substr($cand_mapInfo_Href->{sbjct_info}->[1], $ReplcLen);
                    }
                    # concatenate 'makeup'+'cand' to 'fore'
                    ## still record in 'fore_mapInfo_Href'
                    $fore_mapInfo_Href->{score} += $cand_mapInfo_Href->{score};
                    $fore_mapInfo_Href->{ident} = ($fore_mapInfo_Href->{ident}+$cand_mapInfo_Href->{ident}) / 2;
                    $fore_mapInfo_Href->{map_info} .= $cand_mapInfo_Href->{map_info};
                    $fore_mapInfo_Href->{query_info}->[1] .= $cand_mapInfo_Href->{query_info}->[1];
                    $fore_mapInfo_Href->{query_info}->[2]  = $cand_mapInfo_Href->{query_info}->[2];
                    $fore_mapInfo_Href->{sbjct_info}->[1] .= $cand_mapInfo_Href->{sbjct_info}->[1];
                    $fore_mapInfo_Href->{sbjct_info}->[2]  = $cand_mapInfo_Href->{sbjct_info}->[2];
                    # inform
                    stout_and_sterr  "[INFO]:\tConnect two blast mapping blocks (large tdp-INS) of query_id($query_id):\n"
                                    ."\tfore_mapping interval, idx=$fore_map_idx:\n"
                                    ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_map_QedP;"
                                    .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_map_SedP\n"
                                    ."\ttail_mapping interval, idx=$cand_map_idx:\n"
                                    ."\t Qst:$cand_mapInfo_Href->{query_info}->[0]; Qed:$cand_mapInfo_Href->{query_info}->[2];"
                                    .  " Sst:$cand_mapInfo_Href->{sbjct_info}->[0]; Sed:$cand_mapInfo_Href->{sbjct_info}->[2]\n"
                                    ."\tnew_mapping interval, idx=$fore_map_idx:\n"
                                    ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_mapInfo_Href->{query_info}->[2];"
                                    .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_mapInfo_Href->{sbjct_info}->[2]\n";
                                    # debug, show the new connected mapping block
                                    # ."\t '$fore_mapInfo_Href->{query_info}->[1]'\n"
                                    # ."\t '$fore_mapInfo_Href->{map_info}'\n"
                                    # ."\t '$fore_mapInfo_Href->{sbjct_info}->[1]'\n"
                    # record cand_map_idx as old_discarded
                    $discard_idx{$cand_map_idx} = 1;
                    # reocrd update sign
                    $update_sign = 1;
                    # stop to start the whole new detection loop
                    last;
                }
            }
            # stop to start the whole new detection loop
            if( $update_sign ){
                last;
            }
        }
        # stop when no update operation is done
        unless( $update_sign ){
            last;
        }
    }
    # discard any recorded idx
    if( scalar(keys %discard_idx) != 0 ){
        my @keep_idx;
        for my $idx ( 0 .. scalar(@$MapInfo_Aref)-1 ){
            unless( exists($discard_idx{$idx}) ){
                push @keep_idx, $MapInfo_Aref->[$idx];
            }
        }
        # re-assign new mapping result array-ref to this query_id
        $mapIF_Hf->{$query_id} = \@keep_idx;
    }
}

#--- summary mapIF of given query_id to refpos  ---
sub mapIF2RefPos{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $mapIF_Hf = $parm{mapIF_Hf};
    my $refposHf = $parm{refposHf}; # refpos hash ref
    my $query_id = $parm{query_id};
    my $firstMap = $parm{firstMap} || 0; # only use the first map

    my $MapInfo_Aref = $mapIF_Hf->{$query_id};
    for my $map_idx ( 0 .. scalar(@$MapInfo_Aref)-1 ){
        last if $firstMap && $map_idx != 0; # only first alignment?
        my $QueryMapIF_Hf = $MapInfo_Aref->[$map_idx];
        # information
        my $sbjct_id   =    $QueryMapIF_Hf->{sbjct_id};
        my $sbjct_strd =    $QueryMapIF_Hf->{sbjct_strd};
        my $align_strd =    $QueryMapIF_Hf->{align_strd};
        my $align_info =    $QueryMapIF_Hf->{map_info};
        my @query_info = @{ $QueryMapIF_Hf->{query_info} };
        my @sbjct_info = @{ $QueryMapIF_Hf->{sbjct_info} };
        # make the subject plus orientated
        if($sbjct_strd == -1){ # Minus
            @sbjct_info = reverse @sbjct_info;
            ($sbjct_info[1] = reverse $sbjct_info[1]) =~ tr/acgtACGT/tgcaTGCA/;
            @query_info = reverse @query_info;
            ($query_info[1] = reverse $query_info[1]) =~ tr/acgtACGT/tgcaTGCA/;
            $align_info = reverse $align_info;
            # check
            if($sbjct_info[2] < $sbjct_info[0]){
                cluck_and_exit "<ERROR>\twrong sbjct info:$sbjct_info[2] < $sbjct_info[0] for query_id $query_id\n".Dumper($QueryMapIF_Hf);
            }
        }
        # details
        my @align_info = split //, $align_info;
        my @query_seq  = split //, $query_info[1];
        my @sbjct_seq  = split //, $sbjct_info[1];
        my $sbjct_st = $sbjct_info[0];
        my $sbjct_ed = $sbjct_info[2];
        my $last_sbjct_pos = $sbjct_st - 1;
        # record to refpos
        ## judge each base, and operate according to the mutation type
        my $sjRefPosHf = $refposHf->{$sbjct_id};
        for (my $i=0; $i<=$#align_info; $i++){
            my $map_sign = $align_info[$i];
            my $query_allel = $query_seq[$i];
            my $sbjct_allel = $sbjct_seq[$i];
            # mutation types
            if($map_sign eq '|'){ # same
                my $sbjct_pos = $last_sbjct_pos + 1;
                my $sjRefPosOB = $sjRefPosHf->{$sbjct_pos};
                # accumulate whole depth
                $sjRefPosOB->addDepth(add=>1);
                $sjRefPosOB->addRefDepth(add_fw=>1) if $align_strd > 0;
                $sjRefPosOB->addRefDepth(add_rv=>1) if $align_strd < 0;
                # update last sbjct_pos
                $last_sbjct_pos = $sbjct_pos;
            }
            elsif($query_allel ne '-' && $sbjct_allel ne '-'){ # mismatch
                my $sbjct_pos = $last_sbjct_pos + 1;
                my $sjRefPosOB = $sjRefPosHf->{$sbjct_pos};
                #check
                if($query_allel eq $sbjct_allel){
                    cluck_and_exit "<ERROR>\tShould be a Mismatch: $sbjct_pos\n"
                                         ."\tQ: $query_allel\n"
                                         ."\tS: $sbjct_allel\n";
                }
                # accumulate whole depth
                $sjRefPosOB->addDepth(add=>1);
                # record mutation
                my $mut_id = 'snp,' . uc($query_allel);
                # FORMAT (Sum, Plus, Minus, RefAllelDepth[deprecated])
                $sjRefPosOB->addMut(mut_id=>$mut_id, depthAf=>[0,0,0,0]) unless $sjRefPosOB->has_mutation(mut_id=>$mut_id);
                # accumulate mut depth
                $sjRefPosOB->addMutDepth(mut_id=>$mut_id, add_fw=>1) if $align_strd > 0;
                $sjRefPosOB->addMutDepth(mut_id=>$mut_id, add_rv=>1) if $align_strd < 0;
                # update last sbjct_pos
                $last_sbjct_pos = $sbjct_pos;
            }
            elsif($sbjct_allel eq '-'){ # insertion
                my $sbjct_pos = $last_sbjct_pos; # insertion is assigned to the base before it.
                my $sjRefPosOB = $sjRefPosHf->{$sbjct_pos};
                # find the inserted part from query sequence
                my $j = 1;
                while($sbjct_seq[$i+$j] eq '-'){ # insertion extend
                    $j ++;
                }
                $j --; # adjust
                my $insert_seq = join('',@query_seq[$i .. $i+$j]);
                $i += $j; # jump to the last base of the insertion
                # do not need to accumulate whole depth for this pos 'again', because just added before ('#same' case above)
                # $sjRefPosOB->addDepth(add=>1);
                # decrease one from 'refallel_depth' of this pos, because just added before ('#same' case above)
                $sjRefPosOB->addRefDepth(add_fw=>-1) if $align_strd > 0;
                $sjRefPosOB->addRefDepth(add_rv=>-1) if $align_strd < 0;
                # record mutation
                my $mut_id = 'ins,' . uc($insert_seq);
                # FORMAT (Sum, Plus, Minus, RefAllelDepth[deprecated])
                $sjRefPosOB->addMut(mut_id=>$mut_id, depthAf=>[0,0,0,0]) unless $sjRefPosOB->has_mutation(mut_id=>$mut_id);
                # accumulate mut depth
                $sjRefPosOB->addMutDepth(mut_id=>$mut_id, add_fw=>1) if $align_strd > 0;
                $sjRefPosOB->addMutDepth(mut_id=>$mut_id, add_rv=>1) if $align_strd < 0;
                # update last sbjct_pos
                $last_sbjct_pos = $sbjct_pos;
            }
            elsif($query_allel eq '-'){ # deletion
                my $sbjct_pos = $last_sbjct_pos + 1; # deletion is assigned to the first base of itself.
                my $sjRefPosOB = $sjRefPosHf->{$sbjct_pos};
                my $j = 1;
                while($query_seq[$i+$j] eq '-'){ # insertion extend
                    $j ++;
                }
                $j --; # adjust
                my $delete_seq = join('',@sbjct_seq[$i .. $i+$j]);
                $i += $j; # jump to the last base of the deletion
                # accumulate whole depth, it is just used to measure the reads' supportings of this pos
                $sjRefPosHf->{$_}->addDepth(add=>1) for ($sbjct_pos .. $sbjct_pos+$j);
                # record mutation
                my $mut_id = 'del,' . uc($delete_seq);
                # FORMAT (Sum, Plus, Minus, RefAllelDepth[deprecated])
                $sjRefPosOB->addMut(mut_id=>$mut_id, depthAf=>[0,0,0,0]) unless $sjRefPosOB->has_mutation(mut_id=>$mut_id);
                # accumulate mut depth
                if($align_strd > 0){
                    $sjRefPosOB->addMutDepth(mut_id=>$mut_id, add_fw=>1);
                }
                else{
                    $sjRefPosOB->addMutDepth(mut_id=>$mut_id, add_rv=>1);
                }
                # update last sbjct_pos
                $last_sbjct_pos = $sbjct_pos + $j;
            }
            else{
                cluck_and_exit "<ERROR>\tWrong align info:\n$query_info[1]\n$align_info\n$sbjct_info[1]\n".Dumper($QueryMapIF_Hf);
            }
        }
        # last check
        if($last_sbjct_pos != $sbjct_ed){
            cluck_and_exit "<ERROR>\tThe last_sbjct_pos:$last_sbjct_pos is not the sbjct_ed:$sbjct_ed\n".Dumper($QueryMapIF_Hf);
        }
    }
}

1; ## tell the perl script the successful access of this module.
