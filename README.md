# BioFuse

Util PERL module of general functions and objects applied in bioinformatics software development.

- Author: Wenlong Jia
- Email:  wenlongkxm@gmail.com

## Check.pm
### BioFuse::Check
### VERSION = "0.01"
- check

## Util/Web.pm
### BioFuse::Util::Web
### VERSION = "0.01"
- urlToHtmlText

## Util/Log.pm
### BioFuse::Util::Log
### VERSION = "0.32"
- warn_and_exit
- cluck_and_exit
- stout_and_sterr

## Util/Index.pm
### BioFuse::Util::Index
### VERSION = "0.05"
- Pos2Idx
- IndexRegion
- FindOverlapIdxRegion

## Util/Interval.pm
### BioFuse::Util::Interval
### VERSION = "0.34"
- Get_Two_Seg_Olen
- merge
- intersect
- exclude
- arrange_region_with_clip
- deal_circular_extended_part
- clip_region_edge

## Util/Sort.pm
### BioFuse::Util::Sort
### VERSION = "0.01"
- sortByStrAndSubNum

## Util/FileHeader.pm
### BioFuse::Util::FileHeader
### VERSION = "0.01"
- getHeaderTag
- lineInfoToHash

## Util/Random.pm
### BioFuse::Util::Random
### VERSION = "0.31"
- GetRandBool
- PickRandAele

## Util/GZfile.pm
### BioFuse::Util::GZfile
### VERSION = "0.08"
- Try_GZ_Read
- is_idx_bgz
- Try_GZ_Write

## Util/Array.pm
### BioFuse::Util::Array
### VERSION = "0.33"
- binarySearch

## Util/String.pm
### BioFuse::Util::String
### VERSION = "0.33"
- getStrRepUnit
- getStrUnitRepeatTime
- getRegexRegion

## Util/Sys.pm
### BioFuse::Util::Sys
### VERSION = "0.34"
- file_exist
- trible_run_for_success
- check_java_version
- reset_folder

## Dist/DistStat.pm
### BioFuse::Dist::DistStat
### VERSION = "0.02"
- get_value_mean
- get_trimmed_mean
- engineer_Ntimes_SD_evaluation

## GetPath.pm
### BioFuse::GetPath
### VERSION = "0.01"
- GetPath

## BioInfo/Depth/MixToPureT.pm
### BioFuse::BioInfo::Depth::MixToPureT
### VERSION = "0.02"
- get_Tcell_GMpart
- get_Tcell_GMratio
- get_Ncell_GMpart
- get_Ncell_GMratio
- get_ObjCountOfPureTumorCell
- get_ObjSingleCNdepthInMixed

## BioInfo/GeneAnno/GTFtoGenePSL.pm
### BioFuse::BioInfo::GeneAnno::GTFtoGenePSL
### VERSION = "0.83"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- GTFtoGenePSL

## BioInfo/GeneAnno/PSLtoBED.pm
### BioFuse::BioInfo::GeneAnno::PSLtoBED
### VERSION = "0.03"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- PSLtoBED
- output_region_of_transOB
- record_protein_coding_regions
- record_exon_regions

## BioInfo/GeneAnno/PSLtoFASTA.pm
### BioFuse::BioInfo::GeneAnno::PSLtoFASTA
### VERSION = "0.85"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- PSLtoFASTA

## BioInfo/GeneAnno/GTFtoTransPSL.pm
### BioFuse::BioInfo::GeneAnno::GTFtoTransPSL
### VERSION = "0.83"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- GTFtoTransPSL

## BioInfo/Objects/SeqData/PairEnd_OB.pm
### BioFuse::BioInfo::Objects::SeqData::PairEnd_OB
### VERSION = "0.09"
- new
- load_reads_OB
- pid
- peIdx
- rOB_Af
- sorted_rOB_Af
- tryDiscardAlign
- test_need_RefSeg
- test_pair_RefSeg
- onlyKeep_need_RefSeg
- makePrimeAlignment
- discardAbnormalSP
- printSAM

## BioInfo/Objects/SeqData/HicReads_OB.pm
### BioFuse::BioInfo::Objects::SeqData::HicReads_OB
### VERSION = "0.03"
- new
- load_AlignJudge
- load_SuppHaplo
- del_SuppHaplo
- onlyKeep_SuppHaplo
- get_AlignJudge
- get_SuppHaploHref
- get_SuppHaploStr
- has_SuppHaplo
- is_fromUnPhasedRegRand
- addHapIDtoOptfd
- recover_SuppHaploAttr

## BioInfo/Objects/SeqData/Bam_OB.pm
### BioFuse::BioInfo::Objects::SeqData::Bam_OB
### VERSION = "0.16"
- new
- verify_bam
- verify_index
- delete_file
- addTool
- filepath
- tag
- tissue
- rgOB_Hf
- header_Af
- isNsort
- toNsort
- isCsort
- toCsort
- toIndex
- toFixmate
- merge
- toMarkdup
- start_read
- start_write
- stop_write
- write
- load_reads_for_ReadsGroup
- rg_count_need_reads_ForIns
- pick_rgOB
- add_rgOB
- get_regionCovStat
- get_region_depth
- delete_regionDepthFile
- get_region_alt_vcf_gz
- get_pos_marker_stat
- get_allele_marker_stat
- smartBam_PEread

## BioInfo/Objects/SeqData/HicPairEnd_OB.pm
### BioFuse::BioInfo::Objects::SeqData::HicPairEnd_OB
### VERSION = "0.08"
- new
- get_rEndWholeAlignJudge
- get_rEndWholeSuppHaplo
- isInValidPair
- testLinkRefSeg
- dEndSameHapJudge
- sEndSoloHapJudge
- sEndInterHapJudge
- dEndInterHapJudge
- addHapIDtoReadsOptfd

## BioInfo/Objects/SeqData/FastQ_OB.pm
### BioFuse::BioInfo::Objects::SeqData::FastQ_OB
### VERSION = "0.01"
- new
- filepath
- tag
- rgOB_Hf
- start_read
- start_write
- stop_write
- write
- add_rgOB

## BioInfo/Objects/SeqData/Reads_OB.pm
### BioFuse::BioInfo::Objects::SeqData::Reads_OB
### VERSION = "0.16"
- new
- pid
- endNO
- mseg
- mpos
- mapQ
- rlen
- mReadLen
- mRefLen
- lenFromCigar
- barcode
- add_BarcToPid
- optfd_str
- add_str_to_optfd
- optfd_has_regex
- judgeAlign
- is_proper_map
- is_fw_map
- is_rv_map
- is_unmap
- is_2ndmap
- free_2ndmap
- is_suppmap
- free_suppmap
- is_dup
- is_mltmap
- is_good_cigar
- is_softclip
- is_hardclip
- is_clip
- has_MDtag
- is_closeAlign
- foreClipLen
- hindClipLen
- biClipLen
- tlen_FixTlen_with_S
- digestMDtag
- fuseCigarMD
- getNearAltDist
- get_pos_allele
- printSAM
- printFQ
- find_rgOB
- add_RGrIDprefToPid
- update_rgOB_maxRlen
- update_rid_RGprefix
- get_pos_allele_v1_BaseOnCigar

## BioInfo/Objects/SeqData/ReadsGroup_OB.pm
### BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB
### VERSION = "0.05"
- new
- RG_ID
- LB_ID
- RG_NO
- set_RG_NO
- rID_prefix
- set_rID_prefix
- maxRlen_Hf
- update_maxRlen
- affix
- addAffix
- bam
- RG_para
- load_reads_for_ins_evalue
- evalue_ins
- test_3p_overlap
- file_prefix
- report_format
- write_report
- load_report
- write_insDistLog

## BioInfo/Objects/GeneAnno/PSL_OB.pm
### BioFuse::BioInfo::Objects::GeneAnno::PSL_OB
### VERSION = "0.05"
- new
- load_GeneOrTrans_from_PSL
- get_objPOOLHf
- extract_GeneOrTrans_seq
- output_exon_seq

## BioInfo/Objects/GeneAnno/Trans_OB.pm
### BioFuse::BioInfo::Objects::GeneAnno::Trans_OB
### VERSION = "0.02"
- new
- get_ENSid
- get_ref_seg
- get_strand
- get_use_name
- get_ori_name
- get_gene_use_name
- get_biotype
- get_exon_region
- get_exon_sumLen
- get_CDS_region
- get_UTR_region
- get_find_seq_mark
- mark_find_seq

## BioInfo/Objects/GeneAnno/GTF_OB.pm
### BioFuse::BioInfo::Objects::GeneAnno::GTF_OB
### VERSION = "0.08"
- new
- load_GTF
- read_refseg_transform
- refine_GTF_info
- refine_names
- add_refseg_cytoband
- create_gene_PSL
- create_trans_PSL
- mark_abnormal_Start_codon
- check_Start_codon_Seq

## BioInfo/Objects/GeneAnno/Gene_OB.pm
### BioFuse::BioInfo::Objects::GeneAnno::Gene_OB
### VERSION = "0.02"
- new
- get_ENSid
- get_ref_seg
- get_strand
- get_use_name
- get_ori_name
- get_biotype
- get_exon_region
- get_find_seq_mark
- mark_find_seq

## BioInfo/Objects/GeneAnno/GTF_transOB.pm
### BioFuse::BioInfo::Objects::GeneAnno::GTF_transOB
### VERSION = "0.01"
- new
- load_exon_region
- load_CDS_region
- load_codon_info
- add_cytoband_info
- get_ENSid
- get_ref_seg
- get_strand
- get_biotype
- get_trans_use_name
- get_trans_ori_name
- get_edge
- get_startCodonPos
- get_stopCodonPos
- get_exon_Aref
- get_CDS_Aref
- update_biotype
- note_abnor_codon_num
- refine_exon
- refined_use_name
- get_trans_psl_line

## BioInfo/Objects/GeneAnno/GTF_geneOB.pm
### BioFuse::BioInfo::Objects::GeneAnno::GTF_geneOB
### VERSION = "0.07"
- new
- load_gtf_info
- add_cytoband_info
- get_ENSid
- get_ref_seg
- get_strand
- get_gene_use_name
- get_gene_ori_name
- get_trans_use_name
- get_trans_ori_name
- get_transOB_Aref
- refine_gtf_info
- refined_use_name
- get_gene_psl_line

## BioInfo/Objects/GeneAnno/GTF_lineOB.pm
### BioFuse::BioInfo::Objects::GeneAnno::GTF_lineOB
### VERSION = "0.01"
- new
- get_refSeg
- get_strand
- get_regionType
- get_stpos
- get_edpos
- get_ENSid
- get_name
- get_source
- get_version
- get_biotype
- get_transLocusTag
- get_transNote
- get_transProduct
- get_transExonNO
- get_proteinENSid
- update_refSeg
- match_gtf_source

## BioInfo/Objects/Segment/RefSeg_OB.pm
### BioFuse::BioInfo::Objects::Segment::RefSeg_OB
### VERSION = "0.01"
- new
- id
- set_length
- length
- set_note
- note

## BioInfo/Objects/Allele/AlleleOnReads_OB.pm
### BioFuse::BioInfo::Objects::Allele::AlleleOnReads_OB
### VERSION = "0.03"
- new
- setMiss
- loadInfo
- get_chr
- get_pos
- get_refBase
- get_rOBmLen
- get_type
- get_score
- get_qual
- get_rEdgeDist5
- get_rEdgeDist3
- get_nAltDist5
- get_nAltDist3
- get_alleleSeq
- get_delSize
- get_deloffset

## BioInfo/Objects/Allele/PhasedMut_OB.pm
### BioFuse::BioInfo::Objects::Allele::PhasedMut_OB
### VERSION = "0.10"
- new
- load_posIdx
- load_ref_allele
- load_haplotype_allele
- load_judgeReadEdgeDist
- get_infoSummary
- get_chr
- get_pos
- get_refAllele
- get_alleleTypeToHapAref
- get_hapInfo
- get_phMutNO
- get_judgeReadEdgeDist
- is_hom
- has_alleleType
- has_ref
- has_snv
- has_del
- has_ins
- has_indel
- allele2haploID
- release_memory

## BioInfo/Objects/Sample/Sample_OB.pm
### BioFuse::BioInfo::Objects::Sample::Sample_OB
### VERSION = "0.02"
- new
- id
- setType
- type
- setPatient
- patient
- setTissue
- tissue
- set_note
- note
- addBam
- bam
- addFQ

## BioInfo/Objects/Sample/Patient_OB.pm
### BioFuse::BioInfo::Objects::Sample::Patient_OB
### VERSION = "0.01"
- new
- set_note
- addSample

## BioInfo/FASTA/GetNonNBed.pm
### BioFuse::BioInfo::FASTA::GetNonNBed
### VERSION = "0.02"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- GetNonNBedfromFasta
- get_nonN_region_from_segseq

## BioInfo/FASTA/GetMutFA.pm
### BioFuse::BioInfo::FASTA::GetMutFA
### VERSION = "0.04"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- GetMutFasta
- load_mutations
- filer_mutations
- merge_mutations
- alertDiffCN
- get_mut_refseq
- add_mut_on_chr
- add_single_mut

## BioInfo/Codon.pm
### BioFuse::BioInfo::Codon
### VERSION = "0.31"
- Load_Codon

## BioInfo/CytoBand.pm
### BioFuse::BioInfo::CytoBand
### VERSION = "0.01"
- load_cytoband
- get_cytoband

## BioInfo/Depth.pm
### BioFuse::BioInfo::Depth
### VERSION = "0.05"
- get_windowSmoDepth
- deal_window_depth_info
- get_ctrl_copyR
- get_givenWinItvalMeanDepth

## BioInfo/Alignment/BWA.pm
### BioFuse::BioInfo::Alignment::BWA
### VERSION = "0.01"
- PEfqToSortBam

## BioInfo/Position.pm
### BioFuse::BioInfo::Position
### VERSION = "0.01"
- load_region_for_randPos
- get_random_pos

## BioInfo/Quality.pm
### BioFuse::BioInfo::Quality
### VERSION = "0.31"
- baseQ_char2score

## BioInfo/BED.pm
### BioFuse::BioInfo::BED
### VERSION = "0.02"
- read_bed_file

## BioInfo/FASTA.pm
### BioFuse::BioInfo::FASTA
### VERSION = "0.33"
- read_fasta_file
- write_fasta_file
- BWA_index_fasta
- Faidx_Dict_fasta
- read_Fai

## LoadOn.pm
### BioFuse::LoadOn
### VERSION = "0.50"
- load_variants_dict

## RunFunc.pm
### BioFuse::RunFunc
### VERSION = "0.51"
- options_alert_and_run
- func_run
- load_functions
- return_HELP_INFO
- HelpInfo_ExtractOpt
- para_alert

## Stat/SpearmanCorr.pm
### BioFuse::Stat::SpearmanCorr
### VERSION = "0.01"
- get_spearman_corr

## Stat/FisherTest.pm
### BioFuse::Stat::FisherTest
### VERSION = "0.01"
- fisher_P
- log_sum

## Stat/ConfInt.pm
### BioFuse::Stat::ConfInt
### VERSION = "0.01"
- CIL2Z

## Stat/MultiTest.pm
### BioFuse::Stat::MultiTest
### VERSION = "0.01"
- new
- load_test
- get_all_test
- get_FDR_sig_test
- get_p_adjust_method
- p_adjust
- set_P_sig_under_FDR

## Stat/ChiSquareTest/FourFoldTable.pm
### BioFuse::Stat::ChiSquareTest::FourFoldTable
### VERSION = "0.03"
- new
- get_odds_ratio
- get_risk_ratio
- get_ratio_diff
- get_theroy_value # a, b, c, d
- get_P_value
- get_chi_square
- get_method

## Stat/PearsonCorr.pm
### BioFuse::Stat::PearsonCorr
### VERSION = "0.01"
- get_pearson_corr

## Visual/SVG_Util/RadSysEle.pm
### BioFuse::Visual::SVG_Util::RadSysEle
### VERSION = "0.29"
- draw_circle_seg
- draw_a_sector

## Visual/SVG_Util/Color.pm
### BioFuse::Visual::SVG_Util::Color
### VERSION = "0.01"

## Visual/SVG_Util/SVGWork.pm
### BioFuse::Visual::SVG_Util::SVGWork
### VERSION = "0.01"
- initialize_SVG_obj
- output_SVG_file

## Visual/SVG_Util/Font.pm
### BioFuse::Visual::SVG_Util::Font
### VERSION = "0.30"
- confirm_transform_ratio
- show_text_in_line
- show_text_on_arc
- get_size_of_text_to_show
- decode_char_size_symbol

## Visual/SVG_Util/RectSysEle.pm
### BioFuse::Visual::SVG_Util::RectSysEle
### VERSION = "0.36"
- draw_a_parallelogram
- draw_a_triangle
- draw_a_arrow
- draw_a_ellipse

## Visual/SVG_Util/RadSys.pm
### BioFuse::Visual::SVG_Util::RadSys
### VERSION = "0.16"
- get_coordinate_on_circle
- normalize_radian
- draw_an_arc

## Visual/Objects/Axis.pm
### BioFuse::Visual::Objects::Axis
### VERSION = "0.02"
- new
- add_resol
- set_label
- set_stub
- add_stub
- set_tic
- set_resolItvGap
- validate_resol
- update_axisLen
- get_origX
- get_origY
- get_axisLen
- get_resol
- get_stub
- get_headAng
- get_headRad
- get_origValue
- valueToAxisDist
- valueToSVGcoord
- extendCoord
- draw_quick
- draw_axisBodyLine
- draw_stub
- draw_label

## Visual/Objects/BiAxis.pm
### BioFuse::Visual::Objects::BiAxis
### VERSION = "0.01"
- new
- calc_orig
- get_origX
- get_origY
- get_axis
- valueToSVGcoord
- draw

## Visual/Objects/Histogram.pm
### BioFuse::Visual::Objects::Histogram
### VERSION = "0.01"
- new
- load_data
- set_attr
- draw

## Visual/Objects/GradColor.pm
### BioFuse::Visual::Objects::GradColor
### VERSION = "0.01"
- new
- set_axis
- get_stRGB
- get_edRGB
- get_stValue
- get_edValue
- valueToRGB
- draw_quick
- draw_body

## BioFuse.pm
### BioFuse::BioFuse

