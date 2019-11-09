package BioFuse::BioFuse;

use strict;
use warnings;

# basic
use BioFuse::LoadOn;
# functions
use BioFuse::BioInfo::GeneAnno::GTFtoGenePSL;
use BioFuse::BioInfo::GeneAnno::GTFtoTransPSL;
use BioFuse::BioInfo::GeneAnno::PSLtoFASTA;
use BioFuse::BioInfo::GeneAnno::PSLtoBED;
use BioFuse::BioInfo::FASTA::GetNonNBed;
use BioFuse::BioInfo::FASTA::GetMutFA;
use BioFuse::BioInfo::Alignment::BamQC;

#--- 
1; ## tell the perl script the successful access of this module.
