\name{DECIPHER-package}
\alias{DECIPHER-package}
\alias{DECIPHER}
\docType{package}
\title{
Tools for curating, analyzing, and manipulating biological sequences
}
\description{
DECIPHER is a software toolset that can be used for deciphering and managing biological sequences efficiently using the R statistical programming language.  The program is designed to be used with non-destructive workflows for importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.
}
\details{
\tabular{ll}{
Package: \tab DECIPHER\cr
Type: \tab Package\cr
Depends: \tab R (>= 3.5.0), Biostrings (>= 2.59.1), stats, parallel\cr
Imports: \tab methods, DBI, S4Vectors, IRanges, XVector\cr
Suggests: \tab RSQLite (>= 1.1)\cr
LinkingTo: \tab Biostrings, S4Vectors, IRanges, XVector\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
AA_REDUCED               Reduced amino acid alphabets
Add2DB                   Add Data to a Database
AdjustAlignment          Improve An Existing Alignment By Adjusting Gap
                         Placements
AlignDB                  Align Two Sets of Aligned Sequences in a Sequence
                         Database
AlignPairs               Align pairs of sequences
AlignProfiles            Align Two Sets of Aligned Sequences
AlignSeqs                Align a Set of Unaligned Sequences
AlignSynteny             Pairwise Aligns Syntenic Blocks
AlignTranslation         Align Sequences By Their Amino Acid Translation
AmplifyDNA               Simulate Amplification of DNA by PCR
Array2Matrix             Create a Matrix Representation of a Microarray
BLOSUM                   BLOSUM Amino Acid Substitution Matrices
BrowseDB                 View a Database Table in a Web Browser
BrowseSeqs               View Sequences in a Web Browser
CalculateEfficiencyArray Predict the Hybridization Efficiency of
                         Probe/Target Sequence Pairs
CalculateEfficiencyFISH  Predict Thermodynamic Parameters of Probe/Target
                         Sequence Pairs
CalculateEfficiencyPCR   Predict Amplification Efficiency of Primer Sequences
Clusterize               Cluster Sequences By Distance
Codec                    Compression/Decompression of Character Vectors
ConsensusSequence        Create a Consensus Sequence
Cophenetic               Compute cophenetic distances on dendrogram objects
CorrectFrameshifts       Corrects Frameshift Errors In Protein Coding
                         Sequences
CreateChimeras           Create Artificial Chimeras
DB2Seqs                  Export Database Sequences to a FASTA or FASTQ File
deltaGrules              Free Energy of Hybridization of Probe/Target
                         Quadruplets
deltaHrules              Change in Enthalpy of Hybridization of DNA/DNA
                         Quadruplets in Solution
deltaHrulesRNA           Change in Enthalpy of Hybridization of RNA/RNA
                         Quadruplets in Solution
deltaSrules              Change in Entropy of Hybridization of DNA/DNA
                         Quadruplets in Solution
deltaSrulesRNA           Change in Entropy of Hybridization of RNA/RNA
                         Quadruplets in Solution
DesignArray              Design a Set of DNA Microarray Probes for Detecting
                         Sequences
DesignPrimers            Design Primers Targeting a Specific Group of
                         Sequences
DesignProbes             Design FISH Probes Targeting a Specific Group of
                         Sequences
DesignSignatures         Design PCR Primers for Amplifying Group-Specific
                         Signatures
DetectRepeats            Detect Repeats in a Sequence
DigestDNA                Simulate Restriction Digestion of DNA
Disambiguate             Expand Ambiguities into All Permutations of a
                         DNAStringSet
DistanceMatrix           Calculate the Distance Between Sequences
ExtractGenes             Extract Predicted Genes from a Genome
FindChimeras             Find Chimeras in a Sequence Database
FindGenes                Find Genes in a Genome
FindNonCoding            Find Non-Coding RNAs in a Genome
FindSynteny              Finds Synteny in a Sequence Database
FormGroups               Forms Groups By Rank
Genes-class              Genes objects and accessors
HEC_MI                   Mutual Information for Protein Secondary Structure
                         Prediction
IdConsensus              Create Consensus Sequences by Groups
IdentifyByRank           Identify By Taxonomic Rank
IdLengths                Determine the Number of Characters in Each Sequence
                         of Each Sequence
IdTaxa                   Assign Sequences a Taxonomic Classification
IndexSeqs                Build an inverted index
InvertedIndex-class      InvertedIndex objects
LearnNonCoding           Learn a Non-Coding RNA Model
LearnTaxa                Train a Classifier for Assigning Taxonomy
MapCharacters            Map Changes in Ancestral Character States
MaskAlignment            Mask Highly Variable Regions of An Alignment
MeltDNA                  Simulate Melting of DNA
MIQS                     MIQS Amino Acid Substitution Matrix
MMLSUM                   MMLSUM Amino Acid Substitution Matrices
MODELS                   Available Models of Sequence Evolution
NNLS                     Sequential Coordinate-wise Algorithm for the
                         Non-negative Least Squares Problem
NonCoding                NonCoding Models for Common Non-Coding RNA Families
NonCoding-class          NonCoding Objects and Methods
OrientNucleotides        Orient Nucleotide Sequences
PAM                      PAM Amino Acid Substitution Matrices
PFASUM                   PFASUM Amino Acid Substitution Matrices
PredictDBN               Predict RNA Secondary Structure in Dot-Bracket
                         Notation
PredictHEC               Predict Protein Secondary Structure as Helix,
                         Beta-Sheet, or Coil
Read Dendrogram          Read a Dendrogram from a Newick Formatted File
RemoveGaps               Remove Gap Characters in Sequences
RESTRICTION_ENZYMES      Common Restriction Enzyme's Cut Sites
ScoreAlignment           Score a Multiple Sequence Alignment
SearchDB                 Obtain Specific Sequences from a Database
SearchIndex              Search an inverted index
Seqs2DB                  Add Sequences from Text File to Database
StaggerAlignment         Produce a Staggered Alignment
Synteny-class            Synteny blocks and hits
Taxa-class               Taxa training and testing objects
TerminalChar             Determine the Number of Terminal Characters
TileSeqs                 Form a Set of Tiles for Each Group of Sequences
TrainingSet_16S          Training Set for Classification of 16S rRNA Gene
                         Sequences
TreeLine                 Construct a Phylogenetic Tree
TrimDNA                  Trims DNA Sequences to the High Quality Region
                         Between Patterns
WriteDendrogram          Write a Dendrogram to Newick Format
WriteGenes               Write Genes to a File
}
}
\author{
Erik Wright

Maintainer: Erik Wright <eswright@pitt.edu>
}
\keyword{ package }
