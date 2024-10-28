\name{MODELS}
\alias{MODELS}

\title{
Available Models of Sequence Evolution
}
\description{
The \code{MODELS} contains the models of sequence evolution that can be used by \code{Treeline}.
}
\usage{
MODELS
}
\details{
\code{MODELS} is a list of two elements: a character vector of (eight) nucleotide models and a character vector of (37) protein models.  All \code{MODELS} are time reversible.

Nucleotide models are described in order of increasing number of parameters as follows:

\bold{JC69 (Jukes and Cantor, 1969)}
The simplest substitution model that assumes equal base frequencies (1/4) and equal mutation rates.

\bold{K80 (Kimura, 1980)}
Assumes equal base frequencies, but distinguishes between the rate of transitions and transversions.

\bold{T92 (Tamura, 1992)}
In addition to distinguishing between transitions and transversions, a parameter is added to represent G+C content bias.

\bold{F81 (Felsenstein, 1981)}
Assumes equal mutation rates, but allows all bases to have different frequencies.

\bold{HKY85 (Hasegawa, Kishino and Yano, 1985)}
Distinguishes transitions from transversions and allows bases to have different frequencies.

\bold{TN93 (Tamura and Nei, 1993)}
Allows for unequal base frequencies and distinguishes between transversions and the two possible types of transitions (i.e., A <-> G & C <-> T).

\bold{SYM (Zharkikh, 1994)}
Equal base frequencies but all substitution rates are free parameters.

\bold{GTR (Tavare, 1986)}
The general time reversible model allowing for unequal base frequencies and substitution rates.

Protein models are described in the following publications:

AB (Mirsky, 2015), BLOSUM62 (Henikoff, 1992), cpREV (Adachi, 2000), cpREV64 (Zhong, 2010), Dayhoff (Dayhoff, 1978), DCMut-Dayhoff (Kosiol, 2005), DCMut-JTT (Kosiol, 2005), DEN (Le, 2018), FLAVI (Le, 2020), FLU (Dang, 2010), gcpREV (Cox, 2013), HIVb (Nickle, 2007), HIVw (Nickle, 2007), JTT (Jones, 1992), LG (Le, 2008), MtArt (Abascal, 2007), mtDeu (Le, 2017), mtInv (Le, 2017), mtMam (Yang, 1998), mtMet (Le, 2017), mtOrt (Chang, 2020), mtREV (Adachi, 1996), mtVer (Le, 2017), MtZoa (Rota-Stabelli, 2009), PMB (Veerassamy, 2003), Q.bird (Minh, 2021), Q.insect (Minh, 2021), Q.LG (Minh, 2021), Q.mammal (Minh, 2021), Q.pfam (Minh, 2021), Q.plant (Minh, 2021), Q.yeast (Minh, 2021), rtREV (Dimmic, 2002), stmtREV (Liu, 2014), VT (Muller, 2000), WAG (Whelan, 2001), WAGstar (Whelan, 2001).

\bold{+G (Yang, 1993)}
Specifying any model+G4 adds a single parameter to any of the above models to relax the assumption of equal rates among sites in the sequence.  The single parameter specifies the shape of the Gamma Distribution.  The continuous distribution is represented with 2-10 discrete rates and their respective probabilities as determined by equal bins or the Laguerre Quadrature method (Felsenstein, 2001).  For example, specifying a model+G8 would represent the continuous Gamma Distribution with eight rates and their associated probabilities.

\bold{+F}
Specifying any model+F uses empirical frequencies rather than optimized state frequencies.  This is only applicable for models having state frequencies with free parameters.
}
\references{
Abascal, F., Posada, D., and Zardoya, R. (2007) Molecular Biology and Evolution, \bold{24}, 1-5.

Adachi, J. and Hasegawa, M. (1996) Journal of Molecular Evolution, \bold{42}, 459-468. 

Adachi, J., Waddell, P., Martin, W., and Hasegawa, M. (2000) Journal of Molecular Evolution, \bold{50}, 348-358. 

Chang, H., Nie, Y., Zhang, N., Zhang, X., Sun, H., Mao, Y., Qiu, Z., and Huang, Y. (2020) BMC Ecology and Evolution, \bold{20}, 57.

Cox, C. and Foster, P. (2013) Molecular Phylogenetics and Evolution, \bold{68}, 218-220.

Dang, C., Le, S., Gascuel, O., and Le, V. (2010) BMC Evolutionary Biology, \bold{10}, 99. 

Dayhoff, M., Schwartz, R., and Orcutt, B. (1978) Atlas of Protein Sequence and Structure, National Biomedical Research Foundation, Washington DC, \bold{5}, 345-352.

Dimmic, M., Rest, J., Mindell, D., and Goldstein, R. (2002) Journal of Molecular Evolution, \bold{55}, 65-73. 

Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum likelihood approach. Journal of Molecular Evolution, \bold{17(6)}, 368-376.

Felsenstein, J. (2001) Taking Variation of Evolutionary Rates Between Sites into Account in Inferring Phylogenies. Journal of molecular evolution, \bold{53(4-5)}, 447-455.

Hasegawa, M., Kishino H., Yano T. (1985) Dating of human-ape splitting by a molecular clock of mitochondrial DNA. Journal of Molecular Evolution, \bold{22(2)}, 160-174.

Henikoff, S. and Henikoff, J. (1992) Proceedings of the National Academy of Sciences of the USA, \bold{89}, 10915-10919.

Jones, D., Taylor, W., and Thornton, J. (1992) Computer Applications in the Biosciences, \bold{8}, 275-282. 

Jukes, T. and Cantor C. (1969) Evolution of Protein Molecules. New York: Academic Press. pp. 21-132.

Kimura, M. (1980) A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution, \bold{16(2)}, 111-120.

Kosiol, C. and Goldman, N. (2005) Molecular Biology and Evolution, \bold{22}, 193-199. 

Le, S. and Gascuel, O. (2008) Molecular Biology and Evolution, \bold{25}, 1307-1320. 

Le, T., Dang, C., and Le, S. (2018) Proceedings of 10th International Conference on Knowledge and Systems Engineering (KSE 2018), Ho Chi Minh City, Vietnam, 242-246.

Le, T., and Vinh, L. (2020) Journal of Molecular Evolution, \bold{88}, 445-452.

Le, V., Dang, C., and Le, S. (2017) BMC Evolutionary Biology, \bold{17}, 136.

Liu, Y., Cox, C., Wang, W., and Goffinet, B. (2014) Systematic Biology, \bold{63}, 862-878. 

Minh, B., Dang, C., Le, S., and Lanfear, R. (2021) Systematic Biology, syab010.

Mirsky, A., Kazandjian, L., and Anisimova, M. (2015) Molecular Biology and Evolution, \bold{32}, 806-819.

Muller, T. and Vingron, M. (2000) Journal of Computational Biology, \bold{7}, 761-776. 

Nickle, D., Heath, L., Jensen, M., Gilbert P., and Mullins, J., Kosakovsky Pond SL (2007) PLoS ONE, \bold{2}, e503. 

Rota-Stabelli, O., Yang, Z., and Telford, M. (2009) Molecular Phylogenetics and Evolution, \bold{52}, 268-272. 

Tamura, K. (1992) Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C content biases. Molecular Biology and Evolution, \bold{9(4)}, 678-687.

Tamura, K. and Nei M. (1993) Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Molecular Biology and Evolution, \bold{10(3)}, 512-526.

Tavare, S. (1986) ``Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences.'' Lectures on Mathematics in the Life Sciences, \bold{17}: 57-86.

Veerassamy, S., Smith, A., and Tillier, E. (2003) Journal of Computational Biology, \bold{10}, 997-1010. 

Whelan, S. and Goldman, N. (2001) Molecular Biology and Evolution, \bold{18}, 691-699. 

Yang, Z., Nielsen, R., and Hasegawa, M. (1998) Molecular Biology and Evolution, \bold{15}, 1600-1611. 

Yang, Z. (1993) Maximum-likelihood estimation of phylogeny from DNA sequences when substitution rates differ over sites. Molecular Biology and Evolution, \bold{10(6)}, 1396-1401.

Zharkikh, A. (1994) Estimation of evolutionary distances between nucleotide sequences. Journal of Molecular Evolution, \bold{39}, 315-329.

Zhong, B., Yonezawa, T., Zhong, Y., and Hasegawa, M. (2010) Molecular Biology and Evolution, \bold{27}, 2855-2863.
}
\seealso{
\code{\link{Treeline}}
}
\examples{
str(MODELS)
}
\keyword{data}
