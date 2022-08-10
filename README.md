WGCNA Script for: Monocyte, Neutrophil and Whole Blood Transcriptome
Dynamics Following Ischemic Stroke
================
07/28/2022

-----

**Monocyte, Neutrophil and Whole Blood Transcriptome Dynamics Following
Ischemic Stroke**  
Paulina Carmona-Mora, Bodie Knepp, Glen C Jickling, Xinhua Zhan, Marisa
Hakoupian, Heather Hull, Noor Alomar, Hajar Amini, Frank R Sharp,
Boryana Stamova, Bradley P Ander  
medRxiv 2022.03.03.22271866; doi:
<https://doi.org/10.1101/2022.03.03.22271866>

# Overview

This script was used to generate the Monocyte (MON), Neutrophil (NEU),
and Whole Blood (WB) Weighted Gene Co-Expression Network Analyses
(WGCNA) Networks for Carmona-Mora et al’s “Monocyte, Neutrophil and
Whole Blood Transcriptome Dynamics Following Ischemic Stroke”
publication. This study analyzed changes in the dyamics of the
peripheral blood transcriptome of human Ishcemic Stroke patients. If
this script is used, please cite the above paper.

# Running the Script

Study analyses were run using Microsoft R Open 4.0.2

By default, the script has the required parameters to recreate the MON
Network. Commented in are the parameters required to generate the NEU
and WB networks (only two places need modification: beta1 and kCut).
Various output files will also need name changes
(“DynamicsOfIS\_MON\_Network” to “DynamicsOfIS\_NEU\_Network”, etc).
Input data file name will need to be modified.

Testing for module relationships to Diagnosis and other clinical
parameters was conducted in another program - Partek Genomics Suite®.
Spearman Correlations and Kruskal-Wallis tests were used to determine
continuous and categorical parameters’ association with the
ModuleEigengeneValues file (noted in the script).

# References

### WGCNA

Weighted Gene Co-Expression Network Analysis (WGCNA) was first described
in the following publications:

\-Langfelder P, Horvath S (2008). “WGCNA: an R package for weighted
correlation network analysis.” BMC Bioinformatics, 559.
<https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559>.  
\-Langfelder P, Horvath S (2012). “Fast R Functions for Robust
Correlations and Hierarchical Clustering.” Journal of Statistical
Software, 46(11), 1-17. <https://www.jstatsoft.org/v46/i11/>.

### CarmonaMoraEtAl\_TranscriptomeDynamicsAfterStroke\_WGCNA.R

This script is based on and modified from Jeremy Miller’s “Meta-analyses
of data from two (or more) microarray data sets” tutorial. The tutorial
and related files can be found at:
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/>.
Miller’s tutorial is based on:

\-Miller JA, Horvath S, Geschwind DH. (2010) Divergence of human and
mouse brain transcriptome highlights Alzheimer disease pathways. Proc
Natl Acad Sci U S A. 2010 Jul 13;107(28):12698-703.

*Additionally, the following resources were utilized in creating this
script:*

Horvath, S. Weighted Network Analysis. Applications in Genomics and
Systems Biology. Book. 2011.
<https://link.springer.com/book/10.1007/978-1-4419-8819-5>

Steve Horvath’s Tutorial “Weighted Gene Co-Expression Network Analysis
(WGCNA) R Tutorial”
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ASPMgene/>,
which was based on:  
\-Horvath S, Zhang B, Carlson M, Lu KV, Zhu S, Felciano RM, Laurance MF,
Zhao W, Shu, Q, Lee Y, Scheck AC, Liau LM, Wu H, Geschwind DH, Febbo PG,
Kornblum HI, Cloughesy TF, Nelson SF, Mischel PS (2006) “Analysis of
Oncogenic Signaling Networks in Glioblastoma Identifies ASPM as a Novel
Molecular Target”, PNAS | November 14, 2006 | vol. 103 | no. 46 |
17402-17407

Langfelder, P. Signed vs. Unsigned Topological Overlap Matrix. Technical
Report. 2013.
<https://www.researchgate.net/file.PostFileLoader.html?id=57bdeaad40485404eb0753d4&assetKey=AS%3A398680193552384%401472064173254>

Langfelder and Horvath’s Tutorial “Network analysis of liver expression
data from female mice: finding modules related to body weight”  
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/>  
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf>  
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf>  
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf>  
based on:  
\-Ghazalpour A, Doss S, Zhang B, Wang S, Plaisier C, et al. (2006)
Integrating Genetic and Network Analysis to Characterize Genes Related
to Mouse Weight. PLOS Genetics 2(8): e130.
<https://doi.org/10.1371/journal.pgen.0020130>

WGCNA FAQ:
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html>

### TutorialFunctions.R

This script was taken from Jeremy Miller’s “Meta-analyses of data from
two (or more) microarray data sets” tutorial with a minor addition to
automatically output Hub Gene lists for each module (defined as the top
5% most interconnected genes in each module). The tutorial and related
files can be found at:
<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/>.
Miller’s tutorial is based on:

\-Miller JA, Horvath S, Geschwind DH. (2010) Divergence of human and
mouse brain transcriptome highlights Alzheimer disease pathways. Proc
Natl Acad Sci U S A. 2010 Jul 13;107(28):12698-703.
