\name{stat.wilcox}
\alias{stat.wilcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analysis: Analysis of pooled CRISPR screening data using a Wilcoxon Test
}
\description{
__Wilcox__

Within this approach, the read counts of all sgRNAs in one dataset are first normalized by the function set in the MIACCS file. By default, normalization is done by read count division with the dataset median.  
Then, the fold change of each population of sgRNAs for a gene is tested against the population of either the non-targeting controls or randomly picked sgRNAs, as defined by the random picks option within the MIACCS file, using a two-sided Mann-Whitney-U test. P-values are corrected for multiple testing using FDR.  

}
\usage{
stat.wilcox(untreated.list=list(NULL, NULL),treated.list=list(NULL, NULL),
namecolumn=1, fullmatchcolumn=2,normalize=TRUE,norm.fun=median,
extractpattern=expression("^(.+?)_.+"), controls=NULL, control.picks=300, sorting=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{untreated.list}{
A list of data.frames of untreated, control samples. e.g. list(df.control1, df.control2)
}
  \item{treated.list}{
A list of data.frames of treated samples. e.g. list(df.treated1, df.treated2)
}
  \item{namecolumn}{
In which the target names are located, e.g. namecolumn=1 for the first columns.
}
  \item{fullmatchcolumn}{
Column, in which readcounts are located, e.g. fullmatchcolumn=2 for the second column.
}
  \item{normalize}{
Datasets can be normalized by norm.fun if normalize=TRUE.
}
  \item{norm.fun}{
The function used to normalize the datasets if normalize=TRUE. By default, normalization is done using the dataset median, but any other function e.g. mean, can be used in principle.
}
 \item{extractpattern}{
%%     ~~Describe \code{extractpattern} here~~
Regular Expression, used to extract the gene name from the sgRNA name. Please make sure that the gene name extracted is accesible by putting its regular expression in brackets (). The default value expression("^(.+?)_.+") will look for the gene name (.+?) in front of the separator _ and any character afterwards .+ e.g. gene1_anything .
}
  \item{controls}{
DSS requires a set of non-targeting sgRNAs (negative controls) within the datasets. You can specify the arbitrary gene name for these controls using controls="arbitrary.gene.name.of.controls".
}
  \item{sorting}{
Analysis output is by default sorted by gene name (sorting=FALSE). If desired, the output table can be sorted according to the p-value of the genes (sorting=TRUE).
}
\item{control.picks}{
If no non-targeting controls are present or set, wilcox will pick a randum number of sgRNAs from the data set as the alternative population. This is only used if `controls=NULL`.  
*Default* 300  
*Values* numeric  
}
}

\value{
stat.wilcox return a data.frame, which can be visualized by plot.hitident.
The data.frame has the following format:

\tabular{lcccc}{

  \tab untreated\tab  treated\tab foldchange\tab   p.value \cr
AAK1 \tab 2.061346\tab 3.007924\tab   1.351672\tab 0.2966311\cr
AATK \tab 3.413357\tab 5.129985\tab   1.398695\tab 0.1146190\cr
ABI1 \tab 2.997385\tab 4.384881\tab   1.418959\tab 0.1437962\cr
ABL1 \tab 2.269906\tab 2.874087\tab   1.211499\tab 0.3681327\cr
ABL2 \tab 2.519391\tab 4.539583\tab   1.732575\tab 0.6335575\cr
}

For each gene, the foldchange as well as the p-value, derived by the Mann-Whitney U test against the non-targeting controls, are listed.
}

\author{
Jan Winter
}
\note{
%%  ~~further notes~~
none
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(caRpools)

data.wilcox = stat.wilcox(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2,
  normalize=TRUE, norm.fun=median, sorting=FALSE, controls="random",
  control.picks=NULL)
  
knitr::kable(data.wilcox[1:10,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
