\name{stat.DESeq}
\alias{stat.DESeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analysis: DESeq2 Analysis of pooled CRISPR NGS data
}
\description{
For the DESeq2 analysis implementation, the read counts of all sgRNAs for a given gene are first summed up to increase the available read count.  
Then, DESeq2 analysis is perfomed, which includes the estimation of size-factors, the variance stabilization using a parametric fit and a Wald-Test for differnece in log2 fold changes between the untreated and treated data.  
More information about this can be found in _Love et al._  
[Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](http://www.ncbi.nlm.nih.gov/pubmed/25516281)  
_Genome Biology_ 2014  
}
\usage{
stat.DESeq(untreated.list,treated.list,namecolumn=1, fullmatchcolumn=2,
agg.function=sum, extractpattern=expression("^(.+?)_.+"), sorting=FALSE,
sgRNA.pval = 0.01, filename.deseq="data", fitType="parametric", p.adjust="holm")
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
  \item{agg.function}{
Function used to aggregate gene data from individual sgRNA data. By default, agg.function=mean, but it can be any other function e.g. sum or median.
}
  \item{extractpattern}{
%%     ~~Describe \code{extractpattern} here~~
Regular Expression, used to extract the gene name from the sgRNA name. Please make sure that the gene name extracted is accesible by putting its regular expression in brackets (). The default value expression("^(.+?)_.+") will look for the gene name (.+?) in front of the separator _ and any character afterwards .+ e.g. gene1_anything .

}
  \item{sorting}{
Defines whether the final output is sorted by the calculated p-value. By default, sorting=FALSE will return a table sorted by gene name.
}
\item{sgRNA.pval}{
p-value threshold to count significant sgRNAs for each gene.  
*Default* 0.001  
*Value* (numeric)  
}
\item{filename.deseq}{
Filename of raw DESeq2 data output.  
*Default* "data"  
*Values* (character)  
}
\item{fitType}{
See `?DESeq2`.  
*Default* "parametric"  
*Values* "parametric", "local" "mean"  
}
\item{p.adjust}{
Method to adjust p-value for multiple testing. See `?DEseq2`.  
*Default* "holm"
*Values* see `?DESeq2`  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
stat.DESeq returns a formal class that contains gene names including the calculated p-value. The returned class can be visualized using carpools.hitident (see ?carpools.hitident).
The output is formatted as follows:

log2 fold change (MAP): condition untreated vs treated\cr
Wald test p-value: condition untreated vs treated\cr
DataFrame with 813 rows and 6 columns

\tabular{lcccccc}{
     \tab  baseMean \tab log2FoldChange   \tab  lfcSE   \tab    stat  \tab   pvalue    \tab  padj \cr
AAK1 \tab  73.90565  \tab  -0.23319491 \tab 0.2927459 \tab-0.7965779 \tab 0.42569619\tab 0.7018234\cr
AATK \tab 159.43350  \tab  -0.11312924 \tab 0.2740927 \tab-0.4127408 \tab 0.67979655\tab 0.8514905\cr
ABI1 \tab 131.03013  \tab  -0.09915855 \tab 0.2693971 \tab-0.3680758 \tab 0.71281670\tab 0.8691949\cr
ABL1 \tab  77.51711  \tab   0.07837768 \tab 0.3155477  \tab 0.2483862 \tab 0.80383562\tab 0.9114121\cr
ABL2 \tab 119.22621  \tab  -0.49412039 \tab 0.2846396 \tab-1.7359507 \tab 0.08257254\tab 0.3128525\cr
...  \tab       ...    \tab       ...    \tab   ...  \tab      ... \tab       ...  \tab     ...\cr
}
}

\author{
Jan Winter, DESEq2 was developed by the Wolfgang Huber lab (EMBL, Heidelberg)
}
\note{
%%  ~~further notes~~
none
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

data(caRpools)
data.deseq = stat.DESeq(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1,TREAT2), namecolumn=1,
  fullmatchcolumn=2, extractpattern=expression("^(.+?)(_.+)"),
  sorting=FALSE, filename.deseq = "ANALYSIS-DESeq2-sgRNA.tab",
  fitType="parametric")
  
knitr::kable(data.deseq$genes[1:10,])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
