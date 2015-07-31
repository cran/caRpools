\name{carpools.hit.overview}
\alias{carpools.hit.overview}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analysis: Analysis of pooled CRISPR screening data using a Wilcoxon Test
}
\description{
Candidate genes from all methods can be plotted in an overview to identify overlapping signficant candidate genes using `carpools.hit.overview`.  
}
\usage{
carpools.hit.overview(wilcox=NULL, deseq=NULL, mageck=NULL, cutoff.deseq = 0.001,
cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL,
plot.genes="overlapping", type="all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{wilcox}{
Data output from `stat.wilcox`.  
*Default* NULL  
*Values* Data output from `stat.wilcox`.  
}
\item{deseq}{
Data output from `stat.deseq`.  
*Default* NULL  
*Values* Data output from `stat.deseq`.  
}
\item{mageck}{
Data output from `stat.mageck`.  
*Default* NULL  
*Values* Data output from `stat.mageck`.  
}
\item{cutoff.deseq}{
P-Value threshold used to determine significance.  
*Default* 0.001  
*Values* numeric  
}
\item{cutoff.wilcox}{
P-Value threshold used to determine significance.  
*Default* 0.001  
*Values* numeric  
}
\item{cutoff.mageck}{
P-Value threshold used to determine significance.  
*Default* 0.001  
*Values* numeric  
}
\item{cutoff.override}{
Shall the p-value threshold be ignored? If this is TRUE, the top percentage gene of `cutoff.hits` is used instead.  
*Default* FALSE  
*Values* TRUE, FALSE  
}
\item{cutoff.hits}{
The percentatge of top genes being used if `cutoff.override=TRUE`.  
*Default** NULL  
*Values* numeric  
}
\item{plot.genes}{
Defines what kind of data is used. By default, overlapping genes are highlighted in red color.  
*Default* "overlapping"  
*Values* "overlapping"  
}
\item{type}{
Defines whether all genes are plotted or only those being enriched or depleted.  
*Default* "all"  
*Values* "all", "enriched", "depleted"  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Returns a generic plot.
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


data.deseq = stat.DESeq(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1,TREAT2), namecolumn=1,
  fullmatchcolumn=2, extractpattern=expression("^(.+?)(_.+)"),
  sorting=FALSE, filename.deseq = "ANALYSIS-DESeq2-sgRNA.tab",
  fitType="parametric")
  
data.mageck = stat.mageck(untreated.list = list(CONTROL1, CONTROL2),
treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2,
norm.fun="median", extractpattern=expression("^(.+?)(_.+)"),
mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr", filename = "TEST" , fdr.pval = 0.05)

carpools.hit.overview(wilcox=data.wilcox, deseq=data.deseq, mageck=data.mageck,
    cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05,
    cutoff.override=FALSE, cutoff.hits=NULL, plot.genes="overlapping", type="enriched")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
