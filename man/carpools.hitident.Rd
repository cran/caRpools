\name{carpools.hitident}
\alias{carpools.hitident}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of hit analysis performed by Wilcox, DESeq2 and MAGeCK
}
\description{
The output from `stat.wilcox`, `stat.DEseq` and `stat.mageck` can be visualized with `carpools.hitident`.  
In this case, log2 fold changes are plotted against the gene names for all methods as well as the number of significant sgRNAs for data analyzed with DESeq2 or MAGeCK.  
}
\usage{
carpools.hitident(data, type="deseq2", title="DESeq2 plot", print.names=FALSE,
cutoff=c(0,0,0,0), inches=0.1, offsetplot=1.2, plot.p=0.01, sgRNA.top=1, separate=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
Output data from either `stat.wilcox`, `stat.DEseq` or `stat.mageck`.  
*Default* empty  
*Values* Output from either `stat.wilcox`, `stat.DEseq` or `stat.mageck`. 
}
\item{type}{
Which type of analysis method was used?  
*Default* deseq2  
*Values* "wilcox", "deseq2", "mageck"  
}
\item{title}{
Title of the plot.  
*Default* "DESeq2 plot"  
*Values* (character)  
}
\item{print.names}{
Shall the names of significant or top candidates being plotted?  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
\item{cutoff}{
A vector containing plotting cutoffs if `print.names=TRUE`.  c("top enriched", "top depleted", "most sgRNA enriched", "most sgRNA depleted").  
*Default* c(0,0,0,0)  
*Values* Vector of length 4 (numeric)  
}
\item{inches}{
see `?par`.  
*Default* 0.1  
*Values* (numeric)  
}
\item{offsetplot}{
Multiplication factor for stretching the plotting area to get a better plot experience.  
*Default* 1.2  
*Values* > 1 (numeric)  
}
\item{plot.p}{
Which p-value shall be plotted and used for visualization?  
*Default* 0.05  
*Values* (numeric)  
}
\item{sgRNA.top}{
For sgRNA plots, this indicates how many genes will be labeled (the top X genes).  
*Default* 1  
*Values* (numeric, integer)  
}
\item{separate}{
Gene that showed enrichment can be plotted separately from those that have shown a depletion for better overview, works only for wilcox.  
*Default* FALSE  
*Values* TRUE, FALSE  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
carpools.hitident returns a generic plot, which can be passed on to any device.
}

\author{
Jan Winter
}
\note{
see ?plot for detailed plotting information.
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
mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr",
filename = "TEST" , fdr.pval = 0.05)

mageck.result = carpools.hitident(data.mageck, type="mageck",
title="MAGeCK", inches=0.1, print.names=TRUE, plot.p=0.05, offsetplot=1.2, sgRNA.top=1)

wilcox.result = carpools.hitident(data.wilcox, type="wilcox",
title="Wilcox", inches=0.1, print.names=TRUE, plot.p=0.05, offsetplot=1.2, sgRNA.top=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Visualization }
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
