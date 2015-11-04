\name{carpools.waterfall.pval}
\alias{carpools.waterfall.pval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of p-value distribution
}
\description{
Each of the analysis methods returns an adjusted p-value (corrected for multiple testing) as well as a fold change (Wilcox, DESeq2) or gene rank (MAGeCK).  
Therefore the -log10 p-value can be plotted against the gene names with `carpools.waterfall.pval`: 
}
\usage{
carpools.waterfall.pval (type=NULL,dataset=NULL, pval=0.05, mageck.type="pos", log=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{type}{
This indicates which kind of analysis method was used for p-value calculation.  
*Default* NULL  
*Values* "mageck", "deseq2", "wilcox"  
}
\item{dataset}{
Result from either `stat.wilcox`, `stat.DEseq` or `stat.mageck`.  
*Default* NULL  
*Values* Result from either `stat.wilcox`, `stat.DEseq` or `stat.mageck`  
}
\item{pval}{
The significance value set for the analysis which is to be plotted.  
*Default* 0.05  
*Values* numeric  
}
\item{mageck.type}{
Only for plotting p-value calculate by MAGeCK. Indicates whether enriched ("pos") or depleted ("neg") genes are used.  
*Default* "pos"  
*Values* "pos", "neg"  
}
\item{log}{
-log10 of the p-valiues is plotted if set to TRUE.  
*Default* TRUE  
*Values* TRRUE, FALSE (boolean)  
}
}
\value{

Return a generic plot.
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

data.mageck = stat.mageck(untreated.list = list(CONTROL1, CONTROL2),
treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2,
norm.fun="median", extractpattern=expression("^(.+?)(_.+)"),
mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr",
filename = "TEST" , fdr.pval = 0.05)

carpools.waterfall.pval(type="mageck", dataset=data.mageck, pval=0.05, log=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
