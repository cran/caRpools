\name{compare.analysis}
\alias{compare.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exporting Hit Candidate Gene Information
}
\description{
Although the candidate lists of each analysis method can be saved separately, caRpools offer a comparative approach, which creates tables that include the information from all analysis methods at once for a faster overview.  

This is done using the function `compare.analysis`, which offers not only ouptu for Venn Diagrams, but also for tables.  
}
\usage{
compare.analysis(wilcox=NULL, deseq=NULL, mageck=NULL, type="enriched",
cutoff.deseq = NULL, cutoff.wilcox = NULL, cutoff.mageck = NULL,
cutoff.override=TRUE, cutoff.hits=5, output="list",
sort.by=c("mageck","fdr","fdr"), plot.method=c("wilcox","mageck", "deseq"),
plot.feature=c("pval","fdr","pval"), pch=16)
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
\item{type}{
Either enriched or depleted.
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
\item{output}{
Three different types of output can be generated: A list with all genes including the information from `stat.wilcox`, `stat.DEseq` and `stat.mageck`, a sorted ranked output, a venn diagram compatible output and 3D scatterplot.  
*Default* "list"  
*Values* "list", "rank", "venn", "3dplot"  
}
\item{sort.by}{
This indicates the sorting for `type="rank" and type="list"` and is a vector. By default, data is sorted by the FDR of MAGeCK. needs to be a vector.  
*Default* c("mageck","fdr","fdr")  
*Values* c("mageck","fdr","fdr"), c("mageck","fdr","rank"), c("mageck","fdr","rank"), c("wilcox","pval","pval"), c("wilcox","pval","genes"), c("deseq","pval","pval"), c("deseq","pval","genes")  
}
\item{plot.method}{
Used only if `type="3dplot"`. This indicates what is plotted at the x, y and z-axis and thus needs to be a vector of length 3.  
*Default* c("wilcox","mageck", "deseq"), plots wilcox on X-axis, mageck on y-axis and deseq on z-axis  
*Values* c("wilcox","mageck", "deseq") or any other combination  
}
\item{plot.feature}{
If `type="3dplot"`, this indicates the type of data plotted on each axis of the 3d plot. This can only be set according to the features available of the method used to be plotted as indicated in `plot.method`.  
*Default* c("pval","fdr", "pval")  which uses the p-value of wilcox, the fdr or MAGeCK and p-value of DESeq2.  
*Values* c("pval","fdr", "pval"), or ANY combination according to `plot.method`  
}
\item{pch}{
The type of point used in the plot. See `?par()`.  
*Default* 16  
*Values* Any number describing the point, e.g. 16 (numeric)  
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Returns a table with information.
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
mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr",
filename = "TEST" , fdr.pval = 0.05)

# Perform the comparison
data.analysis.enriched = compare.analysis(wilcox=data.wilcox,
    deseq=data.deseq, mageck=data.mageck, type="enriched",
    cutoff.override = FALSE, cutoff.hits=NULL, output="list",
    sort.by=c("mageck","fdr","rank"))
## Write to a file
xlsx::write.xlsx(data.analysis.enriched,
    file="COMPARE-HITS.xls",
    sheetName="Enriched")
# Print to console
knitr::kable(data.analysis.enriched[1:10,c(2:7)])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
