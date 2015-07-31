\name{carpools.hit.scatter}
\alias{carpools.hit.scatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot: Plotting Scatters for hit candidate genes for all provided sampled
}
\description{
As described before, scatter plots can be generated for all datasets. `carpools.hit.scatter` serves as a wrapper for `carpools.read.count.vs` and allows faster plotting for individual candidate genes or all overlapping candidate genes. It generated a pairs plot with the representation of all provided samples and highlights the candidate gene. 
}
\usage{
carpools.hit.scatter(wilcox=NULL, deseq=NULL, mageck=NULL, dataset, dataset.names = NULL,
namecolumn=1, fullmatchcolumn=2, title="Read Count", xlab="Readcount Dataset1",
ylab="Readcount Dataset2", labelgenes=NULL, labelcolor="orange",
extractpattern=expression("^(.+?)_.+"),
plotline=TRUE, normalize=TRUE, norm.function=median, offsetplot=1.2,
center=FALSE, aggregated=FALSE, type="enriched",
cutoff.deseq = 0.001, cutoff.wilcox = 0.05,
cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL,
plot.genes="overlapping", pch=16, col = rgb(0, 0, 0, alpha = 0.65))
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
\item{dataset}{
A list of data frames of read-count data as created by load.file().  
*Default* none  
*Values* A list of data frames  
}
\item{namecolumn}{
In which column are the sgRNA identifiers?  
*Default* 1  
*Values* column number (numeric)  
}
\item{fullmatchcolumn}{
In which column are the read counts?  
*Default* 2  
*Values* column number (numeric)  
}
\item{dataset.names}{
A list of names that must be according to the list of data sets given in *dataset*.  
*Default* NULL  
*Value* NULL or list of data names (list)  
}
\item{norm.function}{
The mathematical function to normalize data. By default, the median is used.  
*Default* median  
*Values* Any mathematical function of R (function)  
}
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
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
\item{labelgenes}{
For which gene shall the sgRNA effects being plotted? This expects a gene identifier or a vector of gene identifiers.  
*Default* NULL  
*Values* A gene identifier or vector of gene identifiers (character)  
}
\item{xlab}{
Label of X-Axis, only if `pairs=FALSE`   
*Default* "X-Axis"  
*Values* "Label of X-Axis" (character)  
}
\item{ylab}{
Label of Y-Axism only if `pairs=FALSE`  
*Default* "Y-Axis"  
*Values* "Label of Y-Axis" (character)  
}
\item{pch}{
The type of point used in the plot. See `?par()`.  
*Default* 16  
*Values* Any number describing the point, e.g. 16 (numeric)  
}
\item{col}{
The color of the plotted data. Can be any R color or RGB object. See ?rgb() for further information.  
*Default* rgb(0, 0, 0, alpha = 0.65)  
*Values* Any R color name or RGB color object (character OR color object)  
}
\item{plotline}{
You can draw additional lines indicating a fold change of 0, 2, 4.  
*Default* TRUE  
*Values** TRUE, FALSE (boolean)  
}
\item{normalize}{
Whether you would like to normalize read-counts first. Recommended if not done already.  
*Default* TRUE  
*Values* TRUE, FALSE (boolean)  
}
\item{offsetplot}{
Offetplot is used to stretch the x- and y-axis for nicer graphs. This will extend plotting area by offsetplot.  
*Default* 1.2 (Plotting area is streched to 1.2 times)  
*Values* any number (numeric)  
}
\item{center}{
If you like you can center your data within the plot.  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
\item{aggregated}{
If you want to highlight genes, set this to true if you provide already aggregated gene read count instead of sgRNA read counts.  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
\item{labelcolor}{
Color to highlight genes stated in `labelgenes`.  
*Default* "organge"  
*Values* Any R color or RGB color object.  
}
\item{title}{
Title of the plot.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Return generic plots.
See ?plot and ?pairs.

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

#Single Gene
plothitsscatter.enriched = carpools.hit.scatter(wilcox=data.wilcox,
deseq=data.deseq, mageck=data.mageck, dataset=list(TREAT1, TREAT2, CONTROL1, CONTROL2),
dataset.names = c(d.TREAT1, d.TREAT2, d.CONTROL1, d.CONTROL2),
namecolumn=1, fullmatchcolumn=2, title="Title", labelgenes="CASP8",
labelcolor="orange", extractpattern=expression("^(.+?)(_.+)"),
normalize=TRUE, norm.function=median, offsetplot=1.2, center=FALSE,
aggregated=FALSE, type="enriched", cutoff.deseq = 0.001,
cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE,
cutoff.hits=NULL,  pch=16)

#Overlapping candidate genes

plothitsscatter.enriched = carpools.hit.scatter(wilcox=data.wilcox,
deseq=data.deseq, mageck=data.mageck, dataset=list(TREAT1, TREAT2, CONTROL1, CONTROL2),
dataset.names = c(d.TREAT1, d.TREAT2, d.CONTROL1, d.CONTROL2), namecolumn=1,
fullmatchcolumn=2, title="Title", labelgenes=NULL, labelcolor="orange",
extractpattern=expression("^(.+?)(_.+)"), normalize=TRUE, norm.function=median,
offsetplot=1.2, center=FALSE, aggregated=FALSE, type="enriched",
cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05,
cutoff.override=FALSE, cutoff.hits=NULL,  pch=16)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
\keyword{Visualization}