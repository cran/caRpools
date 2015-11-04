\name{carpools.read.count.vs}
\alias{carpools.read.count.vs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QC: Scatterplots of Read-Counts
}
\description{
CaRpools also allows you to compare the readcount for different samples using `carpools.read.count.vs`.  
By this, you can easily compare the screen and replicate performance as well as highlighting your non-targeting or positive controls. Moreover, you can highlight any gene as well.  
For details regarding all arguments and option see `?carpools.read.count.vs`.  
}
\usage{
carpools.read.count.vs(dataset, namecolumn=1, fullmatchcolumn=2, title="Read Count",
dataset.names = NULL, xlab="Readcount Dataset1", ylab="Readcount Dataset2", xlim=NULL,
ylim=NULL, pch=16, col = rgb(0, 0, 0, alpha = 0.65), labelgenes=NULL, labelcolor="red",
extractpattern=expression("^(.+?)_.+"), plotline=TRUE, normalize=TRUE,
norm.function=median, offsetplot=1.2, center=FALSE, aggregated=FALSE,
pairs=FALSE, type=NULL, plot.identify=FALSE, plot.log=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
\item{title}{
The title of the plot.  
*Default* "Read Count"  
*Values* "Any title" (character)  
}
\item{dataset.names}{
A list of names that must be according to the list of data sets given in *dataset*.  
*Default* NULL  
*Value* NULL or list of data names (list)  
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
\item{xlim}{
You can define the x-axis range being plotted, e.g. `c(0,1)`.  
*Default* empty  
*Values* empty or a vector with the lower and upper limit.
}
\item{ylim}{
You can define the y-axis range being plotted, e.g. `c(0,1)`.  
*Default* empty  
*Values* empty or a vector with the lower and upper limit.
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
\item{labelgenes}{
You can highlight certain genes within the plot. This expects a gene identifier or a fector of gene identifiers.  
*Default* NULL  
*Values* A gene identifier or vector of gene identifiers (character)  
}
\item{labelcolor}{
Color to highlight genes stated in `labelgenes`.  
*Default* "organge"  
*Values* Any R color or RGB color object.  
}
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
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
\item{norm.function}{
The mathematical function to normalize data if `normalize=TRUE`. By default, the median is used.  
*Default* median  
*Values* Any mathematical function of R (function)  
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
\item{pairs}{
In the case of plotting all four data sets at once, you can use a pairs plot for easier overview (see `?pairs()`).  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
\item{type}{
This indicates whether you would like to color all highlighted genes in either red ("enriched") or blue ("depleted") color according to the standrds in caRpools for plotting enriched or depleted genes after analysis.  
*Default* NULL  
*Values* NULL, "enriched", "depleted"  
}
\item{plot.identify}{
You can ask R to let you identify genes by clikcing on the dots in the graph. This only works if `pairs=FALSE`.  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
\item{plot.log}{
If all plots are created using log-transformed data.  
*Default* TRUE  
*Values* TRUE, FALSE (boolean)  
}
}
\details{
For generic plot arguments, see ?plot.
}
\value{
plot.read.count.vs returns a basic plot.
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

carpools.read.count.vs(dataset=list(TREAT1,CONTROL1),
dataset.names = c(d.TREAT1, d.CONTROL1),
  pairs=FALSE, namecolumn=1, fullmatchcolumn=2, title="", pch=16,
  normalize=TRUE, norm.function=median, labelgenes="random", labelcolor="blue",
  center=FALSE, aggregated=FALSE)
  
carpools.read.count.vs(dataset=list(TREAT1, TREAT2, CONTROL1, CONTROL2),
  dataset.names = c(d.TREAT1, d.TREAT2, d.CONTROL1, d.CONTROL2),
  pairs=TRUE, namecolumn=1, fullmatchcolumn=2, title="", pch=16,
  normalize=TRUE, norm.function=median,
  labelgenes="random", labelcolor="blue", center=FALSE, aggregated=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~compare }
\keyword{ ~readcount }% __ONLY ONE__ keyword per line
