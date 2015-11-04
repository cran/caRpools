\name{carpools.read.distribution}
\alias{carpools.read.distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QC: Plot Readcount Distribution
}
\description{
A distribution for NGS data readcount can be created by `carpools.read.distribution` to visualize how the data set is distributed.  
This allows to check for data skewness and to estimate the overall assay quality.  
For further details see `?carpools.read.distribution`.  
}
\usage{
carpools.read.distribution(dataset,namecolumn=1, fullmatchcolumn=2, breaks="",
title="Title", xlab="X-Axis", ylab="Y-Axis",statistics=TRUE,
col=rgb(0, 0, 0, alpha = 0.65), extractpattern=expression("^(.+?)_.+"),
plotgene=NULL, type="distribution", logscale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dataset}{
Data frame of read-count data as created by load.file().  
*Default* none  
*Values* A data frame  
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
\item{breaks}{
Histogramm breaks see `?hist`. By default, will be calculated according to the dataset length.  
*Default* NULL
*Values* (numeric)  
}
\item{title}{
Main title of plot  
*Default* "Title"  
*Values* "The title you want" (character)  
}
\item{xlab}{
Label of X-Axis  
*Default* "X-Axis"  
*Values* "Label of X-Axis" (character)  
}
\item{ylab}{
Label of Y-Axis  
*Default* "Y-Axis"  
*Values* "Label of Y-Axis" (character)  
}
\item{statistics}{
Whether basic stattistics will be shown in the plot.  
*Default* TRUE  
*Values* TRUE, FALSE (boolean)  
}
\item{col}{
The color of the plotted data. Can be any R color or RGB object. See ?rgb() for further information.  
*Default* rgb(0, 0, 0, alpha = 0.65)  
*Values* Any R color name or RGB color object (character OR color object)  
}
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
}
\item{plotgene}{
You can only plot the read count distribution of sgRNAs belonging to a certain gene, which is given to the function via plotgene.  
*Default* NULL  
*Value* NULL or gene identifier (character)  
}
\item{type}{
You can plot either the read count distribution either as a normal histogram, or a box-and-whisker plot.  
*Default* "distribution"  
*Values* "distribution" to plot a histogram, or "whisker" to plot a whisker plot (character)  
}
\item{logscale}{
Indicates whether the read-count is plotted in a logarithmic scale.  
*Default* TRUE  
*Values* TRUE, FALSE (boolean)  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
plot.read.distribution return a generic plot, that can be passed on to any device.
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

carpools.read.distribution(CONTROL1, fullmatchcolumn=2,breaks=200,
  title=d.CONTROL1, xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE) 
  
carpools.read.distribution(CONTROL1, fullmatchcolumn=2,breaks=200,
  title=d.CONTROL1, xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE,
  type="whisker") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Distribution }
\keyword{ ~Reads }% __ONLY ONE__ keyword per line
