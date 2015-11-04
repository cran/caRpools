\name{carpools.read.depth}
\alias{carpools.read.depth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QC: Plot Sequencing Read Depth
}
\description{
You can also visualize the read depth of genes per sgRNA in order to check for sufficient sequencing depth using `carpools.read.depth`.  
For further details see `?carpools.read.depth`.  
You can either plot single dat samples or all four data samples at once.  
}
\usage{
carpools.read.depth(datasets, namecolumn=1, fullmatchcolumn=2, dataset.names=NULL,
extractpattern=expression("^(.+?)_.+"), col=rgb(0, 0, 0, alpha = 0.65), xlab="Genes",
ylab="Read Count per sgRNA", statistics=TRUE, labelgenes = NULL,
controls.target = controls.target,
controls.nontarget=controls.nontarget, labelcolor="orange", waterfall=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{datasets}{
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
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
}
\item{col}{
The color of the plotted data. Can be any R color or RGB object. See ?rgb() for further information.  
*Default* rgb(0, 0, 0, alpha = 0.65)  
*Values* Any R color name or RGB color object (character OR color object)  
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
\item{controls.target}{
Highlights the positive control in red color.  
*Default* NULL  
*Value* Gene Identifier (character)  
}
\item{controls.nontarget}{
Highlights the non-targeting control in blue color.  
*Default* "random"  
*Value* Gene Identifier (character)  
}
\item{waterfall}{
You can either plot the read depth sorted by gene identifier (FALSE, default) or according to the read depth.  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  s
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
notes
}
\value{
plot.read.depth returns a generic plot.
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

carpools.read.depth(datasets = list(CONTROL1), namecolumn=1 ,fullmatchcolumn=2,
  dataset.names=list(d.CONTROL1), extractpattern=expression("^(.+?)_.+"),
  xlab="Genes", ylab="Read Count per sgRNA",statistics=TRUE, labelgenes = NULL,
  controls.target = "CASP8", controls.nontarget="random", waterfall=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~readdepth }
\keyword{ ~qc }% __ONLY ONE__ keyword per line
