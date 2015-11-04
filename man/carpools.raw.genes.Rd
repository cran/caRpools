\name{carpools.raw.genes}
\alias{carpools.raw.genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting sgRNA phenotype effects of a given gene
}
\description{
CaRpools also allows you to visualize the phenotypic effects of sgRNA belonging to the same gene via `carpools.raw.genes`. This includes plotting of sgRNA foldchanges, z-score, z-ratios or read-counts. Moreover, `type="vioplot"` will present fold change data in comparison to the whole dataset and controls.
}
\usage{
carpools.raw.genes(untreated.list,treated.list, genes=NULL, namecolumn=1,
fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"),
do.plot=TRUE, log=FALSE, put.names=FALSE, type="foldchange", controls.target= NULL,
controls.nontarget=NULL, sort=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{untreated.list}{
A list of untreated sample data frames of read-count data as created by load.file().  
*Default* none  
*Values* A list of data frames of the untreated samples 
}
\item{treated.list}{
A list of treated sample data frames of read-count data as created by load.file().  
*Default* none  
*Values* A list of data frames of the treated samples 
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
\item{norm.function}{
The mathematical function to normalize data if `normalize=TRUE`. By default, the median is used.  
*Default* median  
*Values* Any mathematical function of R (function)  
}
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
}
\item{do.plot}{
Whether a plot is drawn or only tabular output is returned.  
*Default* TRUE  
*Values* TRUE, FALSE (boolean)  
}
\item{log}{
Plot in log-scale?  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
\item{put.names}{
Do you want the sgRNA identifiers to be plotted?  
*Default* FALSE  
*Values* TRUE, FALSE  
}
\item{type}{
Provides different types. "foldchange" for log2 foldchange, "readcount" for read-count, "z-score" for Z-scores, "z-ratio" for a Z-ratio or "vioplot" for a log2 FC of sgRNA effects.  
*Default* "foldchange"  
*Values* "foldchange", "readcount", "z-score", "z-ratio", "vioplot"  
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
\item{sort}{
This leads to output sorted by foldchange or z-ratio instead of names.  
*Default* TRUE  
*Values* TRUE, FALSE  
}
\item{genes}{
For which gene shall the sgRNA effect plots being generated?
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Return either generic plots or tables.
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

# Foldchange
p1 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="foldchange" )

# Z-Ratio
p2 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="z-ratio" )

# Read Count
p3 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="readcount" )

# Violine plot
p4 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="vioplot" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
\keyword{Read-count}