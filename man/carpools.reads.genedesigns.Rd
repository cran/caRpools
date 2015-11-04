\name{carpools.reads.genedesigns}
\alias{carpools.reads.genedesigns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QC: Plot representation of sgRNAs per gene
}
\description{
Since in most cases several sgRNAs are used to target a gene, the information how many sgRNAs are present in the data for each gene is of interest to make sure the number of sgRNAs present is still sufficient. Typically, only few sgRNAs should get "lost" during the screening procedure, so that the full sgRNA coverage is maintained throughout the assay. The only exception would be drop-out screens with a stringent setup.  
The representation of sgRNAs per gene can be plotted using `carpools.reads.genedesigns`.  
For further details see `?carpools.reads.genedesigns`.  
}
\usage{
carpools.reads.genedesigns(dataset, namecolumn=1, fullmatchcolumn=2, title="Read Count",
xlab="Percentage of sgRNAs present", ylab="Number of Genes", agg.function=sum,
extractpattern=expression("^(.+?)_.+"), col = rgb(0, 0, 0, alpha = 0.65))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dataset}{
A data frame of read-count data as created by load.file().  
*Default* none  
*Values* Adata frame  
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
\item{agg.function}{
The function to aggregate sgRNA read-count.
*Default* sum  
*Values* any mathematical function (function)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
carpools.reads.genedesigns returns a generic plot.
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

control1.readspergene = carpools.reads.genedesigns(CONTROL1, namecolumn=1, fullmatchcolumn=2,
title=paste("sgRNA Represenation:", d.CONTROL1, sep=" "),
xlab="Percentage of sgRNAs present", ylab="# of Genes")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~coverage }
\keyword{ ~sgRNA }% __ONLY ONE__ keyword per line
