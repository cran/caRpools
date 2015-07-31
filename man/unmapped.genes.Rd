\name{unmapped.genes}
\alias{unmapped.genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sgRNAs without reads
}
\description{
CaRpools also provides you with the number of missing sgRNA, that means sgRNAs without a single read during NGS.
If you want to know WHICH sgRNAs dropped out for a given gene, please consider using `genes` as an optional argument with the gene identifier of interest.
}
\usage{
unmapped.genes(data, namecolumn=1, fullmatchcolumn=2,
genes=NULL, extractpattern=expression("^(.+?)_.+"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
A data.frame as created by `load.file`.  
*Default* empty  
*Values* read-count data.frame  
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
\item{genes}{
If you want to know how many sgRNAs are not present for a single gene, set `genes` to your gene identifier of interest.  
*Default* NULL  
*Values* gene identifier (character)  
}
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
}
}

\value{
Tabular output with number of missing sgRNAs for each gene or the name of the missing sgRNA if genes!=NULL.
}

\author{
Jan Winter
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(caRpools)
U1.unmapped = unmapped.genes(data=CONTROL1, namecolumn=1,
fullmatchcolumn=2, genes=NULL, extractpattern=expression("^(.+?)_.+"))

knitr::kable(U1.unmapped)

U1.unmapped = unmapped.genes(data=CONTROL1, namecolumn=1,
fullmatchcolumn=2, genes="random", extractpattern=expression("^(.+?)_.+"))

knitr::kable(U1.unmapped)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
