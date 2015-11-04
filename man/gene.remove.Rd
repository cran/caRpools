\name{gene.remove}
\alias{gene.remove}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove gene information from sgRNA data.frame
}
\description{
This function is used to remove genes/gene information from a data.frame containing pooled CRISPR screen data.
It is meant to exclude genes from the analysis and removes all entries belonging to a gene from the sgRNA data.frame.
}
\usage{
gene.remove(data, namecolumn = 1, toremove = NULL,
extractpattern = expression("^(.+?)_.+") )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data.frame with sgRNA readcounts. Must have one column with sgRNA names and one column with readcounts. Please note that the data must be formatted in a way, that gene names are included within the sgRNA name and can be extracted using the extractpattern expression.
e.g. GENE_sgRNA1 -> GENE as gene name, _ as the separator and sgRNA1 as the sgRNA identifier.
}
  \item{namecolumn}{
%%     ~~Describe \code{namecolumn} here~~
integer, indicates in which column the names are stored
}

  \item{toremove}{
Vector of gene names that will be removed from sgRNA dataset. The gene name must be included in the sgRNA names in order to be extracted using the pattern defined in extractpattern. e.g. c=("gene1","gene2")
}
  \item{extractpattern}{
Regular Expression, used to extract the gene name from the sgRNA name. Please make sure that the gene name extracted is accesible by putting its regular expression in brackets (). The default value expression("^(.+?)_.+") will look for the gene name (.+?) in front of the separator _ and any character afterwards .+ e.g. gene1_anything .
}
}
\details{
In a table with

\tabular{rr}{
    DesignID \tab fullmatch   \cr
    AAK1_104_0 \tab 0         \cr
       AAK1_105_0 \tab  197        \cr
       AAK1_106_0 \tab  271         \cr
      AAK1_107_0  \tab   1         \cr
       AAK1_108_0 \tab    0         \cr
  }
  
calling gene.remove(data.frame, toremove="AAK1", extractpattern = expression("^(.+?)_.+")) will remove all entries shown above, since AAK1 is the gene name, separated by an undescore _ from the sgRNA identifier.

}
\value{
gene.remove returns a data.frame that has the same column dimensions as the input data.frame, however all rows in which toremove=gene is present, are deleted.
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
gene.remove(CONTROL1, toremove="AAK1", extractpattern = expression("^(.+?)_.+"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Gene }% __ONLY ONE__ keyword per line