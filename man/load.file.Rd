\name{load.file}
\alias{load.file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load sgRNA NGS Data especially for caRpools
}
\description{
This function is a parser of read.table to load sgRNA NGS data into a data.frame
}
\usage{
load.file(filename, header = TRUE, sep = "\t", comment.char="", type=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The filename of the NGS dataset file.
}
  \item{header}{
Specifies whether a header is present in the file or not.
}
  \item{sep}{
Specifies how data is separated column-wise. See ?read.table for further information.
}
\item{comment.char}{
comment.char see ?read.table
}
\item{type}{
Type of data being loaded. Bu default NULL, which loads tabular data.
Other values: xlsx for MIACCS file and fastalib to read the library reference fasta file
}
}
\details{
See ?read.table for further information.
}
\value{
load.file returns a data.frame.
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
#data.frame = load.file("sgRNA.txt", header= TRUE, sep="\t")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Loading Data }
\keyword{ ~read-Count }% __ONLY ONE__ keyword per line
