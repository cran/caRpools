\name{stats.data}
\alias{stats.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating data set statistics
}
\description{
General statistics for a given dataset can be obtained by `stats.data`.
}
\usage{
stats.data(dataset, namecolumn = 1, fullmatchcolumn = 2,
extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = NA,
controls.target = NULL, controls.nontarget = "random", type="stats")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dataset}{
Data frame of read-count object.
*Default* none
*Values* data frame as created by `load.file()`
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
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)
}
\item{readcount.unmapped.total}{
Number of raw NGS reads, only used if `type="mapping`.
*Default* NA
*Values* Number of raw reads (integer)
}
\item{controls.target}{
If `type="controls"`, this is the gene identifier of the positive control.
*Default* NULL
*Value* Gene Identifier (character)
}
\item{controls.nontarget}{
If `type="controls"`, this is the gene identifier of the non-targeting control.
*Default* "random"
*Value* Gene Identifier (character)
}
\item{type}{
Which type os statistic will be generated.
*Default* "stats" 
*Values* "stats" will generate short statistics like median and mean for the data set,
"mapping" will generate an overview of how many reads are present,
"datatset" is used to generate in-depth statistics for each gene of a dataset,
"controls" is used for in-depth statistics of the controls.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Returns different tabular outputs.
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
U1.stats = stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
                      extractpattern=expression("^(.+?)_.+"), type="stats")

knitr::kable(stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
  extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = 1786217, type="mapping"))
  
knitr::kable(stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
  extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = 1786217,
  type="stats"))
  
knitr::kable(stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
  extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = 1786217,
  type="dataset")[1:10,1:5])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
