\name{aggregatetogenes}
\alias{aggregatetogenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregates pooled CRISPR screen sgRNA data to gene data
}
\description{
Aggregate all sgRNA data from pooled CRISPR screens to their corresponding gene level.
}
\usage{
aggregatetogenes(data.frame, namecolumn = 1, countcolumn = 2,
agg.function = sum, extractpattern = expression("^(.+?)_.+"), type="aggregate")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.frame}{
%%     ~~Describe \code{data.frame} here~~
data.frame with sgRNA readcounts. Must have one column with sgRNA names and one column with readcounts. Please note that the data must be formatted in a way, that gene names are included within the sgRNA name and can be extracted using the extractpattern expression.
e.g. GENE_sgRNA1 -> GENE as gene name, _ as the separator and sgRNA1 as the sgRNA identifier.
}
  \item{namecolumn}{
%%     ~~Describe \code{namecolumn} here~~
integer, indicates in which column the names are stored
}
  \item{countcolumn}{
%%     ~~Describe \code{countcolumn} here~~
integer, indicates in which column the readcount are stored
}
  \item{agg.function}{
%%     ~~Describe \code{agg.function} here~~
expression, the function to be used for aggregating data. Since for sgRNAs, aggregating data to the corresponding gene, sum will be right function in this case. Other possibilities include any other mathematical function R is capable of, e.g. median, mean.
}
  \item{extractpattern}{
%%     ~~Describe \code{extractpattern} here~~
Regular Expression, used to extract the gene name from the sgRNA name. Please make sure that the gene name extracted is accesible by putting its regular expression in brackets (). The default value expression("^(.+?)_.+") will look for the gene name (.+?) in front of the separator _ and any character afterwards .+ e.g. gene1_anything .
}
\item{type}{
CaRpools can either aggregate the data frame (`type = "annotate"`) or annotate the gene identifiers only as an additional column (`type = "annotate"`).  
*Default* "aggregate"  
*Values* "aggregate", "annotate" 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
aggregatetogenes can be used after load.file() to create quality control plots for aggregated gene data instead of single sgRNA data.

Before:

\tabular{rr}{
    DesignID \tab fullmatch   \cr
    AAK1_104_0 \tab 0         \cr
       AAK1_105_0 \tab  197        \cr
       AAK1_106_0 \tab  271         \cr
      AAK1_107_0  \tab   1         \cr
       AAK1_108_0 \tab    0         \cr
  }


Afterwards:

\tabular{rr}{
DesignID \tab fullmatch \cr
  AAK1   \tab    880        \cr
 AATK    \tab  2105         \cr
 ABI1    \tab  1610         \cr

}

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data.frame is returned with namecolumn (which no includes only gene names) and all readcount information aggregated by the agg.function.
}

\author{
%%  ~~who you are~~
Jan Winter
}
\note{
none
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(caRpools)

CONTROL1.g=aggregatetogenes(data.frame = CONTROL1, agg.function=sum,
                            extractpattern = expression("^(.+?)(_.+)"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CRISPR }

