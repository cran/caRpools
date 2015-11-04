\name{stat.mageck}
\alias{stat.mageck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analysis: Analysis of pooled CRISPR screening data using a MAGeCK
}
\description{
CaRpools also uses MAGeCK to look for enriched or depleted genes within your screening data. Please note that MAGeCK needs to be installed correctly, this can be tested by `check.caRpools`.  

Within this approach, the read counts of all sgRNAs in one dataset are first normalized by the function set in the MIACCS file. By default, normalization is done by read count division with the dataset median.  
Then, the fold change of each population of sgRNAs for a gene is tested against the population of either the non-targeting controls or randomly picked sgRNAs, as defined by the random picks option within the MIACCS file, using a two-sided Mann-Whitney-U test. P-values are corrected for multiple testing using FDR.  
}
\usage{
stat.mageck(untreated.list, treated.list, namecolumn=1, fullmatchcolumn=2,
norm.fun=median, extractpattern=expression("^(.+?)_.+"), mageckfolder=NULL,
sort.criteria="neg", adjust.method="fdr", filename=NULL, fdr.pval=0.05)
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
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
}
\item{sort.criteria}{
MAGeCK argument *--sort-criteria*  
*Default* "neg"  
*Values* see MAGeCK documentation  
}
\item{mageckfolder}{
Folder for MAGeCK raw data output (internally used).  
*Default* NULL  
*Value* (character) 
}
\item{filename}{
Filename of raw MAGeCK data output.  
*Default* "data"  
*Values* (character)    
}
\item{adjust.method}{
Method to adjust p-value for multiple testing. See MAGeCK documentation.  
*Default* "fdr"
*Values* see MAGeCK documentation  
}
\item{fdr.pval}{
FDR used for correction.
*Default* 0.05  
*Values* (numeric)  
}
\item{norm.fun}{
The mathematical function to normalize data. By default, the median is used.  
*Default* median  
*Values* Any mathematical function of R (function)  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
stat.mageck retrieves a list of two data frames.
One with gene information, the other with sgRNA information.
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
data.mageck = stat.mageck(untreated.list = list(CONTROL1, CONTROL2),
treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2,
norm.fun="median", extractpattern=expression("^(.+?)(_.+)"),
mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr",
filename = "TEST" , fdr.pval = 0.05)

knitr::kable(data.mageck$genes[1:10,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
