\name{carpools.sgrna.table}
\alias{carpools.sgrna.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Table Output of sgRNA effect and Target Sequence
}
\description{
Since there is more than just one single sgRNA targeting your gene of interest, you can user caRpools to plot different sgRNA phenotype effects, e.g. the fold change or z-ratio, as desribed before in `carpools.raw.genes`.  
In addition to that, caRpools also generated a tabular view which includes the log2 fold change as well as the target sequence, so that the user can directly pick the target sequence of the sgRNA he or shew wants.  

**This function, `carpools.sgrna.table` is best combined with `carpools.raw.genes` to give a fast overview of the sgRNA performance.**  
}
\usage{
carpools.sgrna.table (wilcox=NULL, deseq=NULL, mageck=NULL, dataset=NULL,
dataset.names = NULL, namecolumn=1, fullmatchcolumn=2, norm.function=median,
extractpattern=expression("^(.+?)_.+"), type="enriched", cutoff.deseq = 0.05,
cutoff.wilcox = 0.05, cutoff.mageck = 0.05,
cutoff.override=FALSE, plot.genes="overlapping", cutoff.hits=NULL, sgrna.file=NULL,
labelgenes=NULL, write=FALSE, datapath=getwd(), analysis.name="Screen")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{wilcox}{
Data output from `stat.wilcox`.  
*Default* NULL  
*Values* Data output from `stat.wilcox`.  
}
\item{deseq}{
Data output from `stat.deseq`.  
*Default* NULL  
*Values* Data output from `stat.deseq`.  
}
\item{mageck}{
Data output from `stat.mageck`.  
*Default* NULL  
*Values* Data output from `stat.mageck`.  
}
\item{cutoff.deseq}{
P-Value threshold used to determine significance.  
*Default* 0.001  
*Values* numeric  
}
\item{cutoff.wilcox}{
P-Value threshold used to determine significance.  
*Default* 0.001  
*Values* numeric  
}
\item{cutoff.mageck}{
P-Value threshold used to determine significance.  
*Default* 0.001  
*Values* numeric  
}
\item{dataset}{
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
\item{norm.function}{
The mathematical function to normalize data. By default, the median is used.  
*Default* median  
*Values* Any mathematical function of R (function)  
}
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
}
\item{cutoff.override}{
Shall the p-value threshold be ignored? If this is TRUE, the top percentage gene of `cutoff.hits` is used instead.  
*Default* FALSE  
*Values* TRUE, FALSE  
}
\item{cutoff.hits}{
The percentatge of top genes being used if `cutoff.override=TRUE`.  
*Default** NULL  
*Values* numeric  
}
\item{plot.genes}{
Defines what kind of data is used. By default, overlapping genes are highlighted in red color.  
*Default* "overlapping"  
*Values* "overlapping"  
}
\item{type}{
Defines whether all genes are plotted or only those being enriched or depleted.  
*Default* "all"  
*Values* "all", "enriched", "depleted"  
}
\item{sgrna.file}{
This is the library reference file loaded via `load.file` providing the sgRNA target sequence.  
*Default* NULL  
*Values* object from `load.file`  
}
\item{labelgenes}{
For which gene shal the sgRNA effects being generated? This expects a gene identifier or a fector of gene identifiers.  
*Default* NULL  
*Values* A gene identifier or vector of gene identifiers (character)  
}
\item{write}{
If you want to write directly to a file, this must be TRUE. Leave FALSE if you want the function to return a table.  
*Default* FALSE  
*Values* TRUE, FALSE  
}
\item{datapath}{
If `write=TRUE`, this is the directory the file is written.  
*Default* getwd()  
*Values* absolute path  
}
\item{analysis.name}{
The name of the file if `write=TRUE`  
*Default* "Screen"  
*Values* any file name (character)  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Output is a table or file (if write=TRUE).

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

data.wilcox = stat.wilcox(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2,
  normalize=TRUE, norm.fun=median, sorting=FALSE, controls="random",
  control.picks=NULL)
  
data.deseq = stat.DESeq(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1,TREAT2), namecolumn=1,
  fullmatchcolumn=2, extractpattern=expression("^(.+?)(_.+)"),
  sorting=FALSE, filename.deseq = "ANALYSIS-DESeq2-sgRNA.tab",
  fitType="parametric")
  
data.mageck = stat.mageck(untreated.list = list(CONTROL1, CONTROL2),
treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2,
norm.fun="median", extractpattern=expression("^(.+?)(_.+)"),
mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr",
filename = "TEST" , fdr.pval = 0.05)

sgrnas.en.table = carpools.sgrna.table(wilcox=data.wilcox, deseq=data.deseq,
    mageck=data.mageck, dataset=list(CONTROL1, CONTROL2, TREAT1, TREAT2),
    dataset.names = c(d.CONTROL1, d.CONTROL2, d.TREAT1, d.TREAT2), namecolumn=1,
    fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)(_.+)"),
    type="enriched", labelgenes="CASP8", cutoff.deseq = 0.001, cutoff.wilcox=0.05,
    cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL, sgrna.file = libFILE,
    write=FALSE)
    
knitr::kable(sgrnas.en.table)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~sgRNA }% __ONLY ONE__ keyword per line
