\name{carpools.hit.sgrna}
\alias{carpools.hit.sgrna}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting sgRNA effects for all candidate genes or single genes
}
\description{
Since there is more than just one single sgRNA targeting your gene of interest, you can user caRpools to plot different sgRNA phenotype effects, e.g. the fold change or z-ratio, as desribed before in `carpools.raw.genes`.  
A set of plots can be generated with `carpools.hit.sgrna`, which serves as a wrapper for `carpools.raw.genes`.  
By default, a foldchange plot as well as a violine plot are generated.
}
\usage{
carpools.hit.sgrna(wilcox=NULL, deseq=NULL, mageck=NULL, dataset=NULL,
dataset.names = NULL, namecolumn=1, fullmatchcolumn=2,
norm.function=median, extractpattern=expression("^(.+?)_.+"),
put.names=TRUE, type="enriched", labelgenes=NULL, cutoff.deseq = 0.05,
cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE,
plot.genes="overlapping", cutoff.hits=NULL,
plot.type=NULL, controls.target=NULL, controls.nontarget=NULL)
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
\item{labelgenes}{
For which gene shall the sgRNA effects being plotted? This expects a gene identifier or a vector of gene identifiers. If NULL, plots will be generated for all overlapping hit candidate genes.  
*Default* NULL  
*Values* A gene identifier or vector of gene identifiers (character)  
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
\item{put.names}{
Do you want the sgRNA identifiers to be plotted?  
*Default* FALSE  
*Values* TRUE, FALSE
}
\item{plot.type}{
WHich kind of plot is to be drawn? If NULL, foldchange and violine plots are generated.  
*Default* NULL  
*Values* NULL, "foldchange", "z-score", "z-ratio", "vioplot"    
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Return generic plots according to `type`.

By default, a foldchange plot as well as a violine plot are generated representing log2 fold changes of single sgRNAs.

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


sgrnas.en = carpools.hit.sgrna(wilcox=data.wilcox, deseq=data.deseq,
    mageck=data.mageck, dataset=list(CONTROL1, CONTROL2, TREAT1, TREAT2),
    dataset.names = c(d.CONTROL1, d.CONTROL2, d.TREAT1, d.TREAT2), namecolumn=1,
    fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)(_.+)"),
    put.names=TRUE, type="enriched", labelgenes="CASP8", plot.type=NULL, 
    cutoff.deseq = 0.001, cutoff.wilcox=0.05, cutoff.mageck = 0.05,
    cutoff.override=FALSE, cutoff.hits=NULL, controls.target="CASP8", controls.nontarget="random")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
