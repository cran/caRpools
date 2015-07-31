\name{get.gene.info}
\alias{get.gene.info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieving Gene Annotation and Gene Identifier Conversion from BiomaRt
}
\description{
It is also possible to either enrich the screening dataset file with additional information provided by the biomaRt interface.  
For example, gene identifiers can be changed from EnsemblIDs to official gene symbols are Gene Ontology terms can be added to the dataset.  
This can be done using `get.gene.info`, which serves as a wrapper for the **biomaRt** package with its load of options and possibilities (more information see `?biomaRt`).  

You can convert any gene identifier which is included in your sgRNA identifer to e.g. EnsemblID or HGNC Gene Symbol using caRpools.  
**Please note that Internet Access is required for biomaRt.**  
For further information about biomaRt conversion, please see the [biomaRt Manual](www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf).  
}
\usage{
get.gene.info(data, namecolumn=1, extractpattern=expression("^(.+?)(_.+)"),
database="ensembl", dataset="hsapiens_gene_ensembl", filters="ensembl_gene_id",
attributes = c("hgnc_symbol"), return.val = "dataset", controls=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
Data frame that contains read-count data.
*Default* none
*Values* data.frame containing read-count data (data.frame)  
}
\item{namecolumn}{
In which column are the sgRNA identifiers?  
*Default* 1  
*Values* column number (numeric)  
}
\item{extractpattern}{
PERL regular expression that is used to retrieve the gene identifier from the overall sgRNA identifier.  
e.g. in **AAK1_107_0** it will extract **AAK1**, since this is the gene identifier beloning to this sgRNA identifier. **Please see: Read-Count Data Files**  
*Default* expression("^(.+?)(_.+)"), will work for most available libraries.  
*Values* PERL regular expression with parenthesis indicating the gene identifier (expression)  
}
\item{database}{
BiomaRt database to be used. See `?listMarts()` or biomaRt documentation.  
*Default* "ensembl", is using the ensembl database  
*Values* Any biomaRt database (character)  
}
\item{dataset}{
The biomaRt dataset to be used. For *homo sapiens*, *hsapiens_gene_ensembl* is recommended. See `?listDatasets` or biomaRt documentation.  
*Default* "hsapiens_gene_ensembl"  
*Values* Any biomaRt dataset (character)  
}
\item{filters}{
The input filter information to retrieve biomaRt annotation, usually is the type of gene identifier used in the read-count files, e.g. "ensemble_gene_id". see `?listFilters`  
*Default* "ensembl_gene_id"  
*Values* Any biomaRt filter (character)  
}
\item{attributes}{
The output attribute to retrieve from biomaRt, usually the annotations that need to be fetched, e.g. "hgnc_symbol". see `?listAttributes`  
*Default* "hgnc_symbol"  
*Values* Any biomaRt attribute (character)  
}
\item{return.val}{
The type of object that is returned. For whole dataset, e.g. conversion of gene identifiers, use "dataset".  
*Default* "dataset"  
*Values* "dataset" (will give back the same data frame, but with exchanged gene identifiers), "info" (will return a data frame with all attributes fetched for genes, is used to annotate gene with additional information)  
}
\item{controls}{
Is set to TRUE if `data` is not a data frame, but a vector.  
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Return either a data.frame with converted gene identifier or a data frame with annotations.
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
#CONTROL1.replaced = get.gene.info(CONTROL1, namecolumn=1,
#extractpattern=expression("^(.+?)(_.+)"), database="ensembl",
#dataset="hsapiens_gene_ensembl", filters="hgnc_symbol",
#attributes = c("ensembl_gene_id"), return.val = "dataset")

#knitr::kable(CONTROL1.replaced[1:10,])

#CONTROL1.replaced.info = get.gene.info(CONTROL1, namecolumn=1,
#extractpattern=expression("^(.+?)(_.+)"), database="ensembl",
#dataset="hsapiens_gene_ensembl", filters="hgnc_symbol",
#attributes = c("ensembl_gene_id","description"), return.val = "info")

#knitr::kable(CONTROL1.replaced.info[1:10,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
