\name{check.caRpools}
\alias{check.caRpools}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test caRpools installation and dependent software
}
\description{
You can verify that the MIACCS.xls file as well as the used template file and all necessary scripts are found by calling `check.caRpools()`.  
CaRpools also uses MAGeCK to look for enriched or depleted genes within your screening data. Please note that MAGeCK needs to be installed correctly, this can be tested by `check.caRpools`.  
}
\usage{
check.caRpools(packages=TRUE, files=TRUE, mageck=TRUE, bowtie2=TRUE,
pandoc=TRUE, skip.updates=TRUE, template=NULL, scripts=TRUE, miaccs="MIACCS.xls")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{packages}{
if TRUE, packages will be checked using load.packages()
}
\item{files}{
If TRUE, MIACCS as well as data and scripts folder will be checked in addition to CRISPR-mapping.pl and CRISPR-extract.pl.
}
\item{mageck}{
If TRUE, mageck installation is checked.
}
\item{bowtie2}{
if TRUE, bowtie2 installation is checked.
}
\item{pandoc}{
if TRUE, pandoc installation is checked.
}
\item{skip.updates}{
if TRUE, updates are skipped during package check.
}
\item{template}{
Rmd template file name to use.
}
\item{scripts}{
if TRUE, checks for perl scripts CRISPR-mapping and CRISPR-extract.pl.
}
\item{miaccs}{
Filename of MIACCS file. Will be checked for proper loading.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
This function does not return any value.
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
#check.caRpools()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
