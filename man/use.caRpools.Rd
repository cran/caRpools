\name{use.caRpools}
\alias{use.caRpools}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Starting caRpools eport generation from R console
}
\description{
Moreover, caRpools report generation can also be initiated without R-studio installation, so that this can be done via R command line even on remote computers.  
In this case, caRpools report generation can be started via `use.caRpools` with additional parameters, which are described below. 
}
\usage{
use.caRpools(type=NULL, file="CaRpools-extended-PDF.Rmd",
miaccs="MIACCS.xls", check=TRUE, work.dir=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{type}{
*Description* If you provide a custom Rmd template that can generate both, PDF and HTML reports you can indicate which version you want to generate.
*Default* NULL
*Values* "PDF", "HTML"
}
\item{file}{
*Description* The file name of your custom Rmd template file (with extension).
*Default* "CaRpools-extended-PDF.Rmd"
*Values* filename as character
}
\item{miaccs}{
*Description* The filename of your MIACCS file. 
*Default* "MIACCS.xls"
*Values* filename as character
}
\item{check}{
*Description* Indicates whether caRpools will check for correct installation and file access.
*Default* TRUE
*Values* TRUE or FALSE (boolean)
}
\item{work.dir}{
*Description* You can provide the absolute path to the working directory in which all files are placed (e.g. the MIACCS.xls and Rmd template).
*Default* NULL
*Values* absolute path (character) or NULL if standard R working directory is used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Start caRpools report generation, so no direct return value is generated.
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
#use.caRpools(check=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Report }% __ONLY ONE__ keyword per line
