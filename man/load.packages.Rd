\name{load.packages}
\alias{load.packages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loading and Installing packages used for caRpools
}
\description{
This function is used to check for presence of all packages and install them if not.
}
\usage{
load.packages(noupdate=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{noupdate}{
Indicates whether packages will NOT be updated, by default TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Is only used to check R packages
}
\value{
load.packages does not give any return value, however it will give you errors if something is wrong.
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
load.packages()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Analysis }% __ONLY ONE__ keyword per line
