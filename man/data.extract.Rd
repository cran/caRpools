\name{data.extract}
\alias{data.extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting sgRNA information from NGS FASTQ files to create read-count files for caRpools Analysis
}
\description{
CaRpools offers two ways of providing CRISPR/Cas9 screening data.  
Either raw **read-count files** are directly used as described before, or read-count files are generated from NGS FASTQ files by extracting the 20 nt target sequence, mapping it against a reference library and extracting the read-count information for each sgRNA identifier.

In a first step, NGS FASTQ data is extracted and mapped against a reference library file using bowtie2.
}
\usage{
data.extract(scriptpath=NULL, datapath=NULL, fastqfile=NULL, extract = FALSE,
pattern = "default", machinepattern = "default", createindex = FALSE,
referencefile = NULL, mapping = FALSE, reversecomplement = FALSE,
threads = 1, bowtieparams = "", sensitivity = "very-sensitive-local", match = "perfect")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{scriptpath}{
Absolute path of the folder that contains `CRISPR-extract.pl` and `CRISPR-mapping.pl`  
*Default* NULL  
*Values* absolute path (character)  
}
\item{datapath}{
Absolute path of the folder that contains the data files (e.g. file.FASTQ)  
*Default* NULL  
*Values* absolute path (character)  
}
\item{fastqfile}{
Filename of FASTQ file WITHOUT .fastq extension  
*Default* NULL  
*Values* filename (character)  
}
\item{extract}{
Whether CRISPR-extract.pl is used to extract the 20 nt target sequence from the NGS reads using `pattern` 
*Default* FALSE  
*Values* TRUE, FALSE (boolean)  
}
\item{pattern}{
PERL regular Expression to extract 20 nt target sequence from NGS reads. Please see *extract pattern* in this manual for more information.  
*Default* Regular Expression (character)  
}
\item{machinepattern}{
Maschine ID of your Sequencing maschine. Used ot identify the read id.  
}
\item{createindex}{
Do you want caRpools to generate a bowtie2 index? Only necessary if `mapping=TRUE`.
*Default* FALSE  
*Values* TRUE, FALSE  
}
\item{referencefile}{
Filename of the library reference FASTA file, without extension. Is the same as bowtie2 file, if `createindex=TRUE`.  
}
\item{mapping}{
Indicates whether FASTQ files need to be mapped against `referencefile`/`bowtie2file`. FALSE by default.  
*Default* FALSE  
*Values* TRUE, FALSE  
}
\item{reversecomplement}{
Is the NGS sequence in reverse complement order?  
*Default* FALSE  
*Values* TRUE, FALSE  
}
\item{threads}{
How many threads can bowtie2 use for mapping? Only used if `mapping=TRUE`. Usually cores of CPU.  
*Default* 2  
*Values* any integer  
}
\item{bowtieparams}{
If you want to pass additional parameters to bowtie2.  
}
\item{sensitivity}{
You can djust the sensitivity of bowtie2 using this parameter. By default, bowtie2 is used in a very-sensitive-local setting. More information about different sensitivy parameters can be found at the [bowtie2 options](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#options).  
*Default* "very-sensitive-local"  
*Other options: very-fast, fast, sensitive, very-fast-local, fast-local, sensitive-local*  
}
\item{match}{
After bowtie2 mapping, the aligment is converted into read count files *filename_extracted-design.txt* and *filename_extracted-genes.txt*.  
You can indiciate how well the alignment must be in order to be used for generating the read count for each sgRNA.  
By default, this is set to *perfect*, which only employs a mapped read if the full 20 nt from the sequencing match perfectly to the sgRNA found in your library reference. The following options can be used:  

* __perfect__ - Read is used of all 20 nt from the sequencing are matching the target sequence given in the library reference
* __high__ - Read is used if at least 18 nt (starting from the PAM) are matching the target sequence in the reference
* __seed__ - Read is used if at least 14 nt (starting from the PAM) are a perfect match against the target sequence in the reference  
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
none
}
\value{
Returns file name for load.file().
Generated additional read-count files.
}

\author{
Jan Winter
}
\note{
%%  ~~further notes~~
Needs bowtie2 and PERL working. use check.caRpools() first.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(caRpools)
# fileCONTROL1 = data.extract(scriptpath="path.to.scripts",
# datapath="path.to.FASTQ", fastqfile="filename1", extract=TRUE,
# seq.pattern, maschine.pattern, createindex=TRUE,
# bowtie2file=filename.lib.reference, referencefile="filename.lib.reference", 
# mapping=TRUE, reversecomplement=FALSE, threads, bowtieparams,
#sensitivity="very-sensitive-local",match="perfect")  
# Now we can load the generated Read-Count file directly!
#CONTROL1 = load.file(paste(datapath, fileCONTROL1, sep="/")) # Untreated sample 1 loaded

# Don't forget the library reference
# libFILE = load.file( paste(datapath, paste(referencefile,".fasta",sep=""), sep="/"),
# header = FALSE, type="fastalib")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~NGS }% __ONLY ONE__ keyword per line
