load.file = function(filename, header= TRUE, sep="\t", comment.char="", type = NULL)
{
  
  
  if(identical(type,"fastalib"))
  {
    # Read fasta reference file for sgRNA sequence list at the end
    # fasta file has      <IDENTIFIER
    # next line           SEQUENCE
    if("seqinr" %in% rownames(installed.packages()) == FALSE) {install.packages("seqinr")}
    #library("seqinr")
    file.raw <- seqinr::read.fasta(file = filename, 
                           seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)

    
    #str(attributes(file.raw))
    # make as data frame
    file = data.frame(
      ID = attributes(file.raw),
      sequence = unlist(file.raw),
      stringsAsFactors=FALSE)
  }
  
  # read XLSX file
  else if(identical(type,"xlsx"))
  {
    #library("xlsx")
    file = xlsx::read.xlsx (filename, sheetName="MIACCS", header=FALSE, stringsAsFactors=FALSE)
    
    # only take those with information
    file = file[which(file[,1] !=""),]
    
    #set identifiers to as rownames
    rownames(file) = file[,1]
  }
  
  
  else
  {
    file = read.table(filename,header=header, sep=sep, comment.char=comment.char)
  }
  
  return(file)
}