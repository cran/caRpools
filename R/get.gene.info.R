get.gene.info = function(data, namecolumn=1, extractpattern=expression("^(.+?)(_.+)"), database="ensembl", dataset="hsapiens_gene_ensembl", filters="ensembl_gene_id", attributes = c("hgnc_symbol"), return.val = "dataset", controls=FALSE)
{
#   if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
#     source("http://bioconductor.org/biocLite.R")
#     biocLite("biomaRt",suppressUpdates=update , ask=FALSE)
#   }
#requireNamespace(biomaRt)
  
if(!is.null(data) && nrow(data)>=1)
{
    
  
  # we apply this to an sgRNA dataset, not single genes
  if(return.val=="dataset")
  {
    # check for attributes and filters to be only one, attributes must be also the filter
    if (length(filters)!=1) stop("None or more than one single filter selected")
    
    # if dataset is control
    if(identical(controls, TRUE))
    {
      # get gene names from controls
      
      gene.names = data[,namecolumn]
    }
    else
    {
      # get gene names from design file
      gene.names = sub(extractpattern,"\\1",data[,namecolumn],perl=TRUE)
    }
    
    
    #start biomaRt interface and check if biomaRt is available
    handling = biomaRt::useMart(database)
    if(!exists("handling"))
    {stop("biomaRt connection is not working. This can be a connectivity issue (e.g. proxy settings, internet connection) or the biomaRt service is currently not avaible. \n
          You can skip any data annotation by setting it to FALSE in the MIACCS file.")}
    handling = biomaRt::useDataset(dataset,mart=handling)
    gene.info = biomaRt::getBM(
      filters=filters,
      attributes= c(filters,attributes),
      values= gene.names,
      mart = handling)
    
    if(nrow(gene.info) >= 1)
    {
      data$replace = as.character(gene.names)
      
      #print(data$replace)
      #gene.info.ext <-gene.info
      #print(gene.info)
      for(i in 1:nrow(gene.info))
      {
        
        #print(gene.info[i,])
        #add replace information to dataset
        if(gene.info[i,2] != "")
        {
          # if dataset is control
          if(identical(controls, TRUE))
          {
            data[data[,namecolumn] == gene.info[i,1],"replace"] = gene.info[i,2]
          }
          else
          {
            data[sub(extractpattern,"\\1",data[,namecolumn],perl=TRUE) == gene.info[i,1],"replace"] = gene.info[i,2]
          }
          
          # test2 <<- data[sub(extractpattern,"\\1",data[,namecolumn],perl=TRUE) == gene.info[i,1],"replace"]
        }
        else
        {
          #do nothing
        }
      }
      
      #Overwrite gene identifier by what is listed as attribute
      # if empty, give them the original name back
      if(identical(controls,TRUE))
      {
        data[,namecolumn] = as.factor(data[,"replace"])
        #data[,namecolumn] = as.factor(data[,"replace"])
      }
      else
      {
        data[,namecolumn] = apply(data,1, function(x){
          
          as.factor(paste(as.character(x["replace"]),sub(extractpattern,"\\2",as.character(x[namecolumn]),perl=TRUE),sep=""))
          
          
        }
        )
          #as.factor(paste(data[,"replace"],sub(extractpattern,"\\2",data[,namecolumn],perl=TRUE),sep=""))
        #data[,namecolumn] = as.factor(paste(data[,"replace"],sub(extractpattern,"\\2",data[,namecolumn],perl=TRUE),sep=""))
      }
      
      data$replace=NULL
      return(data)
      
    }
    else
    {
      return(as.data.frame(data))
    }
    
     
  }
  if(return.val=="info")
  {
    # We can enrich gene information based on what is provided for biomaRt
    # get gene names from design file
    # if dataset is control
    if(identical(controls, TRUE))
    {
      # get gene names from controls
      gene.names = data[,namecolumn]
    }
    else
    {
      # get gene names from design file
      gene.names = sub(extractpattern,"\\1",data[,namecolumn],perl=TRUE)
    }
    #start biomaRt interface
    handling = biomaRt::useMart(database)
    handling = biomaRt::useDataset(dataset,mart=handling)
    gene.info = biomaRt::getBM(
      filters=filters,
      attributes= c(filters,attributes),
      values= gene.names,
      mart= handling)
    
    cols = ncol(gene.info)
    gene = aggregate.data.frame(gene.names,by=list(gene.names), function(x) return(x[1]))
    gene$Group.1=NULL
    data.return=data.frame(
      gene = as.character(gene$x),
      stringsAsFactors=FALSE)
    
    
    data.return[,colnames(gene.info)] = NA
    for(m in 1:cols)
    {
      for(i in 1:nrow(gene.info))
      {
        
        data.return[data.return[,namecolumn]== gene.info[i,1],m] = gene.info[i,m]
       
      }
      #data.return = cbind.data.frame(data.return, gene.info[,1:cols])
    }
  
    colnames(data.return) = colnames(gene.info)
    data.return[,is.na(colnames(data.return))] = NULL
    return(data.return)
    
  }
  
}
  
}