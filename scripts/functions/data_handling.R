#Split data based on the index
#Return a tibble with 
##Type : training/test
##name : index
##value : sequence name
##label : label sequence (0/1)
ManualSplitIndex <- function(Index=NULL,percTrain = 0.8,sampling=T){
  Index <- Index %>% enframe()
  if(sampling){
    Index <- Index[sample(1:nrow(Index)),]
  }
  Index <- Index %>% initial_split(prop = percTrain)
  
  return(
    list("training" = training(Index),
         "test" = testing(Index)) %>% bind_rows(.id = "Type") 
  )
}

#Transform a character chain into a vector of integer values
DNAToNumerical <- function(x){
  maxicall <- x %>% width %>% as.vector() %>% max()
  x <- x %>% str_replace_all(c("N"="5","T"="4","G"="3","C"="2","A"="1"))%>% str_split("")
  x %>% lapply(function(x){as.numeric(x)}) %>% do.call("rbind",.)
}
# Function to convert sentence (like DNA) into matrix of size len*dictionnary
convertOneHot<-function(matin,tabv = NULL){
  if(is.null(tabv)){
    tab=table(matin)
    tabv=as.numeric(names(tab))
  }
  listMat=list()
  for(i in 1:length(tabv)){
    mat=matrix(0,nrow(matin),ncol(matin))
    mat[matin==tabv[i]]=1
    listMat[[i]]=mat
  }
  arrayout=array(unlist(listMat), dim = c(nrow(listMat[[1]]), ncol(listMat[[1]]), length(listMat)))
  return(arrayout)
}
# Function to convert a string char into kmers
Kmerisation <- function(sentence,kmer_size = 3){
  res <- lapply(1:kmer_size,
         function(i){
           str_extract_all(substring(sentence,first = i)
                           ,"...",
                           simplify=F) %>%
             unlist()
         }
  ) %>% transpose() %>% unlist()
  return(res)
}
#USe motifmatchR to get motifs occurancy over a DNAStringSet, given a PFMMatrix Set
ToMotif <- function(MyDNAStringSet=NULL,PFMatrixList=NULL){
  motifmatchr::matchMotifs(PFMatrixList,MyDNAStringSet,out="scores",p.cutoff=5e-5) %>% 
    motifmatchr::motifCounts()
}

#Get score from bigwig data based on GRanges
getScoreBW <- function (one.w, x) 
{
  require(magrittr)
  lapply(split(x, droplevels(seqnames(x))), function(zz) {
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- IRanges::Views(cov, start = start(zz), end = end(zz)) %>% 
      as.matrix()
  }) %>% do.call("rbind",.)
}

#For the hyper parameter tuning part
# list runs witin the specified runs_dir
#Because AUC is not directly computed, we can extracted it using regex on the script output
extract_AUC <- function(output){
  output %>% str_extract("roc_auc\\s+binary.+\n") %>% str_extract("[0-9]\\.[0-9]+") %>% as.numeric()
}
extract_acc <- function(output){
  output %>% str_extract("accuracy\\s+binary.+\n") %>% str_extract("[0-9]\\.[0-9]+") %>% as.numeric()
}
#Given a directory of runs, extract informations, and AUC (extract_AUC) and return a tibble
get_run_val <- function(my_dir){
  ls_runs(runs_dir = my_dir) %>% 
    as_tibble() %>% 
    mutate(AUC = map_dbl(output, extract_AUC)) %>% 
    mutate(acc_test = map_dbl(output, extract_acc))
}

