################################################################################
# Requires the libraries "doMC", "foreach", "itertools", "xgboost" (0.4-4),    #
# "caret" and "plyr".                                                          #
# If not already installed in R, you can install those by typing:              #
# install.packages(c("doMC", "foreach", "itertools", "plyr", "caret"))         #
# Specific steps are needed to install the xgboost version 0.4-4, as detailed  #
# in the Readme file.                                                          #
# The default fastread==TRUE option in readfasta requires "data.table" library.#
################################################################################
#setwd("./lib")
#source("bitcompile.R")
#rm(list=ls())
#setwd("..")
setwd(dir = "/home/disleg/Documents/DeepG4/DeepG4/Quadron")
load("Quadron.lib")
setwd(dir = "/media/ElissarDisk/Elissar/Projects/G4pred/Data/Fasta")
print("NOTE: Loading Quadron core...", quote=FALSE)


Quadron(FastaFile= "test.fasta", 
        OutFile  = "out.txt",
        nCPU     = 4,
        SeqPartitionBy = 1000000)

#file.remove("Quadron.lib")
#rm(list=ls())
################################################################################


