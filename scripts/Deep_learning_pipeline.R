# Matthieu Genais
# created on friday, 23/05/2020
# CBI LBCMCP Legube Team, Toulouse
# This code is the pipeline used to create all output files from a given DeepG4 model in H5 file format
#it needs to launch tomtom analyses and matrix clustering

#this pipeline does the following steps
#- step1: uses activation_finder() to recover all positive activation of each sequence
#- step2: uses KernelToPWM() to get summarised PWM for each kernel in DeepG4 model
#- step3: uses PFMToTrimPFM() to trim each PWM (kernel Motif)
####run Matrix clustering tool (put the output file in the result file created by the script) http://rsat-tagc.univ-mrs.fr/matrix-clustering_form.cgi
#linkage : median, cor:0.6, Ncor:0.6
#- step4: correcting matrixClustering jaspar file format error in the (matrix-clustering_cluster_root_motifs.tf).
#The code changes actg from lower case to upper case (important for R and tomtom recognition of transfac format)
####run tomtom analysis on kernel KernelTrim and cluster http://meme-suite.org/tools/tomtom
#uses result/interpret_DeepG4_features/model/database/jasparTrimedSeuil0.9.jaspar
#uses result/interpret_DeepG4_features/model/database/MOTIF_FORMAT_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42
#uses result/interpret_DeepG4_features/model/database/MOTIF_FORMAT_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42
#uses result/interpret_DeepG4_features/model/database/rootMotifsJaspar
#uses also JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt database for the comparison
#- step5: random forest analysis (variable importance computation for each motif)

# setwd("/media/mourad/disk750Gb/G4pred_Matthieu/")

#required packages
require(PWMEnrich)
require(ranger)
require(yardstick)
require(cowplot)
require(gridExtra)
require(ggseqlogo)
# require(JASPAR2020)
source("/media/ElissarDisk/Elissar/Projects/DeepG4/G4pred_Matthieu/scriptsR/Deep_learning_pipeline_functions.R")
source("/media/ElissarDisk/Elissar/Projects/DeepG4/G4pred_Matthieu/scriptsR/Variable_importance_functions.R")
source("/media/ElissarDisk/Elissar/Projects/DeepG4/G4pred_Matthieu/scriptsR/Motif_Trimming_function.R")
'%ni%' <- Negate('%in%')
###Find activation zones (motifs)
activation_finder <- function(pathOneHot,pathModel,layerIndex=c(2),cores,PositiveSequence = TRUE,ATAC=F,file=NULL){

    ##One hot data is load
    OneHotDat= readRDS(pathOneHot)

    ###get model features from a h5 file
    model = load_model_hdf5(filepath = pathModel)
    #set kernel size
    kernel_size=lapply(layerIndex,function(u){
        dim(get_layer(model, index = u)$weights[[1]])[[1]]
    })
    #create the model for the convolutional output (the layerName is called)
    if(ATAC){
        inputs <- model$input[[1]]

    }else{
        inputs <- model$input
    }
    if(length(layerIndex)>1){
        outputs <- list(
            get_layer(model, index = layerIndex[1])$output,
            get_layer(model, index = layerIndex[2])$output
        )
    }else{
        outputs <- list(
            get_layer(model, index = layerIndex[1])$output
        )
    }


    model <- keras_model(inputs = inputs,
                         outputs = outputs)
    #Bin sequences
    # bin =c(seq(from = 1,to =dim(OneHotDat$test$x)[1],by =round(dim(OneHotDat$test$x)[1]/100)),dim(OneHotDat$test$x)[1]+1)

    ########################################################
    if(PositiveSequence == TRUE){

        #bind all positive sequences from RDS file:
        OneHotDat = abind(OneHotDat$train$x[OneHotDat$train$y==1,,],OneHotDat$test$x[OneHotDat$test$y==1,,],along = 1)

    }else{
        #bind all negative sequences from RDS file:
        OneHotDat = abind(OneHotDat$train$x[OneHotDat$train$y==0,,],OneHotDat$test$x[OneHotDat$test$y==0,,],along = 1)
    }

    batchSize =100
    #create batch dataset in a tensor object
    test_dataset <- tensor_slices_dataset(OneHotDat) %>% dataset_batch(.,batch_size = ceiling(dim(OneHotDat)[1]/batchSize))


    #create iteration
    iter = make_iterator_one_shot(test_dataset)

    #create final list
    final_motif_predict <- list()

    #starting batch index
    batch_index = 0

    #starting id sequence
    seq_id_start <- 1

    for(i in 1:batchSize){

        #start from batch 1 to the end
        batch_index = batch_index + 1

        #get the data id from the batch
        batch <-  iterator_get_next(iter)

        #get the data activation for this batch
        preds <- as.array(predict(model,batch,steps = 1))
        if(length(layerIndex) ==1){
            preds <- list(preds)
        }
        print(batch_index)

        lapply(1:length(layerIndex),function(j){

            #get the maximum activation for each sequence and kernel
            layer_max=apply(preds[[j]],c(1,3),max)

            #get index of all those maximum activations
            layer_whichmax=apply(preds[[j]],c(1,3),which.max)

            #for each activation >0 do
            pbmcapply::pbmclapply(1:ncol(layer_max),function(i){
                #ncol(conv1D_max)

                #get index of all non zero activations
                seqidx=which(layer_max[,i]>0) ## motif found in the sequence

                activationValue = layer_max[layer_max[,i]>0,i]

                #get the start index of each activation index
                startidx=layer_whichmax[seqidx,i] # idx start motif

                #get the end of each activation index
                endidx =layer_whichmax[seqidx,i]+(kernel_size[[j]]-1)

                #create a tibble containing all data
                tibble(seq = (seqidx+seq_id_start - 1),start = startidx,end = endidx,ActivationValue=activationValue, kernel_id = i)
            },mc.style ="ETA",mc.cores = 8) %>% bind_rows() %>% drop_na()
        }) %>% bind_rows(.id = "Layer") %>% write_tsv(file,append = T)

        #reset the id_start after each batch
        seq_id_start <- seq_id_start + nrow(preds[[1]])
    }

    return(read_tsv(file,col_names = c("Layer","seq","start","end","ActivationValue","kernel_id")))

    #bind all batch together
    # motifs <- final_motif_predict %>% bind_rows()
    #
    # return(motifs)
}

setwd("/home/disleg/Documents/Vincent/DeepG4/scripts/ATACseq_18012020/PIPELINE")
# models=list.files("results/runs",recursive = T,pattern="best_model.h5")
# names(models) <- dirname(models) %>% str_replace("/","_")
# models <- models[!str_detect(names(models),"03_03")]
# layerIndex <- list(
#     c(2),
#     c(2),
#     c(2)
# )
#
# names(layerIndex) <- names(models)
#UPDATE AU 28/04/21 -> we focus only on DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021/2021-04-19T16-32-14Z/best_model.h5
#
# models <- "DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021/2021-04-19T16-32-14Z/best_model.h5"
# names(models) <- dirname(models) %>% str_replace("/","_")
# layerIndex <- list(
#     c(2)
# )
#
# names(layerIndex) <- names(models)
#Update 12/07/21

# models <- "DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021/2021-07-02T08-48-55Z/best_model.h5"
models <- "DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021/2021-07-06T07-59-11Z/best_model.h5"
names(models) <- dirname(models) %>% str_replace("/","_")
layerIndex <- list(
    c(2)
)

names(layerIndex) <- names(models)


#UPDATE AU 02/03/21 -> we focus now only on DeepG4_ATAC_rescale_BW_sampling_01_03_2021_2021-03-02T08-51-20Z and DeepG4_classic_rescale_BW_sampling_01_03_2021_2021-03-02T08-41-25Z

# models <- models[c("DeepG4_ATAC_rescale_BW_sampling_01_03_2021_2021-03-02T08-51-20Z","DeepG4_classic_rescale_BW_sampling_01_03_2021_2021-03-02T08-41-25Z")]
# layerIndex <- layerIndex[c("DeepG4_ATAC_rescale_BW_sampling_01_03_2021_2021-03-02T08-51-20Z","DeepG4_classic_rescale_BW_sampling_01_03_2021_2021-03-02T08-41-25Z")]

for(model_name in names(models)){
    stranded <- str_detect(model_name,"stranded")
    message(model_name)
    if(stranded){
        ##OneHotPath
        pathOneHot = "data/rds/oneHot/Peaks_BG4_G4seq_GSE76688_HaCaT_hg19_201b_stranded_Ctrl_gkmSVM_0.8_42_rescale_BW_sampling_OneHot_train_test.rds"
        ##fastaPath
        fastaPath= "data/rds/Fasta/Peaks_BG4_G4seq_GSE76688_HaCaT_hg19_201b_stranded_Ctrl_gkmSVM_0.8_42_rescale_BW_sampling_Sequence_train_test.rds"


    }else{
        ##OneHotPath
        pathOneHot = "data/rds/oneHot/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_rescale_BW_sampling_OneHot_train_test.rds"
        ##fastaPath
        fastaPath= "data/rds/Fasta/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_rescale_BW_sampling_Sequence_train_test.rds"

    }

    ##modelPath
    # pathModel = paste0("data/model/",model,"/best_model.h5")
    # pathModel = str_c("results/runs/",models[[model_name]])
    pathModel = str_c("DVP_github/DeepG4ToolsComparison/rds/runs/",models[[model_name]])




    if(!dir.exists(paste0("results/interpret_DeepG4_features/",model_name,"/"))){
        dir.create(paste0("results/interpret_DeepG4_features/",model_name,"/"),recursive=T)
        dir.create(paste0("results/interpret_DeepG4_features/",model_name,"/raw/"),recursive=T)
        dir.create(paste0("results/interpret_DeepG4_features/",model_name,"/database/"),recursive=T)
        dir.create(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/"),recursive=T)
        dir.create(paste0("results/interpret_DeepG4_features/",model_name,"/plots/"),recursive=T)

    }

    fileoutkernel <- paste0("results/interpret_DeepG4_features/",model_name,"/raw/KernelActivation_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.tsv.gz")
    if(!file.exists(fileoutkernel)){
        #calling activation_finder function which allow us to get all activation sites for each kernel
        kernel_activation = activation_finder(pathOneHot = pathOneHot,
                                              pathModel = pathModel,
                                              layerIndex = layerIndex[[model_name]],
                                              core=20,
                                              PositiveSequence = T,
                                              ATAC = str_detect(model_name,"ATAC|dot"),
                                              file = fileoutkernel
        )
        # saveRDS(kernel_activation,file = )

    }else{
        kernel_activation <- read_tsv(fileoutkernel,col_names = c("Layer","seq","start","end","ActivationValue","kernel_id"))

    }


    #kernel list to PWM
    if(!file.exists(paste0("results/interpret_DeepG4_features/",model_name,"/raw/PFCWM_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.RDS"))){
        PWM_PFM_PCM_list= KernelToPWM(dataFrameMotif = kernel_activation,pathSequenceTrain_test =fastaPath,cores=10)
        saveRDS(object=PWM_PFM_PCM_list,file=paste0("results/interpret_DeepG4_features/",model_name,"/raw/PFCWM_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.RDS"))

    }else{
        PWM_PFM_PCM_list = readRDS(paste0("results/interpret_DeepG4_features/",model_name,"/raw/PFCWM_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.RDS"))

    }


    #output
    JasparOutput(path =paste0("results/interpret_DeepG4_features/",model_name,"/database/MOTIF_FORMAT_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42"),PWM_PFM_PCM_list$PCM)
    #create the MEME input
    memeOutput(paste0("results/interpret_DeepG4_features/",model_name,"/database/MOTIF_FORMAT_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.meme"),PWM_PFM_PCM_list$PFM)



    ###Triming Motifs
    seuil=0.9
    TrimedMatrixList = PFMToTrimPFM(PWM_PFM_PCM_list,seuil = seuil)

    #output trimmed matrix
    JasparOutput(path = paste0("results/interpret_DeepG4_features/",model_name,"/database/jasparTrimedSeuil",seuil),TrimedMatrixList$PCM)
    memeOutput(path = paste0("results/interpret_DeepG4_features/",model_name,"/database/jasparTrimedSeuil",seuil),PFM_list = TrimedMatrixList$PFM,fromJaspar = T)

}



###ADD BY VINCENT : PLOT PWM FOR EACH KERNEL
require(ggseqlogo)
for(model_name in names(models)){




    cc <- readJASPARMatrix(str_c("results/interpret_DeepG4_features/",model_name,"/database/MOTIF_FORMAT_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.jaspar"))

    cc.table <- names(cc) %>% enframe() %>% mutate(group = ntile(name,round(length(cc)/10)))
    pdf(str_c("results/interpret_DeepG4_features/",model_name,"/plots/MOTIF_FORMAT_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.pdf"),height=16,width=18)
    walk(split(cc.table,cc.table$group),function(zz){
        p.pcm <- lapply(cc[zz$value],function(x){as.matrix(x)})
        p.pcm <- lapply(names(p.pcm),function(i){ggseqlogo(p.pcm[[i]],ncol=1) + ggtitle(i)})
        p.pcm <- cowplot::plot_grid(plotlist = p.pcm)
        print(p.pcm)
    })
    dev.off()



    cc <- readJASPARMatrix(str_c("results/interpret_DeepG4_features/",model_name,"/database/jasparTrimedSeuil0.9.jaspar"))

    cc.table <- names(cc) %>% enframe() %>% mutate(group = ntile(name,round(length(cc)/10)))
    pdf(str_c("results/interpret_DeepG4_features/",model_name,"/plots/jasparTrimedSeuil0.9.pdf"),height=16,width=18)
    walk(split(cc.table,cc.table$group),function(zz){
        p.pcm <- lapply(cc[zz$value],function(x){as.matrix(x)})
        p.pcm <- lapply(names(p.pcm),function(i){ggseqlogo(p.pcm[[i]],ncol=1) + ggtitle(i)})
        p.pcm <- cowplot::plot_grid(plotlist = p.pcm)
        print(p.pcm)
    })
    dev.off()
}

#########Please run matrix clustering tool from the trimed motifs now



for(model_name in  names(models)){
    message(model_name)

    PWM_PFM_PCM_list = readRDS(paste0("results/interpret_DeepG4_features/",model_name,"/raw/PFCWM_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.RDS"))


    file = list.files(paste0("results/interpret_DeepG4_features/",model_name,"/matrix-clustering_archive/"))

    #correct error in matrix cluster file for Meme output, a c t g have to be A C T G in maj
    text = readLines(paste0("results/interpret_DeepG4_features/",model_name,"/matrix-clustering_archive/",file[[1]],"/matrix-clustering_cluster_root_motifs.tf"))
    text[text== "P0           a         c         g         t"] = "P0           A         C         G         T"
    write_lines(text,path = paste0("results/interpret_DeepG4_features/",model_name,"/matrix-clustering_archive/",file[[1]],"/matrix-clustering_cluster_root_motifs.tf"))

    rootMotifs = readMotifs(paste0("results/interpret_DeepG4_features/",model_name,"/matrix-clustering_archive/",file[[1]],"/matrix-clustering_cluster_root_motifs.tf"),remove.acc = T)



    rootMotifs = pbmcapply::pbmclapply(1:length(rootMotifs),function(i){
        PFMatrix(ID=names(rootMotifs[i]),name =names(rootMotifs[i]) ,matrixClass = "PCM",profileMatrix = rootMotifs[[i]])
    }) %>% c()

    rootMotifs=do.call(TFBSTools::PFMatrixList,rootMotifs)

    rootMotifsPFM = TFBSTools::toPWM(rootMotifs,type= 'prob')


    JasparOutput(PCM_list = rootMotifs,path = paste0("results/interpret_DeepG4_features/",model_name,"/database/rootMotifsJaspar"))
    memeOutput(path = paste0("results/interpret_DeepG4_features/",model_name,"/database/rootMotifsMeme"),PFM_list =rootMotifsPFM,fromJaspar = T )



    matrixList=TFBSTools::Matrix(rootMotifsPFM)
    #get Matrices from PWM list


}



##You have to compute tomtom on Kernel, KernelTrim and Cluster to go trough the rest of the pipeline
for(model_name in names(models)){
    PWM_PFM_PCM_list = readRDS(paste0("results/interpret_DeepG4_features/",model_name,"/raw/PFCWM_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.RDS"))
    stranded <- str_detect(model_name,"stranded")
    #not trimmed
    # kernelMatrix = PWM_PFM_PCM_list$PCM
    #
    # #trimmed Matrix
    seuil=0.9
    # TrimedMatrixList = PFMToTrimPFM(PWM_PFM_PCM_list,seuil = seuil)
    # Trimedmatrix = TrimedMatrixList$PCM

    # for(nameMot in c("Kernel","KernelTrim","Cluster")){
    for(nameMot in c("Cluster")){
        # for(nameMot in c("Kernel","KernelTrim")){

        if(!dir.exists(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot))){
            dir.create(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot))
        }

        if(nameMot=="Kernel"){
            MotifsList = readJASPARMatrix(paste0("results/interpret_DeepG4_features/",model_name,"/database/MOTIF_FORMAT_Positive_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.jaspar"))
            tomtomTrim = read_tsv(paste0("results/interpret_DeepG4_features/",model_name,"/tomtom",nameMot,"/tomtom.tsv"))

            # MotifsList =kernelMatrix
        }
        if(nameMot=="KernelTrim"){
            MotifsList = readJASPARMatrix(paste0("results/interpret_DeepG4_features/",model_name,"/database/jasparTrimedSeuil",seuil,".jaspar"))
            tomtomTrim = read_tsv(paste0("results/interpret_DeepG4_features/",model_name,"/tomtom",nameMot,"/tomtom.tsv"))

            # MotifsList =Trimedmatrix
        }
        if(nameMot=="Cluster"){
            MotifsList = readJASPARMatrix(paste0("results/interpret_DeepG4_features/",model_name,"/database/rootMotifsJaspar.jaspar"))
            tomtomTrim = read_tsv(paste0("results/interpret_DeepG4_features/",model_name,"/tomtomCluster/tomtom.tsv"))

            # MotifsList =rootMotifs
        }

        #First, plot IC by motif and in general (by VINCENT)

        ICM_summary <- MotifsList %>% lapply(as.matrix) %>% map(TFBSTools::toICM) %>% map(reshape2::melt) %>% bind_rows(.id = "kernel")
        p1 <- ICM_summary %>% ggplot(aes(x=value)) + geom_histogram(bins=500) +ggtitle("global IC distribution")
        p2 <- ICM_summary %>% group_by(kernel) %>% summarise(meanIC = mean(value)) %>% ggplot(aes(x=meanIC)) + geom_histogram(bins=50) +ggtitle("average IC by motifs")
        pdf(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/ICM_distribution_",nameMot,".pdf"),height=10,width=16)
        table <- summary(ICM_summary$value) %>% enframe() %>% tableGrob
        h <- grobHeight(table)
        title <- textGrob("IC Distribution by position", y=unit(0.5,"npc") + h,
                          vjust=0, gp=gpar(fontsize=20))

        gt <- gTree(children=gList(table, title))
        grid.draw(gt)
        grid.newpage()
        table <- ICM_summary %>% group_by(kernel) %>% summarise(meanvalue = mean(value)) %>% pull(meanvalue) %>% summary() %>% enframe() %>% tableGrob
        h <- grobHeight(table)
        title <- textGrob("IC distribution by motif (mean)", y=unit(0.5,"npc") + h,
                          vjust=0, gp=gpar(fontsize=20))

        gt <- gTree(children=gList(table, title))
        grid.draw(gt)
        print(p1)
        print(p2)
        dev.off()



        #get fasta to count motifs
        # fasta=readDNAStringSet("pipeline/results/201/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b/fasta/merged/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_merged.Fa")
        if(stranded){
            fasta_set= "data/rds/Fasta/Peaks_BG4_G4seq_GSE76688_HaCaT_hg19_201b_stranded_Ctrl_gkmSVM_0.8_42_rescale_BW_sampling_Sequence_train_test.rds" %>% readRDS()




        }else{

            fasta_set= "data/rds/Fasta/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_rescale_BW_sampling_Sequence_train_test.rds" %>% readRDS()


        }
        y <- c(fasta_set$train$y,fasta_set$test$y)
        fasta = c(fasta_set$train$x,fasta_set$test$x)
        train_test <- c(rep(1,length(fasta_set$train$y)),rep(0,length(fasta_set$test$y)))
        # y <- fasta_set$train$y
        # fasta = fasta_set$train$x
        message("ix")
        #matchMotifs function
        ix = motifmatchr::matchMotifs(MotifsList,fasta,out = "scores")

        #count number of motifs in each sequence
        mat = motifmatchr::motifCounts(ix)

        #change into dataframe
        dataKernel = as.data.frame(as.matrix(mat))
        colnames(dataKernel) = ID(MotifsList)


        #add id column
        # rownames(dataKernel)=names(fasta)
        # dataKernel$id = names(fasta)

        #create peak column
        dataKernel=dataKernel %>%
            mutate(id = names(fasta)) %>%
            # mutate(peak = if_else(grepl("pos",rownames(dataKernel)),1,0))
            mutate(peak = y)


        # fastaPos=readDNAStringSet("pipeline/data/samples/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b.Fa")
        fastaPos=fasta[y == 1]
        #matchMotifs function
        ix2 = motifmatchr::matchMotifs(MotifsList,fastaPos,out = "scores")

        #count number of motifs in each sequence
        mat2 = motifmatchr::motifCounts(ix2)

        #change into dataframe
        dataKernel2 = as.data.frame(as.matrix(mat2))
        colnames(dataKernel2) = ID(MotifsList)


        #add id column
        # rownames(dataKernel2)=names(fastaPos)
        # dataKernel2$id = names(fastaPos)

        #create peak column
        dataKernel2=dataKernel2 %>%
            mutate(id = names(fastaPos)) %>%
            mutate(peak = 1)


        #abundance kernel
        abundanceKernel = colSums(dataKernel2 %>% select(-id,-peak))

        # Shuffle
        # set.seed(123)
        # percTrain<-0.7
        # idxS <- sample(1:nrow(dataKernel))
        # dataKernel <- dataKernel[idxS,]

        #
        label <- dataKernel[,"peak"]
        dataKernel = dataKernel %>% dplyr::select(-peak)
        row.names(dataKernel)=dataKernel$id
        dataKernel = dataKernel %>% dplyr::select(-id)
        dataKernel <- dataKernel %>% as("Matrix")




        # Train
        x_train=dataKernel[train_test == 1,]
        y_train=label[train_test == 1]

        # Test
        x_test=dataKernel[train_test == 0,]
        y_test=label[train_test == 0]



        # Random forests
        dataRF=cbind(y_train,x_train)
        if(file.exists(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/RF_",nameMot))){
            RF = readRDS(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/RF_",nameMot))
        }else{
            RF=ranger(dependent.variable.name="y_train",data=dataRF,importance="impurity")
            saveRDS(RF,paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/RF_",nameMot))

        }


        # Predict
        ypred_prob= predict(RF,data=x_test)$predictions %>% as.data.frame()
        treshold <- 0.5
        ypred=ifelse(ypred_prob[,1] < treshold,0,1)
        cc <- evaluate_model(factor(y_test,levels = c(1,0)),ypred_prob,factor(ypred,levels = c(1,0)))



        impAll <- RF %>% print_variable_imp_rf(replace.name = NULL,top = RF$num.independent.variables)

        df=tibble(name = impAll$data$name,
                  importanceVariable=impAll$data$value,
                  abundanceKernel=abundanceKernel[match(impAll$data$name,names(abundanceKernel))])
        #idxDF=data.frame(name=KernelWithoutMotif,idx=idxKernelWithoutMotif)

        #df = left_join(df,idxDF,by=c("name"="name"))
        df_best <- df %>% filter(impAll$data$value>100)




        plot = df %>% ggplot(aes(x=abundanceKernel,y=importanceVariable,label=name))+geom_point()+
            scale_x_continuous(trans = "log10")+
            xlab(paste0("abundance_",nameMot))+
            ggrepel::geom_label_repel(data= df_best)


        p.pcm <- lapply(MotifsList[df_best$name],function(x){as.matrix(x)})
        names(p.pcm) <- str_remove(names(p.pcm),"_number_of.+")
        p.pcm <- ggseqlogo(p.pcm,ncol=1) + theme_void(base_size=10)
        pf <- cowplot::align_plots(p.pcm,plot,align = "hv",axis="bt")


        pdf(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/RF_importance_variable_",nameMot,"_abundance.pdf"),height=14,width=16)
        cowplot::plot_grid(pf[[2]],pf[[1]],ncol=2,rel_widths = c(0.85,0.15)) %>% print()
        dev.off()



        ########all other plots



        imp <- RF %>% print_variable_imp_rf(top = 25)
        nb_max <- ifelse(length(imp$data$name)<20,length(imp$data$name),20)
        pdf(width = 20,paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/ranger_",nameMot,".pdf"))
        cc[[1]] %>% gridExtra::grid.table()
        cc[[2]] %>% print()
        cc[[3]] %>% print()
        imp %>% print()

        lapply(MotifsList[paste0(imp$data$name[1:nb_max])],function(x){as.matrix(x)}) %>% ggseqlogo::ggseqlogo(ncol=5) %>% print()

        dev.off()

        ####For The Rapport Figure
        #std_coefs <- importance(RF) %>% enframe()
        #std_coefs <- std_coefs %>%
        #  mutate(value = abs(value)) %>%
        #  arrange(desc(value)) %>%
        #  mutate(name = forcats::fct_reorder(name,value)) %>%  dplyr::slice(1:25)
        #imp=std_coefs %>% ggplot(aes(x=name,y=value,fill=value)) +
        #  geom_bar(stat="identity") + scale_fill_viridis_c() +
        #  coord_flip()+
        #theme(legend.position = "none") +
        #  ggtitle(str_c("Top",25,"importance variables",sep=" ")) +
        #  xlab("Variables") +
        #  ylab("Importance")+
        #  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        #        strip.background =element_rect(colour = "black",size = 0,fill = "white"),
        #        strip.text.x = element_text(colour = "black", face = "bold",size = 20),
        #        strip.text.y = element_text(colour = "black", face = "bold",size = 12),
        #        axis.text.y.left = element_text(colour = "black",size =15),
        #        axis.title.y = element_text(size=20,face = "bold"),
        #        axis.title.x = element_text(size=20,face = "bold"),
        #        axis.text.x  = element_text(colour = "black",size =12),legend.position = "none",
        #        title = element_text(colour = "black",size =20))

        #ggsave(filename =paste0("results/interpret_DeepG4_features/",model,"/randomForest/tomtom",nameMot,"/Imp_",nameMot,".pdf"),
        #       plot=imp,height = 25,width = 7,units = "cm")


        #get SIgnificant Match to jaspar database
        MotifSignif=tomtomTrim %>%
            group_by(Query_ID) %>%
            arrange((`q-value`),.by_group=T) %>%
            dplyr::filter(`q-value`<0.05) %>%
            na.omit()

        #object not in Jaspar
        NotInMotifs = MotifsList[names(MotifsList)%ni%unique(MotifSignif$Target_ID)] %>% lapply(as.matrix)

        impNotMotif=imp$data$name[imp$data$name %in% names(NotInMotifs)] %>% paste()

        if(!is_empty(impNotMotif)){
            pdf(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/",nameMot,"NotInTomTom.pdf"))
            NotInMotifs[impNotMotif] %>% ggseqlogo::ggseqlogo() %>% print()

            dev.off()
        }


        #object  in Jaspar
        imp <- RF %>% print_variable_imp_rf(top = 150)

        signi =MotifSignif[MotifSignif$Target_ID%in% imp$data$name,]


        #if condition nameMot is "cluster", get the name of cluster for a  given kernel
        if(nameMot=="Cluster"){
            top1 =MotifSignif %>% group_by(Target_ID) %>% arrange((`q-value`),.by_group=T) %>% dplyr::slice(1)

            top1= top1 %>% group_by(Query_ID) %>% arrange((`q-value`),.by_group=T) %>% dplyr::slice(1)


            ##information cluster
            file = list.files(paste0("results/interpret_DeepG4_features/",model_name,"/matrix-clustering_archive/"))[[1]]


            files =list.files(paste0("results/interpret_DeepG4_features/",model_name,"/matrix-clustering_archive/",file,"/matrix-clustering_data/"))


            for(x in 1:(length(files[str_detect(files,"cluster_")]))){

                kernel = universalmotif::read_transfac(paste0("results/interpret_DeepG4_features/",model_name,"/matrix-clustering_archive//",file[[1]],"/matrix-clustering_data/",files[str_detect(files,"cluster_")][[x]]))
                print(x)

                if(x==1){
                    if(length(kernel) ==1){
                        name = kernel["altname"] %>% unlist()

                        df = tibble("ClusterID"=files[str_detect(files,"cluster_")][[x]],"Kernel" = name)
                    }else{
                        df=lapply(1:length(kernel),function(i){
                            name = kernel[[i]]["altname"] %>% unlist()

                            tibble("ClusterID"=files[str_detect(files,"cluster_")][[x]],"Kernel" = name)
                        }) %>% bind_rows()
                    }
                }else{
                    if(length(kernel) ==1){
                        name = kernel["altname"] %>% unlist()

                        df = bind_rows(df,tibble("ClusterID"=files[str_detect(files,"cluster_")][[x]],"Kernel" = name))
                    }else{
                        for(z in 1:length(kernel)){
                            name = kernel[[z]]["altname"] %>% unlist()
                            df = bind_rows(df,tibble("ClusterID"=files[str_detect(files,"cluster_")][[x]],"Kernel" = name))
                        }
                    }
                }

            }

            NameCluster = lapply(str_split(df$ClusterID,pattern = "_"), function(x){
                paste0(x[[1]],"_",x[[2]])
            })

            df$ClusterID=unlist(NameCluster)

            #get length of each cluster
            ClusterWeight = df %>% group_by(ClusterID) %>% summarise(n=dplyr::n())


            #plot Comparison between Jaspar database and Cluster rootmotif (with informaton of number kernel in motif)
            # opts <- list()
            # opts[["species"]] <- 9606 #human
            # opts[["all_versions"]] <- TRUE
            # MatrixList <- getMatrixSet(JASPAR2020@db, opts)
            MatrixList <- readJASPARMatrix("data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

            motifs <- lapply(MotifsList,function(x){as.matrix(x)})

            pdf(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/",nameMot,"Motif.pdf"))
            for( candidate in imp$data$name){
                print(candidate)
                plot1 = motifs[[candidate]] %>% ggseqlogo::ggseqlogo()+
                    ggtitle(paste0(candidate,"_",ClusterWeight$n[ClusterWeight$ClusterID==candidate],"_Kernel"))+
                    theme(panel.background = element_rect(fill = "white", colour = "white"),
                          strip.background =element_rect(colour = "white",size = 0,fill = "white"),
                          strip.text.x = element_text(colour = "black", face = "bold",size = 25),
                          strip.text.y = element_text(colour = "black", face = "bold",size = 25),
                          axis.text.y = element_text(colour = "black",size =25),
                          axis.title.y = element_text(size=25,face = "bold"),
                          axis.title.x = element_text(size=25,face = "bold"),
                          axis.text.x  = element_text(colour = "black",size =20),legend.position = "none",
                          title = element_text(colour = "black",size =20))
                if(sum(ID(MatrixList)%in%top1$Query_ID[top1$Target_ID==candidate])!=0 && sum(names(motifs)%in%candidate)!=0){

                    plot2 = TFBSTools::Matrix(MatrixList[ID(MatrixList)%in%top1$Query_ID[top1$Target_ID==candidate]]) %>% ggseqlogo::ggseqlogo()+
                        ggtitle(paste0(ID(MatrixList[ID(MatrixList)%in%top1$Query_ID[top1$Target_ID==candidate]]),"_",name(MatrixList[ID(MatrixList)%in%top1$Query_ID[top1$Target_ID==candidate]])))+
                        theme(panel.background = element_rect(fill = "white", colour = "white"),
                              strip.background =element_rect(colour = "white",size = 0,fill = "white"),
                              strip.text.x = element_text(colour = "black", face = "bold",size = 25),
                              strip.text.y = element_text(colour = "black", face = "bold",size = 25),
                              axis.text.y = element_text(colour = "black",size =25),
                              axis.title.y = element_text(size=25,face = "bold"),
                              axis.title.x = element_text(size=25,face = "bold"),
                              axis.text.x  = element_text(colour = "black",size =20),legend.position = "none",
                              title = element_text(colour = "black",size =20))
                    gridExtra::grid.arrange(plot1, plot2)
                }else{
                    gridExtra::grid.arrange(plot1, plot1)
                }
            }
            dev.off()
            motifs = readMotifs(file = paste0("results/interpret_DeepG4_features/",model_name,"/database/rootMotifsJaspar.jaspar"))
            imp <- RF %>% print_variable_imp_rf(top = 150)


            motifsCandidate = motifs[paste(imp$data$name)]


            for(candidate in names(motifsCandidate)){
                if(candidate%in%top1$Target_ID){
                    names(motifsCandidate)[names(motifsCandidate)==candidate]=paste0(names(motifsCandidate[candidate]),"_",
                                                                                     ID(MatrixList[ID(MatrixList)==top1$Query_ID[top1$Target_ID==candidate]]),"_",
                                                                                     name(MatrixList[ID(MatrixList)==top1$Query_ID[top1$Target_ID==candidate]]))

                }
            }


            #Plot Figure 4
            # motifsCandidate$cluster_18_MA0139.1_CTCF
            # plotfig4= c(motifsCandidate[1:10],motifsCandidate["cluster_98"],motifsCandidate["cluster_18_MA0139.1_CTCF"]) %>% ggseqlogo::ggseqlogo()+
            #   theme(panel.background = element_rect(fill = "white", colour = "white"),
            #         strip.background =element_rect(colour = "white",size = 0,fill = "white"),
            #         strip.text.x = element_text(colour = "black", face = "bold",size = 25),
            #         strip.text.y = element_text(colour = "black", face = "bold",size = 25),
            #         axis.text.y = element_text(colour = "black",size =25),
            #         axis.title.y = element_text(size=25,face = "bold"),
            #         axis.title.x = element_text(size=25,face = "bold"),
            #         axis.text.x  = element_text(colour = "black",size =20),legend.position = "none",
            #         title = element_text(colour = "black",size =10,face = "bold")) %>% print()
            #
            #
            #
            # ggsave(filename =paste0("results/interpret_DeepG4_features/",model,"/randomForest/tomtom",nameMot,"/",nameMot,"MotifFig4.pdf"),
            #        plot=plotfig4,height = 15,width = 15,units = "cm",device = "pdf")
            #
            pdf(paste0("results/interpret_DeepG4_features/",model_name,"/randomForest/tomtom",nameMot,"/",nameMot,"NotMotif.pdf"))
            for( candidate in names(motifs)[names(motifs)%ni%top1$Target_ID]){
                if(sum(ID(MotifsList)%in%top1$Query_ID[top1$Target_ID==candidate])==0 && sum(names(motifs)%in%candidate)!=0){
                    plot1 = motifs[candidate] %>% ggseqlogo::ggseqlogo() %>% print()
                }

            }
            dev.off()

        }


    }


}

