require(tidyverse)
setwd("/home/disleg/Documents/Vincent/DeepG4/scripts/ATACseq_18012020/PIPELINE/results/AUC_AUPR/AUPR/")
namefile <- list.files(recursive=T,pattern=".tsv") 
names(namefile) <- namefile

namefile %>% map(read_tsv) -> cc
cc <- cc %>% bind_rows(.id = "Model") %>% mutate(Model = str_remove(Model,"_aupr.tsv")) %>% 
  separate(Model,into = c("file","Model"),sep="/")
cc <- cc %>% group_by(file,Model) %>% arrange(desc(AUC)) %>% dplyr::slice(1) %>% ungroup()
cc <- cc %>% 
  mutate(GSE = str_extract(file,"GSE[0-9]+|breastCancer")) %>% 
  mutate(stranded = ifelse(str_detect(file,"stranded"),"stranded","unstranded")) %>% 
  mutate(cell_line = str_extract(file,"HeLaS3|H1975|A549|293T|K562|HEKnp|HaCaT|breastCancer")) %>% 
  dplyr::select(-file) %>% unite(file,c(GSE,cell_line,stranded)) %>% 
  mutate(Model = str_replace(Model,"_rescale.+","_same_hyperparams"))

cc2 <- cc%>% 
  group_by(Model) %>% 
  arrange(desc(AUC)) %>% mutate(Classement = 1:dplyr::n()) %>% 
  ungroup() %>% 
  mutate(my_text = str_c(Classement," (",round(AUC,3),")")) %>% 
  mutate(file =as.factor(file)) %>% filter(!str_detect(file,"_stranded") )


cctext <- cc2  %>% group_by(Model) %>% summarise(meanscore = mean(AUC),AUC = sum(AUC))


p2 <- cc2 %>% 
  ggplot(aes(x=fct_reorder(file,AUC,mean),y=AUC,fill=Model)) +
  geom_bar(alpha=1,stat="identity",col="black") +
  theme_classic(base_size=22) +
  coord_flip() +
  theme(legend.text=element_text(size=20),
        legend.direction = "vertical",
        legend.box = "horizontal",
        # legend.position = c(1,0.025),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.justification = c(0, 0)) +
  geom_text(aes(label=my_text),size=6, position=position_stack(0.45)) +
  scale_fill_brewer(palette = "Set1") +
  ylab("sum(AUC)") + xlab("Models")  +
  geom_label(data=cctext,aes(x=Model,y=AUC,label = round(meanscore,3)),fill="white",position = position_dodge(width = 1),hjust = -0.1,size=5) 

