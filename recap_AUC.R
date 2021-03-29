require(tidyverse)
#peaks BG4G4seq
fres <- read_tsv("results/recap_AUC_4.tsv")
cc <- fres %>% filter(TypeExp %in% c("Peaks_G4P_G4seq","TestSet_Peaks_BG4_G4seq","Peaks_BG4_G4seq","Peaks_qG4")) %>% 
  unite(Files,TypeExp,Exp,cell_line,Ctrl) %>%
  mutate(Files = factor(Files,levels = rev(c("TestSet_Peaks_BG4_G4seq_GSE76688_HaCaT_Ctrl_gkmSVM",
                                             "Peaks_BG4_G4seq_GSE99205_HaCaT_Ctrl_gkmSVM",
                                             "Peaks_BG4_G4seq_GSE76688_HEKnp_Ctrl_gkmSVM",
                                             "Peaks_BG4_G4seq_GSE107690_K562_Ctrl_gkmSVM",
                                             "Peaks_qG4_breast-cancer-PDTX_qG4_Ctrl_gkmSVM")))) %>% 
  filter(!Tool %in% c("qparse_mean","qparse_max","quadparser_mean","quadparser_sum","DeepG4_G4seqpmBG4","DeepG4_qG4G4seq","G4hunter_max","G4hunter_mean","G4CatchAll_max","G4CatchAll_mean","gqrs_mapper_max","gqrs_mapper_sum")) %>%
  filter(!str_detect(Tool,"DeepG4Scan")) %>%
  mutate(Tool = str_remove(Tool,"_tsv|_sum|_mean|max")) %>% 
  group_by(Files) %>% 
  arrange(desc(AUC)) %>% mutate(Classement = 1:dplyr::n()) %>% 
  ungroup() %>% 
  mutate(my_text = str_c(Classement," (",round(AUC,3),")"))
cctext <- cc %>% group_by(Tool) %>% summarise(meanscore = mean(AUC),AUC = sum(AUC))
cc%>% 
  mutate(Files =as.factor(Files)) %>% 
  ggplot(aes(x=fct_reorder(Tool,AUC,mean),y=AUC,fill=Files)) +
  geom_bar(alpha=0.25,stat="identity",col="black") +
  theme_classic(base_size=22) +
  theme(legend.text=element_text(size=20),
        legend.direction = "vertical",
        legend.box = "horizontal",
        # legend.position = c(1,0.025),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.justification = c(0, 0)) +
  coord_flip() +
  geom_text(aes(label=my_text),size=6, position=position_stack(0.45)) +
  scale_size_manual(values= c(4,4.8))+
  ylab("sum(AUC)") + xlab("Tools") +ylim(c(0,5)) +
  geom_label(data=cctext,aes(x=Tool,y=AUC,label = round(meanscore,3)),fill="white",position = position_dodge(width = 1),hjust = -0.1,size=5) +
  scale_fill_brewer(palette = "Set1")-> p4
print(p4)



cc <- fres %>% 
  # filter(!Tool %in% c("ATACDeepG4_ATACtuningOH5","ATACDeepG4_ATAC","DeepG4Scan_qG4G4seq","DeepG4_qG4G4seq","DeepG4_BG4G4seq","DeepG4Scan_BG4G4seq","ATACDeepG4_classictuningOH5","ATACDeepG4_classic","DeepG4_G4seqpmBG4","DeepG4Scan_G4seqpmBG4")) %>% 
  filter(!Tool %in% c("ATACDeepG4_ATACtuningOH5","ATACDeepG4_ATAC","DeepG4Scan_qG4G4seq","DeepG4_qG4G4seq","DeepG4_BG4G4seq","DeepG4Scan_BG4G4seq","ATACDeepG4_classic_old","ATACDeepG4_classic","DeepG4_G4seqpmBG4","DeepG4Scan_G4seqpmBG4")) %>% 
  filter(!Tool %in% c("qparse_mean","qparse_max","quadparser_mean","quadparser_sum","DeepG4_G4seqpmBG4","DeepG4_qG4G4seq","G4hunter_max","G4hunter_mean","G4CatchAll_max","G4CatchAll_mean","gqrs_mapper_max","gqrs_mapper_sum")) %>%
  filter(TypeExp %in% c("Peaks_G4P_G4seq","TestSet_Peaks_BG4_G4seq","Peaks_BG4_G4seq","Peaks_qG4","TestSet_Peaks_G4seqpm_BG4")) %>%
  # unite(Files,TypeExp,Exp,cell_line,Ctrl) %>%
  unite(Files,Exp,cell_line) %>%
  group_by(Files) %>%
  arrange(desc(AUC)) %>% mutate(Classement = 1:dplyr::n()) %>%
  ungroup() %>%
  mutate(my_text = str_c(Classement," (",round(AUC,3),")")) %>%
  mutate(Tool = fct_reorder(Tool,AUC,mean))
cctext <- cc %>% group_by(Tool) %>% summarise(meanscore = mean(AUC),AUC = sum(AUC))

cc%>%
  mutate(Files =factor(Files,levels = rev(c("GSE76688_HaCaT","GSE99205_HaCaT","GSE107690_K562","GSE76688_HEKnp","GSE133379_293T","GSE133379_A549","GSE133379_H1975","GSE133379_HeLaS3","breast-cancer-PDTX_qG4")))) %>%
  ggplot(aes(x=fct_reorder(Tool,AUC,mean),y=AUC,fill=Files)) +
  geom_bar(alpha=0.25,stat="identity",col="black") +
  theme_classic(base_size=22) +
  theme(legend.text=element_text(size=20),
        legend.direction = "vertical",
        legend.box = "horizontal",
        # legend.position = c(1,0.025),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.justification = c(0, 0)) +
  coord_flip(clip="off") +
  geom_text(aes(label=my_text), position=position_stack(0.45), size=4) +
  ylab("sum(AUC)") + xlab("Tools") +
  geom_label(data=cctext,aes(x=Tool,y=AUC,label = round(meanscore,3)),fill="white",position = position_dodge(width = 1),hjust = -0.1,size=5) +
  scale_fill_brewer(palette = "Set1") + ggtitle("AUC on (BG4|qG4)-G4Seq datasets") -> p4
print(p4)
pdf("recap_AUC.pdf",height=8,width=1)
print(p4 + theme(legend.position ="none"))
dev.off()
