library("BSgenome.Hsapiens.UCSC.hg19")
library("rtracklayer")
library("dplyr")
library("Biostrings")
#https://www.pnas.org/content/pnas/suppl/2019/09/20/1909047116.DCSupplemental/pnas.1909047116.sapp.pdf
hg19 = BSgenome.Hsapiens.UCSC.hg19
motif_MYC <-DNAString("GGGTGGGGAGGGTGGG")
my_fa <- "/home/rochevin/Téléchargements/gkw006_Supplementary_Data/cc.Fa" %>% readDNAStringSet()
#### On enregistre le nom des differentes séquences (chromosomes) contenu dans le génome
chromosomes = seqnames(hg19)[1:23]

res <- chromosomes %>% lapply(function(chromosome){
  chr_hg19 = hg19[[chromosome]]
  mclapply(my_fa,function(my_seq){
    cc = matchPattern(my_seq, chr_hg19, max.mismatch=0)
    cc2 = matchPattern(reverseComplement(my_seq), chr_hg19, max.mismatch=0)
    tibble(chromosome,start = c(start(cc),start(cc2)),end = c(end(cc),end(cc2)))
  },mc.cores=10) %>% setNames(names(my_fa)) %>% bind_rows(.id = "id") 
  
}) %>% bind_rows()


motif_MYC <- DNAString("GGGTGGGGAGGGTGGG")
res_MYC <- mclapply(chromosomes,function(chromosome){
  chr_hg19 = hg19[[chromosome]]
  cc = matchPattern(motif_MYC, chr_hg19, max.mismatch=0)
  cc2 = matchPattern(reverseComplement(motif_MYC), chr_hg19, max.mismatch=0)
  tibble(chromosome,start = c(start(cc),start(cc2)),end = c(end(cc),end(cc2)))
},mc.cores=10) %>% mutate(id ="MYC_GGGTGGGGAGGGTGGG") %>% mutate(value = length(motif_MYC))

#rtemove not g4
G4_list <- my_fa[1:298] %>% lapply(length) %>% unlist() %>% enframe()

res.G4 <- res %>% filter(id %in% G4_list$name) %>% left_join(G4_list,by = c("id"="name"))

#Get G4 in oncogenes

targets <- c("ENSG00000136997","ENSG00000164362","ENSG00000141510","ENSG00000213281","ENSG00000133703","ENSG00000157404")

data(genesymbol, package = "biovizBase")

stretch_me_up <- . %>% anchor_start() %>% stretch(3000) %>% anchor_end() %>% stretch(3000)


mes_genes <- genesymbol %>% filter(ensembl_id %in% targets) %>% stretch_me_up()

res.G4.f <- rbind(res.G4,res_MYC) %>% rename(seqnames = chromosome) %>% as_granges() %>% filter_by_overlaps(mes_genes)
res.G4.f %>% as_tibble() %>% write_tsv("motifs_G4Hunter_papier_overlap_oncogenes_hg19.bed",col_names=F)
