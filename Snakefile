# Author Vincent ROCHER
# Snakemake pipeline who launch every tool used as comparison against DeepG4
# Each tool is configured in a rule and each rule generate an output, usually the score for each sequence in input fasta
# Each experiment has to be configured in the array "EXPERIMENTS" when you put the positive set (G4 sequences) and negative set (control sequences)



OUT = "results/" # Where the output are written
IN = "fasta/" # Where the fasta are stored

# Experiment to be tested by our pipeline
EXPERIMENTS = {
"Peaks_G4seqpm_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4seqpm_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4seqpm_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_hg19_201b"},
"Peaks_G4seqpm_BG4_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4seqpm_BG4_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4seqpm_BG4_HaCaT_GSE76688_hg19_201b"},
"Peaks_G4seqpm_BG4_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4seqpm_BG4_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4seqpm_BG4_HaCaT_GSE99205_hg19_201b"},
"Peaks_G4seqpm_BG4_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4seqpm_BG4_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4seqpm_BG4_HEKnp_GSE76688_hg19_201b"},
"Peaks_G4seqpm_BG4_K562_GSE107690_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4seqpm_BG4_K562_GSE107690_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4seqpm_BG4_K562_GSE107690_hg19_201b"},

"Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b"},
"Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b"},
"Peaks_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_BG4_G4seq_K562_GSE107690_hg19_201b"},
"Peaks_G4seq_BG4_HaCaT_GSE76688_hg19_Ctrl_gkmSVM_201b":{"CTRL":"Peaks_G4seq_BG4_HaCaT_GSE76688_hg19_Ctrl_gkmSVM_201b","EXP":"Peaks_G4seq_BG4_HaCaT_GSE76688_hg19_201b"},
"Peaks_G4seq_BG4_HaCaT_GSE99205_hg19_Ctrl_gkmSVM_201b":{"CTRL":"Peaks_G4seq_BG4_HaCaT_GSE99205_hg19_Ctrl_gkmSVM_201b","EXP":"Peaks_G4seq_BG4_HaCaT_GSE99205_hg19_201b"},
"Peaks_G4seq_BG4_HEKnp_GSE76688_hg19_Ctrl_gkmSVM_201b":{"CTRL":"Peaks_G4seq_BG4_HEKnp_GSE76688_hg19_Ctrl_gkmSVM_201b","EXP":"Peaks_G4seq_BG4_HEKnp_GSE76688_hg19_201b"},
"Peaks_G4seq_BG4_K562_GSE107690_hg19_Ctrl_gkmSVM_201b":{"CTRL":"Peaks_G4seq_BG4_K562_GSE107690_hg19_Ctrl_gkmSVM_201b","EXP":"Peaks_G4seq_BG4_K562_GSE107690_hg19_201b"},
"Peaks_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_253b_Ctrl_gkmSVM":{"CTRL":"Peaks_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_253b_Ctrl_gkmSVM","EXP":"Peaks_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_253b"},
"Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4seq_qG4_breastCancer_qG4-ChIP-seq-of-breast-cancer-PDTX_hg19_201b"},


"Promoters_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_253b_Ctrl":{"CTRL":"Promoters_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_253b_Ctrl","EXP":"Promoters_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_253b"},
"Promoters_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl","EXP":"Promoters_BG4_G4seq_HaCaT_GSE76688_hg19_201b"},
"Promoters_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl","EXP":"Promoters_BG4_G4seq_HaCaT_GSE99205_hg19_201b"},
"Promoters_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl","EXP":"Promoters_BG4_G4seq_HEKnp_GSE76688_hg19_201b"},
"Promoters_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl","EXP":"Promoters_BG4_G4seq_K562_GSE107690_hg19_201b"},



"Promoters_BG4_G4seq_HaCaT_GSE76688_hg19_500b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_HaCaT_GSE76688_hg19_500b_Ctrl","EXP":"Promoters_BG4_G4seq_HaCaT_GSE76688_hg19_500b"},
"Promoters_BG4_G4seq_HaCaT_GSE99205_hg19_500b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_HaCaT_GSE99205_hg19_500b_Ctrl","EXP":"Promoters_BG4_G4seq_HaCaT_GSE99205_hg19_500b"},
"Promoters_BG4_G4seq_HEKnp_GSE76688_hg19_500b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_HEKnp_GSE76688_hg19_500b_Ctrl","EXP":"Promoters_BG4_G4seq_HEKnp_GSE76688_hg19_500b"},
"Promoters_BG4_G4seq_K562_GSE107690_hg19_500b_Ctrl":{"CTRL":"Promoters_BG4_G4seq_K562_GSE107690_hg19_500b_Ctrl","EXP":"Promoters_BG4_G4seq_K562_GSE107690_hg19_500b"},
"Promoters_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_500b_Ctrl":{"CTRL":"Promoters_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_500b_Ctrl","EXP":"Promoters_qG4-ChIP-seq-of-breast-cancer-PDTX_breastCancer_G4seq_hg19_500b"},
"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Ctrl_gkmSVM":{"CTRL":"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Ctrl_gkmSVM","EXP":"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42"}
}

#A function build to get the positive set from the dictionnary EXPERIMENTS
def GetExp(wildcards):
  return(IN+EXPERIMENTS[wildcards.sample]["EXP"]+".Fa")
  
#A function build to get the negative set from the dictionnary EXPERIMENTS
def GetCtrl(wildcards):
  return(IN+EXPERIMENTS[wildcards.sample]["CTRL"]+".Fa")

#Scripts locations


G4iMGrinder = "scripts/Snakemake/G4iMGrinder.R"

include: "rules/fasta.smk"
include: "rules/G4hunter.smk"
include: "rules/penguinn.smk"
include: "rules/pqsfinder.smk"
include: "rules/quadparser.smk"
include: "rules/quadron.smk"
include: "rules/gqrs_mapper.smk"
include: "rules/G4CatchAll.smk"
include: "rules/G4detector.smk"
include: "rules/qparse.smk"
include: "rules/DeepG4.smk"
include: "rules/AUC.smk"


rule all:
  input:
    expand(OUT+"{sample}/AUC/{sample}.{tools}.tsv",sample=EXPERIMENTS.keys(),tools = ["DeepG4_G4seqpmBG4","DeepG4","DeepG4_scan","G4detector_tsv","G4detector_retrained_tsv","penguinn","penguinn_retrained","qparse_mean","qparse_sum","qparse_max","quadron_score","quadron_retrained","pqsfinder","quadparser_max","quadparser_mean","quadparser_sum","G4CatchAll_sum","G4CatchAll_mean","G4CatchAll_max","G4hunter_max","G4hunter_mean","G4hunter_sum","G4hunter_retrained"]) # ,"gqrs_mapper_max","gqrs_mapper_sum","gqrs_mapper_mean"
  output:
    expand(OUT+"recap_AUC_{nb}.tsv",nb=[1,2,3,4])
  conda : "envs/AUC.yaml"  
  script:
    "scripts/Snakemake/report_AUC.R"

rule G4iMGrinder:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.G4iMGrinder"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4iMGrinder.benchmark"
  threads:30
  conda:
    "envs/Rv4.yaml"
  script:
    G4iMGrinder







