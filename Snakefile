# Author Vincent ROCHER
# Snakemake pipeline who launch every tool used as comparison against DeepG4
# Each tool is configured in a rule and each rule generate an output, usually the score for each sequence in input fasta
# Each experiment has to be configured in the array "EXPERIMENTS" when you put the positive set (G4 sequences) and negative set (control sequences)



OUT = "results/" # Where the output are written
OUTFASTA  = "fasta/" # Where the fasta/accessibility are stored for each experiment
INBED = "bed/" # Where the fasta are stored
INATAC = "bigwig/" # Where the accessibility files are stored (in bigwig format)
# Experiment to be tested by our pipeline
EXPERIMENTS = {
"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM":{"CTRL":"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM","EXP":"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b"},
"Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b"},
"Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b"},
"Peaks_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_BG4_G4seq_K562_GSE107690_hg19_201b"},
"Peaks_G4P_G4seq_GSE133379_293T_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4P_G4seq_GSE133379_293T_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4P_G4seq_GSE133379_293T_hg19_201b"},
"Peaks_G4P_G4seq_GSE133379_A549_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4P_G4seq_GSE133379_A549_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4P_G4seq_GSE133379_A549_hg19_201b"},
"Peaks_G4P_G4seq_GSE133379_H1975_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4P_G4seq_GSE133379_H1975_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4P_G4seq_GSE133379_H1975_hg19_201b"},
"Peaks_G4P_G4seq_GSE133379_HeLaS3_hg19_201b_Ctrl_gkmSVM":{"CTRL":"Peaks_G4P_G4seq_GSE133379_HeLaS3_hg19_201b_Ctrl_gkmSVM","EXP":"Peaks_G4P_G4seq_GSE133379_HeLaS3_hg19_201b"}
}



ATACFILE = {
	"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM":["ATAC_entinostat_mean.bw"],
	"Peaks_BG4_G4seq_HaCaT_GSE99205_hg19_201b_Ctrl_gkmSVM":["ATAC_entinostat_mean.bw"],
	"Peaks_BG4_G4seq_HEKnp_GSE76688_hg19_201b_Ctrl_gkmSVM":[
		"rhh113_HEKnp_ATAC_701_517_24022015_normalized.bw",
		"rhh114_HEKnp_ATAC_702_502_24022015_normalized.bw",
		"rhh115_HEKnp_ATAC_703_503_24022015_normalized.bw",
		"rhh116_HEKnp_ATAC_701_517_27032015_normalized.bw",
		"rhh117_HEKnp_ATAC_702_502_27032015_normalized.bw",
		"rhh118_HEKnp_ATAC_704_504_27032015_normalized.bw"
	],
	"Peaks_BG4_G4seq_K562_GSE107690_hg19_201b_Ctrl_gkmSVM":["GSM4133303_YET96_ATAC_K652.bw"],
	"Peaks_G4P_G4seq_GSE133379_293T_hg19_201b_Ctrl_gkmSVM":["ENCFF716SFD.bw"],
	"Peaks_G4P_G4seq_GSE133379_A549_hg19_201b_Ctrl_gkmSVM":["ENCFF180FXV.bw"],
	"Peaks_G4P_G4seq_GSE133379_H1975_hg19_201b_Ctrl_gkmSVM":[
		"GSM4217852_WT-rep1-ATAC.bw",
		"GSM4217853_WT-rep2-ATAC.bw"

	],
	"Peaks_G4P_G4seq_GSE133379_HeLaS3_hg19_201b_Ctrl_gkmSVM":["SRX2370816.bw"]
}


#A function built to get the positive set from the dictionnary EXPERIMENTS
def GetExp(wildcards):
  return(IN+EXPERIMENTS[wildcards.sample]["EXP"]+".Fa")

#A function built to get the negative set from the dictionnary EXPERIMENTS
def GetCtrl(wildcards):
  return(IN+EXPERIMENTS[wildcards.sample]["CTRL"]+".Fa")
#SAME BUT WITH BED
#A function built to get the positive set from the dictionnary EXPERIMENTS
def GetExpBed(wildcards):
  return(INBED+EXPERIMENTS[wildcards.sample]["EXP"]+".bed")

#A function built to get the negative set from the dictionnary EXPERIMENTS
def GetCtrlBed(wildcards):
  return(INBED+EXPERIMENTS[wildcards.sample]["CTRL"]+".bed")

#A function built to get the atac-seq/dnase-seq file for each set
def GetATAC(wildcards):
	return([INATAC+x for x in ATACFILE[wildcards.sample]])


#Scripts locations


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
include: "rules/NewDeepG4.smk"
include: "rules/AUC.smk"


rule all:
  input:
    expand(OUT+"{sample}/{metrics}/{sample}.{tools}.tsv",metrics = ["other_metrics","AUC"],
    sample=EXPERIMENTS.keys(),
    tools = ["ATACDeepG4_classictuningOH5","ATACDeepG4_ATACnormBG","ATACDeepG4_ATACnorm900BG","G4detector_tsv","G4detector_retrained_tsv","penguinn","penguinn_retrained"]),
    expand(OUT+"{sample}/AUC/{sample}.{tools}.tsv",metrics = ["other_metrics","AUC"],sample=EXPERIMENTS.keys(),
    # tools = ["qparse_mean","qparse_sum","qparse_max","quadron_score","quadron_retrained","pqsfinder","quadparser_max","quadparser_mean","quadparser_sum","G4CatchAll_sum","G4CatchAll_mean","G4CatchAll_max","G4hunter_max","G4hunter_mean","G4hunter_sum","G4hunterRF","gqrs_mapper_max","gqrs_mapper_sum","gqrs_mapper_mean"]) #
    tools = ["quadron_score","quadron_retrained"]) #
  output:
    expand(OUT+"recap_AUC_{nb}.tsv",nb=[1,2,3,4])
  conda : "envs/AUC.yaml"
  script:
    "scripts/Snakemake/report_AUC.R"





