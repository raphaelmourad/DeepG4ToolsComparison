G4Hunter = "../scripts/Snakemake/G4Hunter.R" # for G4Hunter rule
G4Hunter_tsv = "scripts/Snakemake/G4Hunter_to_tsv.R" # for G4Hunter_to_tsv rule
G4Hunter_retrained_script = "scripts/Snakemake/G4Hunter_retrained.R" # for G4Hunter_retrained rule

#rule G4Hunter 
#this rule compute a Fasta file containing positive sequences control sequences to create G4hunter score
#this rule input one fasta file per run 
#this rule run the script "G4Hunter.R" in the script file
#this rule output a tsv file containing all sequence's score from the fasta file 
rule G4Hunter:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.G4hunter"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4hunter.benchmark"
  threads:30
  conda: "../envs/G4hunter.yaml"
  params:
    source="scripts/imports/G4Hunter/seekG4hunt.r",
    start_threshold = 1,
    end_threshold = 2,
    pas = 0.1
  script:
    G4Hunter

#rule G4Hunter_tsv_format
#take as input the output of G4Hunter rule and output one score per sequence using specific treshold
rule G4Hunter_tsv_format:
  input:
    OUT+"{sample}/predict/{sample}.G4hunter"
  output:
    OUT+"{sample}/predict/{sample}.G4hunter_{calc,mean|max|sum}"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4Hunter_{calc}.benchmark"
  conda: "../envs/G4hunter.yaml"
  params:
    calc="{calc}",
    score_selec="1.5"
  shell:
    "Rscript "+G4Hunter_tsv+" {input} {output} {params.calc} {params.score_selec}"

#rule G4HunterRF_retrained
#take as input the output of G4Hunter rule and output one score per sequence
#Use a random forest implementation for each treshold of G4Hunter output
rule G4HunterRF_retrained:
  input:
    OUT+"{sample}/predict/{sample}.G4hunter"
  output:
    score = OUT+"{sample}/predict/{sample}.G4hunter_retrained"
  params:
    calc="max",
    model_path = "rds/G4Hunter_retrained_ranger_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train.rds"
  threads:30
  conda: "../envs/G4hunter.yaml"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4HunterRF_retrained.benchmark"
  shell:
    "Rscript "+G4Hunter_retrained_script+" {params.model_path} {input} {output.score} {params.calc}"
  
