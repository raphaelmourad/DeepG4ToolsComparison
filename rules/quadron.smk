generate_quadron_dat = "scripts/Snakemake/generate_quadron_dat.R" # for quadron rule
quadron_retrained_script = "scripts/Snakemake/quadron_retrained_score.R" # for quadron_retrained rule

#rule quadron 
#this rule compute a Fasta file containing positive sequences control sequences to create Quadron score
#this rule input one fasta file per run 
#this rule run the script "generate_quadron_dat.R" in the script file
#this rule output a tsv.gz file containing all sequence's scores from the fasta file 
rule quadron_score:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    quadron_table = OUT+"{sample}/predict/{sample}.quadron_score"#tsv.gz
  params:
    type = "score",
    source="scripts/Snakemake/Quadron_add_Matthieu.R",
    Quadron_lib = "scripts/imports/Quadron"
  threads:30
  conda:
    "../envs/quadron.yaml"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_score_quadron.benchmark"
  shell:
    "Rscript "+generate_quadron_dat+" {params.type} {params.Quadron_lib} {params.source} {threads} {input.fas} {output.quadron_table}"

#rule quadron_features
#Same as quadron rule but used to extract features instead of extracting the score
#used by quadron_retrained
rule quadron_features:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    quadron_table = OUT+"{sample}/predict/{sample}.quadron_features"#tsv.gz
  params:
    type = "features",
    source="scripts/Snakemake/Quadron_add_Matthieu.R",
    Quadron_lib = "scripts/imports/Quadron"
  threads:30
  conda:
    "../envs/quadron.yaml"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_features_quadron.benchmark"
  shell:
    "Rscript "+generate_quadron_dat+" {params.type} {params.Quadron_lib} {params.source} {threads} {input.fas} {output.quadron_table}"

#rule quadron_retrained
#take quadron_features output as input
#generate score using the retrained quadron algorithm (xgboost)
rule quadron_retrained:
  input:
    quadron_table = OUT+"{sample}/predict/{sample}.quadron_features"
  output:
    score = OUT+"{sample}/predict/{sample}.quadron_retrained"
  params:
    model_path = "rds/Quadron_model_retrained_xgboost_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Sequence_train.model"
  threads:30
  conda:
    "../envs/quadron_retrained.yaml"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_retrained_quadron.benchmark"
  shell:
    "Rscript "+quadron_retrained_script+" {params.model_path} {input.quadron_table} {output.score}"
