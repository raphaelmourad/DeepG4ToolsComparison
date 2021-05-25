G4detector_script = "scripts/imports/G4detector/code/G4detector.py"
G4detector_to_tsv = "scripts/Snakemake/G4detector_to_tsv.R"
rule G4detector:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.G4detector"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4detector.benchmark"
  threads:30
  params:
    model = "scripts/imports/G4detector/models/random/model_rand_K_PDS.h5"
  conda:
    "../envs/G4detector.yaml"
  shell:
    "python "+G4detector_script+" test {input.fas} {params.model} {output}"


rule rule G4detector_tsv_format:
  input:
    fas=OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa",
    table=OUT+"{sample}/predict/{sample}.G4detector"
  output:
    OUT+"{sample}/predict/{sample}.G4detector_tsv"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4detector_tsv.benchmark"
  conda:
    "../envs/fasta.yaml"
  shell:
    "Rscript "+G4detector_to_tsv+" {input.table} {input.fas} {output}"

rule G4detector_retrained:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.G4detector_retrained"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4detector_retrained.benchmark"
  threads:30
  params:
    model = "rds/G4detector_retrained_model_rand_K_PDS.h5"
  conda:
    "../envs/G4detector.yaml"
  shell:
    "python "+G4detector_script+" test {input.fas} {params.model} {output}"



rule rule G4detector_retrained_tsv_format:
  input:
    fas=OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa",
    table=OUT+"{sample}/predict/{sample}.G4detector_retrained"
  output:
    OUT+"{sample}/predict/{sample}.G4detector_retrained_tsv"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4detector_retrained_tsv.benchmark"
  conda:
    "../envs/fasta.yaml"
  shell:
    "Rscript "+G4detector_to_tsv+" {input.table} {input.fas} {output}"
