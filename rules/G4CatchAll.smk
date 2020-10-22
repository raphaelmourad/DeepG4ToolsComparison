G4CatchAll_script = "scripts/imports/G4Catchall/G4Catchall.py"
G4CatchAll_tsv = "scripts/Snakemake/G4CatchAll_to_tsv.R"
rule G4CatchAll:
  input:
    fas=OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    OUT+"{sample}/predict/{sample}.G4CatchAll"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4CatchAll.benchmark"
  conda: "../envs/G4CatchAll.yaml"
  shell:
    "python2 "+G4CatchAll_script+" --fasta {input.fas} > {output}"


rule G4CatchAll_tsv_format:
  input:
    fas=OUT+"{sample}/fasta/merged/{sample}_merged.Fa",
    table=OUT+"{sample}/predict/{sample}.G4CatchAll"
  output:
    OUT+"{sample}/predict/{sample}.G4CatchAll_{calc}"
  params:
    calc="{calc}"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_G4CatchAll_{calc}.benchmark"
  conda: "../envs/fasta.yaml"
  shell:
    "Rscript "+G4CatchAll_tsv+" {input.table} {input.fas} {output} {params.calc}"
