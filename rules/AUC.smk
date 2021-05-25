AUC_script = "scripts/Snakemake/AUC_script.R"
rule AUC:
  input:
    OUT+"{sample}/predict/{sample}.{tool}"
  output:
    pdf=OUT+"{sample}/AUC/{sample}.{tool}.tsv"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_{tool}_AUC.benchmark"
  conda : "../envs/AUC.yaml"
  shell:
    "Rscript "+AUC_script+" {input} {output.pdf}"

supp_metrics_script = "scripts/Snakemake/supp_metrics.R"
rule other_metrics:
  input:
    OUT+"{sample}/predict/{sample}.{tool}"
  output:
    pdf=OUT+"{sample}/other_metrics/{sample}.{tool}.tsv"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_{tool}_others_metrics.benchmark"
  conda : "../envs/AUPR.yaml"
  shell:
    "Rscript "+supp_metrics_script+" {input} {output.pdf}"
