penguinn_script = "scripts/imports/penguinn/penguinn.py"
# rule penguinn_old:
#   input:
#     fas = OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa"
#   output:
#     out = OUT+"{sample}/predict/{sample}.penguinn"#TSV
#   benchmark:
#     OUT+"{sample}/benchmarks/{sample}_penguinn.benchmark"
#   threads:30
#   params:
#     model = "scripts/imports/penguinn/Models/model_1_1.h5"
#   conda:
#     "../envs/penguinn.yaml"
#   shell:
#     "python "+penguinn_script+" --input {input.fas} --output {output} --model {params.model}"
# 
# 
# rule penguinn_retrained_old:
#   input:
#     fas = OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa"
#   output:
#     out = OUT+"{sample}/predict/{sample}.penguinn_retrained"#TSV
#   benchmark:
#     OUT+"{sample}/benchmarks/{sample}_penguinn_retrained.benchmark"
#   threads:30
#   params:
#     model = "rds/penguinn_retrained_best_model_1_1.h5"
#   conda:
#     "../envs/penguinn.yaml"
#   shell:
#     "python "+penguinn_script+" --input {input.fas} --output {output} --model {params.model}"

penguinn_script_R = "scripts/Snakemake/penguinn_custom_implementation.R"
rule penguinn:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.penguinn"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_penguinn.benchmark"
  threads:30
  params:
    model = "scripts/imports/penguinn/Models/model_1_1.h5"
  conda:
    "../envs/penguinn.yaml"
  shell:
    "Rscript "+penguinn_script_R+" {input.fas} {output} {params.model}"

rule penguinn_retrained:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_trimmed.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.penguinn_retrained"#TSV
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_penguinn_retrained.benchmark"
  threads:30
  params:
    model = "rds/penguinn_retrained_best_model_1_1.h5"
  conda:
    "../envs/penguinn.yaml"
  shell:
    "Rscript "+penguinn_script_R+" {input.fas} {output} {params.model}"
