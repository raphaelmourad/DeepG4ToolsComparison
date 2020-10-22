DeepG4_script = "scripts/Snakemake/DeepG4.R"
rule DeepG4:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.DeepG4"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_DeepG4.benchmark"
  params:
    model = "scripts/imports/DeepG4/inst/extdata/model.hdf5"
  threads:30
  conda:
    "../envs/DeepG4.yaml"
  shell:
    "Rscript "+DeepG4_script+" {input.fas} {output} {params.model}"

DeepG4_scan_script = "scripts/Snakemake/DeepG4_scan.R"
rule DeepG4_scan:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.DeepG4_scan"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_DeepG4_scan.benchmark"
  params:
    model = "scripts/imports/DeepG4/inst/extdata/model.hdf5"
  threads:30
  conda:
    "../envs/DeepG4.yaml"
  shell:
    "Rscript "+DeepG4_scan_script+" {input.fas} {output} {params.model} {threads}"
