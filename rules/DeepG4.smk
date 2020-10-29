DeepG4_script = "scripts/Snakemake/DeepG4.R"

DeepG4Models = {
  "BG4G4seq":"scripts/imports/DeepG4/inst/extdata/model.hdf5",
  "G4seqpmBG4":"rds/DeepG4_retrained_G4seqpm_BG4_1_1.h5",
  "qG4G4seq":"rds/DeepG4_retrained_qG4_G4seq.h5"
}
def GetModelDeepG4(wildcards):
  return(DeepG4Models[wildcards.model])

rule DeepG4:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    out = OUT+"{sample}/predict/{sample}.DeepG4_{model}"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_DeepG4_{model}.benchmark"
  params:
    model = GetModelDeepG4
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
    out = OUT+"{sample}/predict/{sample}.DeepG4Scan_{model}"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_DeepG4Scan_{model}.benchmark"
  params:
    model = GetModelDeepG4
  threads:30
  conda:
    "../envs/DeepG4.yaml"
  shell:
    "Rscript "+DeepG4_scan_script+" {input.fas} {output} {params.model} {threads}"
