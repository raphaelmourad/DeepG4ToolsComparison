

ModelDeepG4OH4 = {
"ATAC":"rds/runs/DeepG4_ATAC_rescale_BW_sampling_03_03_2021/2021-03-03T15-49-25Z/best_model.h5",
"classic":"rds/runs/DeepG4_classic_rescale_BW_sampling_03_03_2021/2021-03-03T15-52-06Z/best_model.h5",
"ATACtuningOH5":"rds/runs/DeepG4_ATAC_rescale_BW_sampling_02_03_2021/2021-03-02T16-01-34Z/best_model.h5",
"classictuningOH5":"rds/runs/DeepG4_classic_rescale_BW_sampling_02_03_2021/2021-03-02T16-17-28Z/best_model.h5",
"classic_old":"rds/runs/DeepG4_classic_rescale_BW_sampling_old_hypparams/2021-03-03T11-24-25Z/best_model.h5"


}

def GetModelDeepG4OH4(wildcards):
  return(ModelDeepG4OH4[wildcards.model])

DeepG4OH4_script="scripts/Snakemake/DeepG4OH4_ATAC.R"

rule DeepG4OH4:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa",
    atac_merged = OUT+"{sample}/fasta/merged/{sample}_atac_merged.rds"
  output:
    out = OUT+"{sample}/predict/{sample}.ATACDeepG4_{model}"#TSV		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_ATACDeepG4__{model}.benchmark"
  params:
    model = GetModelDeepG4OH4
  threads:30
  conda:
    "../envs/DeepG4.yaml"
  shell:
    "Rscript "+DeepG4OH4_script+" {input.fas} {input.atac_merged} {output} {params.model}"

