
ModelDeepG4OH4 = {
#"ATACnormBG":"rds/runs/DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021/2021-04-19T16-32-14Z/best_model.h5",
"ATACnormBG":"rds/runs/DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021/2021-07-02T08-48-55Z/best_model.h5",
"ATACnorm900BG":"rds/runs/DeepG4_ATAC_rescale_BW_by_bg_5kb_seuil_2_19_04_2021/2021-07-06T07-59-11Z/best_model.h5",
"classictuningOH5":"rds/runs/DeepG4_classic_rescale_BW_sampling_02_03_2021/2021-03-02T16-17-28Z/best_model.h5"
}

def GetModelDeepG4OH4(wildcards):
  return(ModelDeepG4OH4[wildcards.model])

DeepG4OH4_script="scripts/Snakemake/DeepG4OH4_ATAC.R"

rule DeepG4OH4:
  input:
    fas = OUTFASTA+"{sample}/{sample}_merged.Fa",
    atac_merged = OUTFASTA+"{sample}/{sample}_atac_merged.tsv"
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
