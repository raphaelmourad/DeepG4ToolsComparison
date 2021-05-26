generate_fasta_atac_from_bed = "../scripts/Snakemake/generate_fasta_atac_from_bed.R"# for merge_fasta rule
generate_fasta_from_bed = "../scripts/Snakemake/generate_fasta_from_bed.R"# for merge_fasta rule
merging_ctrl_pos = "../scripts/Snakemake/merging_fasta.R"# for merge_fasta rule
subseq_fasta = "scripts/Snakemake/subset_seq.R" # for subseq_fasta rule
#rule merge fasta 
#this rule compute a Fasta file containing positive sequences (here G4seq) and another Fasta file containing control sequences (here 3 controls)
#this rule input two fasta files per run (one positive and one control)
#this rule run the script "merging_fasta.R" in the script file
#this rule output as many merged_fasta_files as control files 
# rule merge_fasta:
#   input:
#     fasta_pos = GetExp,
#     fasta_ctrl = GetCtrl
#   output:
#     fasta_merged = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
#   benchmark:
#     OUT+"{sample}/benchmarks/{sample}_merged.benchmark"
#   conda: "../envs/fasta.yaml"
#   script:
#     merging_ctrl_pos

#rule subset fasta
#this rule take the output of merge_fasta and output the subset sequences of size "sub"
#also used to generate one line fasta to be used by G4detector
rule subseq_fasta:
  input:
    fas = OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    fasta_trimmed = OUTFASTA+"{sample}/{sample}_trimmed.Fa"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_trimmed.benchmark"
  params:
    sub = 200
  conda: "../envs/fasta.yaml"
  shell:
    "Rscript "+subseq_fasta+" {params.sub} {input.fas} {output.fasta_trimmed}"



rule fastafrom_bed:
  input:
    bed_pos = GetExpBed,
    bed_ctrl = GetCtrlBed
  output:
    fasta_merged = OUTFASTA+"{sample}/{sample}_merged.Fa"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_fasta_bed_merged.benchmark"
  conda: "../envs/fasta.yaml"
  script:
    generate_fasta_from_bed


# rule atac_from_bed:
#   input:
#     bed_pos = GetExpBed,
#     bed_ctrl = GetCtrlBed,
#     atac_data = GetATAC
#   output:
#     atac_merged = OUT+"{sample}/fasta/merged/{sample}_atac_merged.rds",
#   params:random_regions = "bed/random_region_for_scaling_min_max.bed"
#   benchmark:
#     OUT+"{sample}/benchmarks/{sample}_atac_bed_merged.benchmark"
#   conda: "../envs/fasta.yaml"
#   script:
#     generate_fasta_atac_from_bed


rule atac_from_bed:
  input:
    bed_pos = GetExpBed,
    bed_ctrl = GetCtrlBed,
    atac_data = GetATAC
  output:
    atac_merged = OUTFASTA+"{sample}/{sample}_atac_merged.tsv"
  params:
    random_regions = "bed/random_region_for_scaling_min_max.bed",
	  seuil=2,
	  window=5000
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_atac_bed_merged_bgnorm.benchmark"
  conda: "../envs/fasta.yaml"
  script:
    "../scripts/Snakemake/generate_fasta_atac_from_bed_bgnorm.R"
