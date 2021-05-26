pqsfinder = "../scripts/Snakemake/pqsfinder.R"  # for pqsfinder rule

#rule pqsfinder 
#this rule compute a Fasta file containing positive sequences control sequences to create pqs score
#this rule input one fasta file per run 
#this rule run the script "pqsfinder.R" in the script file
#this rule output a tsv file containing all sequence's score from the fasta file 
rule pqsfinder:
  input:
    fas=OUTFASTA+"{sample}/{sample}_merged.Fa"
  output:
    tsv=OUT+"{sample}/predict/{sample}.pqsfinder"#tsv
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_pqs.benchmark"
  conda:
    "../envs/pqsfinder.yaml"
  script:
    pqsfinder

