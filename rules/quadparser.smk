quadparser_tsv = "scripts/Snakemake/quadparser_to_tsv.R" # for quadparset_to_tsv rule
#rule quadparser_command 
#this rule compute merged fasta files to compute quadparser
#this rule input one merge fasta file per run 
#this rule run the shell line "g4predict intra -M -f {input} -b {output}" in the tools file
#this rule output a bed file containing the output of quadparser 
rule quadparser_command:
  input:
    fas = OUTFASTA+"{sample}/{sample}_merged.Fa"
  output:
    quad_parser_table = OUT+"{sample}/predict/{sample}.quadparser"#bed		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_quadparser.benchmark"
  shell:
    "g4predict intra -M -f {input} -b {output}"

#rule quadparser_tsv_format 
#this rule take as input the output of quadparser_command rule and export it at tsv format, with 'params.calc' score for each sequence
rule quadparser_tsv_format:
  input:
    fas = OUTFASTA+"{sample}/{sample}_merged.Fa",
    table=OUT+"{sample}/predict/{sample}.quadparser"
  output:
    OUT+"{sample}/predict/{sample}.quadparser_{calc}"
  params:
    calc="{calc}"
  conda: "../envs/fasta.yaml"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_quadparser_{calc}.benchmark"
  shell:
    "Rscript "+quadparser_tsv+" {input.table} {input.fas} {output} {params.calc}"

