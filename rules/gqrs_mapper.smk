GQRS_Mapper_script = "scripts/imports/qgrs-cpp/qgrs" # for gqrs_mapper rule
GQRS_Mappper_tsv = "scripts/Snakemake/qgrs_to_tsv.R" # for GQRS_Mappper_tsv_format rule
#rule GQRS_Mapper
#for each sequence of {input.fas} echo in shell and compute G4 score(s) using GQRS_Mapper_script
rule GQRS_Mapper:
  input:
    fas=OUT+"{sample}/fasta/merged/{sample}_merged.Fa"
  output:
    OUT+"{sample}/predict/{sample}.gqrs_mapper"#bed		
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_gqrs_mapper.benchmark"
  run:
    from Bio import SeqIO
    with open(output[0],"w") as outfile:
      with open(input.fas, "rU") as handle:
          for record in SeqIO.parse(handle, "fasta"):
              for line in shell("echo "+str(record.seq)+" | "+GQRS_Mapper_script+" -csv -notitle",iterable=True):
                if len(line) == 0:
                  outfile.write(record.id+",,,,,,,,,\n")
                else:
                  outfile.write(record.id+","+line+"\n")
      
#rule GQRS_Mappper_tsv_format
#take as input GQRS_Mapper output and generate a {params.calc} of scores by sequence
rule GQRS_Mappper_tsv_format:
  input:
    table=OUT+"{sample}/predict/{sample}.gqrs_mapper"
  output:
    OUT+"{sample}/predict/{sample}.gqrs_mapper_{calc}"
  params:
    calc="{calc}"
  conda: "../envs/gqrs_mapper.yaml"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_gqrs_mapper_{calc}.benchmark"
  shell:
    "Rscript "+GQRS_Mappper_tsv+" {input.table} {output} {params.calc}"

