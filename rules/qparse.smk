
qparse_script = "scripts/imports/qparse/QPARSE_1.10.py"
qparse_tsv = "scripts/Snakemake/qparse_to_tsv.R"

rule qparse:
  input:
    fas=OUTFASTA+"{sample}/{sample}_merged.Fa"
  output:
    OUT+"{sample}/predict/{sample}.qparse"
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_qparse.benchmark"
  conda:"../envs/G4CatchAll.yaml"
  shell:
    "python2 "+qparse_script+" -i {input.fas} -o {output}"



rule qparse_tsv_python:
  input:
    qparse_file=OUT+"{sample}/predict/{sample}.qparse"
  output:
    temp(OUT+"{sample}/predict/{sample}.qparse_tsv_inter")
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_qparse_tsv_python.benchmark"
  run:
    fasta = {}
    with open(output[0],"w") as outfile:
      with open(input.qparse_file, "r") as qparsefile:
          for line in qparsefile:
            line = line.strip()
            if not line:
              continue
            if line.startswith("#"):
              continue
            if line.startswith(">"):
              active_sequence_name = line[1:]
              if active_sequence_name not in fasta:
                  fasta[active_sequence_name] = []
              continue
            sequence = line
            fasta[active_sequence_name].append(sequence)
      for i,v in fasta.items():
         outfile.write("\n".join([i+"\t"+j for j in v])+"\n")

rule qparse_tsv_R:
  input:
    fas=OUTFASTA+"{sample}/{sample}_merged.Fa",
    qparse_file=OUT+"{sample}/predict/{sample}.qparse_tsv_inter"
  output:
    OUT+"{sample}/predict/{sample}.qparse_{calc}"	
  benchmark:
    OUT+"{sample}/benchmarks/{sample}_qparse_{calc}_R.benchmark"
  params:
    calc="{calc}"
  conda: "../envs/quadron.yaml"
  shell:
    "Rscript "+qparse_tsv+" {input.fas} {input.qparse_file} {output} {params.calc}"
