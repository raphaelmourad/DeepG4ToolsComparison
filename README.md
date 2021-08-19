
<!-- README.md is generated from README.Rmd. Please edit that file -->

## DeepG4ToolsComparison: A snakemake pipeline to run and compare G4 DNA prediction tools with DeepG4

The predictions for differents tissues and cancer with DeepG4 is
available [here](prediction_genome/bed/predictions)

The code to generate the precision/recall curve is available
[here](prediction_genome).

### Overview

It’s based on [Snakemake](https://snakemake.readthedocs.io/en/stable/)
to manage the workflow and [Docker](https://www.docker.com/) to isolate
the application and run it with the appropriate tool versions.

### Installation

#### Clone the repository :

``` bash
git clone https://github.com/morphos30/DeepG4ToolsComparison.git
cd DeepG4ToolsComparison
```

#### Install the docker image and run it :

``` bash
docker build . -t morphos30/g4docker -f Dockerfile/Dockerfile
docker run -it -v $(pwd):/DeepG4ToolsComparison morphos30/g4docker /bin/bash
```

Where `$(pwd)` is the working directory of `DeepG4ToolsComparison` on
your computer.

#### Launch the pipeline :

``` bash
cd /DeepG4ToolsComparison
snakemake --use-conda -j 30
```

You have to set the option `--use-conda` in order to install and run
each tool in its proper environment.

### Workflow specifications

#### Input

-   DNA sequences into bed format, split into positive set and negative
    set, written into the bed directory.

**Note :** if you want add a new dataset, edit the `Snakefile` file and
add the bed files in the dictionnary `EXPERIMENTS`, without the `.bed`
extension. Example :

`TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Ctrl_gkmSVM.bed`
`TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42.bed`

``` python
EXPERIMENTS = {
  "TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Ctrl_gkmSVM":{"CTRL":"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42_Ctrl_gkmSVM","EXP":"TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_0.8_42"}
}
```

Where `CTRL` is the negative set and `EXP` is the positive set.

-   DNA Accessibility (ATAC-seq/DNAse-seq/MNase-seq) in bigwig format or
    directly the averaged value for each sequence in a `one-column tsv`
    file.

``` python
ATACFILE = {
    "TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM":["ATAC_entinostat_mean.bw"]
}
```

or `one-column tsv` file in
`fasta/{Experiment_name}/{Experiment_name}_atac_merged.tsv`. Example :

`fasta/TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM/TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_atac_merged.tsv`

    head TestSet_Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM_atac_merged.tsv 
    0.01628741641898675
    0.028752257447422012
    0.028878783223623482
    0.055516399884055316
    0.02825982069785745
    0.03582923041809851
    0.023904436394151577
    0.07724288611280039
    0.01740800116454673
    0.05779605688479145

#### Rulegraph :

<img src="rulegraph.svg" width="150%" />

#### Workflow output for each tools :

| Outputs                      | Tools      | Methods                                                       |
|:-----------------------------|:-----------|:--------------------------------------------------------------|
| ATACDeepG4\_ATACnormBG       | ATACDeepG4 | DeepG4 using accessibily (DeepG4 in paper)                    |
| ATACDeepG4\_classictuningOH5 | ATACDeepG4 | DeepG4 without accessibility (DeepG4\* in paper)              |
| penguinn\_retrained          | penguinn   | penguinn using custom model trained on **BG4G4seq** dataset   |
| penguinn                     | penguinn   | penguinn using default model                                  |
| G4detector\_retrained        | G4detector | G4detector using custom model trained on **BG4G4seq** dataset |
| G4detector                   | G4detector | G4detector using default model                                |
| quadron\_retrained           | quadron    | quadron using custom model trained on **BG4G4seq** dataset    |
| quadron\_score               | quadron    | quadron using default model                                   |
