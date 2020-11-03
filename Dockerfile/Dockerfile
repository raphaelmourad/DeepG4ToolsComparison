FROM rocker/tidyverse:4.0.3
RUN R -e "BiocManager::install(c('biomaRt', 'Biostrings'))"
COPY scripts/imports/G4iMGrinder /G4iMGrinder
RUN R -e "devtools::install_local('/G4iMGrinder/',upgrade ='never')"
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && /bin/bash miniconda.sh -b
ENV PATH /opt/conda/bin:$PATH
RUN conda install -c conda-forge biopython --yes
RUN conda install -c conda-forge mamba --yes && mamba create --yes -c conda-forge -c bioconda -n snakemake snakemake
RUN echo "conda activate snakemake" >> ~/.bashrc

COPY scripts/imports/g4predict /g4predict

RUN cd /g4predict && python setup.py build && python setup.py install
RUN cd /root/
RUN apt-get update && apt-get install --yes libbz2-dev liblzma-dev libgsl-dev
RUN R -e "BiocManager::install(c('PWMEnrich','JASPAR2020','BSgenome','rtracklayer','Rsamtools','Rsamtools','Rhtslib','TFBSTools','motifmatchr'))"
RUN R -e "install.packages(c('ranger','pROC'))"
