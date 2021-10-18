FROM conda/miniconda2@sha256:d30ea858b9d3425127c2e835d39653cdca74368af832d65748cf912f6c1da67b

#FROM conda/miniconda3@sha256:2141bfc4b60cc5eb332402207c80d884daa72fcefa48b50f3ccadd934d1f3d03

MAINTAINER loipf

#RUN apt-get update 
#RUN conda update --all

RUN apt-get update && apt install -y wget

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


### bioinformatic tools
### strelka needs python 2.7 - conflict with multiqc (replace with older version)
RUN conda install \
		gatk4==4.2.2.0 \
		lofreq==2.1.5 \
		vardict==2019.06.04 \
		strelka==2.9.10 \
		samtools==1.11.0 

# RUN pip install multiqc

#RUN conda create --name strelka_py27 #python=2.7
#RUN conda run -n strelka_py27 --no-capture-output install -c bioconda strelka==2.9.10
### multiqc==1.11

# SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]



### download SomaticCombiner
RUN wget https://github.com/mingyi-wang/somatic-combiner/blob/V1.02/example/somaticCombiner.jar
RUN mv somaticCombiner.jar /usr/src/somaticCombiner.jar


### download dbSNP - UNZIP MAYBE ??
# RUN wget https://ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.gz
# RUN wget https://ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.gz.tbi


# /bin/sh -c cd /opt && wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 && tar -xvf strelka-2.9.10.centos6_x86_64.tar.bz2 && rm strelka-2.9.10.centos6_x86_64.tar.bz2 && ln -s strelka-2.9.10.centos6_x86_64 strelka


# install https://github.com/nygenome/Conpair
# add java



