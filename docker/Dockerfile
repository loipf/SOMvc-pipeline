FROM conda/miniconda2@sha256:d30ea858b9d3425127c2e835d39653cdca74368af832d65748cf912f6c1da67b

MAINTAINER loipf

RUN apt-get update && apt install -y wget unzip

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


### bioinformatic tools
### strelka needs python 2.7 - conflict with multiqc (replace with older version)
RUN conda install \
		gatk4==4.2.2.0 \
		lofreq==2.1.5 \
		strelka==2.9.10 \
		manta==1.6.0 \
		samtools==1.12.0 \
		vardict==2019.06.04
RUN conda install -y numpy==1.16.5 scipy==1.2.1
RUN conda install r-base==3.2.2 r-essentials


### download VarDictJava - faster than perl conda version - DOES NOT WORK BUT NEED /bin scripts
RUN wget https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.2/VarDict-1.8.2.tar
RUN tar -xf VarDict-1.8.2.tar -C /usr/src/
RUN rm VarDict-1.8.2.tar


### download SomaticCombiner
RUN wget https://github.com/mingyi-wang/somatic-combiner/blob/V1.02/example/somaticCombiner.jar
RUN mv somaticCombiner.jar /usr/src/somaticCombiner.jar


### install Conpair (master from Apr 3 2020) https://github.com/nygenome/Conpair
RUN wget https://github.com/nygenome/Conpair/archive/144be0aecd97c289b095230c09c11e12665a1dae.zip -O conpair.zip
RUN unzip conpair.zip
RUN mv Conpair-* /usr/src/conpair







