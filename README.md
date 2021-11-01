# SOMvc pipeline

a variant calling pipeline from matched normal-tumor `.bam` files to `.vcf` files using NeuSomatic ensembl approach.


---
### set up pipeline


before running, you have to set up the attached Docker images, DNAvc-pipeline docker image can take some time:
```sh
docker build -t somvc-pipeline https://raw.githubusercontent.com/loipf/SOMvc-pipeline/master/docker/Dockerfile
docker build -t dnavc-pipeline https://raw.githubusercontent.com/loipf/DNAvc-pipeline/master/docker/Dockerfile

```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker somvc-pipeline` argument.



if you already have an reference genome, you can specify it with the nextflow argument `--reference_genome [genome.fa.gz]` or in the `main.nf` file. otherwise run the following script to download it into `data` folder:
```sh
bash scripts/data_acquisition.sh
```



---
### run mapping pipeline

it can be run locally with downloaded github-repo and edited `nextflow.config` file with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/DNAseq-pipeline -r main --project_dir /path/to/folder --reads_dir /path/to/samples --num_threads 10 --adapter_3_seq_file adapter_3.fasta --adapter_5_seq_file adapter_5.fasta --reference_genome genome.fasta --bed_file cds.bed -with-docker dnaseq-pipeline
```
for this execution to work properly, you have to be in the current project directory.


optional extendable with:
```sh
-resume
-with-report report_SOMvc-pipeline
-with-timeline timeline_SOMvc-pipeline
-w work_dir
```

by default, all output will be saved into the `data` folder of the current directory

input `.bam` files must be sorted with marked duplicates and index created in the same directory.

the sample matching `.csv` file must be constructed like following (with same column names, filenames per sample must be different!):
```sh
sample_id,normal_file,tumor_file
sample_1,/path/to/sample_1.bam,/path/to/sample_1_tum.bam
sample_2,/path/to/sample_2.bam,/path/to/sample_2_tum.bam
sample_3,/path/to/sample_3.bam,/path/to/sample_3_tum.bam
```








