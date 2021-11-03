# SOMvc pipeline

a variant calling pipeline from matched normal-tumor `.bam` files to `.vcf` files using somaticcombiner consensus approach utilizing lofreq, mutect2, strelka and vardict.


---
### set up pipeline


before running, you have to set up the attached Docker images, DNAvc-pipeline docker image (for variant interpretation) can take some time:
```sh
docker build -t somvc-pipeline https://raw.githubusercontent.com/loipf/SOMvc-pipeline/master/docker/Dockerfile
docker build -t dnavc-pipeline https://raw.githubusercontent.com/loipf/DNAvc-pipeline/master/docker/Dockerfile

```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker somvc-pipeline` argument.


no mapping-refinement is needed since VC-callers do so internally.

if you already have an reference genome and bed file, you can specify it with the nextflow argument `--reference_genome [genome.fa.gz]` and `--bed_file [bed_file_filtered.bed]` (only chromosomes which are available in the genome). otherwise run the following script to download it into `data` folder:
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
nextflow run loipf/SOMvc-pipeline --project_dir /path/to/folder --sample_match_file normal_tumor_pairs.csv --reference_genome genome.fa.gz -resume --bed_file bed_file_filtered.bed --num_threads 10 -with-docker somvc-pipeline

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

input `.bam` files must be sorted with marked duplicates (and index created in the same directory).

the sample matching `normal_tumor_pairs.csv` file must be constructed like following (with same column names, filenames per sample must be different!):
```sh
sample_id,normal_file,tumor_file
sample_1,/path/to/sample_1.bam,/path/to/sample_1_tum.bam
sample_2,/path/to/sample_2.bam,/path/to/sample_2_tum.bam
sample_3,/path/to/sample_3.bam,/path/to/sample_3_tum.bam
```








