

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$launchDir/data"



process INDEX_REFERENCE { 
	publishDir "$params.data_dir", mode: "copy"

	input:
		path reference_genome

	output:
		tuple path("*.fa"), path("*.fa.fai"), emit: reference_genome

	shell:
	'''
	reference_name=!{reference_genome}
	reference_name=${reference_name%.*}  # removes last file extension .gz

	gunzip -c !{reference_genome} > $reference_name 
	samtools faidx $reference_name -o $reference_name.fai
	'''
}

// data aquisition dbsnp
//https://ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.gz


process SOMATICSEQ_CALLING { 
	container "somvc-pipeline:latest"
	tag "$sample_id"
	publishDir "$params.data_dir/somaticseq", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		//tuple val(sample_id), tuple(normal_file, tumor_file)
		tuple val(sample_id), path(normal_file), path(tumor_file) 
		path reference_genome

	output:
		//path "${sample_id}.bam.bai", emit: reads_mapped_index
		path "*", emit: all


	script:
	"""
	cp -r /usr/src/neusomatic/ensemble_docker_pipelines .
	chmod -R 777 ensemble_docker_pipelines
	
	bash ensemble_docker_pipelines/prepare_callers_scripts.sh --normal-bam $normal_file --tumor-bam $tumor_file --human-reference $reference_genome --output-dir varcaller_commands --splits 10 --mutect2 --somaticsniper --vardict --varscan2 --strelka --wrapper




	"""
}


process SOMVC_LOFREQ { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/lofreq", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		//path "${sample_id}.bam.bai", emit: reads_mapped_index
		path "*.vcf", emit: lofreq_vcf

	shell:
	'''
	### https://csb5.github.io/lofreq/commands/#somatic

	lofreq viterbi -f !{reference_genome} !{normal_file} | samtools sort -@ !{num_threads} -o normal_file_viterbi.bam 
	lofreq indelqual --dindel -f !{reference_genome} -o normal_file_indelqual.bam normal_file_viterbi.bam

	lofreq viterbi -f !{reference_genome} !{tumor_file} | samtools sort -@ !{num_threads} -o tumor_file_viterbi.bam 
	lofreq indelqual --dindel -f !{reference_genome} -o tumor_file_indelqual.bam tumor_file_viterbi.bam

	lofreq somatic -n normal_file_viterbi.bam -t tumor_file_viterbi.bam -f !{reference_genome} --threads !{num_threads} -o lofreq_ [-d dbsnp.vcf.gz] -l !{bed_file} --call-indels
	
	## maybe: bedtools sort -i Homo_sapiens.GRCh38.cds.all.bed


#Pre-processing:
#lofreq indelqual --dindel -f Ref GRCh37.67.fasta -o sampleX q.bam sampleX.bam
#samtools index -b sampleX q.bam > sampleX q.bai
#Variant calling:
#lofreq call --call-indels -f Ref GRCh37.67.fasta -o sampleX.vcf sampleX q.bam -s
#-S dbSNP.polymorphisms.vcf.gz



	'''
}


process SOMVC_MUTECT2 { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/mutect2", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file) 
		path reference_genome
		val num_threads

	output:
		path "*.vcf", emit: mutect2_vcf

	shell:
	'''
	gatk Mutect2 -R !{reference_genome} -I !{normal_file} -I !{tumor_file} --native-pair-hmm-threads !{num_threads} -O mutect2_unfiltered.vcf
	gatk FilterMutectCalls -R !{reference_genome} -V mutect2_unfiltered.vcf -O mutect2_filtered.vcf

	'''
}



process SOMVC_STRELKA { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/strelka", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file) 
		path reference_genome
		val num_threads

	output:
		path "*", emit: strelka_vcf


	shell:
	'''
	configManta.py --tumorBam !{tumor_file} --normalBam !{normal_file} --referenceFasta !{reference_genome} --runDir manta_dir
	manta_dir/runWorkflow.py -m local -j !{num_threads}	

	configureStrelkaSomaticWorkflow.py --tumorBam !{tumor_file} --normalBam !{normal_file} --referenceFasta !{reference_genome} --exome --runDir strelka_dir --indelCandidates manta_dir/results/variants/candidateSmallIndels.vcf.gz
	strelka_dir/runWorkflow.py -m local -j !{num_threads}
	'''
}


process SOMVC_VARDICT { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/vardict", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file) 
		path reference_genome
		val num_threads

	output:
		path "*", emit: vardict_vcf


	shell:
	'''
	/usr/src/VarDict-1.8.2/bin/VarDict -G !{reference_genome} -k 1 -b !{tumor_file}|!{normal_file} -Q 5 -th !{num_threads}

#	 /usr/src/VarDict-1.8.2/bin/var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M vars.txt ?


	gatk VariantFiltration -R !{reference_genome} -V input.vcf.gz -O output.vcf.gz --filterName "bcbio_advised" \
			--filterExpression "((AF*DP < 6) && ((MQ < 55.0 && NM > 1.0) || (MQ < 60.0 && NM > 2.0) || (DP < 10) || (QUAL < 45)))" 

	# https://github.com/bcbio/bcbio-nextgen/blob/5cfc02b5974d19908702fa21e6d2f7a50455b44c/bcbio/variation/vardict.py#L248

# var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M” and 
# filtered with “((AF*DP < 6) && ((MQ < 55.0 && NM > 1.0) || (MQ < 60.0 && NM > 2.0) || (DP < 10) || (QUAL < 45)))” using GATK VariantFiltration module according to the setting used by Bcbio.

# check
#AF THR="0.01"
#vardict -C -G Ref GRCh37.67.fasta -f $AF THR -N sampleX -b sampleX.bam -h -c 1 -S 2
#-E 3 -g 4 target.bed > sampleX.txt


	'''
}



process TEST { 
	tag "$sample_id"
	publishDir "$params.data_dir/somaticseq", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		//tuple val(sample_id), tuple(normal_file, tumor_file)
		tuple val(sample_id), path(normal_file), path(tumor_file) 
		path reference_genome

	output:
		//path "${sample_id}.bam.bai", emit: reads_mapped_index
		path "*", emit: all


	shell:
	'''
	# cd ./neusomatic/ensemble_docker_pipelines

	echo !{sample_id}
	echo !{normal_file}

	#bash prepare_callers_scripts.sh --normal-bam !{mapping_files[0]} --tumor-bam /home/stefan/Documents/umcg/somatic_test/cancer_reads_mapped/DRR023395/DRR023395.bam --human-reference /home/stefan/Documents/umcg/dna-seq_pipeline/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-dir /home/stefan/Documents/umcg/somatic_test/output --splits 10 --mutect2 --somaticsniper --vardict --varscan2 --strelka --wrapper

	#cp -r /bin/neusomatic/ensemble_docker_pipelines .
	#ls -alFh ensemble_docker_pipelines
	#chmod 777 ensemble_docker_pipeline	
	#ls -alFh ensemble_docker_pipelines


	#bash ensembl_docker_pipelines/prepare_callers_scripts.sh --normal-bam !{normal_file}


	#NEUSOMATIC_PREPRO --normal-bam
	
	#docker cp -r somvc-pipeline:latest ./neusomatic/ensemble_docker_pipelines .
	#cd ./neusomatic/ensemble_docker_pipelines

	#echo !{sample_id}
	#echo !{normal_file}


	#neusomatic/ensemble_docker_pipelines/prepare_callers_scripts.sh --normal-bam !{normal_file}

	'''
}






