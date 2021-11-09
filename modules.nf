
//TODO
// dont public index_reference

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$launchDir/data"



process INDEX_REFERENCE { 
	publishDir "$params.data_dir", mode: "copy"

	input:
		path reference_genome
		path bed_file

	output:
		tuple path("*.fa"), path("*.fa.fai"), path("*.dict"), emit: reference_genome
		tuple path("*.bed_sorted.gz"), path("*.bed_sorted.gz.tbi"), emit: bed_file

	shell:
	'''
	reference_name=!{reference_genome}
	reference_name=${reference_name%.*}  # removes last file extension .gz

	gunzip -c !{reference_genome} > $reference_name 
	samtools faidx $reference_name -o $reference_name.fai
	samtools dict $reference_name -o ${reference_name%.*}.dict  # remove .fa so name is only .dict

	bedtools sort -i !{bed_file} | bgzip -c > !{bed_file}_sorted.gz
	tabix --zero-based -b 2 -e 3 !{bed_file}_sorted.gz
	'''
}



process SOMVC_LOFREQ { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/lofreq", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file), path(normal_file_index), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		path "*"
		tuple path ("lofreq_somatic_final.snvs.vcf.gz"), path("lofreq_somatic_final.snvs.vcf.gz.tbi"), emit: lofreq_snvs_vcf
		tuple path ("lofreq_somatic_final.indels_vt.vcf.gz"), path("lofreq_somatic_final.indels_vt.vcf.gz.tbi"), emit: lofreq_indel_vcf
		val $sample_id, emit: lofreq_sample_id

	shell:
	'''
	### https://csb5.github.io/lofreq/commands/#somatic

	gunzip -c !{bed_file[0]} > bed_file_unzipped.bed   # vardict need unzipped
	

	lofreq viterbi -f !{reference_genome[0]} !{normal_file} | samtools sort -@ !{num_threads} -o normal_file_viterbi.bam -
	lofreq indelqual --dindel -f !{reference_genome[0]} -o normal_file_indelqual.bam normal_file_viterbi.bam
	samtools index -b -@ !{num_threads} normal_file_viterbi.bam

	lofreq viterbi -f !{reference_genome[0]} !{tumor_file} | samtools sort -@ !{num_threads} -o tumor_file_viterbi.bam -
	lofreq indelqual --dindel -f !{reference_genome[0]} -o tumor_file_indelqual.bam tumor_file_viterbi.bam
	samtools index -b -@ !{num_threads} tumor_file_viterbi.bam

	lofreq somatic -n normal_file_viterbi.bam -t tumor_file_viterbi.bam -f !{reference_genome[0]} --threads !{num_threads} -o lofreq_ -l bed_file_unzipped.bed --call-indels
	
	### vt normalization for indels
	vt normalize lofreq_somatic_final.indels.vcf.gz -r !{reference_genome[0]} -o lofreq_somatic_final.indels_vt.vcf.gz
	tabix -p vcf lofreq_somatic_final.indels_vt.vcf.gz

	'''
}


process SOMVC_MUTECT2 { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/mutect2", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file), path(normal_file_index), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		path "*"
		tuple path ("mutect2_filtered_vt.vcf.gz"), path("mutect2_filtered_vt.vcf.gz.tbi"), emit: mutect2_vcf

	shell:
	'''
	gunzip -c !{bed_file[0]} > bed_file_unzipped.bed   # mutect2 need unzipped
	
	gatk Mutect2 -R !{reference_genome[0]} -I !{normal_file} -I !{tumor_file} --native-pair-hmm-threads !{num_threads} --intervals bed_file_unzipped.bed -O mutect2_unfiltered.vcf
	gatk FilterMutectCalls -R !{reference_genome[0]} -V mutect2_unfiltered.vcf -O mutect2_filtered.vcf

	bgzip -c mutect2_filtered.vcf > mutect2_filtered.vcf.gz
	vt normalize mutect2_filtered.vcf.gz -r !{reference_genome[0]} -o mutect2_filtered_vt.vcf.gz
	tabix -p vcf mutect2_filtered_vt.vcf.gz
	'''
}



process SOMVC_STRELKA { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/strelka", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file), path(normal_file_index), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		path "*"
		tuple path ("strelka_dir/results/variants/somatic.snvs.vcf.gz"), path("strelka_dir/results/variants/somatic.snvs.vcf.gz.tbi"), emit: strelka_snv_vcf
		tuple path ("somatic.indels_vt.vcf.gz"), path("somatic.indels_vt.vcf.gz.tbi"), emit: strelka_indel_vcf

	shell:
	'''
	configManta.py --tumorBam !{tumor_file} --normalBam !{normal_file} --referenceFasta !{reference_genome[0]} --runDir manta_dir
	manta_dir/runWorkflow.py -m local -j !{num_threads}	

	configureStrelkaSomaticWorkflow.py --tumorBam !{tumor_file} --normalBam !{normal_file} --referenceFasta !{reference_genome[0]} --exome --runDir strelka_dir --indelCandidates manta_dir/results/variants/candidateSmallIndels.vcf.gz --callRegions !{bed_file[0]}
	strelka_dir/runWorkflow.py -m local -j !{num_threads}

	### vt normalization for indels
	vt normalize strelka_dir/results/variants/somatic.indels.vcf.gz -r !{reference_genome[0]} -o somatic.indels_vt.vcf.gz
	tabix -p vcf somatic.indels_vt.vcf.gz
	'''
}


process SOMVC_VARDICT { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/vardict", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file), path(normal_file_index), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		path "*"
		path "*", emit: vardict_snv_vcf


	shell: 
	'''
	gunzip -c !{bed_file[0]} > bed_file_unzipped.bed   # vardict need unzipped

	/usr/src/VarDict-1.8.2/bin/VarDict -G !{reference_genome[0]} -k 1 -b "!{tumor_file}|!{normal_file}" -Q 5 -z 1 -c 1 -S 2 -E 3 -g 4 -th !{num_threads} bed_file_unzipped.bed | /usr/src/VarDict-1.8.2/bin/testsomatic.R | /usr/src/VarDict-1.8.2/bin/var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M -N "tumor_sample|normal_sample" > vardict_output.vcf

	# https://github.com/bcbio/bcbio-nextgen/blob/5cfc02b5974d19908702fa21e6d2f7a50455b44c/bcbio/variation/vardict.py#L248
	gatk VariantFiltration -R !{reference_genome[0]} -V vardict_output.vcf -O vardict_output_filtered.vcf --filterName "bcbio_advised" \
			--filterExpression "((AF*DP < 6) && ((MQ < 55.0 && NM > 1.0) || (MQ < 60.0 && NM > 2.0) || (DP < 10) || (QUAL < 45)))" 

	
### TODO DELETE
### weird C writing error introducing empty bytes - deleted with tr - ONLY WITH JAVA version
#	/usr/src/VarDict-1.8.2/bin/VarDict -G !{reference_genome[0]} -k 1 -b "!{tumor_file}|!{normal_file}" -Q 5 -th !{num_threads} | tr -d '\000' | /usr/src/VarDict-1.8.2/bin/testsomatic.R | /usr/src/VarDict-1.8.2/bin/var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M -N "tumor_sample|normal_sample" > vardict_output.vcf

# /usr/src/VarDict-1.8.2/bin/VarDict -G !{reference_genome[0]} -k 1 -b "!{tumor_file}|!{normal_file}" -Q 5 -th !{num_threads} -z 1 -c 1 -S 2 -E 3 -g 4 !{bed_file[0]} | tr -d '\000' | /usr/src/VarDict-1.8.2/bin/testsomatic.R | /usr/src/VarDict-1.8.2/bin/var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M -N "tumor_sample|normal_sample" > vardict_output.vcf

#VarDict-1.8.2/bin/VarDict -G Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -b "/home/stefanloipfinger/Documents/somvc_pipeline_test/baseline_reads_mapped/IMPACT_MBC1_001_BKDN190631707-1A_HTKWGDSXX_L2/IMPACT_MBC1_001_BKDN190631707-1A_HTKWGDSXX_L2.bam|/home/stefanloipfinger/Documents/somvc_pipeline_test/tumor_reads_mapped/IMPACT_MBC1_002_BKDN190631708-1A_HTKWGDSXX_L4/IMPACT_MBC1_002_BKDN190631708-1A_HTKWGDSXX_L4.bam" -Q 5 -th 1 -R Homo_sapiens_short.GRCh38.cds.all.bed -z 1 -c 1 -S 2 -E 3 -g 4| tr -d '\000' | VarDict-1.8.2/bin/testsomatic.R | VarDict-1.8.2/bin/var2vcf_paired.pl -P 0.9 -f 0.8 -M -N "tumur_sample|normal_sample" > vardict_output.vcf

# -R chr7:50000-200000

#VarDict-1.8.2/bin/VarDict -G Homo_sapiens.GRCh38.dna.primary_assembly.fa -k 1 -b "/home/stefanloipfinger/Documents/somvc_pipeline_test/tumor_reads_mapped/IMPACT_MBC1_002_BKDN190631708-1A_HTKWGDSXX_L4/IMPACT_MBC1_002_chr17.bam|/home/stefanloipfinger/Documents/somvc_pipeline_test/baseline_reads_mapped/IMPACT_MBC1_001_BKDN190631707-1A_HTKWGDSXX_L2/IMPACT_MBC1_001_chr17.bam" -Q 5 -z 1 -c 1 -S 2 -E 3 -g 4 -th 10 -R chr17:7500005-7599990 > vardict_output_all.vcf


# check
#AF THR="0.01"
#vardict -C -G Ref GRCh37.67.fasta -f $AF THR -N sampleX -b sampleX.bam -h -c 1 -S 2
#-E 3 -g 4 target.bed > sampleX.txt


	'''
}


process SOMATIC_COMBINER { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/somatic_combiner", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		val sample_id
		path lofreq_indel_vcf
		path lofreq_snv_vcf
		path mutect2_vcf
		path strelka_indel_vcf
		path strelka_snv_vcf
		path vardict_vcf

	output:
		//path "somatic_combiner_vcf_renamed.vcf", emit: somatic_combiner_vcf
		tuple val($sample_id), path ("somatic_combiner_vcf_renamed.vcf"), emit: somatic_combiner_vcf


	shell:
	'''
	java -jar /usr/src/somaticCombiner.jar -L !{lofreq_indel_vcf} -l !{lofreq_snv_vcf} -M !{mutect2_vcf} -s ${strelka_snv} -S !{strelka_indel_vcf} -D !{vardict_vcf} -o somatic_combiner_vcf.vcf

	printf '%s\n' !{sample_id}_tumor !{sample_id}_normal > sample_names.txt 
	bcftools reheader --samples sample_names.txt -o somatic_combiner_vcf_renamed.vcf
	
	'''
}




process CONPAIR_CONTAMINATION { 
	tag "$sample_id"
	publishDir "$params.data_dir/conpair", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(normal_file), path(tumor_file), path(normal_file_index), path(tumor_file_index) 
		path reference_genome

	output:
		path "*", emit: conpair_info

	shell:
	'''

	MARKER_FILE_BED="/usr/src/conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
	sed -e 's/chr//g' $MARKER_FILE_BED > markers_GRCh38_snv_formatted.bed  ### rename chromosomes

	MARKER_FILE_TXT="/usr/src/conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
	sed -e 's/chr//g' $MARKER_FILE_TXT > markers_GRCh38_snv_formatted.txt  ### rename chromosomes


	run_gatk_pileup_for_sample.py -B !{tumor_file} -O tumor_pileup --reference !{reference_genome[0]} --conpair_dir /usr/src/conpair/ --markers markers_GRCh38_snv_formatted.bed
	run_gatk_pileup_for_sample.py -B !{normal_file} -O normal_pileup --reference !{reference_genome[0]} --conpair_dir /usr/src/conpair/ --markers markers_GRCh38_snv_formatted.bed

	verify_concordance.py -T tumor_pileup -N normal_pileup --markers markers_GRCh38_snv_formatted.txt --outfile concordance_stats.txt

	estimate_tumor_normal_contamination.py -T tumor_pileup -N normal_pileup --markers markers_GRCh38_snv_formatted.txt --outfile contamination_stats.txt

	'''
}



process VARIANT_CALLING_STATS { 
	container "dnavc-pipeline:latest"
	tag "$sample_id"
	publishDir "$params.data_dir/variants_vcf", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(vcf_file)
		val num_threads

	output:
		path "${sample_id}_vcfstats.txt", emit: vcf_stats

	shell:
	'''
	bcftools stats -f ADJ_PASS --threads !{num_threads} !{vcf_file} > !{sample_id}_vcfstats.txt
	### bcftools merge !!!
	'''
}


process MERGE_VCF { 
	container "dnavc-pipeline:latest"
	tag "$sample_id"
	publishDir "$params.data_dir/variants_vcf/_all", mode: "copy", overwrite: false

	input:
		path(vcf_files)
		val num_threads

	output:
		path "all_vcf_merged.vcf", emit: vcf_all

	shell:
	'''
	bcftools merge -o all_vcf_merged.vcf --threads !{num_threads} !{vcf_files}
	'''
}




process VARIANT_ANNOTATION { 
	container "dnavc-pipeline:latest"
	publishDir "$params.data_dir/variants_vcf/_all", mode: "copy", overwrite: false

	input:
		path vcf_all
		val num_threads

	output:
		path "*"

	shell:
	'''
	oc run -l hg38 -t csv -x --mp !{num_threads} !{vcf_all}
	'''
}


process MULTIQC_VCF { 
	container "dnavc-pipeline:latest"
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files
		path conpair_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o variants_vcf .
	'''
}

















