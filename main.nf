
/* 
 * SOMATIC-VC PIPELINE 
 * for paired normal-tumor mapped reads
 */


/* 
 * import modules 
 */

nextflow.enable.dsl=2

include { 
	INDEX_REFERENCE;
	SOMVC_LOFREQ;
	SOMVC_MUTECT2;
	SOMVC_STRELKA;
	SOMVC_VARDICT;
	SOMATIC_COMBINER;
	CONPAIR_CONTAMINATION;
	VARIANT_CALLING_STATS;
	MULTIQC_VCF;
	VARIANT_ANNOTATION;
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev_samples = -1

params.project_dir	= "$projectDir"
params.sample_match_file	= "$params.project_dir/normal_tumor_pairs.csv"
params.data_dir		= "$params.project_dir/data"
params.scripts_dir	= "$params.project_dir/scripts"


/*
 * other parameters
 */

params.num_threads		= 3
params.reference_genome	= "$params.project_dir/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
params.bed_file	= "$params.project_dir/data/Homo_sapiens.GRCh38.cds.all.bed"



log.info """\
SOMATIC-VC PIPELINE
===================================================
sample_match_file		: $params.sample_match_file
data_dir		: $params.data_dir
reference_genome	: $params.reference_genome
bed_file		: $params.bed_file

===================================================

"""



/* 
 * main pipeline logic
 */
workflow {

	channel_sample_match = Channel
			.fromPath(params.sample_match_file)
			.splitCsv(header:true)
			.map{ row -> tuple(row.sample_id, row.normal_file, row.tumor_file ) }
			//.map{ row -> tuple(row.sample_id, tuple(row.normal_file, row.tumor_file)) }
			.ifEmpty { error "cannot read in sample_match_file correctly: ${params.sample_match_file}" }
			.take( params.dev_samples )  // only consider a few files for debugging
	channel_sample_match.view()


	INDEX_REFERENCE(params.reference_genome)

	SOMATICSEQ_CALLING(channel_sample_match, INDEX_REFERENCE.out.reference_genome)


	//VARIANT_CALLING(channel_reads_mapped, params.num_threads, INDEX_REFERENCE.out.reference_genome)
	//VARIANT_CALLING_STATS(VARIANT_CALLING.out.vcf, params.num_threads)
	//MULTIQC_VCF(VARIANT_CALLING_STATS.out.vcf_stats.collect())
	
	//VARIANT_MERGING(VARIANT_CALLING.out.global_vcf.collect(), params.num_threads)
	//GLNEXUS_BCF_TO_VCF(VARIANT_MERGING.out.glnexus_bcf, params.num_threads)

	//VARIANT_ANNOTATION(GLNEXUS_BCF_TO_VCF.out.pvcf_file, params.num_threads)

}




workflow.onComplete { 
	println ( workflow.success ? "\ndone! check the quality reports in --> $params.data_dir/quality_reports\n" : "oops .. something went wrong" ) } 


// ~/tools/nextflow run main.nf --project_dir ../somatic_test/ --reference_genome /home/stefan/Documents/umcg/somatic_test/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz






