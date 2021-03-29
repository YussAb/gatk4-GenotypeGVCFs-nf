#!/usr/bin/env nextflow

//ref
//genomicdbimport
//https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport
//joint calling
//https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants
//variant filtering
//https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering


params.help = null

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  gatk4-GenotypeGVCFs v1: Exact Joint Genotyping GATK4 Best Practices         "
log.info "-------------------------------------------------------------------------"
log.info "yussline_vda"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/gatk4-GenotypeGVCFs-nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input                         VCF FILES                 All cohort gVCF files (between quotes)"
    log.info "--output_dir                    OUTPUT FOLDER             Output for VCF file"
    log.info "--cohort                        STRING                    Cohort name"
    log.info "--ref_fasta                     FASTA FILE                Reference FASTA file"
    log.info "--gatk_exec                     BIN PATH                  Full path to GATK4 executable"
    log.info "--dbsnp                         VCF FILE                  dbSNP VCF file"
    log.info "--mills                         VCF FILE                  Mills and 1000G gold standard indels VCF file"
    log.info "--axiom                         VCF FILE                  Axiom Exome Plus genotypes all populations poly VCF file"
    log.info "--hapmap                        VCF FILE                  hapmap VCF file"
    log.info "--omni                          VCF FILE                  1000G omni VCF file"
    log.info "--onekg                         VCF FILE                  1000G phase1 snps high confidence VCF file"
    exit 1
}

//
// Parameters Init
params.interval = "/hpcshare/genomics/references/gatk_bundle/resources/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
//
params.sample = ""
//
params.input         = null
params.outdir        = "results"
params.cohort        = "cohort_db" 
params.ref_fasta     = "/hpcshare/genomics/references/gatk_bundle/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
params.gatk_exec     = "gatk"
params.dbsnp         = "/hpcshare/genomics/references/gatk_bundle/resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
params.mills         = "/hpcshare/genomics/references/gatk_bundle/resources/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.axiom         = "/hpcshare/genomics/references/gatk_bundle/resources/resources_broad_hg38_v0_Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
params.hapmap        = "/hpcshare/genomics/references/gatk_bundle/resources/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
params.omni          = "/hpcshare/genomics/references/gatk_bundle/resources/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
params.onekg         = "/hpcshare/genomics/references/gatk_bundle/resources/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"

//
// Parse Input Parameters
//
gvcf_ch = Channel
			.fromPath(params.input)

gvcf_idx_ch = Channel
			.fromPath(params.input)
			.map { file -> file+".tbi" }

			
GATK                              = params.gatk_exec
ref                               = file(params.ref_fasta)
dbsnp_resource_vcf                = file(params.dbsnp)
mills_resource_vcf                = file(params.mills)
axiomPoly_resource_vcf            = file(params.axiom)
hapmap_resource_vcf               = file(params.hapmap)
omni_resource_vcf                 = file(params.omni)
one_thousand_genomes_resource_vcf = file(params.onekg)



//
process GenomicsDBImport {

    
    publishDir "${params.outdir}/genomicdb", mode:'copy'

    input:
    file (gvcf) from gvcf_ch.collect()
    file (gvcf_idx) from gvcf_idx_ch.collect()

	  output:
    file ("${params.cohort}") into gendb_ch
	
    script:
  	"""
  	${GATK} GenomicsDBImport --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/tmp" \
  	${gvcf.collect { "-V $it " }.join()} \
          -L $params.interval  \
        	--genomicsdb-workspace-path ${params.cohort} \
          -R ${params.ref_fasta}
  	"""
}	

//
// Process launching GenotypeGVCFs on the previously created genDB, per chromosome
//
process GenotypeGVCFs {

    publishDir "${params.outdir}/genotype_gvcf", mode:'copy' , pattern: '*.{vcf,idx}'
  
    input:
	  file (workspace) from gendb_ch
   	file genome from ref

	output:
    set file("${params.sample}.vcf"), file("${params.sample}.vcf.idx") into (vcf_ch,vcf_snv_ch, vcf_sid_ch, vcf_recal_ch)
    file "${genome}.fai" into faidx_sid_ch,faidx_snv_ch
	  file "${genome.baseName}.dict" into dict_sid_ch,dict_snv_ch

    script:
  	"""
    samtools faidx ${genome}

    java -jar /apps/picard/2.17.11/picard-2.17.11.jar \
    CreateSequenceDictionary \
    R=${genome} \
    O=${genome.baseName}.dict

    WORKSPACE=\$( basename ${workspace} )

    ${GATK} --java-options "-Xmx5g -Xms5g" \
     GenotypeGVCFs \
     -R ${genome} \
     -O ${params.sample}.vcf \
     -D ${dbsnp_resource_vcf} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --use-new-qual-calculator \
     -V gendb://\$WORKSPACE \
     -L $params.interval

	"""
}	


//
// Process SID recalibration
//
process SID_VariantRecalibrator {

    input:
	set file (vcf), file (vcfidx) from vcf_sid_ch
    file genome from ref
    file faidx from faidx_sid_ch
    file dict from dict_sid_ch

	output:
    set file("${params.sample}.sid.recal"),file("${params.sample}.sid.recal.idx"),file("${params.sample}.sid.tranches") into sid_recal_ch

    script:
	"""
    ${GATK} --java-options "-Xmx24g -Xms24g" \
      VariantRecalibrator \
      -R ${genome} \
      -V ${vcf} \
      --output ${params.sample}.sid.recal \
      --tranches-file ${params.sample}.sid.tranches \
      --trust-all-polymorphic \
      -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
      -mode INDEL \
      --max-gaussians 4 \
      --resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \
      --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf}
	
	"""
}	



//
// Process SNV recalibration
//
process SNV_VariantRecalibrator {

    input:
	set file (vcf), file (vcfidx) from vcf_snv_ch
    file genome from ref
    file faidx from faidx_snv_ch
    file dict from dict_snv_ch

	output:
    set file("${params.sample}.snv.recal"),file("${params.sample}.snv.recal.idx"),file("${params.sample}.snv.tranches") into snv_recal_ch

    script:
	"""
    ${GATK} --java-options "-Xmx24g -Xms24g" \
      VariantRecalibrator \
      -R ${genome} \
      -V ${vcf} \
      --output ${params.sample}.snv.recal \
      --tranches-file ${params.sample}.snv.tranches \
      --trust-all-polymorphic \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
      -mode SNP \
      --max-gaussians 6 \
      --resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
      --resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
      --resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf}
	
	"""
}	



//
// Process Apply SNV and SID recalibrations
//
process ApplyRecalibration {

  	cpus 1 
  	memory '7 GB'
  	time '12h'
	
	tag "${params.cohort}"

	//publishDir params.output_dir, mode: 'copy'
  publishDir "${params.outdir}/vqsr", mode:'copy' 
 
    input:
	set file (input_vcf), file (input_vcf_idx) from vcf_recal_ch
	set file (indels_recalibration), file (indels_recalibration_idx), file (indels_tranches) from sid_recal_ch
	set file (snps_recalibration), file (snps_recalibration_idx), file (snps_tranches) from snv_recal_ch

	output:
    set file("${params.sample}.recalibrated.vcf"),file("${params.sample}.recalibrated.vcf.idx") into vcf_final_ch

    script:
	"""
    ${GATK} --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level 99.0 \
      --exclude-filtered \
      --create-output-variant-index true \
      -mode INDEL

    ${GATK} --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O ${params.sample}.recalibrated.vcf \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ${snps_recalibration} \
      --tranches-file ${snps_tranches} \
      --truth-sensitivity-filter-level 99.5 \
      --exclude-filtered \
      --create-output-variant-index true \
      -mode SNP
		
	"""
}	



/*
////JUNK_CODE

// ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
// than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
excess_het_threshold = 54.69

/*
//
// Process Hard Filtering on ExcessHet, per chromosome
//
process HardFilter {

	cpus 1
	memory '24 GB'
	time '12h'
	
	tag { chr }

    input:
	set chr, file (vcf), file (vcfidx) from vcf_ch

	output:
    file("${params.cohort}.${chr}.filtered.vcf") into (vcf_hf_ch)
    file("${params.cohort}.${chr}.filtered.vcf.idx") into (vcf_idx_hf_ch)

    script:
	"""
	${GATK} --java-options "-Xmx3g -Xms3g" \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -V ${vcf} \
      -O ${params.cohort}.${chr}.markfiltered.vcf

	${GATK} --java-options "-Xmx3g -Xms3g" \
      SelectVariants \
      --exclude-filtered \
      -V ${params.cohort}.${chr}.markfiltered.vcf \
      -O ${params.cohort}.${chr}.filtered.vcf

	"""
}	


process GatherVcfs {


    input:
   	set file (vcf), file (vcfidx) from vcf_ch.collect()
    
	  output:
    set file("${params.sample}.vcf"), file("${params.sample}.vcf.idx") into (vcf_snv_ch, vcf_sid_ch, vcf_recal_ch)

    // WARNING : complicated channel extraction! 
    // GATK GatherVcfs only accepts as input VCF in the chromosomical order. Nextflow/Groovy list are not sorted. The following command does :
    // 1 : look for all VCF with "chr[0-9]*" in the filename (\d+ means 1 or + digits)
    // 2 : Tokenize the filenames with "." as the separator, keep the 2nd item (indexed [1]) "chr[0-9]*"
    // 3 : Take from the 3rd character till the end of the string "chr[0-9]*", ie the chromosome number
    // 4 : Cast it from a string to an integer (to force a numerical sort)
    // 5 : Sort 
    // 6 : Add chrX and chrY to the list

    script:
	"""
	${GATK} --java-options "-Xmx3g -Xms3g" GatherVcfs \
      ${vcf.findAll{ it=~/chr\d+/ }.collect().sort{ it.name.tokenize('.')[1].substring(3).toInteger() }.plus(vcf.find{ it=~/chrX/ }).plus(vcf.find{ it=~/chrY/ }).collect{ "--INPUT $it " }.join() } \
      --OUTPUT ${params.sample}.vcf
	"""
}	

*/
