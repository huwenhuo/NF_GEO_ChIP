nextflow.enable.dsl=2

params.genome = 'GRCh38'
params.csv_file = "samples.csv"
params.reference_genome = params.genomes[params.genome].fasta
params.index_bowtie2      = params.genomes[params.genome].bowtie2
params.reference_gtf    = params.genomes[params.genome].gtf
params.adapters = "adapters.fa"

Channel
    .fromPath(params.csv_file)
    .splitCsv(header: true, sep: '\t')
    .map { row -> 
        return tuple(row.sample_name, row.fastq)
    }
    .set { sample_run_ch }

process TRIM_ADAPTERS {                                                         
    tag "Trim adapters for ${sample_name}"                                      
                                                                                
    input:                                                                      
    tuple val(sample_name), 
    path(fastq1) 
                                                                                
    output:                                                                     
    tuple val(sample_name),                                                     
    path("${sample_name}_1.trimmed.fastq.gz"),                                  
                                                                                
    script:                                                                     
    """                                                                         
    singularity exec ~/images/trimmomatic_latest.sif trimmomatic SE -threads 4 ${fastq1} ${sample_name}_1.trimmed.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """                                                                         
}                                                                               

process ALIGN_READS {                                                           
    tag "Align ${sample_name} to genome"                                        
                                                                                
    input:                                                                      
    tuple val(sample_name),                                                     
    path(trimmed1)
                                                                                
    output:                                                                     
    tuple val(sample_name), 
    path("${sample_name}.sorted.bam"),                   
    path("${sample_name}.sorted.bam.bai")                   
                                                                                
    script:                                                                     
    """                                                                         
    singularity exec ${params.img_bowtie2} bowtie2 -p 9 -x ${params.index_bowtie2} -U ${trimmed1} | singularity exec ${params.img_samtools} samtools view -Sb - | singularity exec ${params.img_samtools} samtools sort -@ 5 -o ${sample_name}.sorted.bam

    singularity exec ${params.img_samtools} samtools index ${sample_name}.sorted.bam

    """                                                                         
}                                                                               

process REMOVE_DUPLICATES {                                                     
    tag "Remove duplicates for ${sample_name}"                                  
                                                                                
    input:                                                                      
    tuple val(sample_name), 
    path(bam_file),
    path(bai_file)                                      
                                                                                
    output:                                                                     
    tuple val(sample_name), 
    path("${sample_name}.dedup.bam"),
    path("${sample_name}.dedup.bam.bai")                    
                                                                                
    script:                                                                     
    """                                                                         
    singularity exec ${params.img_picard} picard MarkDuplicates I=${bam_file} O=${sample_name}.dedup.bam M=${sample_name}.metrics.txt REMOVE_DUPLICATES=true                     
    singularity exec ${params.img_samtools} samtools index ${sample_name}.dedup.bam

    """                                                                         
}                                                                               

process CALL_PEAKS {                                                            
    tag "Call peaks for ${sample_name}"                                         
                                                                                
    input:                                                                      
    tuple val(sample_name), 
    path(readtype_file),                                                        
    path(dedup_bam_file),
    path(dedup_bai_file)                                      
                                                                                
    output:                                                                     
    tuple val(sample_name), 
    path(dedup_bam_file),
    path(dedup_bai_file),
    path("${sample_name}_peaks.narrowPeak")             

    publishDir path: "${params.publish_dir}/peaks", mode: 'copy', pattern: "*${sample_name}_peaks.narrowPeak"
                                                                                
    script:                                                                     
    """                                                                         
    singularity exec ${params.img_macs2} macs2 callpeak -t ${dedup_bam_file} -f BAM -g hs --name ${sample_name} --outdir .
    """                                                                         
}                                                                               

process GENERATE_BIGWIG {                                                       
    tag "Convert BAM to BigWig for ${sample_name}"                              
                                                                                
    input:                                                                      
    tuple val(sample_name), 
    path(dedup_bam_file),
    path(dedup_bai_file),
    path(peak_file)                                      
                                                                                
    output:                                                                     
    tuple val(sample_name), path("${sample_name}.bw")                           

    publishDir path: "${params.publish_dir}/bigwig", mode: 'copy', pattern: "*${sample_name}.bw"
                                                                                
    script:                                                                     
    """                                                                         
    singularity exec ${params.img_deeptools} bamCoverage -b ${dedup_bam_file} -o ${sample_name}.bw --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 2013022398
                                                                                
    """                                                                         
}                                                                               

// hg19: 2864785220                                                        
// hg38: 2913022398                                                        
                                                                                

workflow {

    parsed_samples = sample_run_ch

    trimmed_samples = parsed_samples  | TRIM_ADAPTERS                             

    aligned_samples = trimmed_samples | ALIGN_READS                             
                                                                                
    dedup_samples = aligned_samples   | REMOVE_DUPLICATES                         
                                                                                
    peak_samples = dedup_samples      | CALL_PEAKS                                   
                                                                                
    bigwig_samples = peak_samples     | GENERATE_BIGWIG                             

}
