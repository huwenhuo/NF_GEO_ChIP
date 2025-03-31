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
        return tuple(row.gsm_id, row.sample_name)
    }
    .set { sample_run_ch }

process EXTRACT_SRR {
    tag "Extract SRR for ${gsm_id}"

    input:
    tuple val(gsm_id), val(sample_name)

    output:
    tuple val(gsm_id), val(sample_name), 
    path("${gsm_id}_srr_list.txt"),
    path("${gsm_id}_read_type.txt")

    script:
    """
    singularity exec ${params.img_edirect} \\
        esearch -db gds -query ${gsm_id} | \\
        elink -target sra | \\
        efetch -format runinfo > ${gsm_id}_runinfo.txt

    # Extract SRR IDs
    grep ${gsm_id} ${gsm_id}_runinfo.txt | cut -d ',' -f 1 | grep SRR > ${gsm_id}_srr_list.txt

    # Determine read type
    if grep -q ${gsm_id} ${gsm_id}_runinfo.txt | grep -q 'PAIRED'; then
        read_type='pe'
        echo "pe" > ${gsm_id}_read_type.txt
    else
        read_type='se'
        echo "se" > ${gsm_id}_read_type.txt
    fi

    """
    
}


process DOWNLOAD_SRA {                                                          
                                                                                
    tag "Download SRA for ${sample_name}"                                       
                                                                                
    input:                                                                      
    tuple val(gsm_id),                                                          
    val(sample_name),                                                           
    path(srr_file),                                                             
    path(readtype_file)                                                         
                                                                                
    output:                                                                     
    tuple val(sample_name),                                                     
    path(readtype_file),                                                        
    path("${sample_name}_1.fastq"),                                          
    path("${sample_name}_2.fastq")

    script: 
    """                                                                         
    mkdir -p ${sample_name}_fastq                                                  
                                                                                
    read_type=\$(cat "${readtype_file}")                                        
                                                                                
    while read srr_id; do                                                       
        singularity exec ${params.img_sratool} /opt/sratoolkit.3.1.0-ubuntu64/bin/prefetch "\$srr_id"

        if [ "\$read_type" = "pe" ]; then                                     
            singularity exec ${params.img_sratool} /opt/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump --split-files "\$srr_id" -O ${sample_name}_fastq     
        else                                                                    
            singularity exec ${params.img_sratool} /opt/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump "\$srr_id" -O ${sample_name}_fastq                   
            mv ${sample_name}_fastq/\${srr_id}.fastq ${sample_name}_fastq/\${srr_id}_1.fastq 
            touch ${sample_name}_fastq/\${srr_id}_2.fastq 
        fi                                                                      
    done < "${srr_file}"                                                          

    pattern="${sample_name}_fastq/*_1.fastq"
    file_count=\$(ls \$pattern 2>/dev/null | wc -l)
    if [ \$file_count -eq 1 ]; then
        mv \$(ls \$pattern) ${sample_name}_1.fastq
    else
        cat \$(ls \$pattern | sort) > ${sample_name}_1.fastq
    fi

    pattern="${sample_name}_fastq/*_2.fastq"
    file_count=\$(ls \$pattern 2>/dev/null | wc -l)
    if [ \$file_count -eq 1 ]; then
        mv \$(ls \$pattern) ${sample_name}_2.fastq
    else
        cat \$(ls \$pattern | sort) > ${sample_name}_2.fastq
    fi

    """ 
}  

process TRIM_ADAPTERS {                                                         
    tag "Trim adapters for ${sample_name}"                                      
                                                                                
    input:                                                                      
    tuple val(sample_name), 
    path(readtype_file),                                                        
    path(fastq1), 
    path(fastq2)
                                                                                
    output:                                                                     
    tuple val(sample_name),                                                     
    path(readtype_file),                                                        
    path("${sample_name}_1.trimmed.fastq.gz"),                                  
    path("${sample_name}_2.trimmed.fastq.gz")
                                                                                
    script:                                                                     
    """                                                                         
    read_type=\$(cat "${readtype_file}")                                        

    if [ "\${read_type}" == "pe" ]; then                                         
        singularity exec ~/images/trimmomatic_latest.sif trimmomatic PE -threads 4 ${fastq1} ${fastq2} ${sample_name}_1.trimmed.fastq.gz ${sample_name}_1.unpaired.fastq.gz ${sample_name}_2.trimmed.fastq.gz ${sample_name}_2.unpaired.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else                                                                        
        singularity exec ~/images/trimmomatic_latest.sif trimmomatic SE -threads 4 ${fastq1} ${sample_name}_1.trimmed.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        touch ${sample_name}_2.trimmed.fastq.gz 
    fi                                                                          
    """                                                                         
}                                                                               

process ALIGN_READS {                                                           
    tag "Align ${sample_name} to genome"                                        
                                                                                
    input:                                                                      
    tuple val(sample_name),                                                     
    path(readtype_file),                                                        
    path(trimmed1),                                  
    path(trimmed2)                                  
                                                                                
    output:                                                                     
    tuple val(sample_name), 
    path(readtype_file),                                                        
    path("${sample_name}.sorted.bam"),                   
    path("${sample_name}.sorted.bam.bai")                   
                                                                                
    script:                                                                     
    """                                                                         
    read_type=\$(cat "${readtype_file}")                                        

    if [ "\${read_type}" == "pe" ]; then                                         
        singularity exec ${params.img_bowtie2} bowtie2 -p 9 -x ${params.index_bowtie2} -1 ${trimmed1} -2 ${trimmed2} | singularity exec ${params.img_samtools} samtools view -Sb - | singularity exec ${params.img_samtools} samtools sort -@ 5 -o ${sample_name}.sorted.bam
    else                                                                        
        singularity exec ${params.img_bowtie2} bowtie2 -p 9 -x ${params.index_bowtie2} -U ${trimmed1} | singularity exec ${params.img_samtools} samtools view -Sb - | singularity exec ${params.img_samtools} samtools sort -@ 5 -o ${sample_name}.sorted.bam
    fi                                                                          

    singularity exec ${params.img_samtools} samtools index ${sample_name}.sorted.bam

    """                                                                         
}                                                                               

process REMOVE_DUPLICATES {                                                     
    tag "Remove duplicates for ${sample_name}"                                  
                                                                                
    input:                                                                      
    tuple val(sample_name), 
    path(readtype_file),                                                        
    path(bam_file),
    path(bai_file)                                      
                                                                                
    output:                                                                     
    tuple val(sample_name), 
    path(readtype_file),                                                        
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
    path(readtype_file),                                                        
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
    path(readtype_file),                                                        
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

    srr_samples = parsed_samples | EXTRACT_SRR

    fastq_samples = srr_samples | DOWNLOAD_SRA

    trimmed_samples = fastq_samples | TRIM_ADAPTERS                             

    aligned_samples = trimmed_samples | ALIGN_READS                             
                                                                                
    dedup_samples = aligned_samples | REMOVE_DUPLICATES                         
                                                                                
    peak_samples = dedup_samples | CALL_PEAKS                                   
                                                                                
    bigwig_samples = peak_samples | GENERATE_BIGWIG                             

}
