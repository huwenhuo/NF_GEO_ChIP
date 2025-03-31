ChipSeq Analysis Pipeline
=========================

This repository provides a streamlined workflow for performing ChipSeq analysis, starting from a sample sheet with GSM IDs and sample names. The workflow proceeds through several stages, including data retrieval, conversion, and analysis. Below are the steps that the pipeline follows:

1.  **Input Data**  
    The process begins with a sample sheet containing GSM IDs and sample names.
    
2.  **SRR ID Identification**  
    Based on the provided GSM IDs, the pipeline identifies corresponding SRR IDs.
    
3.  **Download BAM Files**  
    The SRR IDs are then used to download the associated BAM files.
    
4.  **Convert BAM to FASTQ**  
    The downloaded BAM files are converted to FASTQ format, ready for further analysis.
    
5.  **ChipSeq Analysis**  
    The pipeline follows the standard ChipSeq analysis steps, including quality control, alignment, and peak calling.
    
6.  **Output Files**  
    After completing the analysis, the pipeline generates the following files:
    
    *   **Peak file**: Contains the identified peaks from the ChipSeq data.
        
    *   **BigWig files**: Provide visualization-ready output for genome browser tools.
        
