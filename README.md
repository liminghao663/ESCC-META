# ESCC-META

This repo contains the code for establishing ESCC-META dataset.

We firstly collected all public avaible ESCC genomes from raw sequcence databases like SRA, cancer genome databases like TCGA and ICGC, and mutational records from published article.
We also processed and checked the available clinical information of each individual, including age, sex, drinking and smoking history, tumor stage, tumor location and tumor grade. 
The patient ID was renamed by pasting their source and original sample name. 
If the raw sequence data were also available, we directly included and re-analysis their raw sequence data, ignoring their published results. 
If the cohort was both involved in published article and genome database, we compared them and used the more complete records. 
If the multiple tumor samples were collected from different time point, only the earliest tumor sample (at diagnosis or before any treatment) were used. 

In our work, we separated the dataset of SRP034680 into SRP072858_WGS (data of WGS part) and SRP072858_WES (data of WES part), the SRP072858 into SRP072858_WGS (data of WGS part) and SRP072858_WES (data of WES part), the SRP099292 into SRP099292_S (single tumor sample per patient) and SRP099292_M (multiple tumor sample per patient). 


The codes in prepare_from_mutational_list is for processing raw sequence data

If the reads data were NCBI SRA format, it was converted to fastq file by SRA-Tools (v2.11). The fastq files were firstly performed quality control by fastp (v0.23) with default parameters. 
The files of different sequence lanes from the same library, or the different SRA reads files from the same sample were combined before mapping. 
The mapping process was performed by BWA (v0.7.1) to hg38.p13 genome. 
The bam files were then deduplicated and applied base quality score recalibration by GATK (v4.1) according to the recommended practice. 
The single nucleotide variants (SNVs) and insertion or deletion mutations (INDELs) were called by Mutect2 in GATK (v4.1). 
The filter criteria included the following three situations. 

1, for WES sequence of one tumor sample with normal control, the coverage should be at least 30 in tumor sample and at least 20 in normal sample, at least three alternative reads in tumor sample to support the variant call, and mutation frequency at least 0.05. 
2, for WGS sequence of one tumor sample with normal control, the coverage should be at least 20 in tumor sample and at least 20 in normal sample, the rest is the same as above. 
3, for WES sequence of multiple tumor samples with one normal control from the same patient the included variant should be detected in at least one tumor sample meeting the above type 2 criteria, additionally the variant base should be identified (at least two alternative support reads) in another tumor sample. 


The codes in prepare_from_SRA is for preparation of mutational records

For mutational records from reported list (such as in MAF format) or database without raw data, they were prepared in the following two steps.  
Firstly, the authenticity of obtained variant list was checked by comparison of the provided reference base to the loci in genome of it used. For example, if the raw mutational record is chr19:63554635, G>T in hg18, the true base of chr19:63554635 in hg18 should be G, if not, this record was suspicious and must be re-examined or exclude. 
Secondly , the verified records from each dataset were prepared to VCF format (Version 4.2) and transformed to hg38, which will also remove a few records because of failure to convert. The converted mutational lists were re-check with hg38 as the first step.


The codes in combined_and_annotation is for the datasets integration,annotation and conversion to MAF file.

The processed data of integrated ESCC-META dataset are available at synapse under the accession code of syn27304838.
