#####################################
## the input_vcf file is prepared from merge_datasets.r
## otherinfo is retained in the annoated file to reserve the patient ID
## Note that the annotation output of would be changed in updated refGene 
## The annotation process is performed in 2022-01-18
####################################
anno_dict=/location_to_annovar/humandb
table_annovar=/location_to_annovar/table_annovar.pl
input_vcf=/ESCC_meta_2022/ESCC_META_all.vcf
out_pre=${input_vcf/.vcf/}

$table_annovar $input_vcf $anno_dict -buildver hg38 -out $out_pre -remove \
-protocol refGene \
-operation g \
--otherinfo -nastring . -vcfinput -polish --thread 35
