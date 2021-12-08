
maf_to_vcf<-function(
  maf_df=tcga_escc_maf_all_vcf,
  ref_allele_col="Reference_Allele",
  alr_allele_col="Tumor_Seq_Allele2",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start_Position",
  genome=BSgenome.Hsapiens.UCSC.hg19
){
  
  require(Biostrings)
  
  maf_df[,chro_id_col]->maf_df$CHROM
  maf_df[,start_col]->maf_df$Start_Position
  
  insert_index<- maf_df[,ref_allele_col]=="-" 
  del_index<- maf_df[,alr_allele_col]=="-" 
  
  maf_df[insert_index,]->insert_df
  maf_df[del_index,]->del_df
  maf_df[(!insert_index)*(!del_index)==1,]->snp_df
  
  
  as.character(
    getSeq(genome, insert_df$CHROM,
           start=insert_df$Start_Position, 
           end=insert_df$Start_Position))->left_bs
  
  insert_df$REF<-left_bs
  insert_df$ALT<-paste(left_bs,insert_df[,alr_allele_col],sep = "")
  insert_df$POS<-insert_df$Start_Position
  
  
  as.character(
    getSeq(genome, del_df$CHROM,
           start=del_df$Start_Position-1, 
           end=del_df$Start_Position-1))->left_bs2
  
  
  del_df$REF<-paste(left_bs2,insert_df[,ref_allele_col],sep = "")
  del_df$ALT<-left_bs2
  del_df$POS<-del_df$Start_Position-1
  
  snp_df$REF<-snp_df[,ref_allele_col]
  snp_df$ALT<-snp_df[,alr_allele_col]
  snp_df$POS<-snp_df$Start_Position
  
  rbind(
    snp_df[,c("CHROM","POS","REF","ALT",sample_id_col)],
    insert_df[,c("CHROM","POS","REF","ALT",sample_id_col)],
    del_df[,c("CHROM","POS","REF","ALT",sample_id_col)]
  )->vcf_df
  
  return(vcf_df)
}


