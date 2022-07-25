
pre_mutation<-function(
  maf_df=tcga_escc_maf_all_vcf,
  ref_allele_col="Reference_Allele",
  alr_allele_col="Tumor_Seq_Allele2",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start_Position",
  genome=BSgenome.Hsapiens.UCSC.hg19,
  blank_symble="-"
  # the blank_symble represent the blank value in ref_allele_col or alr_allele_col some mutational list,might be "-","." or ""
  # In the MAF file, the "-" in ref_allele_col indicates insertion, and the "-" in alr_allele_col indicates deletion 
){
  
  require(Biostrings)
  print(paste(nrow(maf_df),"rows input"))
  
  paste(maf_df[,chro_id_col],
        maf_df[,start_col],
        maf_df[,ref_allele_col],
        maf_df[,alr_allele_col],
        maf_df[,sample_id_col],
        sep = ";")->mut_id0
  
  duplicated(mut_id0)->dup_index
  
  print(paste(sum(dup_index),"rows dup"))
  
  maf_df[!dup_index,]->maf_df
  
  
  maf_df[,chro_id_col]->maf_df$CHROM
  maf_df[,start_col]->maf_df$Start_Position
  maf_df[,start_col]+ str_length(maf_df[,ref_allele_col])-1 ->maf_df$temp_end
  
  
  insert_index<- maf_df[,ref_allele_col]==blank_symble
  maf_df[insert_index,]->insert_df
  
  if(nrow(insert_df)>0){
    as.character(
      getSeq(genome, insert_df$CHROM,
             start=insert_df$Start_Position, 
             end=insert_df$Start_Position))->left_bs
    
    insert_df$REF<-left_bs
    insert_df$ALT<-paste(left_bs,insert_df[,alr_allele_col],sep = "")
    insert_df$POS<-insert_df$Start_Position
    
  }
  
  
###################################################  
  maf_df[!insert_index,]->test_df
  print(paste(nrow(test_df),"rows to test"))
  
  as.character(
    getSeq(genome, test_df$CHROM,
           start=test_df$Start_Position, 
           end=test_df$temp_end))->pos_bs
  
  test_df[,ref_allele_col]!=pos_bs->exclude_index
  print(paste(sum(exclude_index),"rows excluded"))
  
  test_df[exclude_index,]->exclude_df
  test_df[!exclude_index,]->valid_df
  
  valid_df$REF<-valid_df[,ref_allele_col]
  valid_df$ALT<-valid_df[,alr_allele_col]
  valid_df$POS<-valid_df$Start_Position
  
  del_index<- valid_df[,alr_allele_col]==blank_symble
  
  if(sum(del_index)>0){
  
    as.character(
      getSeq(genome, valid_df[del_index,"CHROM"],
             start=valid_df[del_index,"Start_Position"]-1, 
             end=valid_df[del_index,"Start_Position"]-1))->left_bs2
    
    
    valid_df[del_index,"REF"]<-paste(left_bs2,valid_df[del_index,ref_allele_col],sep = "")
    valid_df[del_index,"ALT"]<-left_bs2
    valid_df[del_index,"POS"]<-valid_df[del_index,"Start_Position"]-1
    
  }
  

  
  if(nrow(insert_df)>0){
    rbind(
      valid_df[,c("CHROM","POS","REF","ALT",sample_id_col)],
      insert_df[,c("CHROM","POS","REF","ALT",sample_id_col)]
    )->vcf_df
  }else{
    valid_df[,c("CHROM","POS","REF","ALT",sample_id_col)]->vcf_df
  }
  
  print(paste(nrow(vcf_df),"rows prepared"))
  
  return(list(prepared_df=vcf_df,
              exclude_df=exclude_df
              ))
}


liftOver_convert<-function(
  used_chain=chain19to38,
  raw_df=mutations_hg19,
  genome=BSgenome.Hsapiens.UCSC.hg38
){
  

  liftOver(as(data.frame(
    chr=raw_df$CHROM,
    Start=raw_df$POS,
    End=raw_df$POS), 
    "GRanges"),used_chain)->converted_df0
  
  as.data.frame(converted_df0)->converted_df0
  
  converted_df0$group->valiad_rows
  setdiff(1:nrow(raw_df),valiad_rows)->invaliad_rows
  
  
  print(paste(length(invaliad_rows),"rows out of coordinate"))
  

  rownames(raw_df)<- 1:nrow(raw_df)  
  raw_df[valiad_rows,]->converted_df
  converted_df$POS<-converted_df0$start
  converted_df$CHROM<-converted_df0$seqnames
  
  
  converted_df[,"POS"]+ str_length(converted_df[,"REF"])-1 ->converted_df$temp_end
  
  as.character(
    getSeq(genome, converted_df$CHROM,
           start=converted_df$POS, 
           end=converted_df$temp_end))->pos_bs
  
  converted_df[,"REF"]!=pos_bs->exclude_index
  rownames(converted_df)[exclude_index]->invaliad_rows2
  
  
  raw_df[c(invaliad_rows,invaliad_rows2),]->invaliad_converted_df

  
  print(paste(nrow(invaliad_converted_df),"rows excluded"))
  
  
  converted_df[!exclude_index,
               -which(colnames(converted_df)=="temp_end")
               ]->converted_df2
  
  
  print(paste(nrow(converted_df2),"rows converted"))
  
  return(
    list(
      converted_df=converted_df2,
      invaliad_df=invaliad_converted_df
    )
  )
}
