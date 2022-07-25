

pre_from_vcf<-function(
  study="SRP116657",
  all_vcfs=list.files(path = "E:/escc/29142225/PASS",
                      pattern = ".PASS.vcf",
                      full.names = T),
  
  min_cov_ref=20,
  min_cov=30,
  min_AF=0.05,
  min_alt=3,
  min_cov2=10,
  min_AF2=0.01,
  min_alt2=1,
  min_sample=1
  # the min_alt2,min_AF2 and min_cov2 are need to specified when min_sample>1
){
  
  temp_merge_df<-data.frame()
  
  for(i in 1:length(all_vcfs)){
    
    sample_name<-gsub(".PASS.vcf","",fixed = T,basename(all_vcfs[i]))
    
    print(i)
    print(sample_name)
    
    dat <- readLines(all_vcfs[i])
    skip <- max(grep("##", dat))
    
    
    gsub("##normal_sample=","",dat[grep("##normal_sample=", dat)])->normal_id
    
	
    read.delim(
      all_vcfs[i],
      sep = "\t",
      skip=skip,
      header = T,
      stringsAsFactors=FALSE)->temp_vcf
    
    colnames(temp_vcf)[1]<-"CHROM"
    
    temp_vcf$sample_id<-paste(study,sample_name,sep = ";")
    
    dat[skip+1]->ggg1
    str_split(ggg1,"\t")[[1]][-(1:9)]->all_samples
    setdiff(all_samples,normal_id)->tumor_samples

	
    make.names(all_samples)->all_samples
    make.names(tumor_samples)->tumor_samples
    
    paste(make.names(normal_id),"Coverage",sep=":")->ref_cover_col
    
    for(ii in 1:length(all_samples)){
      cbind(
        get_mut_info(temp_vcf[,all_samples[ii]],
                     all_samples[ii]),
        temp_vcf
      )->temp_vcf
    }
    
    temp_vcf$n_pass_tumor<-0
    temp_vcf$n_pass_tumor2<-0
    
    
    for(iii in 1:length(tumor_samples)){
      
      paste(tumor_samples[iii],"Coverage",sep=":")->tmp_cover_col
      paste(tumor_samples[iii],"AF_cat",sep=":")->tmp_AF_col
      paste(tumor_samples[iii],"ALT",sep=":")->tmp_alt_col
      
	  
      temp_vcf[,tmp_cover_col]>=min_cov &
        temp_vcf[,tmp_AF_col]>=min_AF &
        temp_vcf[,tmp_alt_col]>=min_alt -> tmp_index
      
      temp_vcf[,tmp_cover_col]>=min_cov2 &
        temp_vcf[,tmp_AF_col]>=min_AF2 &
        temp_vcf[,tmp_alt_col]>=min_alt2 -> tmp_index2
      
      
      print(paste("mean AF of raw SNVs",tumor_samples[iii], mean(temp_vcf[,tmp_AF_col],na.rm=T)))
       
      temp_vcf$n_pass_tumor<-  temp_vcf$n_pass_tumor + as.numeric(tmp_index)
      temp_vcf$n_pass_tumor2<-  temp_vcf$n_pass_tumor2 + as.numeric(tmp_index2)
    }
    
	
    if(min_sample==1){
	  #for paired tumor-normal samples
      temp_vcf[which(temp_vcf[,ref_cover_col]>=min_cov_ref &
              temp_vcf$n_pass_tumor>=1),
			  ]->temp_vcf2
			  }else{
	  #for multiplie tumor samples and one normal samples from same patients
      temp_vcf[which(temp_vcf[,ref_cover_col]>=min_cov_ref &
               temp_vcf$n_pass_tumor>=1 &
               temp_vcf$n_pass_tumor2>=min_sample),
				]->temp_vcf2}
    
    
    print(paste("passed rows",nrow(temp_vcf2)))
    
    for(iiii in 1:length(tumor_samples)){
	
	paste(tumor_samples[iiii],"AF_cat",sep=":")->tmp_AF_col
	print(paste("mean AF of passed SNVs",tumor_samples[iiii], mean(temp_vcf2[,tmp_AF_col],na.rm=T)))
	  
	  }
    
    
    rbind(temp_merge_df,
          temp_vcf2[,c("CHROM","POS","REF","ALT","sample_id")]
    )->temp_merge_df
    
  }
  
  return(temp_merge_df)
}


get_mut_info<-function(
  x,pre_label
){
  colsplit(x,":",1:50)->tmmm
  tmmm[,4]->v_deeps
  colsplit(tmmm[,2],",",1:2)->tmmm1
  
  cbind(v_deeps,tmmm1)->tmmm2
  tmmm2$af<-tmmm2[,3]/tmmm2[,1]
  
  colnames(tmmm2)<-paste(pre_label,c("Coverage","REF",'ALT',"AF_cat"),
                         sep=":")
  return(tmmm2)}


#Example
######################################3
#In SRP127593, min_sample=1
#################################

#pre_from_vcf(
#  study="SRP127593",
#  all_vcfs=list.files(path = "SRP127593/PASS_vcf",
#                      pattern = ".PASS.vcf",
#                      full.names = T),
#  
#  min_cov_ref=20,
#  min_cov=30,
#  min_AF=0.05,
#  min_alt=3,
#  min_sample=1)->mutation_hg38

#In SRP099292_M, min_sample=2

#pre_from_vcf(
#  study="SRP099292_M",
#  all_vcfs=list.files(path = "SRP099292_M/PASS_vcf",
#                      pattern = ".PASS.vcf",
#                      full.names = T), 
#  min_cov_ref=20,
#  min_cov=20,
#  min_AF=0.05,
#  min_cov2=10,
#  min_AF2=0.01,
#  min_alt=3,
#  min_alt2=1,
#  min_sample=2)->mutation_hg38

