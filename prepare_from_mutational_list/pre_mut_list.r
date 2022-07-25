



#https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
#https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz


library(dplyr)
library(stringr)
library(reshape2)
library(Biostrings)

library("BSgenome.Hsapiens.UCSC.hg18")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")

chain19to38 <- import.chain("hg19ToHg38.over.chain")
chain18to38 <- import.chain("hg18ToHg38.over.chain")

############################
#PMID32974170(hg19)
#######################

dat <- readLines("PMID32974170/from_pdf.txt")
used_lines<-grep("^Smoker|^Non-user", dat)
dat[used_lines]->used_dat


colsplit(used_dat," ",c(
  "group",
  "Patient",
  "Chr",
  "start",
  "REF",
  "ALT",
  "gene"))->used_dat1

used_dat1$sample_id<-paste("PMID32974170",
                           used_dat1$Patient,
                          sep=";")


pre_mutation(
  maf_df=used_dat1,
  ref_allele_col="REF",
  alr_allele_col="ALT",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="start",
  genome=BSgenome.Hsapiens.UCSC.hg19)->res_ls

write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32974170/mutations_hg19.csv")

liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32974170/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32974170/unconverted.csv")

############################################# 
#PMID25151357(hg18)
#########################################

read.csv(
  file = "PMID25151357/all_mut_hg18.csv",
  header=T,
  stringsAsFactors = F)->raw_df


within(
  raw_df,{
    sample_id<-paste("PMID25151357",Case,sep=";")
    Chr<-NA
    POS<-NA
    REF<-NA
    ALT<-NA}
)->raw_df


colsplit(raw_df$Nucleotide..genomic..hg18.,"_",
         c("chrom",
           "start_end",
           "ref",
           "alt"))->temp_hg18

colsplit(temp_hg18$start_end,"-",c("start","end"))->temp_hg18b


cbind(raw_df,
      temp_hg18,
      temp_hg18b)->raw_df1

raw_df1[raw_df1$ref=="","ref"]<-"-"
raw_df1[raw_df1$alt=="","alt"]<-"-"



pre_mutation(
  maf_df=raw_df1,
  ref_allele_col="ref",
  alr_allele_col="alt",
  sample_id_col="sample_id",
  chro_id_col="chrom",
  start_col="start",
  genome=BSgenome.Hsapiens.UCSC.hg18
)->res_ls

write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID25151357/mutations_hg18.csv")

write.csv(
  res_ls$exclude_df,
  row.names = F,
  quote=FALSE,
  file = "PMID25151357/exclude_row.csv")

liftOver_convert(
  used_chain=chain18to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID25151357/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID25151357/unconverted.csv")
  
#################################
##PMID34413060(hg19)
####################################

read.csv(file = "PMID34413060/mut_df.csv",
         header = T)->mut_df

mut_df$sample_id<-paste("PMID34413060",
                        mut_df$Patient.ID,sep=";")

mut_df$CHRO<-paste("chr",mut_df$Chr,sep="")


library(BSgenome.Hsapiens.UCSC.hg19)


pre_mutation(
  maf_df=mut_df,
  ref_allele_col="Ref",
  alr_allele_col="Obs",
  sample_id_col="sample_id",
  chro_id_col="CHRO",
  start_col="Start",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls

write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34413060/mutations_hg19.csv")


liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34413060/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34413060/unconverted.csv")

#######################################
#PMID30012096(hg19)
######################################


read.csv(
  file = "PMID30012096/all_mut.csv",
  header = T,
  stringsAsFactors = F
)->raw_df


within(
  raw_df,{
    sample_id<-paste("PMID30012096",Sample,sep=";")
    Chr<-paste("chr",Chromosome,sep="")})->raw_df

pre_mutation(
  maf_df=raw_df,
  ref_allele_col="Reference_Allele",
  alr_allele_col="Alternative_Allele",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start_position",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls


write.csv(res_ls$prepared_df,
          row.names = F,
          quote=FALSE,
          file = "PMID30012096/mutations_hg19.csv")

liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID30012096/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID30012096/unconverted.csv")

#########################################
###PMID27058444(hg19)
##########################################

read.csv(
  file = "PMID27058444/all_mut.csv",
  header = T,
  stringsAsFactors = F
)->raw_df


apply(raw_df, 1, 
      function(x){
        gsub(paste(x["Ref.Allele"],"_",sep=""),"",
             x["Tumor.Seq.Allele1_Tumor.Seq.Allele2"])->temp_allele
        if(str_detect(temp_allele,"_")){
          str_split(temp_allele,"_",2)[[1]][1]->temp_allele}
        return(temp_allele)}
)-> raw_df$tumor_allele



within(
  raw_df,{
    sample_id<-paste("PMID27058444",Tumor.Sample.Barcode.,sep=";")
    chr<-paste("chr",Chr,sep = "")
  })->raw_df

subset(
  raw_df,
  chr %in% paste("chr",c(1:22,"X","Y"),sep="")
)->raw_df1


pre_mutation(
  maf_df=raw_df1,
  ref_allele_col="Ref.Allele",
  alr_allele_col="tumor_allele",
  sample_id_col="sample_id",
  chro_id_col="chr",
  start_col="Start.Position",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls

res_ls$exclude_df->raw_dfb
raw_dfb$Start.Position<- raw_dfb$Start.Position +1 


pre_mutation(
  maf_df=raw_dfb,
  ref_allele_col="Ref.Allele",
  alr_allele_col="tumor_allele",
  sample_id_col="sample_id",
  chro_id_col="chr",
  start_col="Start.Position",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls2


write.csv(
  rbind(res_ls$prepared_df,
        res_ls2$prepared_df),
  row.names = F,
  quote=FALSE,
  file = "PMID27058444/mutations_hg19.csv")



liftOver_convert(
  used_chain=chain19to38,
  raw_df=  rbind(res_ls$prepared_df,
                 res_ls2$prepared_df)
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID27058444/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID27058444/unconverted.csv")

##########################
##PMID25839328(hg19)
##############################

read.csv(
  file = "PMID25839328/all_mut.csv",
  header = T,
  stringsAsFactors = F
)->all_mut

within(
    all_mut,{
    sample_id<-paste("PMID25839328",Sample,sep=";")
    ref<-colsplit(Allele.change,">",1:2)[,1]
    alt<-colsplit(Allele.change,">",1:2)[,2]
    Chr=paste("chr",Chrom,sep="")
    }
)->all_mut


pre_mutation(
  maf_df=all_mut,
  ref_allele_col="ref",
  alr_allele_col="alt",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls



write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID25839328/mutations_hg19.csv")



liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID25839328/converted_hg38.csv")


write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID25839328/unconverted.csv")

###########################################
##PMID26873401(hg19)
######################################

read.csv(
  file = "PMID26873401/all_mut.csv",
  header = T,
  stringsAsFactors = F
)->raw_df


within(
  raw_df,{
    paste("PMID26873401",sample.ID,sep = ";")->sample_id
    paste("chr",chromosome,sep="")->Chr }
)->raw_df

head(raw_df)


pre_mutation(
  maf_df=raw_df,
  ref_allele_col="normal.base",
  alr_allele_col="tumor.base",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="coordinate",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls



write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID26873401/mutations_hg19.csv")


liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID26873401/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID26873401/unconverted.csv")

######################################
##PMID32398863(hg19)
###########################################

read.csv(
  file = "PMID32398863/41422_2020_333_MOESM23_ESM.csv",
  header = T,
  stringsAsFactors = F)->x

read.csv("PMID32398863/id_convert.csv",
         header = T,
         stringsAsFactors = F)->id_convert


id_convert$maf_ids<-NA


x %>% filter(Variant_Classification=="IGR") %>%
  group_by(Tumor_Sample_Barcode) %>% tally() ->id_convert2
x %>% filter(Variant_Classification=="Intron") %>%
  group_by(Tumor_Sample_Barcode) %>% tally() ->id_convert3
id_convert3$intergene<-id_convert2$n
id_convert3$Intron<-id_convert3$n 


for(i in 1:nrow(id_convert)){
  
  which(id_convert3$intergene==id_convert[i,"intergenic"] &
        id_convert3$Intron==id_convert[i,"intronic"])->temp_iddd
  
  id_convert2$Tumor_Sample_Barcode[temp_iddd]->id_convert$maf_ids[i]
  
}

which(duplicated(id_convert$id))


id_convert[
  which(id_convert$id=="BDESCC-01-838"),
  "id"]<-c("BDESCC-01-838A",
           "BDESCC-01-838B")

x$sample_id<-NA

for(i in 1:nrow(id_convert)){
  
  x[x$Tumor_Sample_Barcode==id_convert[i,"maf_ids"],
    "sample_id"
    ]<-paste("PMID32398863",id_convert[i,"id"],sep = ";")

}


pre_mutation(
  maf_df=x,
  ref_allele_col="Reference_Allele",
  alr_allele_col="Tumor_Seq_Allele2",
  sample_id_col="sample_id",
  chro_id_col="Chromosome",
  start_col="Start_Position",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls

res_ls$prepared_df->mutations_hg19


write.csv(
  mutations_hg19,
  row.names = F,
  quote=FALSE,
  file = "PMID32398863/mutations_hg19.csv")

liftOver_convert(
  used_chain=chain19to38,
  raw_df=mutations_hg19
)->convert_ls



write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32398863/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32398863/unconverted.csv")
  
################################
##PMID30975989(hg19)
#################################

read.csv(file = "PMID30975989/mut.csv",
         header = T)->mut_df

table(mut_df$X.Trunk.branch.private.)

subset(mut_df,
       X.Trunk.branch.private. %in% c(
         "0_ubiquitous",
         "1_branched"))->mut_df1


mut_df1$sample_id<-paste(
  "PMID30975989",
  gsub("_N","",mut_df1$Normal.sample.name),sep=";")

mut_df1$CHRO<-paste("chr",mut_df1$Chromosome,sep="")

pre_mutation(
  maf_df=mut_df1,
  ref_allele_col="Reference_Allele",
  alr_allele_col="Tumor_Seq_Allele2",
  sample_id_col="sample_id",
  chro_id_col="CHRO",
  start_col="Start_position",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls


write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID30975989/mutations_hg19.csv")


liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID30975989/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID30975989/unconverted.csv")
  
  
###########################################
###PMID34263978(hg19)
######################################

read.csv(file = "PMID34263978/mut_df.csv",
         header = T)->mut_df

mut_df$sample_id<-paste("PMID34263978",
                        mut_df$SampleID,sep=";")


pre_mutation(
  maf_df=mut_df,
  ref_allele_col="Ref",
  alr_allele_col="Alt",
  sample_id_col="sample_id",
  chro_id_col="Chromosome",
  start_col="Coordinate",
  genome=BSgenome.Hsapiens.UCSC.hg19)->res_ls


write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34263978/mutations_hg19.csv")

liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34263978/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34263978/unconverted.csv")

#############################################
####TCGA
############################################

load("tcga_28052061/tcga_escc_maf.rdata")

rbind(
  tcga_escc_maf@data,
  tcga_escc_maf@maf.silent)->tcga_escc_maf_all


subset(
  tcga_escc_maf_all,
  Tumor_Sample_Barcode!="TCGA-V5-A7RC-06A-11D-A403-09"
)->tcga_escc_maf_all

data.frame(tcga_escc_maf_all)->tcga_escc_maf_all

within(
  tcga_escc_maf_all,{
    sample_id<-paste("TCGA",
                     substr(as.character(Tumor_Sample_Barcode),1,12),
                     sep=";")
    Chr<-paste("chr",Chromosome,sep="")}
)->tcga_escc_maf_all_vcf


pre_mutation(
  maf_df=tcga_escc_maf_all_vcf,
  ref_allele_col="Reference_Allele",
  alr_allele_col="Tumor_Seq_Allele2",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start_Position",
  genome=BSgenome.Hsapiens.UCSC.hg19)->res_ls

write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "tcga_28052061/mutations_hg19.csv")

liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "tcga_28052061/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "tcga_28052061/unconverted.csv")

##############################################
##PMID32929369(hg19)
#######################################



read.csv("PMID32929369/thnov10p10606s2.csv",
         header = T,
         stringsAsFactors = F)->x

x[which(x$ref != x$'T'),]->x
x$sample_id<-paste("PMID32929369",x$Sample,sep=";")
x$POS<-x$Pos
x$REF<-x$ref
x$ALT<-x$'T'
x$CHROM<-x$Chr


pre_mutation(
  maf_df=x,
  ref_allele_col="REF",
  alr_allele_col="ALT",
  sample_id_col="sample_id",
  chro_id_col="CHROM",
  start_col="POS",
  genome=BSgenome.Hsapiens.UCSC.hg19)->res_ls


write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32929369/mutations_hg19.csv")


liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32929369/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID32929369/unconverted.csv")
 
##########################################
##PMID31289612(hg19)
############################################


read.csv("PMID31289612/all_mut.csv",
         header = T,
         stringsAsFactors = F)->raw_df

raw_df$sample_id<-paste("PMID31289612",raw_df$ID,sep=";")
raw_df$Chr<-paste("chr",raw_df$Chr.,sep = "")

pre_mutation(
  maf_df=raw_df,
  ref_allele_col="Ref.",
  alr_allele_col="Allele",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Region",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls


write.csv(res_ls$prepared_df,
          row.names = F,
          quote=FALSE,
          file = "PMID31289612/mutations_hg19.csv")


liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID31289612/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID31289612/unconverted.csv")
 
#################################################
##PMID28608921(hg19)
################################################
#41tumor+15LNs,here are all Tumors 

read.csv(
  file = "PMID28608921/all_muts.csv",
  header = T,
  stringsAsFactors = F)->all_snv

within(
  all_snv,
  sample_id<-paste("PMID28608921",sampleID,sep=";")
)->all_snv


pre_mutation(
  maf_df=all_snv,
  ref_allele_col="Ref",
  alr_allele_col="Alt",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls

write.csv(res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID28608921/mutations_hg19.csv")


liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID28608921/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID28608921/unconverted.csv")

#############################################
#PMID28548104(hg19)
########################################

read.csv(
  file = "PMID28548104/all_mut.csv",
  header = T,
  stringsAsFactors = F
)->raw_df


within(
  raw_df,{
    sample_id<-paste("PMID28548104",Sample.ID,sep=";")
    Chr<-paste("chr",Chromosome,sep="")
  })->raw_df


pre_mutation(
  maf_df=raw_df,
  ref_allele_col="Ref.",
  alr_allele_col="Alt.",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls

res_ls$exclude_df->raw_dfb
raw_dfb$Start<- raw_dfb$Start +1 


pre_mutation(
  maf_df=raw_dfb,
  ref_allele_col="Ref.",
  alr_allele_col="Alt.",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="Start",
  genome=BSgenome.Hsapiens.UCSC.hg19
)->res_ls2

rbind(res_ls$prepared_df,
      res_ls2$prepared_df)->all_df
		
		
write.csv(
  all_df, row.names = F,
  quote=FALSE,
  file = "PMID28548104/mutations_hg19.csv")

liftOver_convert(
  used_chain=chain19to38,
  raw_df=  all_df
  )->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID28548104/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID28548104/unconverted.csv")

##########################################
##PMID22877736(hg18)
################################################


read.csv(
  file = "PMID22877736/all_mut.csv",
  header = T,
  stringsAsFactors = F
)->all_df
  
colsplit(all_df$Nucleotide..genomic..,":|>",
         c("chrom",
           "start0",
           "alt")
)->temp2_hg18

gsub(" (Homo)","", fixed = T, temp2_hg18$alt)->temp2_hg18$alt
temp2_hg18$alt[which(temp2_hg18$alt=="")]<-"-"


gsub("^[0-9]+","",temp2_hg18$start0
)->temp2_hg18$ref

gsub("del","",temp2_hg18$ref)->temp2_hg18$ref

str_detect(temp2_hg18$ref,"ins")->insert_index
gsub("ins","",temp2_hg18$ref)[insert_index]->insert_bases

temp2_hg18$ref[insert_index]<-"-"
temp2_hg18$alt[insert_index]<-insert_bases

gsub("[^0-9.-]+","",temp2_hg18$start0)->temp2_hg18$start
temp2_hg18$start<-as.numeric(temp2_hg18$start)

temp2_hg18$Chr <- gsub("g.","",temp2_hg18$chrom)

cbind(all_df,temp2_hg18)->all_df1

within(
  all_df1,{
    sample_id<-paste("PMID22877736",Tumor.......Sample,sep=";")}
)->all_df1

subset(all_df1,
       !Chr %in% c("","chrM","chrMT"))->all_df1

pre_mutation(
  maf_df=all_df1,
  ref_allele_col="ref",
  alr_allele_col="alt",
  sample_id_col="sample_id",
  chro_id_col="Chr",
  start_col="start",
  genome=BSgenome.Hsapiens.UCSC.hg18)->res_ls

write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID22877736/mutations_hg18.csv")

liftOver_convert(
  used_chain=chain18to38,
  raw_df=res_ls$prepared_df)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "PMID22877736/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "PMID22877736/unconverted.csv")

####################################
##ICGC(hg19)
##########################################
#from https://dcc.icgc.org/releases/current/Projects/ESCA-CN


read.table(
  file = "ICGC/simple_somatic_mutation.open.ESCA-CN.tsv",
  sep = "\t",
  header = T,
  stringsAsFactors = F)->ICGC_all


unique(ICGC_all$submitted_sample_id)->all_samples
all_samples[1:88] #PMID24670651
all_samples[229:332] #PMID25839328

all_samples[89:228]->unique_samples

subset(
  ICGC_all, 
  submitted_sample_id %in% unique_samples
)->ICGC_unique

within(
  ICGC_unique,
  {paste(submitted_sample_id,chromosome,
         chromosome_start,chromosome_end,
         reference_genome_allele,
         mutated_to_allele)->mut_id
  })->ICGC_unique

ICGC_unique[!duplicated(ICGC_unique$mut_id),
            ]->ICGC_unique


ICGC_unique$sample_id<-paste("ICGC",ICGC_unique$submitted_sample_id,sep=";")

ICGC_unique$CHRO<-paste("chr",ICGC_unique$chromosome,sep="")

pre_mutation(
  maf_df=ICGC_unique,
  ref_allele_col="reference_genome_allele",
  alr_allele_col="mutated_to_allele",
  sample_id_col="sample_id",
  chro_id_col="CHRO",
  start_col="chromosome_start",
  genome=BSgenome.Hsapiens.UCSC.hg19)->res_ls


write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "ICGC/mutations_hg19.csv")


liftOver_convert(
  used_chain=chain19to38,
  raw_df=res_ls$prepared_df
)->convert_ls

write.csv(
  convert_ls$converted_df,
  row.names = F,
  quote=FALSE,
  file = "ICGC/converted_hg38.csv")

write.csv(
  convert_ls$invaliad_df,
  row.names = F,
  quote=FALSE,
  file = "ICGC/unconverted.csv")

####################################
##PMID34285259(hg38)
######################################


read.csv(file = "PMID34285259/mut.csv",
         header = T
)->mut_df

recode(mut_df$Tumor_Sample_Barcode,
       "E09"="T9",
       "E14"="T14",
       "E15"="T15",
       "E16"="T16",
       "E19"="T19",
       "S1"="T1",
       "S5"="T5",
       "S6"="T6",
       "S7"="T7",
       "S8"="T8")->mut_df$Tumor_Sample_Barcode2

paste(
  "PMID34285259",sep = ";",
  mut_df$Tumor_Sample_Barcode2)->mut_df$sample_id

mut_df$CHROM<-paste("chr",mut_df$Chromosome,sep="")



pre_mutation(
  maf_df=mut_df,
  ref_allele_col="Reference_Allele",
  alr_allele_col="Tumor_Seq_Allele2",
  sample_id_col="sample_id",
  chro_id_col="CHROM",
  start_col="Start_Position",
  genome=BSgenome.Hsapiens.UCSC.hg38)->res_ls

write.csv(
  res_ls$prepared_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34285259/mutations_hg38.csv")

write.csv(
  res_ls$exclude_df,
  row.names = F,
  quote=FALSE,
  file = "PMID34285259/exclude_row.csv")