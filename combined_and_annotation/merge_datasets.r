
read.csv("PMID32398863/converted_hg38.csv")->pmid32398863
read.csv("ICGC/converted_hg38.csv")->ICGC
read.csv("PMID22877736/converted_hg38.csv")->pmid22877736
read.csv("PMID25151357/converted_hg38.csv")->pmid25151357
read.csv("PMID25839328/converted_hg38.csv")->pmid25839328
read.csv("PMID26873401/converted_hg38.csv")->pmid26873401
read.csv("PMID27058444/converted_hg38.csv")->pmid27058444
read.csv("PMID28548104/converted_hg38.csv")->pmid28548104
read.csv("PMID28608921/converted_hg38.csv")->pmid28608921
read.csv("PMID30012096/converted_hg38.csv")->pmid30012096
read.csv("PMID31289612/converted_hg38.csv")->pmid31289612
read.csv("PMID32929369/converted_hg38.csv")->pmid32929369
read.csv("PMID32974170/converted_hg38.csv")->pmid32974170
read.csv("tcga_28052061/converted_hg38.csv")->TCGA_escc
read.csv("PMID34263978/converted_hg38.csv")->pmid34263978
read.csv("PMID34285259/mutations_hg38.csv")->pmid34285259
read.csv("PMID34413060/converted_hg38.csv")->pmid34413060
read.csv("PMID30975989/converted_hg38.csv")->pmid30975989


rbind(
 pmid32398863,
 ICGC,
 pmid22877736,
 pmid25151357,
 pmid25839328,
 pmid26873401,
 pmid27058444,
 pmid28548104,
 pmid28608921,
 pmid30012096,
 pmid31289612,
 pmid32929369,
 pmid32974170,
 TCGA_escc,
 pmid34263978,
 pmid34285259,
 pmid34413060,
 pmid30975989
)->all_articles

read.csv("vcf/SRP034680_WES/mutation_hg38.csv")->SRP034680_WES
read.csv("vcf/SRP034680_WGS/mutation_hg38.csv")->SRP034680_WGS
read.csv("vcf/SRP033394/mutation_hg38.csv")->SRP033394
read.csv("vcf/SRP072858_WES/mutation_hg38.csv")->SRP072858_WES
read.csv("vcf/SRP072858_WGS/mutation_hg38.csv")->SRP072858_WGS
read.csv("vcf/SRP116657/mutation_hg38.csv")->SRP116657
read.csv("vcf/SRP127593/mutation_hg38.csv")->SRP127593
read.csv("vcf/SRP099292/Multi_tumor/mutation_hg38.csv")->SRP099292_M
read.csv("vcf/SRP099292/Single_tumor/mutation_hg38.csv")->SRP099292_S
read.csv("vcf/SRP072112/mutation_hg38.csv")->SRP072112
read.csv("vcf/SRP327447/mutation_hg38.csv")->SRP327447
read.csv("vcf/SRP150544/mutation_hg38.csv")->SRP150544
read.csv("vcf/SRP179388/mutation_hg38.csv")->SRP179388
read.csv("vcf/SRP059537/mutation_hg38.csv")->SRP059537
read.csv("vcf/ECRT/mutation_hg38.csv")->ECRT


rbind(
 SRP034680_WES,
 SRP034680_WGS,
 SRP033394,
 SRP072858_WES,
 SRP072858_WGS,
 SRP116657,
 SRP127593,
 SRP099292_M,
 SRP099292_S,
 SRP072112,
 SRP327447,
 SRP150544,
 SRP179388,
 SRP059537
)->all_SRAs


rbind(
  all_articles,
  all_SRAs,
  ECRT
)->ESCC_META_all


write.csv(
  all_articles,
  row.names = F,
  quote=FALSE,
  file = "all_articles.csv")


write.csv(
  all_SRAs,
  row.names = F,
  quote=FALSE,
  file = "all_SRAs.csv")

write.csv(
  ESCC_META_all,
  row.names = F,
  quote=FALSE,
  file = "ESCC_META_all.csv")

data.frame(
  CHROM=ESCC_META_all$CHROM,
  POS=ESCC_META_all$POS,
  ID=ESCC_META_all$sample_id,
  REF=ESCC_META_all$REF,
  ALT=ESCC_META_all$ALT,
  QUAL=".",
  FILTER=".",
  INFO=".")->ESCC_META_all_vcf

colnames(ESCC_META_all_vcf)[1]<-"#CHROM"


write.table(
  ESCC_META_all_vcf,
  file="ESCC_META_all.vcf",
  #fileEncoding = "ascii",
  eol="\n",
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F)
