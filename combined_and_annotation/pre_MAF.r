

### the "ESCC_META_all.hg38_multianno.txt" is from the results of annotation.sh
## the "Otherinfo6" colunm is patient ID 

maftools::annovarToMaf("ESCC_META_all.hg38_multianno.txt", 
             refBuild="hg38",
             tsbCol="Otherinfo6")->maf_all

data.frame(
  Hugo_Symbol=maf_all$Hugo_Symbol,
  Entrez_Gene_Id=".",
  Center=".",
  NCBI_Build="GRCh38",
  Chromosome=maf_all$Chromosome,
  Start_Position=maf_all$Start_Position,
  End_Position=maf_all$End_Position,
  Strand=".",
  Variant_Classification=maf_all$Variant_Classification,
  Variant_Type=maf_all$Variant_Type,
  Reference_Allele=maf_all$Reference_Allele,
  Tumor_Seq_Allele1=maf_all$Tumor_Seq_Allele2,
  Tumor_Seq_Allele2=maf_all$Tumor_Seq_Allele2,
  ExonicFunc.refGene=maf_all$ExonicFunc.refGene,
  Func.refGene=maf_all$Func.refGene,
  Tumor_Sample_Barcode=maf_all$Tumor_Sample_Barcode
)->maf_all_hg38


write.table(maf_all_hg38,
      file = "maf_all_hg38.maf", 
      sep = "\t", quote = FALSE, 
      col.names = T,
      row.names = FALSE)
