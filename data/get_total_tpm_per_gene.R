library(tidyverse)

#Add names of gene from ENSEMBL
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "cfamiliaris_gene_ensembl",host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
t2g

#load kalisto abundance file 
sample_id <- c("SRR6206911","SRR6206901","SRR6206906","SRR6206921")
kal_dirs<-file.path("../R Training/Meta-venom/Comparative_gene_exp/Kallisto/Canis/HKLS",sample_id)
kal_dirs



get_total_tpm_per_gene<-function(file_path) {
  
  #file = "path_to_kallisto_abundance.tsv"
  dat<-read.delim(file_path) #118489 genes
  
  #remove the number after the decimal in ENSMUST00000xxxx.1
  dat<-dat %>% separate(target_id, c("target_id"),extra = 'drop')
  
  
  #filtered t2g with what is present in kallisto
  t2g_2<-t2g %>% filter(target_id %in% dat$target_id) #118395 genes; missing 95 transcripts from kallisto, likely couldn't be matched with ens gene_id
  
  #filtered kallisto data with what is present in t2g; same number as above.
  dat2<-dat %>% filter(target_id %in% t2g$target_id) #118395 genes
  
  #arrange both the data in same order
  dat2_desc_ar<-dat2 %>% arrange(desc(target_id))
  t2g_2_desc_ar<-t2g_2 %>% arrange(desc(target_id))
  
  #combine data to for data frame with gene_id, and gene_symb
  f_dat<-dat2_desc_ar %>% mutate(ens_gene = t2g_2_desc_ar$ens_gene, gene_sym = t2g_2_desc_ar$ext_gene) %>% arrange(desc(ens_gene))
  #summarise to obtain total tpm per gene
  file_name<-regmatches(file_path,regexpr("SRR\\d*",file_path))
  f_dat %>% group_by(ens_gene) %>% summarise(tpm = sum(tpm),n = n()) %>%
    write.csv(file.path("../R Training/Meta-venom/Comparative_gene_exp/Kallisto/Canis/Res",file_name)) #the counts helps confirm that all the transcripts are included
  
}

#enter the following
sapply(file.path(kal_dirs,"abundance.tsv"), get_total_tpm_per_gene)



concat<-function(file) {
  inpt<-read.csv(file)
  try %>% mutate(tpm2 = inpt$tpm) %>% view
}

concat("../../../../Comparative_gene_exp/Kallisto/Sleuth/Res/Anole/")
#sapply(file.path("./Res",s2c$sample),concat)
#############################################################################################################################################################
dir<-"../R Training/Meta-venom/Comparative_gene_exp/Kallisto/Xeno/Res"

#Xeno_sal_1<-read.csv(file.path(dir,"SRR10246417"))
#Xeno_sal_2<-read.csv(file.path(dir,"SRR10246418"))
#Xeno_sal_3<-read.csv(file.path(dir,"SRR8133136"))
#Xeno_sal_4<-read.csv(file.path(dir,"SRR8133137"))

Xeno_liv_1<-read.csv(file.path(dir,"SRR579561"))
Xeno_liv_2<-read.csv(file.path(dir,"SRR5412274"))
Xeno_liv_3<-read.csv(file.path(dir,"SRR5412276"))
#Mus_liv_4<-read.csv(file.path(dir,"SRR453155"))

Xeno_kid_1<-read.csv(file.path(dir,"SRR5412269"))
Xeno_kid_2<-read.csv(file.path(dir,"SRR5412270"))
Xeno_kid_3<-read.csv(file.path(dir,"SRR5412271"))

Xeno_hr_1<-read.csv(file.path(dir,"SRR579562"))
Xeno_hr_2<-read.csv(file.path(dir,"SRR579563"))
Xeno_hr_3<-read.csv(file.path(dir,"SRR5412268"))
#Xeno_hr_4<-read.csv(file.path(dir,"SRR5412155"))


Xeno_hr_1$X<-NULL
Xeno_hr_2$n<-NULL

Xeno_tissues_tpm<- Xeno_hr_1 %>% 
  mutate(Xeno_hrt_2 = Xeno_hr_2$tpm) %>%
  mutate(Xeno_hrt_3 = Xeno_hr_3$tpm) %>%
  mutate(Xeno_kidn_1 = Xeno_kid_1$tpm) %>% 
  mutate(Xeno_kidn_2 = Xeno_kid_2$tpm) %>% 
  mutate(Xeno_kidn_3 = Xeno_kid_3$tpm) %>% 
  mutate(Xeno_livr_1 = Xeno_liv_1$tpm) %>% 
  mutate(Xeno_livr_2 = Xeno_liv_2$tpm) %>% 
  mutate(Xeno_livr_3 = Xeno_liv_3$tpm) %>%
   write.csv(file.path(dir,"Xeno_New_tissues_tpm.csv"))

                                                       