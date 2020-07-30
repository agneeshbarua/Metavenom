setwd("./Data_sets/")

library(dplyr)

#### Making the mammals ortholog and tpm dataset ####
#get each taxa tpm data
Canis<-read.csv("Canis_geneid_tpm.csv")
Homo<-read.csv("Homo_geneid_tpm.csv")
Mus<-read.csv("Mus_geneid.tpm.csv")
Pm_naja<-read.csv("./Pm_naja_alltissueExpr.csv")

#Get the ortholog list obntained from NCBI eukaryotic genome annotation pipeline
orths<-read.csv("Hum_mus_pm_gal_anole_xen_pan_canis_orthologs.csv")
orths<-orths %>% select(h_genes,m_genes,s_genes,dog_gene)

#Perform filtering step
#Different datasets will have different orthlogs expressed. This filtering step will ensure only the orthologs 
#that have expression data in all taxa will make up the dataset.
#Always start with the largest dataset and descend down.
Mus<-Mus %>% filter(gene_id %in% orths$m_genes)
orths<-orths %>% filter(m_genes %in% Mus$gene_id)

Homo<-Homo %>% filter(gene_ids %in% orths$h_genes)
orths<-orths %>% filter(h_genes %in% Homo$gene_ids)

Canis<-Canis %>% filter(gene_ids %in% orths$dog_gene)
orths<-orths %>% filter(dog_gene %in% Canis$gene_ids)

Pm_naja<-Pm_naja %>% filter(gene %in% orths$s_genes)
orths<-orths %>% filter(s_genes %in% Pm_naja$gene)
orths<-orths %>% arrange(s_genes)

#Final filtering to get one-to-one ortholog tpms
Canis<-Canis %>% dplyr::slice(match(orths$dog_gene,Canis$gene_id)) %>% select(-c(X,ens_gene))
Homo<-Homo %>% dplyr::slice(match(orths$h_genes,Homo$gene_ids)) %>% select(-c(X,ens_gene))
Pm_naja<-Pm_naja %>%dplyr::slice(match(orths$s_genes,Pm_naja$gene)) %>% select(-c(Pm_gene_id,Nn_target_id))
Mus<-Mus %>% dplyr::slice(match(orths$m_genes,Mus$gene_id)) %>% select(-c(X,ens_gene)) 

#make the final data frame
#remove gene_id column after checking
Canis$gene_ids<-NULL
Homo$gene_ids<-NULL
Mus$X.1<-NULL
Mus$gene_id<-NULL

# making salivary gland data
Pm_naja_sg<-Pm_naja[,2:34] %>% select(-c(class))
Mus_sg<-Mus[,1:4]
Homo_sg<-Homo[,1:2]
Canis_sg<-Canis[,1:2]
Mammals_sg<-bind_cols(Pm_naja_sg,Mus_sg,Homo_sg,Canis_sg)

#save final dataset
Mammals_sg %>% write.csv("./Mammals_sg_tpm.csv")

# Making non salivary gland data
Pm_naja_liv<-Pm_naja[,c(2:3,5:34)]
Mus_liv<-Mus[,11:14]
Homo_liv<-Homo[,9:11]
Canis_liv<-Canis[,9:11]

Mammals_liv<-bind_cols(Pm_naja_liv,Mus_liv,Homo_liv,Canis_liv)
Mammals_liv %>% write.csv("./Mammals_liv_Pm_vg_tpm.csv")

Pm_naja_hrt<-Pm_naja[,c(2:3,5:34)]
Mus_hrt<-Mus[,5:7]
Homo_hrt<-Homo[,3:5]
Canis_hrt<-Canis[,3:5]

Mammals_hrt<-bind_cols(Pm_naja_hrt,Mus_hrt,Homo_hrt,Canis_hrt)
Mammals_hrt %>% write.csv("./Mammals_hrt_tpm_Pm_vg_tpm.csv")

Pm_naja_kid<-Pm_naja[,c(2:3,5:34)]
Mus_kid<-Mus[,8:10]
Homo_kid<-Homo[,6:8]
Canis_kid<-Canis[,6:8]

Mammals_kid<-bind_cols(Pm_naja_kid,Mus_kid,Homo_kid,Canis_kid)
Mammals_kid %>% write.csv("./Mammals_kid_Pm_vg_tpm.csv")

Pm_naja %>% write.csv("Pm_naja_mammals_orth_notwritteninfile.csv")
