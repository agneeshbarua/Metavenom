
library(tidyverse)

#add taxa ids
h = "9606"
m = "10090"
s = "103944"
x = "8364"
g = "9031"
a = "28377"
p = "9598"
clf = "9615" 

dat<-read.table("./gene_orthologs.txt")

#filter datasets; gene in snake and human
h_s<-dat %>% filter(V1==h) %>% filter(V4 == c(s))

#filter datasets; gene in mouse and human
h_m<-dat %>% filter(V1==h) %>% filter(V4 == c(m))

#filter datasets; gene in xenopus and human
h_x<-dat %>% filter(V1==h) %>% filter(V4 == c(x))

#filter datasets; gene in gallus and human
h_g<-dat %>% filter(V1==h) %>% filter(V4 == c(g))

#filter datasets; gene in anolis and human
h_a<-dat %>% filter(V1==h) %>% filter(V4 == c(a))

#filter datasets; gene in pan and human
h_p<-dat %>% filter(V1==h) %>% filter(V4 == c(p))

#filter datasets; gene dog and human
h_clf<-dat %>% filter(V1==h) %>% filter(V4 == c(clf))


#common between mouse and snake
o_ms<-h_m %>% filter(V2 %in% h_s$V2)
#common between snake and human; should be same as above
o_sh<-h_s %>% filter(V2 %in% h_m$V2)
#common between xeno and snake
o_xs<-h_x %>% filter(V2 %in% h_s$V2)
#common between gallus and snake
o_gs<-h_g %>% filter(V2 %in% h_s$V2)
#common between anole and snake
o_as<-h_a %>% filter(V2 %in% h_s$V2)
#common between pan and snake
o_ps<-h_p %>% filter(V2 %in% h_s$V2)
#common between god and snake
o_clfs<-h_clf %>% filter(V2 %in% h_s$V2)


h_m_s_data<-o1 %>% mutate(h_genes = o1$V2) %>% mutate(m_genes = o1$V5) %>% 
  mutate(s_genes = o2$V5) %>% select(h_genes,m_genes,s_genes) #make dataframe with hms_orthologs

#start filtering the orthologs from lowest to the highest
gal_filter<- h_m_s_data %>% filter(h_genes%in% o_gs$V2) #filter the hms_orthologs in gallus data
gal_filter_2<- o_gs %>% filter(V2 %in% gal_filter$h_genes)#filter hms_orths in gs orths(should be same as above)
h_m_s_g_data<- gal_filter %>% mutate(g_gene = gal_filter_2$V5)

#Now use second lowest orthoset
anole_filter<-h_m_s_g_data %>% filter(h_genes %in% o_as$V2)
anole_filter_2<- o_as %>% filter(V2 %in% anole_filter$h_genes)#filter hms_orths in gs orths(should be same as above)
h_m_s_g_a_data<- anole_filter %>% mutate(a_gene = anole_filter_2$V5)

#Next
xeno_filter<-h_m_s_g_a_data %>% filter(h_genes %in% o_xs$V2)
xeno_filter_2<- o_xs %>% filter(V2 %in% xeno_filter$h_genes)#filter hms_orths in gs orths(should be same as above)
h_m_s_g_a_x_data<- xeno_filter %>% mutate(x_gene = xeno_filter_2$V5)

#Next
pan_filter<-h_m_s_g_a_x_data %>% filter(h_genes %in% o_ps$V2)
pan_filter_2<- o_ps %>% filter(V2 %in% pan_filter$h_genes)#filter hms_orths in gs orths(should be same as above)
h_m_s_g_a_x_p_data<- pan_filter %>% mutate(p_gene = pan_filter_2$V5)

#Next
dog_filter<-h_m_s_g_a_x_p_data %>% filter(h_genes %in% o_clfs$V2)
dog_filter_2<- o_clfs %>% filter(V2 %in% dog_filter$h_genes)#filter hms_orths in gs orths(should be same as above)
h_m_s_g_a_x_p_clf_data<- dog_filter %>% mutate(dog_gene = dog_filter_2$V5)

h_m_s_g_a_x_p_clf_data %>% write.csv("./Hum_mus_pm_gal_anole_xen_pan_canis_orthologs.csv")
######################################################################################################