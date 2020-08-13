setwd('/Users/kellenxu/Documents/Rutgers/Lab_Work/Project_Fusion_gene'); getwd(); list.files()

df <- readRDS("Table1.rds") # Restore the object
class(df)
library(stringr)
names(df)<-str_replace_all(names(df), c(" " = "_" , "," = "" )) #remove space in column names
names(df)

library(dplyr)
library(tidyr)
gene1 <- df %>% 
  group_by(gene1_RE, MSIstatus) %>% 
  tally() %>%
  pivot_wider(names_from = MSIstatus, values_from = n)

gene1[is.na(gene1)] <- 0               #处理 'NA'
g1_cross <- data.frame(gene1[, -4])



N <- sum(g1_cross$MSI.H) #92
M <- sum(g1_cross$MSI.H) + sum(g1_cross$MSS) #243
k <- sum(g1_cross$MSI.H)
kk <- sum(g1_cross$MSI.H) + sum(g1_cross$MSS)

for (i in g1_cross$gene1_RE) {
  q <- filter(g1_cross, gene1_RE == i)[['MSI.H']] 
  m <- rowSums(filter(g1_cross, gene1_RE == i)[,c('MSI.H', 'MSS')])
  n <- kk - m
  p_val <- round(phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE), digits = 6)
  print(c(i, p_val))
}



MSI_depletion <- function(i,table,K,KK){
  q <- filter(table, gene1_RE == i)[['MSI.H']] 
  m <- rowSums(filter(table, gene1_RE == i)[,c('MSI.H', 'MSS')])
  n <- KK - m
  p_val <- phyper(q, m, n, K, lower.tail = TRUE, log.p = FALSE)
  return(p_val)
}
MSI_depletion("ALK",g1_cross,k,kk)
  
MSI_enrichment <- function(i,table,K,KK){
  q <- filter(table, gene1_RE == i)[['MSI.H']] 
  m <- rowSums(filter(table, gene1_RE == i)[,c('MSI.H', 'MSS')])
  n <- KK - m
  p_val <- phyper(q-1, m, n, K, lower.tail = FALSE, log.p = FALSE)
  return(p_val)
}
MSI_enrichment("ALK",g1_cross,k,kk)


print(c("Gene-1", "p_depletion", "p_enrichment"));
for (i in g1_cross$gene1_RE) {
  p1 <- round(MSI_depletion(i,g1_cross,k,kk), digits = 6)
  p2 <- round(MSI_enrichment(i,g1_cross,k,kk), digits = 6)
  print(sprintf("%-5s| p_deplet %f | p_enrich %f", i, p1, p2))
}


setwd('/Users/kellenxu/Documents/Rutgers/Lab_Work/Project_Fusion_gene'); getwd(); list.files()
js1_bi <- read.csv(file = 'js1_bi.csv')



library(stringr)
df <- readRDS("Table1.rds") # Restore the object
names(df)<-str_replace_all(names(df), c(" " = "_" , "," = "" )) #remove space in column names
A <- nrow(subset(df, MSIstatus=='MSI-H' & Junction_Sequence_2!='-'))
B <- nrow(subset(df, MSIstatus=='MSS' & Junction_Sequence_2!='-'))
C <- nrow(subset(df, MSIstatus=='MSI-H' & Junction_Sequence_2=='-'))
D <- nrow(subset(df, MSIstatus=='MSS' & Junction_Sequence_2=='-'))
testor<- matrix(c(A,B,C,D), nrow = 2, byrow = TRUE)
fisher.test(testor)



# draw box plot
library(tidyr)
library(dplyr)
BPDistance <- read.csv("BreakPointDistance.csv")
tmp <- BPDistance %>% 
  select("deidentifiedSpecimenName","MSIstatus","G1_ΔBP","G2_ΔBP") %>% 
  gather("Junction","distance", 3:4) %>%
  mutate_if(is.character,as.factor)
str(tmp)
library(ggplot2)
ggplot(tmp, aes(x=Junction, y=distance, fill=MSIstatus)) + 
  geom_boxplot() + 
  ylim(0,30)


