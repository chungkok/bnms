## load required libraries
library(hacksig)
library(tidyverse)
library(msigdbr)

## get hallmark geneset list
hallmark_list <- msigdbr(species = "Homo sapiens", category = "H") %>%
  distinct(gs_name, gene_symbol) %>%
  nest(genes = c(gene_symbol)) %>%
  mutate(genes = map(genes, compose(as_vector, unname))) %>%
  deframe()


## get reactome geneset list
reactome_list <- msigdbr(species = "Homo sapiens", subcategory = "REACTOME") %>%
  distinct(gs_name, gene_symbol) %>%
  nest(genes = c(gene_symbol)) %>%
  mutate(genes = map(genes, compose(as_vector, unname))) %>%
  deframe()


## load example data
dat <- read_tsv(file="./data/example_GEP_myeloma_FPKM.txt")

## average duplicate genes If there is)
dat1 <- dat %>%
  group_by(gene_name) %>%
  summarise_all(.funs = mean)


dat2 <- dat1 %>% as.data.frame()
rownames(dat2) <- dat2$gene_name
vdat3 <- dat2[, -1]


## Z score for hallmark
res.z.hallmark <- hack_sig(vdat3, signatures = hallmark_list, method = "zscore")

## Z score for reactome
res.z.reactome <- hack_sig(vdat3, signatures = reactome_list, method = "zscore")


## combined reactome and hallmark pathway scores
vtotal <- cbind(res.z.hallmark, res.z.reactome[, -1])
rownames(vtotal) <- vtotal$sample_id

## restricted to our annotated pathways only
annot1 <- read.csv(file="./database/filtered pathways list.csv")
vtotal1 <- vtotal[, c(1, which(colnames(vtotal) %in% annot1$Pathways))]

head(vtotal1)

