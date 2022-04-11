
#Biologist -Shreenk.Project BF528: Applications in Translational Bioinformatics, Project 2
install.packages("tidyr")

# load the tidyr package
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Set the working directory
#setwd("~/Documents/project-2-dachsund/Biologist.R")

# Part 1: Compare trends of genes in FPKM table with those of Figure 1D in paper using sample data for analysis


# READ FILES
fpkm_matrix <- read.table('~/Documents/project-2-dachsund/fpkm_matrix.csv', header = TRUE)


# Subset the matrix
# Get averages for each AD,P4,P7
fpkm_matrix$Ad <- (fpkm_matrix$Ad_1_FPKM + fpkm_matrix$Ad_2_FPKM)/2
fpkm_matrix$P4 <- (fpkm_matrix$P4_1_FPKM + fpkm_matrix$P4_2_FPKM)/2
fpkm_matrix$P7 <- (fpkm_matrix$P7_1_FPKM + fpkm_matrix$P7_2_FPKM)/2

# AVERAGES WITHOUT NAME COL
fpkm_matrix <- fpkm_matrix[ , !(names(fpkm_matrix) %in% c( "Ad_1_FPKM","Ad_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM"))]

# Mitochondria data
genes <- list('ENSMUSG00000023861'='Mpc1','ENSMUSG00000024997'='Prdx3', 'ENSMUSG00000032047' ='Acat1',  'ENSMUSG00000025465'='Echs1',  'ENSMUSG00000014606'='Slc25a11',  
              'ENSMUSG00000026664'='Phyh')

#subset only those col whr tracking id is in genes
mito <- subset(fpkm_matrix, fpkm_matrix$tracking_id %in% names(genes))

for(id in names(genes)){
  for(row in 1:6){
    if(mito$tracking_id[row] == id){
      mito$gene_name[row] <- genes[id]
    }
  }
}
# plotting graph
#P02 is taken as P0
names(mito)[2] <- "P0" 


mito <- mito[, 2:6]


clean <- mito %>% gather(time, expression, c("P0", "P4", "P7", "Ad"), -gene_name)

df2 <- as.data.frame(lapply(clean, unlist))


df2 <- df2[order(df2[, 1]), ]

df2$time <- factor(df2$time, levels=c("P0", "P4", "P7", "Ad"))


a <- ggplot(df2, aes(x=time, y=expression, group=gene_name)) + geom_point(aes(color=gene_name)) + 
  geom_line(aes(color=gene_name)) + labs(title="Mitochondria",x="", y = "FPKM") +
  theme(legend.title = element_blank())


# Sarcomere data
genes <- list("ENSMUSG00000028273"='Pdlim5','ENSMUSG00000032648'='Pygm', 'ENSMUSG00000028116' ='Myoz2',  'ENSMUSG00000026208'='Des',  'ENSMUSG00000030470'='Csrp3',  
              'ENSMUSG00000007877'='Tcap', 'ENSMUSG00000032060'='Cryab')
sarc <- subset(fpkm_matrix, fpkm_matrix$tracking_id %in% names(genes))
sarc[,'gene_name'] <- NA                


for(id in names(genes)){
  for(row in 1:7){
    if(sarc$tracking_id[row] == id){
      sarc$gene_name[row] <- genes[id]
    }
  }
}
# plot for Sacromere
names(sarc)[2] <- "P0"#only P0_2 is present 
sarc <- sarc[, 2:6]
sarc
clean <- sarc %>% gather(time, expression, c("P0", "P4", "P7", "Ad"), -gene_name)
df2 <- as.data.frame(lapply(clean, unlist))
df2 <- df2[order(df2[, 1]), ]
df2$time <- factor(df2$time, levels=c("P0", "P4", "P7", "Ad"))

b <- ggplot(df2, aes(x=time, y=expression, group=gene_name)) + geom_point(aes(color=gene_name)) + 
  geom_line(aes(color=gene_name)) + labs(title="Sarcomere",x="", y = "FPKM") +
  theme(legend.title = element_blank())



# Cell cycle
genes <- list("ENSMUSG00000029283"='Cdc7','ENSMUSG00000046179'='E2f8', 'ENSMUSG00000069089' ='Cdk7',  'ENSMUSG00000066149'='Cdc26',  'ENSMUSG00000017499'='Cdc6',  
              'ENSMUSG00000027490'='E2f1', 'ENSMUSG00000020687'='Cdc27', 'ENSMUSG00000000028' = 'Cdc45', 'ENSMUSG00000027323' = 'Rad51', 'ENSMUSG00000020897'='Aurkb',
              'ENSMUSG00000024370'='Cdc23', 'ENSMUSG00000022070'='Bora')
cell <- subset(fpkm_matrix, fpkm_matrix$tracking_id %in% names(genes)) 
cell[,'gene_name'] <- NA
for(id in names(genes)){
  for(row in 1:12){
    if(cell$tracking_id[row] == id){
      cell$gene_name[row] <- genes[id]
    }
  }
}
# plot for cell cycle data
names(cell)[2] <- "P0"
cell <- cell[, 2:6]
clean <- cell %>% gather(time, expression, c("P0", "P4", "P7", "Ad"), -gene_name)
df2 <- as.data.frame(lapply(clean, unlist))
df2 <- df2[order(df2[, 1]), ]
df2$time <- factor(df2$time, levels=c("P0", "P4", "P7", "Ad"))

c <- ggplot(df2, aes(x=time, y=expression, group=gene_name)) + geom_point(aes(color=gene_name)) + 
  geom_line(aes(color=gene_name)) + labs(title="Cell cycle",x="", y = "FPKM") +
  theme(legend.title = element_blank())

#Arrange in a single row
ggarrange(b, a, c, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)


# PART 3 HEATMAP-
ad_1 <- read.table('~/Documents/project-2-dachsund/AD_1.fpkm_tracking', header = TRUE)
ad_2 <- read.table('~/Documents/project-2-dachsund/AD_2.fpkm_tracking', header = TRUE)
p0_2 <- read.table('~/Documents/project-2-dachsund/P0_2.fpkm_tracking', header = TRUE)
p4_1 <- read.table('~/Documents/project-2-dachsund/P4_1.fpkm_tracking', header = TRUE)
p4_2 <- read.table('~/Documents/project-2-dachsund/P4_2.fpkm_tracking', header = TRUE)
p7_1 <- read.table('~/Documents/project-2-dachsund/P7_1.fpkm_tracking', header = TRUE)
p7_2 <- read.table('~/Documents/project-2-dachsund/P7_2.fpkm_tracking', header = TRUE)
diff <- read.table('~/Documents/project-2-dachsund/gene_exp.diff', header=TRUE)


# Concatenate the FPKM columns of each of the tracking datasets into a single dataframe
track_id <- as.data.frame(p0_2$tracking_id)

gene <- as.data.frame(p0_2$gene_short_name)
full <- cbind(track_id, gene)

colnames(full) <- c('tracking_id', 'gene_name')


add_cols <- function(df, big, name){
  new <- list()
  for(row in 1:37468){
    same <- match(big$tracking_id[row], df$tracking_id)
    new <- c(new, df$FPKM[same])
  }
  new <- t(as.data.frame(new))
  final <- cbind(big, new)
  colnames(final)[-1] <- name
  return(final)
}
#adding cols 
full <- add_cols(p0_2, full, 'P0')
full <- add_cols(p4_1, full, 'P4_1')
full <- add_cols(p4_2, full, 'P4_2')
full <- add_cols(p7_1, full, 'P7_1')
full <- add_cols(p7_2, full, 'P7_2')
full <- add_cols(ad_1, full, 'Ad_1')
full <- add_cols(ad_2, full, 'Ad_2')
colnames(full) <- c('tracking_id', 'gene_name', 'P0_2', "P4_1", 'P4_2', 'P7_1', 'P7_2', "Ad_1", "Ad_2")

# Get top 500 differentially expressed genes between P4 and P7 (FROM SAMPLE MATRIX DATA CSV )

diff <- subset(diff, (diff$log2.fold_change. != "Inf" & diff$log2.fold_change. != "-Inf"))

genes <- top_n(diff, 500, abs(log2.fold_change.)) # get top 1000 genes of highest absolute log fold change
gene_names <- genes$gene

clean <- subset(full, full$gene_name %in% gene_names) #for matching genes

clean <- as.matrix(clean)

really_clean <- clean[, 3:9]

really_clean <- really_clean[apply(really_clean[,-1], 1, function(x) !all(x<0.0000001)),]

mode(really_clean) <- 'numeric'
rownames(really_clean) <- c()
heatmap(really_clean, scale='row')


###7.2 Check overlap of biological process in cluster and append in a col-
processes=list("intrinsic to membrane","golgi apparatus","plasma membrane","intracellular signaling","protein phosphorylation","enzyme linked signaling","guanyl nucleotide binding","endoplasmic reticulum","mitochondrion","oxidation reduction","sacromere","myofibril","contractile fiber","t band","z-disk","t-tubule","nucleolus","ribosome biogenesis","trna processing","mrna processing","spliceosome","mrna transport","exosome","ribosome","proteasome complex","rna pol ||","chromatin","krueppel associated box")

UPREG <- read.csv(file = '~/Documents/project-2-dachsund/UPREG_final.csv',header=TRUE)

H<-as.data.frame(UPREG)

DWNREG <- read.csv(file = '~/Documents/project-2-dachsund/DWNREG_final.csv',header=TRUE)

#taking subset of enrichment col to find out the gene annotation 
a<-substr(UPREG$Term,12,80) 

b<-substr(UPREG$Term,12,80)

#find the biological process common bw paper and those in cluster files i worked on 
#for upreg csv file final
Process_common<-intersect(processes,a)
listofprocesscommon<-list(Process_common)

#added a col for processes common upreg
H$Added_Column <- listofprocesscommon
H

#for downreg csv file final
G<-as.data.frame(DWNREG)
Process_common<-intersect(processes,b)

#added a col for processes common dwnreg
G$Added_Column <- listofprocesscommon
G







