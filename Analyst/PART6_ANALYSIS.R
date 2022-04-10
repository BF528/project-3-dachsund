#PART6 (ANALYST -SHREENK)
# In this part we take differential expression results from both DESeq2 and limma and we measure and examine concordance between platforms.


setwd("~/Documents/project-3-dachsund/ANALYST.R/Part5_ANALYSIS.R/PART6_ANALYSIS.R")

library(dplyr)
library(ggplot2)
library(ggpubr) 

# Read the data files from part 5 limma results:-
methyl <- read.csv('~/Documents/project-3-dachsund/new_3methylcholanthrene_limma_results.csv',as.is=TRUE)

fluco <- read.csv('~/Documents/project-3-dachsund/new_fluconazole_limma_results.csv', as.is = TRUE)

piri <- read.csv('~/Documents/project-3-dachsund/new_pirinixic_acid_limma_results.csv', as.is=TRUE) 

# Creating map data frame from the reference Affy map-
map <-read.csv('~/Documents/project-3-dachsund/refseq_affy_map.csv', as.is=TRUE)

# Get the probe id from the limma results for each chemical :-
methyl_probe <- methyl$X
fluco_probe <- fluco$X
piri_probe <- piri$X

# Find probe id in affy map, add this to the row in limma results to column "REFSEQID"
#FIRST ANALYSIS(3-methylcholanthrene)
for (row in 1:length(methyl_probe)){
  refseq <- list(map$REFSEQ[which(map$PROBEID == methyl_probe[row])])
  methyl$REFSEQID[row] <- refseq
}
#SECOND ANALYSIS (Fluconazole)
for (row in 1:length(fluco_probe)){
  refseq <- list(map$REFSEQ[which(map$PROBEID == fluco_probe[row])])
  fluco$REFSEQID[row] <- refseq
}
# THIRD ANALYSIS (Pirinixic acid)
for (row in 1:length(piri_probe)){
  refseq <- list(map$REFSEQ[which(map$PROBEID == piri_probe[row])])
  piri$REFSEQID[row] <- refseq
}

# ExpandING rows with multiple refseqid and delete rows with no refseqids
#FIRST ANALYSIS( 3-methylcholanthrene)
maptomethyl <- methyl[FALSE, names(methyl) %in% c("X", "REFSEQID", "logFC", "AveExpr")] 
as.data.frame(maptomethyl)
#SECOND ANALYSIS (fluconazole)
maptofluco <- fluco[FALSE, names(fluco) %in% c("X", "REFSEQID", "logFC", "AveExpr")]  
as.data.frame(maptofluco)
#THIRD ANALYSIS (pirinixic acid)
maptopiri <- piri[FALSE, names(piri) %in% c("X", "REFSEQID", "logFC", "AveExpr")] 
as.data.frame(maptopiri)

# FiLL THE EMPTY DATAFRAME-
# 3-methylcholanthrene
for (row in 1:nrow(methyl)) {
  if (methyl[row, "P.Value"] < 0.05) {
    if(abs(methyl[row, "logFC"]) > 1.5) {
      methyl_probe<- methyl[row,"X"] # Get probe id
      refseq <- methyl[row, "REFSEQID"]
      logfc <- methyl[row, "logFC"]
      aveexpr <- methyl[row, "AveExpr"]
      if (length(refseq[[1]])>0) { 
        for (id in refseq[[1]]) { 
          id_row <- data.frame(methyl_probe, id, logfc, aveexpr) 
          names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") 
          maptomethyl <- rbind(maptomethyl, id_row) 
        }
      }
    }
  }
}
# fluconazole
for (row in 1:nrow(fluco)) {
  if (fluco[row, "P.Value"] < 0.05) {
    if(abs(fluco[row, "logFC"]) > 1.5) {
      fluco_probe<- fluco[row,"X"] 
      refseq <- fluco[row, "REFSEQID"] # Get refseqid
      logfc <- fluco[row, "logFC"]
      aveexpr <- fluco[row, "AveExpr"]
      if (length(refseq[[1]])>0) { # if the entry is not character(0)
        for (id in refseq[[1]]) { # For each refseqid in the list
          id_row <- data.frame(fluco_probe, id, logfc, aveexpr) 
          names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") # Label the row
          maptofluco <- rbind(maptofluco, id_row) 
        }
      }
    }
  }
}
# Pirinixic acid
for (row in 1:nrow(piri)) {
  if (piri[row, "P.Value"] < 0.05) {
    if(abs(methyl[row, "logFC"]) > 1.5) {
      piri_probe<- piri[row,"X"] 
      refseq <- piri[row, "REFSEQID"] 
      logfc <- piri[row, "logFC"]
      aveexpr <- piri[row, "AveExpr"]
      if (length(refseq[[1]])>0) { 
        for (id in refseq[[1]]) { 
          id_row <- data.frame(piri_probe, id, logfc, aveexpr) 
          names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") 
          maptopiri <- rbind(maptopiri, id_row) 
        }
      }
    }
  }
}

# After taking the median of the logfoldchange-

# 3-methylcholanthrene 
methylmap <- maptomethyl %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))
# fluconazole
flucomap <- maptofluco %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))
# pirinixic acid
pirimap <- maptopiri %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))

# Read  DESeq data for each-
d_methyl <- read.csv("~/Documents/project-3-dachsund/ahr_cmc_deseq_results.csv", as.is = TRUE)
d_methyl
d_fluco <- read.csv("~/Documents/project-3-dachsund/carpxr_cornoil_deseq_results.csv", as.is = TRUE)
d_piri <- read.csv("~/Documents/project-3-dachsund/ppara_cmc_deseq_results.csv", as.is=TRUE)

# Filtering DEG where Log2foldchange>1.5 and nominal p-value < 0.05-
d_methyl <- subset(d_methyl, pvalue < 0.05)
d_methyl <- subset(d_methyl, abs(log2FoldChange) >1.5)

d_fluco <- subset(d_fluco, pvalue < 0.05)
d_fluco <- subset(d_fluco, abs(log2FoldChange) >1.5)

d_piri <- subset(d_piri, pvalue < 0.05)
d_piri <- subset(d_piri, abs(log2FoldChange) >1.5)

# Calculate concordance for each chemical
# Obtain counts of the number of differentially expressed genes that match between limma and deseq

# FIRST ANALYSIS (3-methylcholanthrene)-
methyl_count <- 0 
for (row in 1:nrow(d_methyl)) {
  id_deseq <- d_methyl[row, "X"] 
  logfc <- d_methyl[row, "log2FoldChange"] 
  limma_row <- methylmap[which(methylmap$REFSEQID == id_deseq), ] 
  if (nrow(limma_row) != 0) { 
    if(sign(logfc) == sign(limma_row[["logFC"]])) { 
      methyl_count <- methyl_count+1 
    }
  }
}
#SECOND ANALYSIS  (Fluconazole)-
fluco_count <- 0 
for (row in 1:nrow(d_fluco)) {
  id_deseq <- d_fluco[row, "X"] 
  logfc <- d_fluco[row, "log2FoldChange"] 
  limma_row <- flucomap[which(flucomap$REFSEQID == id_deseq), ] 
  if (nrow(limma_row) != 0) { 
    if(sign(logfc) == sign(limma_row[["logFC"]])) { 
      fluco_count <- fluco_count+1 
    }
  }
}
# Pirinixic acid
piri_count <- 0 
for (row in 1:nrow(d_piri)) {
  id_deseq <- d_piri[row, "X"]
  logfc <- d_piri[row, "log2FoldChange"] 
  limma_row <- pirimap[which(pirimap$REFSEQID == id_deseq), ] 
  if (nrow(limma_row) != 0) { 
    if(sign(logfc) == sign(limma_row[["logFC"]])) { 
      piri_count <- piri_count+1 
    }
  }
}

# Calculate concordance for each chemical

# 3-methylcholanthrene
methyl_concordance <- (2 * methyl_count)/(nrow(d_methyl) + nrow(methylmap))
# 0.2303143

# fluconazole
fluco_concordance <- (2 * fluco_count)/(nrow(d_fluco) + nrow(flucomap)) 
# 0.5534543

# pirinixic acid
piri_concordance <- (2 * piri_count)/(nrow(d_piri) + nrow(pirimap)) # 0.5214058

# Make scatter plot with treatment on x-axis and concordance of y-axis (figure 2a)
# make a table of the number of degs and concordance from deseq data
scatter <- data.frame("Concordance"=c(methyl_concordance, fluco_concordance, piri_concordance), 
                      "Treatment"=c(nrow(d_methyl), nrow(d_fluco), nrow(d_piri)))
chems <- c("3-methylcholanthrene", "Fluconazole", "Pirinixic acid")
plot1 <- ggplot(scatter, aes(x=Treatment, y=Concordance)) +
  geom_point(color="black", fill=c('red', 'blue', 'lightgreen'), shape=18, size=5) +
  stat_smooth(method = "lm", linetype="dashed", se=FALSE, color="black") +
  geom_text(label=chems, nudge_x = c(75, 0, -60), nudge_y = c(0, -0.01, 0)) +
  stat_cor(aes(label=paste(..rr.label..))) +
  labs(title="RNA-Seq Analysis", x="DEG NO", y="Concordance")


# Make the same table for limma data
scatter2 <- data.frame("Concordance"=c(methyl_concordance, fluco_concordance, piri_concordance), 
                       "Treatment"=c(nrow(methylmap), nrow(flucomap), nrow(pirimap)))
plot2 <- ggplot(scatter2, aes(x=Treatment, y=Concordance)) +
  geom_point(color="black", fill=c('red', 'blue', 'lightgreen'), shape=21, size=3)+
  geom_text(label=chems, nudge_x = c(12, -8, 5), nudge_y = c(0, 0, 0.007)) +
  stat_smooth(method = "lm", linetype="dashed", se=FALSE, color="black") +
  stat_cor(aes(label=paste(..rr.label..))) +
  labs(title="Microarray Analysis", x="DEG NO", y="Concordance")

# Display scatter plots side by side
ggarrange(plot1, plot2, labels = c("A", "B"), ncol = 2, nrow = 1)

# Using the DESEQ we would be finding out the median of the base mean col-
#for  3-methylcholanthrene-
# genes above median
methyl_med <- median(d_methyl$baseMean)
m_deseq_above <- subset(d_methyl, baseMean > methyl_med)
# genes below median
m_deseq_below <- subset(d_methyl, baseMean < methyl_med) 

# same for fluconazole
fluco_med <- median(d_fluco$baseMean)
f_deseq_above <- subset(d_fluco, baseMean > fluco_med) 
f_deseq_below <- subset(d_fluco, baseMean < fluco_med) 

#same for pirinixic acid
piri_med <- median(d_piri$baseMean)
p_deseq_above <- subset(d_piri, baseMean > piri_med) 
p_deseq_below <- subset(d_piri, baseMean < piri_med)

#  from limma results find median of the average expression-

# for 3-methylcholanthrene
methyl_limmamed <- median(methylmap$AveExpr)
m_limma_above <- subset(methylmap, AveExpr > methyl_limmamed)
m_limma_below <- subset(methylmap, AveExpr < methyl_limmamed)

#for  fluconazole
fluco_limmamed <- median(flucomap$AveExpr)
f_limma_above <- subset(flucomap, AveExpr > fluco_limmamed)
f_limma_below <- subset(flucomap, AveExpr < fluco_limmamed)

# for pirinixic acid
piri_limmamed <- median(pirimap$AveExpr)
p_limma_above <- subset(pirimap, AveExpr > piri_limmamed)
p_limma_below <- subset(pirimap, AveExpr < piri_limmamed)

# Recompute concordance for each of the above and below groups

#1st checmical 3-methylcholanthrene-
mcount <- 0
for (row in 1:nrow(m_deseq_below)) {
  id_deseq <- m_deseq_below[row, "X"]
  logfc <- m_deseq_below[row, "log2FoldChange"]
  below_row <- m_limma_below[which(m_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      mcount <- mcount+1
    }
  }
}
m_below_concord <- (2*mcount)/(nrow(m_deseq_below) + nrow(m_limma_below)) 
# 0.1120797

# fluconazole
fcount <- 0
for (row in 1:nrow(f_deseq_below)) {
  id_deseq <- f_deseq_below[row, "X"]
  logfc <- f_deseq_below[row, "log2FoldChange"]
  below_row <- f_limma_below[which(f_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      fcount <- fcount+1
    }
  }
}
f_below_concord <- (2*fcount)/(nrow(f_deseq_below) + nrow(f_limma_below)) 
# 0.3616364

# pirinixic acid
pcount <- 0
for (row in 1:nrow(p_deseq_below)) {
  id_deseq <- p_deseq_below[row, "X"]
  logfc <- p_deseq_below[row, "log2FoldChange"]
  below_row <- p_limma_below[which(p_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      pcount <- pcount+1
    }
  }
}
p_below_concord <- (2*pcount)/(nrow(p_deseq_below) + nrow(p_limma_below)) 
#value is 0.3516524

# Concordance for above groups-

# 3-methylcholanthrene
macount <- 0
for (row in 1:nrow(m_deseq_above)) {
  id_deseq <- m_deseq_above[row, "X"]
  logfc <- m_deseq_above[row, "log2FoldChange"]
  above_row <- m_limma_above[which(m_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      macount <- macount+1
    }
  }
}
m_above_concord <- (2*macount)/(nrow(m_deseq_above) + nrow(m_limma_above)) 
# value is 0.2515567

# fluconazole
facount <- 0
for (row in 1:nrow(f_deseq_above)) {
  id_deseq <- f_deseq_above[row, "X"]
  logfc <- f_deseq_above[row, "log2FoldChange"]
  above_row <- f_limma_above[which(f_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      facount <- facount+1
    }
  }
}
f_above_concord <- (2*facount)/(nrow(f_deseq_above) + nrow(f_limma_above)) # 0.5152451

# pirinixic acid
pacount <- 0
for (row in 1:nrow(p_deseq_above)) {
  id_deseq <- p_deseq_above[row, "X"]
  logfc <- p_deseq_above[row, "log2FoldChange"]
  above_row <- p_limma_above[which(p_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      pacount <- pacount+1
    }
  }
}
p_above_concord <- (2*pacount)/(nrow(p_deseq_above) + nrow(p_limma_above)) 
# value is 0.482746

# Make a barplot of all of these values
# make into a table
bar <- data.frame("Concordance"=c(m_above_concord, m_below_concord, methyl_concordance,
                                  f_above_concord, f_below_concord, fluco_concordance,
                                  p_above_concord, p_below_concord, piri_concordance), 
                  "Analysis"=c("Above", "Below", "Overall", "Above", "Below", "Overall", "Above", "Below", "Overall"),
                  "Chemical"=c("3-methylcholanthrene","3-methylcholanthrene","3-methylcholanthrene",
                               "Fluconazole","Fluconazole","Fluconazole",
                               "Pirinixic Acid","Pirinixic Acid","Pirinixic Acid"))
# Plot data with a side-by-side bar chart
ggplot(bar, aes(fill=Analysis, x=Chemical, y=Concordance)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(title="Above and Below Median Concordance")