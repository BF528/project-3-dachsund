#Analyst -Shreenk.Project BF528: Applications in Translational Bioinformatics, Project 3

# We took group 6 chemicals for our analysis .It consists of  3-methylcholanthrene, fluconazole, and pirinixic acid.
#We run limma on microarray results of these.

# Load packages and set working directory
library(limma)
setwd("~/Documents/project-3-dachsund")

# Read group_6_mic_info.csv file which has the chemical id , names and control grps
samples <- read.csv('group_6_mic_info.csv',as.is=TRUE)

# Read The full RMA expression matrix :-
rma <- read.table('liver-normalization-rma.txt', sep='\t', as.is=TRUE,header=TRUE,row.names=1)

# Subset RMA expression matrix to just the chemicals in this GROUP6 -
#3 Different analysis 
rma.subset <- rma[paste0('X',samples$array_id)]

# Subset of 3-methylcholanthrene and control -
rma.3methyl <- rma.subset[paste0('X',samples$array_id[samples$chemical=='3-METHYLCHOLANTHRENE' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]

# Subset of fluconazole and control-
rma.fluco <- rma.subset[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]

# Subset of pirinixic acid and control-
rma.piri <- rma.subset[paste0('X',samples$array_id[samples$chemical=='PIRINIXIC_ACID' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]

#LIMMA-
# Creating design matrix for all 3 chemicals to be run by limma-

# First analysis, 3-methylcholanthrene :-
design_3methyl <- model.matrix(
  ~factor(samples$chemical[samples$vehicle=='CMC_.5_%'],
    levels=c('Control','3-METHYLCHOLANTHRENE')
  )
)
colnames(design_3methyl) <- c('Intercept','3-METHYLCHOLANTHRENE') 

# Second analysis, fluconazole:-
design_fluco <- model.matrix(
  ~factor(samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','FLUCONAZOLE')
  )
)
colnames(design_fluco) <- c('Intercept','FLUCONAZOLE') 

# Third analysis, pirinixic acid :-
design_piri <- model.matrix(
  ~factor(samples$chemical[samples$vehicle=='CMC_.5_%'],
    levels=c('Control','PIRINIXIC_ACID')
  )
)
colnames(design_piri) <- c('Intercept','PIRINIXIC_ACID') 


#LIMMA ANALYSIS FOR MICROARRAY DATA-
#First we fit the data to design matrix model.Then we run eBayes.Finally we run Top table to get DEG.
#then i plot histograms -

# IN ORDER TO PLOT HISTOGRAMS AND SCATTERPLOTS-
install.packages("ggplot2")                          # Install and load ggplot2
library("ggplot2")

##First analysis (LIMMA)-
fit_3methyl <- lmFit(rma.3methyl, design_3methyl)
fit_3methyl <- eBayes(fit_3methyl)
t <- topTable(fit_3methyl, coef=2, n=nrow(rma.3methyl), adjust='BH')
sig <- subset(t, adj.P.Val < 0.05)
print(nrow(sig)) 
##There are 58 significant genes
head(sig,10)
write.csv(t,'new_3methylcholanthrene_limma_results.csv') 
p<-ggplot(sig, aes(x=logFC)) + 
  geom_histogram(color="black", fill="white")
p
j<- ggplot(sig, aes(x=logFC, y=adj.P.Val)) +
  geom_point()
j


##Second Analysis-
fit_fluco <- lmFit(rma.fluco, design_fluco)
fit_fluco <- eBayes(fit_fluco)
t_fluco <- topTable(fit_fluco, coef=2, n=nrow(rma.fluco), adjust='BH')
sig <- subset(t_fluco, adj.P.Val < 0.05) 
print(nrow(sig))
head(sig,10)
## There are 8761 significant genes
write.csv(t_fluco,'new_fluconazole_limma_results.csv') 
p<-ggplot(sig, aes(x=logFC)) + 
  geom_histogram(color="black", fill="white", bins = 15)
p
j<- ggplot(sig, aes(x=logFC, y=adj.P.Val)) +
  geom_point(size=1)
j

##Third analysis-
fit_piri <- lmFit(rma.piri, design_piri)
fit_piri <- eBayes(fit_piri)
t_piri <- topTable(fit_piri, coef=2, n=nrow(rma.piri), adjust='BH')
sig <- subset(t_piri, adj.P.Val < 0.05)
print(nrow(sig)) 
##There are 1997 significant genes
head(sig,10)
write.csv(t_piri,'new_pirinixic_acid_limma_results.csv') 
p<-ggplot(sig, aes(x=logFC)) + 
  geom_histogram(color="black", fill="white",binwidth=1)
p
j<- ggplot(sig, aes(x=logFC, y=adj.P.Val)) +
  geom_point() 
j


