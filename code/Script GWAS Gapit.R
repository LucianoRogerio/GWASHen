
###Análise GWAS pelo pacote GAPIT

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GAPIT")

library(BiocManager); library(compiler); library(sommer)
source("http://www.bioconductor.org/biocLite.R") 
source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")




dim(DArT)
head(drgphen)
drgphen$X <- NULL
colnames(drgphen)[colnames(drgphen)=="CLONE"] <- "Taxa"
myY <- drgphen[,-c(4+((0:5)*3))]
myY <- myY[,-c(1+((1:6)*2))]

DArT <- as.data.frame(DArT)
DArT2 <- data.frame(Taxa = rownames(DArT), DArT)
rownames(DArT2) <- c(1:nrow(DArT2))
myGD <- DArT2

PosMrkr <- data.frame(Name = colnames(DArT))
PosMrkr$Chromosome <- as.numeric(gsub ("S([0-9]+)_([0-9]+)", "\\1", PosMrkr[,1]))
PosMrkr$Position <- as.numeric(gsub ("S([0-9]+)_([0-9]+)", "\\2", PosMrkr[,1]))
myGM <- PosMrkr

A <- A.mat(DArT)
A2 <- data.frame(Taxa = rownames(A), A)
rownames(A2) <- NULL
colnames(A2) <- NULL
myKi <- A2


#####GWAS Retencao Foliar

myGAPIT <- GAPIT(
Y=myY[,c("Taxa","RF")],
GD=myGD,
GM=myGM,
PCA.total=5,
Model.selection=TRUE,
model=c("GLM","MLM","CMLM","FarmCPU")
)

#####GWAS Mancha Branca

myGAPIT <- GAPIT(
Y=myY[,c("Taxa","MB")],
GD=myGD,
GM=myGM,
PCA.total=5,
Model.selection=TRUE,
model=c("GLM","MLM","CMLM","FarmCPU")
)


#####GWAS Mancha Parda

myGAPIT <- GAPIT(
Y=myY[,c("Taxa","MP")],
GD=myGD,
GM=myGM,
PCA.total=5,
Model.selection=TRUE,
model=c("GLM","MLM","CMLM","FarmCPU")
)



#####GWAS Queima das Folhas

myGAPIT <- GAPIT(
Y=myY[,c("Taxa","QF")],
GD=myGD,
GM=myGM,
PCA.total=4,
Model.selection=TRUE,
model=c("GLM","MLM","CMLM","FarmCPU")
)


#####GWAS Vigor de Haste

myGAPIT <- GAPIT(
Y=myY[,c("Taxa","VH")],
GD=myGD,
GM=myGM,
PCA.total=5,
Model.selection=TRUE,
model=c("GLM","MLM","CMLM","FarmCPU")
)



#####GWAS Produtividade de Amido

myGAPIT <- GAPIT(
Y=myY[,c("Taxa","DRY")],
GD=myGD,
GM=myGM,
PCA.total=5,
Model.selection=TRUE,
model=c("GLM","MLM","CMLM","FarmCPU")
)


save.image(file = "Dados GWAS GAPIT.RData")
