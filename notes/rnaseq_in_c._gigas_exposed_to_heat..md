# RNAseq in C. gigas exposed to heat.

##Once you have the data on your computer perform analysis on R:

###Cargar libreria necesaria para el análisis

library(DESeq2)

#####Aparece 
Loading required package: 
GenomicRanges
BiocGenerics
parallel

###Leer data en R
data <- read.table("#ubicación o directorio/nombre del archivo.txt"#, header = T, sep = "\t")
rownames(data) <- data$Feature
data <- data[,-1]

#####Aparece 
Build Objects

###Specify which columns are in which groups

deseq2.colData <- data.frame(condition=factor(c(rep("Treated", 3), rep("Control", 3))), 
                             type=factor(rep("single-read", 6)))
rownames(deseq2.colData) <- colnames(data)
deseq2.dds <- DESeqDataSetFromMatrix(countData = data,
                                     colData = deseq2.colData, 
                                     design = ~ condition)
                                     
###Run Analysis
deseq2.dds <- DESeq(deseq2.dds)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]

#####Esto es lo que dice que hace
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

###Visualizar tabla

head(deseq2.res)

###Count number of hits with adjusted p-value less then 0.05

dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])

####Respuesta
[1] 3773    6

In [27]:
%%R
tmp <- deseq2.res

###The main plot
plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="darkgray",
     main="DEG Virus Exposure  (pval <= 0.05)",
     xlab="mean of normalized counts",
     ylab="Log2 Fold Change")
     
###Getting the significant points and plotting them again so they're a different color
tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="red")

### 2 F(old)C(hange) lines
abline(h=c(-1,1), col="blue")

### Crea la tabla para exportarla a excel
write.table(tmp.sig, "/Users/sr320/Desktop/ASI-rna-seq/output/Cgigas-DEGlist", row.names = T)

### Visualizar la tabla
head output/Cgigas-DEGlist 

### Plot de fold change de los genes, loos que estan sobre 1 se expresan positivamente y los que están bajo uno negativamente... eso muestra la diferenciación entre los organismos pre y post calentamiento que son los 2 tratamientos.
![image](/Users/alumnomatlab/Desktop/Rplot Cgigas.pdf)



