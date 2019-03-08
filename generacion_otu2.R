#funciones a utilizar
  ## Primero se identifican los datos de secuenciacion del microbioma a trabajar
#Para su obtencion se utilisaran las funciones de SRAdb
dwldFastq<- function(seqdata) {
library(SRAdb)
# para su funcionamieto se requiere descargar los metadatos de SRA como se muestra a continuacón
# getSRAdbFile(destdir = getwd(), destfile = "SRAmetadb.sqlite.gz", method)
# una vez teniendo los metadatos en el directorio se establecen las siguientes variables
print("Selecciona archivo 'SRAmetadb.sqlite'")
sra_dbname <- file.choose()#asegurate de que este en el mismo directorio
sra_con <- dbConnect( dbDriver("SQLite"), sra_dbname )
#obtencion de los archivos SRA
#getSRAfile( in_acc = c("SRR000648","SRR000657"), sra_con = sra_con, destDir = getwd(), fileType = 'sra' )
#obtencion de los archivos FASTQ (existe la función getFASTQfile, pero esta me funcionó)
getSRAfile( in_acc=seqdata, sra_con, destDir = path, fileType = 'fastq', srcType = 'ftp', makeDirectory = FALSE, method = 'curl', ascpCMD = NULL )
dbDisconnect()}
print("Selecciona direccion a trabajar con los archivos FASTQ")
path <-dirname(file.choose()) # CHANGE ME to the directory containing the fastq files after unzipping.
# setwd(path)
print("Seleccionar algun archivo de la carpeta 'tax' para asignar OTU´s")
taxafolder<-dirname(file.choose())
preg<- menu(c("Sí", "No"), T,"¿Quieres descargar los fastq? ")
if(preg ==1L){
  dwldFastq(readline("Escribe SRA o vector con los SRA´s a descargar"
                     ))
}  
## Agrupacion de datos y generación de ASV´s/OTU´s
library(dada2); packageVersion("dada2")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=c("_R1_001.fastq","_1.fastq"), full.names = TRUE))
fnRs <- sort(list.files(path, pattern=c("_R2_001.fastq","_2.fastq"), full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c( 240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
save(out, filtFs, filtRs, file = "filts.RData")
#Porcentaje de error
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) #Grafico con los errores de secuencicacion
save(errF, errR, file = "errorRate.RData")
#Desreplicación
derepFs <- derepFastq(filtFs, verbose=TRUE)
save(derepFs, sample.names, file = "derep.RData")
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names #Nombrar objeto
names(derepRs) <- sample.names #Nombrar objeto
save(derepFs, derepRs, file = "derep.RData")
#Inferir muestras
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]] #Revisar objeto
#Unir secuencias
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])#Revision
#Construir tabla de secuencias
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #Recordar esta dimension
table(nchar(getSequences(seqtab)))
#Remover quimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #Comparar dimension
sum(seqtab.nochim)/sum(seqtab) # Porcentaje de no quimeras
#Revisar resultados del proceso
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
save(mergers, seqtab, seqtab.nochim, track, file= "seqs.RData")
#Asignar taxonomia
#97% identidad
taxa0.97 <- assignTaxonomy(seqtab.nochim,file.path(getwd(), "/silva_nr_v128_train_set.fa.gz"), multithread=TRUE)
#100% Identidad
taxa1.00 <- addSpecies(taxa0.97, file.path(getwd(), "/silva_species_assignment_v128.fa.gz"))
#Inspeccionar taxa
taxa.print <- taxa1.00 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
save(taxa0.97,taxa1.00, file= "taxas.RData")
#phyloseq
## NO OLVIDES AJUSTAR LA INFORMACIÓN A LOS DATOS OBTENIDOS O DARLE FLEXIBILIDAD AL PROGRAMA
library(phyloseq); packageVersion("phyloseq")
# library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1) #This line may change depending on your data
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa1.00))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
# Generating abundances tables 
library(microbiome)
abunE <- abundances(ps,transform = "compositional")
View(taxa1.00)
rownames(abunE)<-taxa1.00[,5] # family = 5 
save(abunE, file= "Abundances.RData")
write.csv(abunE, "AbundancesMatt2.csv")
# Grafica
library(vegan)
datos_abun<-as.data.frame(t(abunE))
S <- specnumber(datos_abun) # observed number of species
(raremax <- min(rowSums(datos_abun)))
Srare <- rarefy(datos_abun, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(datos_abun, step = 20, sample = raremax, col = "blue", cex = 0.6)

# plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
# # Transform data to proportions as appropriate for Bray-Curtis distances
# ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
# ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
