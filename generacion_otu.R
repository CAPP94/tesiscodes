  ## Primero se identifican los datos de secuenciacion del microbioma a trabajar
#Para su obtencion se utilisaran las funciones de SRAdb
library(SRAdb)
# para su funcionamieto se requiere descargar los metadatos de SRA como se muestra a continuacón
# getSRAdbFile(destdir = getwd(), destfile = "SRAmetadb.sqlite.gz", method)
# una vez teniendo los metadatos en el directorio se establecen las siguientes variables
seqdata<-("SRP128025")#clave/es SRA de los datos de secuenciacion/nes 
sra_dbname <- 'SRAmetadb.sqlite'#asegurate de que este en el mismo directorio
sra_con <- dbConnect( dbDriver("SQLite"), sra_dbname )
#obtencion de los archivos SRA
#getSRAfile( in_acc = c("SRR000648","SRR000657"), sra_con = sra_con, destDir = getwd(), fileType = 'sra' )
#obtencion de los archivos FASTQ (existe la función getFASTQfile, pero esta me funcionó)
getSRAfile( in_acc=seqdata, sra_con, destDir = getwd(), fileType = 'fastq', srcType = 'ftp', makeDirectory = FALSE, method = 'curl', ascpCMD = NULL )

  ## Agrupacion de datos y generación de ASV´s/OTU´s
library(dada2); packageVersion("dada2")
#Confirmar que los archivos estan presentes
dirarch<- setwd(choose.dir()) # modificar si los archivos estan encarpetados
list.files(dirarch) #observar que deben tener terminacion ".fastq.gz"
# Archivos forward contienen un ".._1.." y reverse un "..._2.."
fnFs <- sort(list.files(path =  dirarch, pattern=c("_R1_001.fastq","_1.fastq"), full.names = TRUE))
fnRs <- sort(list.files(path =  dirarch, pattern=c("_R2_001.fastq","_2.fastq"), full.names = TRUE))
# Extraer nombre, asumiendo el siguiente formato: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#Inspeccionar calidad
plotQualityProfile(fnFs[1:4]) # forward
plotQualityProfile(fnRs[1:4])# reverse
#con las graficas se observa la calidad y se establece el punto de corte (eje X) donde se obtenga una calidad arriba de 30 (eje Y)
#Corte y filtrado
# Colocar archivos filtrados en subdirectorio
filtFs <- file.path(path= dirarch, "filtered", paste0(sample.names, "cacti_F_filt.fastq.gz"))#modificar cacti
filtRs <- file.path(path = dirarch, "filtered", paste0(sample.names, "cacti_R_filt.fastq.gz"))#modificar cacti
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,175),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out) #observar que la opcion anterior generó los nuevos archivos ya editados
save(filtFs, filtRs, file = "filts.RData")
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
#Asignar taxonomia
#97% identidad
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
#100% Identidad
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v128.fa.gz")



