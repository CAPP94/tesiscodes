library(minet)
library(parmigene)
workingDir = "C:/Users/Charlie/Documents/UAQ/Probable tesis/data/main";
setwd(workingDir);

#Read in the microbiome data set
AbunData = read.csv("AbundancesMatt2.csv");
# Take a quick look at what is in the data set:
dim(AbunData);
names(AbunData);
microbedataset = as.data.frame(t(AbunData[, -c(1)]));
names(microbedataset) = AbunData$X;
rownames(microbedataset) = names(AbunData)[-c(1)];

gsgMI = WGCNA::goodSamplesGenes(microbedataset, verbose = 3);
gsgMI$allOK
if (!gsgMI$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgMI$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(microbedataset)[!gsgMI$goodGenes], collapse = ", ")));
  if (sum(!gsgMI$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(microbedataset)[!gsgMI$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  microbedataset = microbedataset[gsgMI$goodSamples, gsgMI$goodGenes]
}
dim(microbedataset)

MIMatrix<- build.mim(microbedataset)
aracnenet<- aracne(MIMatrix);
aracneanet<- aracne.a(MIMatrix);
aracnemnet<- aracne.m(MIMatrix);
maxRelevanceNet<- mrnet(MIMatrix);
maxRelevanceBNet<- mrnetb(MIMatrix);
