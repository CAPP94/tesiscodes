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
MIMatrix<- build.mim(microbedataset)
aracnenet<- aracne(MIMatrix);
aracneanet<- aracne.a(MIMatrix);
aracnemnet<- aracne.m(MIMatrix);
maxRelevanceNet<- mrnet(MIMatrix);
maxRelevanceBNet<- mrnetb(MIMatrix);
