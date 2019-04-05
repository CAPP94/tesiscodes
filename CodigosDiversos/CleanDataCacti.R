library(vcdExtra)
library(tidyverse)
library(readr)
workingDir = "C:/Users/Charlie/Documents/UAQ/Probable tesis/data/main";
setwd(workingDir);
# Meta data from SRA project
Metadataraw<-read_csv("SraSRP068631.csv")
MetadataFilt<-Metadataraw[c(1,11,21,25,30)] # Keep only the information you need
grep(paste(c("A","C"),collapse="|"), row.names(datExpr0))
MetadataClean <- MetadataFilt %>%
  # mutate(Place = gsub(c("MG|SF", "BSC|RZ|RHIZO|RENDO|LENDO|LEPI", 
  #                          "SUM|SPR", "mg|or", "16S|ITS2","_","-","-1-"),"", Place ))
  tidyr::extract(., col = SampleName
                 , into = c("Place", "Holobiont_Region", "Season","Specie", 
                            "Sequence_Technique" )
                 , regex ="([[:upper:]]{2})-1-([[:upper:]]+)-([[:upper:]]{3})-
                 ([[a-z]]{2})_([[:alnum:]]+)" )
# Use filter to skip data you dont need 
filter(MetadataClean, Sequence_Technique == "16S" )
