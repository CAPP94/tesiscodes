datExpr00<-datExpr0

datExprA<-datExpr0[1:8,]
datExprA.pre<-datExpr0[1:4,]
datExprA.treat<-datExpr0[5:6,]
datExprA.post<-datExpr0[7:8,]
datExprAC<-datExpr0[c(grep(paste(c("A","C"),collapse="|"), row.names(datExpr0))),]
datExpr0.post<-datExpr0[c(grep(paste(c(4,5),collapse="|"), row.names(datExpr0))
),]
datExpr0.pre<-datExpr0[c(grep(paste(c(1,2),collapse="|"), row.names(datExpr0))
),]
datExpr0.treat<-datExpr0[c(grep(3,row.names(datExpr0))),]
datExprAC.post<-datExprAC[c(grep(paste(c(4,5),collapse="|"), row.names(datExprAC))
),]
datExprAC.pre<-datExprAC[c(grep(paste(c(1,2),collapse="|"), row.names(datExprAC))
),]
datExprAC.treat<-datExprAC[c(grep(3,row.names(datExprAC))),]



datExprAC.post<-datExpr0[c(7,8,17,18),]
datExprAC.pre<-datExpr0[c(1:4,15,16),]
datExprAC.treat<-datExpr0[c(5,6,14),]
datExprV<- c("datExpr0","datExpr0.post","datExpr0.pre","datExpr0.treat","datExprA",
             "datExprA.post","datExprA.pre","datExprA.treat","datExprAC","datExprAC.post",
             "datExprAC.pre","datExprAC.treat")
grep(paste(c(1,2)), row.names(datExpr0))
grep(paste(c(1,2),collapse="|"), row.names(datExpr0))

graph_from_adjacency_matrix()

List1 <- names(datExpr00.post)
List2 <- names(datExpr00.treat)
List3 <- names(datExpr00.pre)
Lists <- list(List1, List2, List3)  #put the word vectors into a list to supply lapply
items <- sort(unique(unlist(Lists)))   #put in alphabetical order
MAAT <- matrix(rep(0, length(items)*length(Lists)), ncol=3)  #make a matrix of 0's
colnames(MAAT) <- paste0("List", 1:3)
rownames(MAAT) <- items
lapply(seq_along(Lists), function(i) {   #fill the matrix
  MAAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
})

MAAT   #look at the results
library(venneuler)
v <- venneuler(MAT)
plot(v)
