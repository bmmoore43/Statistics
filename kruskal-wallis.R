setwd()
## filename
nc= "joiv_ADH_orthogroup_genelist_simple.txt_genepair.txt.rate_mod.txt"
## get output name

## read in table
p1 = read.table(nc, header=T, sep="\t",row.names=1)
newdata <- p1[,c(8,11,13)]
#newdata <- subset(p1, type == 'jaz2') 
#newdata <- subset(p1, Class == 'Slyc-GM_GM_Slyc-GM'|Class == 'Slyc-GM_GM_Slyc-SM'|Class == 'Slyc-GM_SM_Slyc-GM'|Class == 'Slyc-GM_SM_Slyc-SM'|Class == 'Slyc-SM_GM_Slyc-GM'|Class == 'Slyc-SM_GM_Slyc-SM'|Class == 'Slyc-SM_SM_Slyc-GM'|Class == 'Slyc-SM_SM_Slyc-SM')

newdata <- subset(newdata, type1 == 'basal') 
newdata <- subset(newdata, type1 == 'basal_angiosperm')
newdata <- subset(newdata, type1 == 'basal_angiosperm-noncanonical')
newdata <- subset(newdata, type1 == 'gymnosperm')
newdata <- subset(newdata, type2 == "eudicot-canonical"| type2 == 'eudicot-noncanonical'| type2 == 'monocot-canonical'| type2 == 'monocot-noncanonical')
newdata <- subset(newdata, type3 == "TyrA1"| type3 == 'TyrA2'| type3 == 'TyrA3')


## kruskal.test for non-parametric data
library(dunn.test)
nb = paste(c(basename(nc),"_k-s_test1.txt"),collapse='')
#single
df= newdata[,"type3"]
#colnames(p1)
df2= cbind.data.frame(newdata[,1],df)
fit <- kruskal.test(newdata[,1] ~ df, data=df2)
print(fit)
krusk <- data.frame("kruskaltest pvalue",fit$p.value)
print(krusk)
PT = dunn.test(df2$`newdata[, 1]`, g=df2$df,
               method="bh")
PT<- as.data.frame(PT)
summary(PT)
name <- "gymnosperm TyrA"
lapply(name, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))
lapply(krusk, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))
#lapply(fit$p.value, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))
lapply(PT, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))

#create empty dataframe with colnames
#x <- c(colnames(p1[3:1817])) 
#d = data.frame(x, krusk_pval=rep(0,1815), chisq=rep(0,1815), Z=rep(0,1815), P=rep(0,1815), p.adjust=rep(0,1815), comparison=rep(0,1815)) 
#loop
for (i in 2:62) {
  df = newdata[,c(1,i)]
  newdf <- na.omit(df)
  vname <- colnames(newdata[i])
  print(vname)
  fit <- kruskal.test(newdf[,2] ~ Class, data=newdf)
  print(fit)
  
  #krusk <- data.frame("kruskaltest pvalue",fit$p.value)
  krusk <- fit$p.value
  print(krusk)
  #d[i-1,2] <- krusk
  #name<-colnames(newdf[2])
  #name
  if (is.nan(krusk)==TRUE){print(vname)}
  else {
    #write.table(name, file = nb, append= T, sep="\t")
    #lapply(krusk, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))
  PT = dunn.test(newdf[,2], g=newdf$Class,
                 method="bh")
  PT<- as.data.frame(PT)
  print(PT)
  feature <- c(vname)
  kruskv <- c(krusk)
  df <- cbind(feature,kruskv,PT)
  #d[i-1:i+6,2:7] <- df
  write.table(df, file = nb, append= T, sep='\t')
  }
}
#write.table(d, file = nb, append= T, sep="\t")
