setwd()
library(doBy)
## filename
nc= "Slyc_classfile_SMvsDAvsGM_Unk_TcycE2P2_20190531_FINAL.txt_continuous_mat.txt_nodups.txt"
## get output name
nb = paste(c(basename(nc),".stats.txt"),collapse='')
## read in table
p1 = read.table(nc, header=T, sep="\t", na.strings = c("", " ","NA"))#,row.names=1)
p1a <- subset(p1, Class == "Slyc-SM" | Class == "Slyc-GM")

#p1a <- na.omit(p1a)
x <- c(colnames(p1[3:1817])) #number of columns, not including the class column
d = data.frame(x, Class=rep(0,1815), Mean=rep(0,1815), Max=rep(0,1815), Min=rep(0,1815), Median=rep(0,1815), SDev=rep(0,1815)) #create empty dataframe based on the number of columns of your original dataframe - 1
p1a <- subset(p1, Class == "Slyc-GM" )
##loop to do median test on all columns, then append pvalue to dataframe
for (i in 3:1817) { #change number to total number of columns 
  p1a[,i] <- as.numeric(p1a[,i])
  #print(colnames(p1a[i]))
  vname <- p1a[,i]
  
  df = data.frame(p1a$Class, vname)
  #print(df)
  newdf <- na.omit(df)
  #if (nlevels(newdf$p1a.Class) == 2){
  #print(newdf)

  y <-summaryBy(vname ~ p1a.Class, data=newdf, 
            FUN = list(mean, max, min, median, sd))
  #print(y)
  d[i-1,2:7] <- y
 
}
#else {print(nlevels(newdf$p1a.Class))}}
#write pvalue dataframe
nb = paste(c(basename(nc),".GMstats.txt"),collapse='')
write.table(d, file = nb, sep="\t")
warnings()
