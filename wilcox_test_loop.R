setwd()
## filename
nc= "filename.txt"
## get output name
nb = paste(c(basename(nc),".wilcoxtest.txt"),collapse='')
## read in table
p1 = read.table(nc, header=T, sep="\t", na.strings = c("", " ","NA"))#,row.names=1)
#subset data
p1a <- subset(p1, Class == "Slyc-SM" | Class == "Slyc-GM")

#p1a <- na.omit(p1a)
features <- c(colnames(p1[3:1817])) #number of columns, not including the class column
d = data.frame(features, pval=rep(0,1815)) #create empty dataframe based on the number of columns of your original dataframe - 1

##loop to do mann-whitney test on all columns, then append pvalue to dataframe
for (i in 3:1817) { #change number to total number of columns 
  p1a[,i] <- as.numeric(p1a[,i])
  vname <- p1a[,i]
  df = data.frame(p1a$Class, vname)
  #print(df)
  newdf <- na.omit(df)
  #print(newdf)
  y <- wilcox.test(vname ~ p1a.Class, data=newdf, alternative = "two.sided")
  #print(y)
  d[i-1,2] <- y$p.value
  }

#write pvalue dataframe
write.table(d, file = nb, sep="\t")

## one at a time
# newdf <- na.omit(df)
x <- c(colnames(p1[2:2])) #number of columns, not including the class column
d = data.frame( x, pval=rep(0,2)) #create empty dataframe based on the number of columns of your original dataframe - 1
d
p1a$Class <- factor(p1a$Class)
p1a$num_co.local_genes_5<- as.numeric(p1a$num_co.local_genes)
test = data.frame(p1a$Class, p1a$num_co.local_genes)
test.x <- na.omit(test)
y <- wilcox.test(p1a.num_co.local_genes ~ p1a.Class, data=test.x, alternative = "two.sided")
y
y$p.value

d[1,1]<- "num_co.local_genes_5"
d[1,]<- cbind("num_co.local_genes_5",y$p.value)

newdf$Class <- factor(newdf$Class)
newdf$FamilySize<- as.numeric(newdf$FamilySize)
test = data.frame(newdf$Class, newdf$FamilySize)
test.x <- na.omit(test)
y <- wilcox.test(newdf.FamilySize ~ newdf.Class, data=test.x, alternative = "two.sided")
y
y$p.value
d[2,]<- cbind("num_co.local_genes_10",y$p.value)
