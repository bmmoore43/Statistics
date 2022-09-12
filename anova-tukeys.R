setwd()
nc= "wound_response_poscontrol_051720.txt"

nb = paste(c(basename(nc),".tukey.txt"),collapse='')
p1 = read.table(nc, header=T, sep="\t", na.strings = c("", " ","NA"),row.names=1)
p2 <- subset(p1, type == 'ger5')
##one test
df = p2[,c(1,3)]
newdf <- na.omit(df)
fit <- aov(RQ ~ class, data=newdf)
x1 <-summary(fit)
x1
lapply(x1, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))
x <- TukeyHSD(fit)
x
x2<-as.data.frame(x$class)
x2
#x2$name<-colnames(newdf[2])
write.table(x2, file = nb, append= T, sep="\t")

df = p1[,c(1:19)]
newdf <- na.omit(df)
colnames(newdf[2])
fit <- aov(newdf[,2] ~ Class, data=newdf)
x1 <-summary(fit)
lapply(x1, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))
x <- TukeyHSD(fit)
x
x2<-as.data.frame(x$Class)
x2$name<-colnames(newdf[2])
x2
write.table(x2, file = nb, append= T, sep="\t")

for (i in 2:59) {print(i)}

##loop tukeys
for (i in 2:59) {
  df = p1[,c(1,i)]
  newdf <- na.omit(df)
  colnames(newdf[2])
  fit <- aov(newdf[,2] ~ Class, data=newdf)
  x1 <-summary(fit)
  lapply(x1, function(x) write.table(data.frame(x), file = nb, append= T, sep='\t'))
  x <- TukeyHSD(fit)
  x
  x2<-as.data.frame(x$Class)
  x2$name<-colnames(newdf[2])
  x2
  write.table(x2, file = nb, append= T, sep="\t")
}
write.table(x2, file = nb, append= T, sep="\t")


