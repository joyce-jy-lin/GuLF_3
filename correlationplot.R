###### spearmans correlation between metals------------------------
library("Hmisc")
library(ggplot2)
library(reshape2)
setwd("/Users/joycelin/Desktop/Gulf/Aim3")
dat<- read_excel("CE_neuroscaled_quantiles.xlsx")

metaldat <- dat[, c(6:23)]
dat <- subset(metaldat, select = -c(Cd, Co, V, Mo, Sb, Ni))

# rearrange columns grouping essential together and toxic together
metaldat1<- dat[, c(2, 9, 11, 12, 4, 7, 5, 6, 3, 1, 10, 8)]

cormat <- cor(metaldat1, method = c("spearman"), use = "complete.obs")
cormat <- round(cormat, 2)
cormat

melted_cormat <- melt(cormat)
head(melted_cormat)

plot <-ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,1), space = "Lab", 
                       name="Correlation\n") +
  theme(legend.title=element_text(size=11), legend.text = element_text(size=10)) + geom_text(aes(label = value), color = "black", size = 3) + coord_fixed()

plot

png("corrplot.png", width = 6, height = 6, units = 'in', res = 300) 
plot
dev.off() 






get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)

plot<-ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman's\nCorrelation") +
  theme_minimal() + geom_text(aes(label = value), color = "black", size = 2) + coord_fixed()

plot

png("corrplot.png", width = 6, height = 6, units = 'in', res = 300) 
plot
dev.off() 
