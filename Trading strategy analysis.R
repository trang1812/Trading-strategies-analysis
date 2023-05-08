rm(list = ls())
library("matrixStats")
library("matrixcalc")
library("zoo")
#load data from 1960
data = read.csv('Industry17PortfoliosDaily.csv', skip = 9, header =T)
data <- data[data[,1]>19600101,]
data[,2:18]<-data[,2:18]/100
names(data)[1]<-'Date'
View(data)

#change date format
year <- substring(data$Date,1,4)
mont <- substring(data$Date,5,6)
day<- substring(data$Date,7,8)
temp <- paste(year, mont, day, sep = '-')
Date <- as.Date(temp)
industry<-colnames(data)[2:18]

#1. PART 1L on average, do more volatile industries have higher average log returns?
# compute average return, logreturn, volatility for each industry
avgreturn<-NULL
vol<-NULL
avglogreturn<-NULL
for (i in (1:17)){
  avgreturn[i]<-mean(data[,industry[i]])*252
  vol[i]<-sd(data[,industry[i]])*sqrt(252)
  avglogreturn[i]<- mean(log(1 + data[,industry[i]]))*252}
# compute correlation between log return-volatility, return-volatility
corlogreturn<-cor(avglogreturn, vol)
correturn<-cor(avgreturn, vol)
# OLS
lm<-lm(avglogreturn~vol)
x <- seq(min(vol), max(vol), length.out=100)
y <- lm$coefficients[1] + lm$coefficients[2]*x
plot(vol, avglogreturn, xlab='average volatility', ylab='average log
return', cex.lab=0.7)
title('average log return in each industry as a function of average
volatility', cex.main = 0.7)
lines(x,y, col='red', lwd=3)
plot(vol, avgreturn)
abline(lm(avgreturn~vol), col="red")

#2. PART2: COMPARE STRATEGIES
# Compute log return
logreturn=cbind(data$Date)
for (i in (1:17)){
  logreturn<-cbind(logreturn,log(1 + data[,industry[i]]))}
#log return dataframe
lret <- data.frame(logreturn)
colnames(lret) <- colnames(data)

#2.1.Equal-weighted portfolio (each industry with weight 1/17)
lretPorfolio1 <- 1/17*(rowSums(lret[,2:18]))
#plot cumulative log returns
CrPorfoliolret1 <- cumsum(lretPorfolio1)
plot(Date, CrPorfoliolret1, type = 'l', ylab = 'Porfolio log
returns', main='Cumulative log returns for equal-weighted portfolio')
#calculate mean, std, sharp ratio
pfret1 <- 1/17*(rowSums(data[,2:18]))
meanpfret1<-mean(pfret1)*252
stdpfret1<-sd(pfret1)*sqrt(252)
sharpe1<-meanpfret1/stdpfret1

#2.2. Risk parity
#Compute weights
weight<-1/vol
nomalizedw<-weight/sum(weight)
sum(nomalizedw)
#plot cummulative log return
lretPorfolio2<-data.matrix(lret[,2:18])%*%matrix(weight/sum(weight))
CrPorfoliolret2 <- cumsum(lretPorfolio2)
plot(Date, CrPorfoliolret2, type = 'l', ylab = 'Porfolio log
returns', main='Cumulative log returns for Risk parity')
#calculate mean, std, sharpe ratio
pfret2<-data.matrix(data[,2:18])%*%matrix(weight/sum(weight))
meanpfret2<-mean(pfret2)*252
stdpfret2<-sd(pfret2)*sqrt(252)
sharpe2<-meanpfret2/stdpfret2

#2.3.12-month momentum
#compute 12 month momentum matrix
rebal_freq <- 21 #rebalance once a month
MOM<-matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2]) 
#momentum matrix
colnames(MOM) <- colnames(data)
MOM[,1]<-data$Date
View(MOM)
d<-data.matrix(data)
for(i in (250:dim(data)[1])){
  MOM[i,2:18] <- colProds((d+1)[(i-249):i,2:18])-1
}
MOM <- MOM[250:dim(data)[1],]
#create ranking matrix
ranks = matrix(NA, nrow = dim(MOM)[1], ncol = dim(MOM)[2])
colnames(ranks) <- colnames(MOM)
ranks[,1]<-MOM[,1]
for (i in (1:dim(MOM)[1])){
  if ((i-1) %% rebal_freq == 0){ranks[i,2:18] = rank(-MOM[i,2:18])}
}
ranks[,2:18]<-shift.down(ranks[,2:18], rows = 1)
ranks<-ranks[-1,]
ranks<-na.locf(na.locf(ranks), fromLast = TRUE)
View(ranks)
#create weight matrix
for (c in (colnames(ranks)[2:18])){ranks[,c]<-(ranks[,c]<=6)*1/6}
#create equal-weighted portfolio
ret3<-data[251:dim(data)[1],2:18]*ranks[,2:18]
#calculate mean, std, sharpe ratio
pfret3 <- rowSums(ret3)
meanpfret3<-mean(pfret3)*252
stdpfret3<-sd(pfret3)*sqrt(252)
sharpe3<-meanpfret3/stdpfret3
View(ret3)


#PART3: VOLATILITY TARGETING
#compute 12 month momentum matrix
MOM<-matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2]) #momentum
matrix
colnames(MOM) <- colnames(data)
MOM[,1]<-data$Date
d<-data.matrix(data)
for(i in (250:dim(data)[1])){
  MOM[i,2:18] <- colProds((d+1)[(i-249):i,2:18])-1
}
View(MOM)
#create ranking matrix
ranks = matrix(NA, nrow = dim(MOM)[1], ncol = dim(MOM)[2])
colnames(ranks) <- colnames(MOM)
ranks[,1]<-MOM[,1]
for (i in (1:dim(MOM)[1])){
  ranks[i,2:18] = rank(-MOM[i,2:18])
}
View(ranks)
dim(ranks)
ranks[,2:18]<-shift.down(ranks[,2:18], rows = 1)
for (c in (colnames(ranks)[2:18])){ranks[,c]<-(ranks[,c]<=6)}
#calculate percent of capital on risky investment
w<-matrix(rep(1/6,6))
alpha<-NULL
for (i in (251:dim(data)[1])){
  r<-matrix(rep(ranks[i,2:18], 63), ncol=17, nrow=63, byrow=TRUE)
  re63days<-r*data[(i-63):(i-1),2:18]
  re63days<-re63days[,colSums(re63days) != 0]
  covariance<-cov(re63days)
  alpha[i]<-0.2/sqrt(252*t(w)%*%covariance%*%w)
}
#compute return of industries in porfolio
return<-
  alpha[251:15627]*1/6*ranks[251:15627,2:18]*data[251:15627,2:18]
#compute mean, standard deviation, sharpe ratio
pfret <-(rowSums(return))
meanpfret<-mean(pfret)*252
stdpfret<-sd(pfret)*sqrt(252)
sharpe<-meanpfret/stdpfret
#compute log return and plot cummulative log return
lretpf<-
  alpha[251:15627]*1/6*ranks[251:15627,2:18]*lret[251:15627,2:18]
lretpf <-rowSums(lretpf)
#plot cumulative log returns
CrPorfoliolret <- cumsum(lretpf)
plot(Date[251:15627], CrPorfoliolret, type = 'l', ylab = 'Porfolio
log returns', main='Cumulative log returns for Colatility targeting',xlab="Date")

#PART4 Equal-weighted portfolio in the share αt of industries with positive momentum
#12 month mometum matrix
MOM17<-matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2]) #
colnames(MOM17) <- colnames(data)
MOM17[,1]<-data$Date
d<-data.matrix(data)
for(i in (250:dim(data)[1])){
  MOM17[i,2:18] <- colProds((d+1)[(i-249):i,2:18])-1
}
#compute alpha, portfolio return and log return
pfreturn<-NULL
pflogreturn<-NULL
for (i in (250:15626)){
  alpha<-sum(MOM17[i,2:18]>0)/17
  pfreturn[i+1]<-alpha*pfret1[i+1]+(1-alpha)*((1+0.03)^(1/252)-1)
  pflogreturn[i+1]<-alpha*lretPorfolio1[i+1]+(1-alpha)*log((1+0.03)^(1/252))
}
#delete NaN in first 250 rows due to momentum calculation
pfreturn<-pfreturn[251:15627]
pflogreturn<-pflogreturn[251:15627]
#compute mean, standard deviation, sharp ratio
meanpfret4<-mean(pfreturn)*252
stdpfret4<-sd(pfreturn)*sqrt(252)
sharpe4<-meanpfret4/stdpfret4
#plot cumulative log return
CrPorfoliolret4 <- cumsum(pflogreturn)
plot(Date[251:15627], CrPorfoliolret4, type = 'l', ylab = 'Porfolio
log returns', main='Cumulative log returns for an nomalized weighted
portfolio with alpha share on risky asset')