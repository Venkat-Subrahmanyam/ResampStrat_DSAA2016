#DATASETS

#DATASETS 1 - 13
load("tseries.RData")
a <- 1 #1-13
new_d <- tseries[[a]]

#DATASETS 14 - 20
d <- read.csv("data_akbilgic.csv")
a <- 4 #4-10
d[,a] <- as.numeric(gsub(",",".",d[,a]))
d$date <- dmy(d$date,locale="pt_PT")
d <- d[,c(1,a)]
new_d <- as.xts(d[,-1],order.by=d[,1])

#DATASET 21
load("NSW1_19990101_20120901.Rdata")
d <- d[,1:2]
d[,1] <- as.POSIXct(d$SETTLEMENTDATE)
new_d <- as.xts(d[,-1],order.by=d[,1])

#DATASET 22
load("NSW1_19990101_20120901.Rdata")
d <- d[,c(1,3)]
d[,1] <- as.POSIXct(d$SETTLEMENTDATE)
new_d <- as.xts(d[,-1],order.by=d[,1])

#DATASETS 23, 24

load("porto.20130206.20160111.Rdata")
a <- 4 #4,10
dat[,1] <- as.character(dat[,1])
d <- dat[dat[,1]==unique(dat[,1])[a],]
d["Timestamp"] <- paste0(d$Data," ",d$Hora)
d[,6] <- as.POSIXct(d$Timestamp)
d <- d[,5:6]
new_d <- as.xts(d[,-2],order.by=d[,2])


