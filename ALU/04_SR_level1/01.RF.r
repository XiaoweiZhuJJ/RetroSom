#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
retro<-args[2]
masterpath<-args[3]
library('randomForest')
library('PRROC')

### read the saved models ###
rf1.filename <- paste(as.character(masterpath),"/ALU/04_SR_level1/RFIV.rds", sep="")
rf1 <- readRDS(rf1.filename)

filename <- paste(as.character(datpath),"/list.txt", sep="")
subs <- read.table(filename, sep="\t", header=F)
for (i in 1:dim(subs)[1])
   {
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_0/",as.character(subs[i,1]),".sr.ALU.matrix", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/",as.character(subs[i,1]),".sr.ALU.matrix", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone  <- rbind(clone1, clone2)
    subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
    subclone$refpos <- (subclone$refpos1 + subclone$refpos2)/2

    test1 <- subclone
    nx1 <- test1[,c("seg", "map", "depth","map_size", "map_ratio", "end3", "end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$chr), test1$cord1, test1$cord2, as.character(test1$read), pred.RF1))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/ALU/",as.character(subs[i,1]),".sr.pred.T1.txt", sep="")
    write(t(pred1), file=filename, ncol=5, sep="\t")
   }
