#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
retro<-args[2]
set.seed(1)

library('randomForest')
rf1.filename <- paste(as.character(masterpath),"/LINE/04_SR_level1/RFII_1.rds", sep="")
rf1 <- readRDS(rf1.filename)
rf2.filename <- paste(as.character(masterpath),"/LINE/04_SR_level1/RFII_2.rds", sep="")
rf2 <- readRDS(rf2.filename)

filename <- paste(as.character(datpath),"/list.txt", sep="")
subs <- read.table(filename, sep="\t", header=F)
for (i in 1:dim(subs)[1])
   {
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/",as.character(subs[i,1]),".sr.LINE.matrix", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_0/",as.character(subs[i,1]),".sr.LINE.matrix", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone  <- rbind(clone1, clone2)
    subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
    subclone$gap <- subclone$gap/(abs(subclone$refpos1 - subclone$refpos2))

    test1 <- subclone
    nx1 <- test1[,c("seg", "map", "depth","map_size","map_ratio", "short","end3","end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS", "oldSR")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$chr), test1$cord1, test1$cord2, as.character(test1$read), pred.RF1))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".sr.pred.G0.txt", sep="")
    write(t(pred1), file=filename, ncol=5, sep="\t")

    test2 <- subset(subclone, oldSR==1)
    nx2 <- test2[,c("seg", "map", "depth","map_size","map_ratio", "short","end3","end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS")]
    pred.RF2 <- predict(rf2,newdata=nx2,type='prob')[,2]
    pred2 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF2))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".sr.pred.G1.txt", sep="")
    write(t(pred2), file=filename, ncol=5, sep="\t")
   }

