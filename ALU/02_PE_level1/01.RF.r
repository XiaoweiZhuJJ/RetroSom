#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
retro<-args[2]
masterpath<-args[3]
library('randomForest')
library('PRROC')

### read the saved models ###
rf1.filename <- paste(as.character(masterpath),"/ALU/02_PE_level1/RFIII_1.rds", sep="")
rf1 <- readRDS(rf1.filename)
rf2.filename <- paste(as.character(masterpath),"/ALU/02_PE_level1/RFIII_2.rds", sep="")
rf2 <- readRDS(rf2.filename)

filename <- paste(as.character(datpath),"/list.txt", sep="")
subs <- read.table(filename, sep="\t", header=F)
for (i in 1:dim(subs)[1])
   {
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/",as.character(subs[i,1]),".pe.ALU.matrix", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_0/",as.character(subs[i,1]),".pe.ALU.matrix", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone <- rbind(clone1, clone2)
    subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
    subclone$gap <- subclone$gap/(abs(subclone$refpos1 - subclone$refpos2))
    subclone$frag <- abs(subclone$refpos1 - subclone$refpos2) / 151

    test1 <- subset(subclone, anchor_split==1)
    nx1 <- test1[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_direct", "anchor_insert", "anchor_seg", "anchor_map", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","end3", "end5", "trans2","direction", "refpos", "dist", "c1", "upstream", "frag")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$chr), test1$cord1, test1$cord2, as.character(test1$read), pred.RF1))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/ALU/",as.character(subs[i,1]),".pe.pred.A1.txt", sep="")
    write(t(pred1), file=filename, ncol=5, sep="\t")

    test2 <- subset(subclone, anchor_split==0)
    nx2 <- test2[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","end3", "end5", "trans2","direction", "refpos", "dist", "c1", "upstream", "frag")]
    pred.RF2 <- predict(rf2,newdata=nx2,type='prob')[,2]
    pred2 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF2))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/ALU/",as.character(subs[i,1]),".pe.pred.A2.txt", sep="")
    write(t(pred2), file=filename, ncol=5, sep="\t")
   }


