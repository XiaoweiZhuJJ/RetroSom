#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
retro<-args[2]
set.seed(1)

library('randomForest')

### read the saved models ###
rf1.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_1.rds", sep="")
rf1 <- readRDS(rf1.filename)
rf2.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_2.rds", sep="")
rf2 <- readRDS(rf2.filename)
rf3.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_3.rds", sep="")
rf3 <- readRDS(rf3.filename)
rf4.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_4.rds", sep="")
rf4 <- readRDS(rf4.filename)
rf5.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_5.rds", sep="")
rf5 <- readRDS(rf5.filename)
rf6.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_6.rds", sep="")
rf6 <- readRDS(rf6.filename)
rf7.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_7.rds", sep="")
rf7 <- readRDS(rf7.filename)
rf8.filename <- paste(as.character(masterpath),"/LINE/02_PE_level1/RFI_8.rds", sep="")
rf8 <- readRDS(rf8.filename)

filename <- paste(as.character(datpath),"/list.txt", sep="")
subs <- read.table(filename, sep="\t", header=F)
for (i in 1:dim(subs)[1])
   {
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/",as.character(subs[i,1]),".pe.LINE.matrix", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_0/",as.character(subs[i,1]),".pe.LINE.matrix", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone  <- rbind(clone1, clone2)
    subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
    subclone$gap <- subclone$gap/(abs(subclone$refpos1 - subclone$refpos2))

    test1 <- subset(subclone, !is.na(ACAG))
    nx1 <- test1[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "ACAG")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$chr), test1$cord1, test1$cord2, as.character(test1$read), pred.RF1))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G1.txt", sep="")
    write(t(pred1), file=filename, ncol=5, sep="\t")

    test2 <- subset(subclone, !is.na(TAG))
    nx2 <- test2[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "TAG")]
    pred.RF2 <- predict(rf2,newdata=nx2,type='prob')[,2]
    pred2 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF2))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G2.txt", sep="")
    write(t(pred2), file=filename, ncol=5, sep="\t")

    test3 <- subset(subclone, is.na(TAG) & is.na(ACAG) & (ORF==1))
    nx3 <- test3[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ","depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "Xscore")]
    pred.RF3 <- predict(rf3,newdata=nx3,type='prob')[,2]
    pred3 <- data.frame(cbind(as.character(test3$chr), test3$cord1, test3$cord2, as.character(test3$read), pred.RF3))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G3.txt", sep="")
    write(t(pred3), file=filename, ncol=5, sep="\t")

    test4 <- subset(subclone, is.na(TAG) & is.na(ACAG) & (ORF==0) & (refpos < 6015))
    nx4 <- test4[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream")]
    pred.RF4 <- predict(rf4,newdata=nx4,type='prob')[,2]
    pred4 <- data.frame(cbind(as.character(test4$chr), test4$cord1, test4$cord2, as.character(test4$read), pred.RF4))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G4.txt", sep="")
    write(t(pred4), file=filename, ncol=5, sep="\t")

    test5 <- subset(subclone, is.na(TAG) & is.na(ACAG) & !is.na(GC5389))
    nx5 <- test5[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "GC5389")]
    pred.RF5 <- predict(rf5,newdata=nx5,type='prob')[,2]
    pred5 <- data.frame(cbind(as.character(test5$chr), test5$cord1, test5$cord2, as.character(test5$read), pred.RF5))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G5.txt", sep="")
    write(t(pred5), file=filename, ncol=5, sep="\t")

    test6 <- subset(subclone, is.na(TAG) & is.na(ACAG) & !is.na(G5533))
    nx6 <- test6[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ","depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "G5533", "C5533")]
    pred.RF6 <- predict(rf6,newdata=nx6,type='prob')[,2]
    pred6 <- data.frame(cbind(as.character(test6$chr), test6$cord1, test6$cord2, as.character(test6$read), pred.RF6))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G6.txt", sep="")
    write(t(pred6), file=filename, ncol=5, sep="\t")

    test7 <- subset(subclone, is.na(TAG) & is.na(ACAG) & !is.na(AT5710))
    nx7 <- test7[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream", "AT5710")]
    pred.RF7 <- predict(rf7,newdata=nx7,type='prob')[,2]
    pred7 <- data.frame(cbind(as.character(test7$chr), test7$cord1, test7$cord2, as.character(test7$read), pred.RF7))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G7.txt", sep="")
    write(t(pred7), file=filename, ncol=5, sep="\t")

    test8 <- subset(subclone, anchor_split==1)
    nx8 <- test8[,c("seg", "map", "anchor_align", "anchor_len", "anchor_mm", "anchor_direct", "anchor_insert", "anchor_seg", "anchor_map", "anchor_XS", "anchor_AS", "anchor_mapQ", "depth","map_len","short","end3","end5","direction", "refpos", "dist", "c1", "upstream")]
    pred.RF8 <- predict(rf8,newdata=nx8,type='prob')[,2]
    pred8 <- data.frame(cbind(as.character(test8$chr), test8$cord1, test8$cord2, as.character(test8$read), pred.RF8))
    filename <- paste(as.character(datpath),"/",as.character(subs[i,1]),"/retro_v",as.character(retro),"_1/LINE/",as.character(subs[i,1]),".pe.pred.G8.txt", sep="")
    write(t(pred8), file=filename, ncol=5, sep="\t")
   }
