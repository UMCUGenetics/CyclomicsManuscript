# @Date created: 16 May 2019
# @Author: Myrthe Jager
# @Description: Plot FP consensus loop (1-40+ times)
# @Abbreviations: BB = backbone/PJET, I = insert, FP = false positive, FPDEL = FP + deletions
# @Version: 1 July 2020




# ---- GET STARTED ----

#1 Load required packages
library(ggplot2)
library(gridExtra)

#2 Define output directory
dir = "/.../Cyclomics_manuscript/"
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/consensus_loop/",sep = "")
datadir = paste(indir,"RCA/",sep="")
remove(dir)

#3 Functions
#3A See separate R script for general functions: functions.R




# ---- GET REF SEQUENCES ----

#1 BB
sequences_backbones <- read.table(paste(indir,"sequences/backbones.txt",sep=""), header=TRUE,sep="\t")

#2 TP53
sequences_tp53 <- read.table(paste(indir,"sequences/tp53.txt",sep=""), header=TRUE,sep="\t")
#2A Get all amplicons analyzed
all_amplicons_tp53 <- unlist(strsplit(as.character(sequences_tp53[which(sequences_tp53$name != "ampliconS0"),]$name), split = "con"))
all_amplicons_tp53[grep(all_amplicons_tp53, pattern = "mpl")] <- c(",")
all_amplicons_tp53 <- paste(all_amplicons_tp53[-1],collapse="")




# ---- LOAD METADATA ----

#1 All
metadata <- read.delim2(paste(indir,"Metadata.txt",sep=""))

#1B Change amplicons 'all' to ones that are analyzed.
metadata$Amplicon.s. <- as.character(metadata$Amplicon.s.)
metadata[grep(pattern = "ALL",toupper(metadata$Amplicon.s.)),]$Amplicon.s. <- paste(unlist(strsplit(all_amplicons_tp53,split = ",12")),collapse="")

#2 Per flowcell type
metadata.r9 <- metadata[which(metadata$Flowcell.version == "R9"),]
metadata.r10 <- metadata[which(metadata$Flowcell.version == "R10"),]
metadata.flongle <- metadata[which(metadata$Flowcell.version == "Flongle"),]

#3 Remove general metadata
remove(metadata)




# ---- GET SAMPLE INFORMATION, OPEN & PREPARE FILES ----

#1 Get sample info from files using metadata
#1A Get info
sample.info.r9 <- getsampleinfo.percount(metadata = metadata.r9)
sample.info.r10 <- getsampleinfo.percount(metadata = metadata.r10)
sample.info.flongle <- getsampleinfo.percount(metadata = metadata.flongle)
#1B Check if numbers match
nrow(sample.info.r9)/40/2 == nrow(metadata.r9)
nrow(sample.info.r10)/40/2 == nrow(metadata.r10)
nrow(sample.info.flongle)/40/2 == nrow(metadata.flongle)
#1C Remove S0 files from flongle data
sample.info.flongle <- sample.info.flongle[-which(sample.info.flongle$position == "TP53 amplicon S0"),]

#2 Read files
#2A Open
myfiles.r9 = lapply(as.character(sample.info.r9$file), read.delim)
myfiles.r10 = lapply(as.character(sample.info.r10$file), read.delim)
myfiles.flongle = lapply(as.character(sample.info.flongle$file), read.delim)
#2B Add names
names(myfiles.r9) = paste(sample.info.r9$samplename,sample.info.r9$type,sample.info.r9$position,sample.info.r9$analysis,sep="__")
names(myfiles.r10) = paste(sample.info.r10$samplename,sample.info.r10$type,sample.info.r10$position,sample.info.r10$analysis,sep="__")
names(myfiles.flongle) = paste(sample.info.flongle$samplename,sample.info.flongle$type,sample.info.flongle$position,sample.info.flongle$analysis,sep="__")

#3 Prepare files for stat calculations
myfiles.r9 <- preparefiles(filelist = myfiles.r9,sampleinfo = sample.info.r9)
myfiles.r10 <- preparefiles(filelist = myfiles.r10,sampleinfo = sample.info.r10)
myfiles.flongle <- preparefiles(filelist = myfiles.flongle,sampleinfo = sample.info.flongle)




# ---- CALCULATE FP AND Q PER TABLE ----

#1 Get error rates per table
stats.r9 <- getstats.pertable(filelist = myfiles.r9,sampleinfo = sample.info.r9, mincoverage = 100,blackposlist = sample.info.r9$mutation,perinsert = NA)
stats.r10 <- getstats.pertable(filelist = myfiles.r10,sampleinfo = sample.info.r10, mincoverage = 100, blackposlist = sample.info.r10$mutation,perinsert = NA)
stats.flongle <- getstats.pertable(filelist = myfiles.flongle,sampleinfo = sample.info.flongle, mincoverage = 100, blackposlist = sample.info.flongle$mutation,perinsert = NA)

#2 Get error rates per table per insert
#2A r9
stats.perI.r9 <- getstats.pertable(filelist = myfiles.r9,
                                         sampleinfo = sample.info.r9,
                                         mincoverage = 100,
                                         blackposlist = sample.info.r9$mutation,
                                         perinsert = sequences_tp53)
#2B Flongle
stats.perI.flongle <- getstats.pertable(filelist = myfiles.flongle,
                                         sampleinfo = sample.info.flongle,
                                         mincoverage = 100,
                                         blackposlist = sample.info.flongle$mutation,
                                         perinsert = sequences_tp53)
#2C r10
stats.perI.r10 <- getstats.pertable(filelist = myfiles.r10,
                                         sampleinfo = sample.info.r10,
                                         mincoverage = 100,
                                         blackposlist = sample.info.r10$mutation,
                                         perinsert = sequences_tp53)

#3 Get error rates for amplicon 12, except the mutated bases in the triple mutant
#3A First get a blacklist and exclude the triple mutant positions from all files (to make it comparable)
black.e12 <- data.frame(position = sample.info.r9$position,
                        mutation = sample.info.r9$mutation)
black.e12[grep(pattern = "12", black.e12$position),]$mutation <- as.character(unique(
  sample.info.r9[which(sample.info.r9$backbone == "PJET" & sample.info.r9$position != "PJET" & sample.info.r9$mutation != "NA"),]$mutation))
#3B Run
stats.r9.e12 <- getstats.pertable(filelist = myfiles.r9,
                                       sampleinfo = sample.info.r9,
                                       mincoverage = 100,
                                       blackposlist = black.e12$mutation,
                                       perinsert = sequences_tp53)
#3C Select amplicon 12
stats.r9.e12 <- stats.r9.e12[grep(pattern = "12", stats.r9.e12$REF),]
#3D Remove blacklist
remove(black.e12)
#3E Add BBtype as analysis
stats.r9.e12[grep(pattern = "PJET", stats.r9.e12$REF),]$ANALYSIS <- c("PCR-free")
stats.r9.e12[grep(pattern = "BB24", stats.r9.e12$REF),]$ANALYSIS <- c("PCR")
stats.r9.e12[grep(pattern = "BB25", stats.r9.e12$REF),]$ANALYSIS <- c("PCR")
stats.r9.e12.all <- stats.r9.e12
stats.r9.e12 <- stats.r9.e12[which(as.numeric(as.character(stats.r9.e12$TYPE)) >= 10),]




# ---- PLOT 1-FP ----

#1 All: FP

plot1.r9<- plotstats(plotdf = stats.r9,
                      groupon = "REF",title="R9",plotcolumn="1-fp",miny=0.95,maxy=1,col = T)

plot1.r10<- plotstats(plotdf = stats.r10,
                       groupon = "REF",title="R10",plotcolumn="1-fp",miny=0.95,maxy=1,col = T)

plot1.flongle <- plotstats(plotdf = stats.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="1-fp",miny=0.95,maxy=1,col = T)

final_plot_1 <- grid.arrange(
  plot1.r9,
  plot1.r10,plot1.flongle, ncol=3)

ggsave(plot = plot1.r9, filename = paste(outdir,"all/error/r9_error_all_fp.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_1, filename = paste(outdir,"all/error/error_all_fp.pdf",sep=""),width =20, height = 8)

#2 All: FPDEL

plot2.r9<- plotstats(plotdf = stats.r9,
                     groupon = "REF",title="R9",plotcolumn="1-error",miny=0.90,maxy=1,col = T)

plot2.r10<- plotstats(plotdf = stats.r10,
                      groupon = "REF",title="R10",plotcolumn="1-error",miny=0.90,maxy=1,col = T)

plot2.flongle <- plotstats(plotdf = stats.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="1-error",miny=0.90,maxy=1,col = T)

final_plot_2 <- grid.arrange(
  plot2.r9,
  plot2.r10,plot2.flongle, ncol=3)

ggsave(plot = plot2.r9, filename = paste(outdir,"all/error/r9_error_all_fpdel.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_2, filename = paste(outdir,"all/error/error_all_fpdel.pdf",sep=""),width =20, height = 8)

#3 All: FP per I

plot3.r9<- plotstats(plotdf = stats.perI.r9,
                     groupon = "REF",title="R9",plotcolumn="1-fp",miny=0.95,maxy=1,col = T)

plot3.r10<- plotstats(plotdf = stats.perI.r10,
                      groupon = "REF",title="R10",plotcolumn="1-fp",miny=0.95,maxy=1,col = T)

plot3.flongle <- plotstats(plotdf = stats.perI.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="1-fp",miny=0.95,maxy=1,col = T)

final_plot_3 <- grid.arrange(
  plot3.r9,
  plot3.r10,plot3.flongle, ncol=3)


ggsave(plot = plot3.r9, filename = paste(outdir,"all/error/r9_error_all_fp_perI.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_3, filename = paste(outdir,"all/error/error_all_fp_perI.pdf",sep=""),width =20, height = 8)

#4 All: FPDEL per I
plot4.r9<- plotstats(plotdf = stats.perI.r9,
                     groupon = "REF",title="R9",plotcolumn="1-error",miny=0.90,maxy=1,col = T)

plot4.r10<- plotstats(plotdf = stats.perI.r10,
                      groupon = "REF",title="R10",plotcolumn="1-error",miny=0.90,maxy=1,col = T)

plot4.flongle <- plotstats(plotdf = stats.perI.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="1-error",miny=0.90,maxy=1,col = T)

final_plot_4 <- grid.arrange(
  plot4.r9,
  plot4.r10,plot4.flongle, ncol=3)


ggsave(plot = plot4.r9, filename = paste(outdir,"all/error/r9_error_all_fpdel_perI.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_4, filename = paste(outdir,"all/error/error_all_fpdel_perI.pdf",sep=""),width =20, height = 8)

 
#5 Cleanup
remove(plot1.flongle,plot1.r10,plot1.r9,final_plot_1,
       plot2.flongle,plot2.r10,plot2.r9,final_plot_2,
       plot3.flongle,plot3.r10,plot3.r9,final_plot_3,
       plot4.flongle,plot4.r10,plot4.r9,final_plot_4)

#6 BB only: FP and FPDEL
fp.bb.r9 <- plotstats(plotdf = stats.r9[grepl(pattern = "^BB2", stats.r9$SAMPLE),],
                      groupon = "REF",title="False positives",plotcolumn="1-fp",miny=0.970,maxy=1,col = 3) #miny = 0.975
error.bb.r9 <- plotstats(plotdf = stats.r9[grepl(pattern = "^BB2", stats.r9$SAMPLE),],
                         groupon = "REF",title="Errors",plotcolumn="1-error",miny=0.94,maxy=1,col = 3) #miny = 0.95
fp.bb.r10 <- plotstats(plotdf = stats.r10[grepl(pattern = "^BB2", stats.r10$SAMPLE),],
                           groupon = "REF",title="False positives R10",plotcolumn="1-fp",miny=0.95,maxy=1,col = NA) #miny = 0.95
error.bb.r10 <- plotstats(plotdf = stats.r10[grepl(pattern = "^BB2", stats.r10$SAMPLE),],
                              groupon = "REF",title="Errors R10",plotcolumn="1-error",miny=0.925,maxy=1,col = NA) #miny = 0.925
fp.bb.F <- plotstats(plotdf = stats.flongle[which(stats.flongle$REF == "BB25"),],
                            groupon = "REF",title="False positives Flongle",plotcolumn="1-fp",miny=0.95,maxy=1,col = NA) #miny = 0.95
error.bb.F <- plotstats(plotdf = stats.flongle[which(stats.flongle$REF == "BB25"),],
                               groupon = "REF",title="Errors Flongle",plotcolumn="1-error",miny=0.925,maxy=1,col = NA) #miny = 0.925
#7 Save plots
ggsave(plot = fp.bb.r9, filename = paste(outdir,"backbone/error_BB_fp.pdf",sep=""),width =7.25, height = 2.5)
ggsave(plot = error.bb.r9, filename = paste(outdir,"backbone/error_BB_fpdel.pdf",sep=""),width =7.25, height = 2.5)
ggsave(plot = fp.bb.r10, filename = paste(outdir,"backbone/error_BB_fp_r10.pdf",sep=""),width =7.25, height = 2.5)
ggsave(plot = error.bb.r10, filename = paste(outdir,"backbone/error_BB_fpdel_r10.pdf",sep=""),width =7.25, height = 2.5)
ggsave(plot = fp.bb.F, filename = paste(outdir,"backbone/error_BB_fp_flongle.pdf",sep=""),width =7.25, height = 2.5)
ggsave(plot = error.bb.F, filename = paste(outdir,"backbone/error_BB_fpdel_flongle.pdf",sep=""),width =7.25, height = 2.5)

#8 exon 12 only
fp.e12.r9.pcr <- plotstats(plotdf= stats.r9.e12,groupon = "ANALYSIS",title="False positives",plotcolumn="1-fp",miny=0.9975,maxy=1,col = "boxplot")
ggsave(plot = fp.e12.r9.pcr, filename = paste(outdir,"insert/error_E12_fp_pcreffect.pdf",sep=""),width =7.8, height = 2.5)
remove(fp.e12.r9.pcr)
fp.e12.r9.pcr.all <- plotstats(plotdf= stats.r9.e12.all,groupon = "ANALYSIS",title="False positives",plotcolumn="1-fp",miny=0.98,maxy=1,col = "boxplot")
ggsave(plot = fp.e12.r9.pcr.all, filename = paste(outdir,"insert/error_E12_fp_pcreffect_all.pdf",sep=""),width =7.8, height = 2.5)
remove(fp.e12.r9.pcr.all)




# ---- PLOT QSCORE ----

#1 All: FP

plot1.r9<- plotstats(plotdf = stats.r9,
                     groupon = "REF",title="R9",plotcolumn="qscore_fp",miny=0,maxy=50, col = T)

plot1.r10<- plotstats(plotdf = stats.r10,
                      groupon = "REF",title="R10",plotcolumn="qscore_fp",miny=0,maxy=50, col = T)

plot1.flongle <- plotstats(plotdf = stats.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="qscore_fp",miny=0,maxy=50, col = T)

final_plot_1 <- grid.arrange(
  plot1.r9,
  plot1.r10,plot1.flongle, ncol=3)

ggsave(plot = plot1.r9, filename = paste(outdir,"all/qscore/r9_qscore_all_fp.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_1, filename = paste(outdir,"all/qscore/qscore_all_fp.pdf",sep=""),width =20, height = 8)

#2 All: FPDEL

plot2.r9<- plotstats(plotdf = stats.r9,
                     groupon = "REF",title="R9",plotcolumn="qscore_fpdel",miny=0,maxy=50, col = T)

plot2.r10<- plotstats(plotdf = stats.r10,
                      groupon = "REF",title="R10",plotcolumn="qscore_fpdel",miny=0,maxy=50, col = T)

plot2.flongle <- plotstats(plotdf = stats.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="qscore_fpdel",miny=0,maxy=50, col = T)

final_plot_2 <- grid.arrange(
  plot2.r9,
  plot2.r10,plot2.flongle, ncol=3)

ggsave(plot = plot2.r9, filename = paste(outdir,"all/qscore/r9_qscore_all_fpdel.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_2, filename = paste(outdir,"all/qscore/qscore_all_fpdel.pdf",sep=""),width =20, height = 8)

#3 All: FP per I

plot3.r9<- plotstats(plotdf = stats.perI.r9,
                     groupon = "REF",title="R9",plotcolumn="qscore_fp",miny=0,maxy=50, col = T)

plot3.r10<- plotstats(plotdf = stats.perI.r10,
                      groupon = "REF",title="R10",plotcolumn="qscore_fp",miny=0,maxy=50, col = T)

plot3.flongle <- plotstats(plotdf = stats.perI.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="qscore_fp",miny=0,maxy=50, col = T)

final_plot_3 <- grid.arrange(
  plot3.r9,
  plot3.r10,plot3.flongle, ncol=3)


ggsave(plot = plot3.r9, filename = paste(outdir,"all/qscore/r9_qscore_all_fp_perI.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_3, filename = paste(outdir,"all/qscore/qscore_all_fp_perI.pdf",sep=""),width =20, height = 8)

#4 All: FPDEL per I
plot4.r9<- plotstats(plotdf = stats.perI.r9,
                     groupon = "REF",title="R9",plotcolumn="qscore_fpdel",miny=0,maxy=50, col = T)

plot4.r10<- plotstats(plotdf = stats.perI.r10,
                      groupon = "REF",title="R10",plotcolumn="qscore_fpdel",miny=0,maxy=50, col = T)

plot4.flongle <- plotstats(plotdf = stats.perI.flongle,
                           groupon = "REF",title="Flongle",plotcolumn="qscore_fpdel",miny=0,maxy=50, col = T)

final_plot_4 <- grid.arrange(
  plot4.r9,
  plot4.r10,plot4.flongle, ncol=3)

ggsave(plot = plot4.r9, filename = paste(outdir,"all/qscore/r9_qscore_all_fpdel_perI.pdf",sep=""),width =12, height = 8)
ggsave(plot = final_plot_4, filename = paste(outdir,"all/qscore/qscore_all_fpdel_perI.pdf",sep=""),width =20, height = 8)

#5 Cleanup
remove(plot1.flongle,plot1.r10,plot1.r9,final_plot_1,
       plot2.flongle,plot2.r10,plot2.r9,final_plot_2,
       plot3.flongle,plot3.r10,plot3.r9,final_plot_3,
       plot4.flongle,plot4.r10,plot4.r9,final_plot_4)




# ---- PLOT COVERAGE ----

#1 Get max coverage (to get the maxy and to calculate how many 'breaks')
max(stats.r9[grepl(pattern = "^BB2", stats.r9$SAMPLE),]$COV)
max(stats.r10[grepl(pattern = "^BB2", stats.r10$SAMPLE),]$COV)
max(stats.flongle[which(stats.flongle$REF == "BB25"),]$COV)
max(stats.r9.e12$COV)

#2 Plot coverage
cov.bb.r9 <- plotcov(plotdf = stats.r9[grepl(pattern = "^BB2", stats.r9$SAMPLE),],
                     maxy = 600000, pertype = NA, breaks = 6)
cov.bb.r10 <- plotcov(plotdf = stats.r10[grepl(pattern = "^BB2", stats.r10$SAMPLE),],
                      maxy =150000, pertype = NA, breaks = 5)
cov.bb.F <- plotcov(plotdf = stats.flongle[which(stats.flongle$REF == "BB25"),],
                    maxy = 20000, pertype = NA, breaks = 5)

ggsave(plot = cov.bb.r9, filename = paste(outdir,"backbone/coverage_BB.pdf",sep=""),width =6.5, height = 2.5)
ggsave(plot = cov.bb.r10, filename = paste(outdir,"backbone/coverage_BB_r10.pdf",sep=""),width =6.5, height = 2.5)
ggsave(plot = cov.bb.F, filename = paste(outdir,"backbone/coverage_BB_flongle.pdf",sep=""),width =6.5, height = 2.5)

#3 Clean up
remove(cov.bb.F,cov.bb.r10,cov.bb.r9)




# ---- CALCULATIONS FOR MANUSCRIPT ----

# FP rate at 1, 5 and 10 repeats
round(mean(rbind(stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB22"),],
                 stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB24"),],
                 stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB25"),])$fp),digits = 4)
paste("(",paste(round(min(rbind(stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB22"),],
           stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB24"),],
           stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB25"),])$fp),digits = 4),
      round(max(rbind(stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB22"),],
           stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB24"),],
           stats.r9[which(stats.r9$TYPE == 1 & stats.r9$REF == "BB25"),])$fp),digits = 4),sep="-"),")",sep="")

round(mean(rbind(stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB22"),],
                 stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB24"),],
                 stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB25"),])$fp),digits = 4)
paste("(",paste(round(min(rbind(stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB22"),],
          stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB24"),],
          stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB25"),])$fp),digits = 4),
          round(max(rbind(stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB22"),],
          stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB24"),],
          stats.r9[which(stats.r9$TYPE == 5 & stats.r9$REF == "BB25"),])$fp),digits = 4),sep="-"),")",sep="")

round(mean(rbind(stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB22"),],
                 stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB24"),],
                 stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB25"),])$fp),digits = 4)
paste("(",paste(round(min(rbind(stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB22"),],
          stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB24"),],
          stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB25"),])$fp),digits = 4),
          round(max(rbind(stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB22"),],
          stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB24"),],
          stats.r9[which(stats.r9$TYPE == 10 & stats.r9$REF == "BB25"),])$fp),digits = 4),sep="-"),")",sep="")

# PCR-free
temp <- stats.r9.e12.all[which(as.numeric(as.character(stats.r9.e12.all$TYPE)) == 1),]
(mean(temp[which(temp$ANALYSIS == "PCR"),]$fp)) / (mean(temp[which(temp$ANALYSIS == "PCR-free"),]$fp))
remove(temp)
