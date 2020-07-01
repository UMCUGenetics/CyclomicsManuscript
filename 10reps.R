# @Date created: 18 February 2020
# @Author: Myrthe Jager
# @Description: Analyses of quality of files with at least 10 repeats cutoff in pbdagcon
# @Abbreviations: BB = backbone/PJET, I = insert, FP = false positive, FPDEL = FP + deletions
# @Version: 1 July 2020




# ---- GET STARTED ----

#1 Load required packages
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)

#2 Define output directory
dir = "/.../Cyclomics_manuscript/"
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/10reps/",sep = "")
datadir = paste(indir,"RCA/",sep ="")
remove(dir)

#3 Functions
#3A See separate R script for general functions: functions.R




# ---- GET REFERENCE SEQUENCES ----

#1 BB
sequences_backbones <- read.table(paste(indir,"sequences/backbones.txt",sep=""), header=TRUE,sep="\t")

#2 TP53
sequences_tp53 <- read.table(paste(indir,"sequences/tp53.txt",sep=""), header=TRUE,sep="\t")
#2A Get all amplicons analyzed
all_amplicons_tp53 <- unlist(strsplit(as.character(sequences_tp53[which(sequences_tp53$name != "ampliconS0"),]$name), split = "con"))
all_amplicons_tp53[grep(all_amplicons_tp53, pattern = "mpl")] <- c(",")
all_amplicons_tp53 <- paste(all_amplicons_tp53[-1],collapse="")




#---- LOAD METADATA ----

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




# ---- LOAD AND PREPARE FILES ----

#1 Get sample info from files using metadata
#1A Get info
sample.info.r9.all <- rbind(getsampleinfo.10reps(metadata = metadata.r9), getsampleinfo.forrev(metadata = metadata.r9))
sample.info.r10.all <- rbind(getsampleinfo.10reps(metadata = metadata.r10), getsampleinfo.forrev(metadata = metadata.r10))
sample.info.flongle.all <- rbind(getsampleinfo.10reps(metadata = metadata.flongle), getsampleinfo.forrev(metadata = metadata.flongle))
#1B Check if numbers match (sample info = ALL&FOR&REV for BB&I, so divide by 2*3 = 6)
nrow(sample.info.r9.all)/6 == nrow(metadata.r9)
nrow(sample.info.r10.all)/6 == nrow(metadata.r10)
nrow(sample.info.flongle.all)/6 == nrow(metadata.flongle)
#1C Remove S0 files from flongle data
sample.info.flongle.all <- sample.info.flongle.all[-which(sample.info.flongle.all$position == "TP53 amplicon S0"),]

#2 Read files
myfiles.r9.all = lapply(as.character(sample.info.r9.all$file), read.delim)
myfiles.r10.all = lapply(as.character(sample.info.r10.all$file), read.delim)
myfiles.flongle.all = lapply(as.character(sample.info.flongle.all$file), read.delim)

#3 Change names of files 
names(myfiles.r9.all) = paste(sample.info.r9.all$samplename,sample.info.r9.all$type,sample.info.r9.all$position,sample.info.r9.all$analysis,sep="__")
names(myfiles.r10.all) = paste(sample.info.r10.all$samplename,sample.info.r10.all$type,sample.info.r10.all$position,sample.info.r10.all$analysis,sep="__")
names(myfiles.flongle.all) = paste(sample.info.flongle.all$samplename,sample.info.flongle.all$type,sample.info.flongle.all$position,sample.info.flongle.all$analysis,sep="__")

#4 Prepare files for stat calculations
myfiles.r9.all <- preparefiles(filelist = myfiles.r9.all,sampleinfo = sample.info.r9.all)
myfiles.r10.all <- preparefiles(filelist = myfiles.r10.all,sampleinfo = sample.info.r10.all)
myfiles.flongle.all <- preparefiles(filelist = myfiles.flongle.all,sampleinfo = sample.info.flongle.all)

#5 Get stats per table
stats.r9.all <- getstats.pertable(filelist = myfiles.r9.all,sampleinfo = sample.info.r9.all, mincoverage = 100,blackposlist = sample.info.r9.all$mutation,perinsert = NA)




# ---- CALCULATE STATS PER POSITION----

#1 Get error rates and add them to the files, per position
myfiles.r9.all <- getstats.perline(filelist = myfiles.r9.all, mincoverage = 100, fakeerror = 0.1, blackposlist = sample.info.r9.all$mutation)
myfiles.r10.all <- getstats.perline(filelist =myfiles.r10.all, mincoverage = 100, fakeerror = 0.1, blackposlist = sample.info.r10.all$mutation)
myfiles.flongle.all <- getstats.perline(filelist =myfiles.flongle.all, mincoverage = 100, fakeerror = 0.1, blackposlist = sample.info.flongle.all$mutation)

#2 Get stats per position for R9
statlist.r9.all <- c(list(
  #2A BB22
  df.stats.r9.all.22 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB22" & sample.info.r9.all$analysis =="ALL")], stat = "mean"),
  df.stats.r9.fwd.22 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB22" & sample.info.r9.all$analysis =="FOR")], stat = "mean"),
  df.stats.r9.rev.22 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB22" & sample.info.r9.all$analysis =="REV")], stat = "mean"),
  #2B BB24
  df.stats.r9.all.24 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB24" & sample.info.r9.all$analysis =="ALL")], stat = "mean"),
  df.stats.r9.fwd.24 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB24" & sample.info.r9.all$analysis =="FOR")], stat = "mean"),
  df.stats.r9.rev.24 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB24" & sample.info.r9.all$analysis =="REV")], stat = "mean"),
  #2C BB25
  df.stats.r9.all.25 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB25" & sample.info.r9.all$analysis =="ALL")], stat = "mean"),
  df.stats.r9.fwd.25 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB25" & sample.info.r9.all$analysis =="FOR")], stat = "mean"),
  df.stats.r9.rev.25 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$position =="BB25" & sample.info.r9.all$analysis =="REV")], stat = "mean"),
  #2D insert: amplicon 12
  df.stats.r9.all.e12 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$analysis == "ALL" & grepl(pattern = "12",sample.info.r9.all$position))], stat = "mean"),
  df.stats.r9.fwd.e12 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$analysis == "FOR" & grepl(pattern = "12",sample.info.r9.all$position))], stat = "mean"),
  df.stats.r9.rev.e12 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$analysis == "REV" & grepl(pattern = "12",sample.info.r9.all$position))], stat = "mean"),
  #2E insert: amplicon 19
  df.stats.r9.all.e19 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$analysis == "ALL" & grepl(pattern = "19",sample.info.r9.all$position))], stat = "mean"),
  df.stats.r9.fwd.e19 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$analysis == "FOR" & grepl(pattern = "19",sample.info.r9.all$position))], stat = "mean"),
  df.stats.r9.rev.e19 = statsperpos(filelist = myfiles.r9.all[which(sample.info.r9.all$analysis == "REV" & grepl(pattern = "19",sample.info.r9.all$position))], stat = "mean")
))

#3 Get stats per position for R10
statlist.r10.all <- c(list(
  #3A BB24
  df.stats.r10.all.24 = statsperpos(filelist = myfiles.r10.all[which(sample.info.r10.all$position =="BB24" & sample.info.r10.all$analysis =="ALL")], stat = "mean"),
  df.stats.r10.fwd.24 = statsperpos(filelist = myfiles.r10.all[which(sample.info.r10.all$position =="BB24" & sample.info.r10.all$analysis =="FOR")], stat = "mean"),
  df.stats.r10.rev.24 = statsperpos(filelist = myfiles.r10.all[which(sample.info.r10.all$position =="BB24" & sample.info.r10.all$analysis =="REV")], stat = "mean"),
  #3B insert: amplicon 12
  df.stats.r10.all.e12 = statsperpos(filelist = myfiles.r10.all[which(sample.info.r10.all$analysis == "ALL" & grepl(pattern = "12",sample.info.r10.all$position))], stat = "mean"),
  df.stats.r10.fwd.e12 = statsperpos(filelist = myfiles.r10.all[which(sample.info.r10.all$analysis == "FOR" & grepl(pattern = "12",sample.info.r10.all$position))], stat = "mean"),
  df.stats.r10.rev.e12 = statsperpos(filelist = myfiles.r10.all[which(sample.info.r10.all$analysis == "REV" & grepl(pattern = "12",sample.info.r10.all$position))], stat = "mean")
))

#4 Get stats per position for flongle
statlist.flongle.all <- c(list(
  #4A BB25
  df.stats.flongle.all.25 = statsperpos(filelist = myfiles.flongle.all[which(sample.info.flongle.all$position =="BB25" & sample.info.flongle.all$analysis =="ALL")], stat = "mean"),
  df.stats.flongle.fwd.25 = statsperpos(filelist = myfiles.flongle.all[which(sample.info.flongle.all$position =="BB25" & sample.info.flongle.all$analysis =="FOR")], stat = "mean"),
  df.stats.flongle.rev.25 = statsperpos(filelist = myfiles.flongle.all[which(sample.info.flongle.all$position =="BB25" & sample.info.flongle.all$analysis =="REV")], stat = "mean"),
  #4B insert: amplicon 19
  df.stats.flongle.all.e19 = statsperpos(filelist = myfiles.flongle.all[which(sample.info.flongle.all$analysis == "ALL" & grepl(pattern = "19",sample.info.flongle.all$position))], stat = "mean"),
  df.stats.flongle.fwd.e19 = statsperpos(filelist = myfiles.flongle.all[which(sample.info.flongle.all$analysis == "FOR" & grepl(pattern = "19",sample.info.flongle.all$position))], stat = "mean"),
  df.stats.flongle.rev.e19 = statsperpos(filelist = myfiles.flongle.all[which(sample.info.flongle.all$analysis == "REV" & grepl(pattern = "19",sample.info.flongle.all$position))], stat = "mean")
))




# ---- FOR REV CORRECTION ----

#1 Get coordinates for Amplicon 12 and 19
coords.e12 <- unlist(strsplit(split = "-", tail(unlist(strsplit(as.character(sequences_tp53[which(sequences_tp53$name == "amplicon12"),]$genomic.coordinates),split=":")),n=1)))
coords.e19 <- unlist(strsplit(split = "-", tail(unlist(strsplit(as.character(sequences_tp53[which(sequences_tp53$name == "amplicon19"),]$genomic.coordinates),split=":")),n=1)))

#2 Join ALL, FWD and REV tables
statlist.all <- c(list(
  df.stats.r9.22 = combine.stats.fwdrev(all = statlist.r9.all[[grep(pattern = "df.stats.r9.all.22",names(statlist.r9.all))]], 
                                         fwd = statlist.r9.all[[grep(pattern = "df.stats.r9.fwd.22",names(statlist.r9.all))]], 
                                         rev = statlist.r9.all[[grep(pattern = "df.stats.r9.rev.22",names(statlist.r9.all))]],
                                         coords =NA, selectcol = "mean_error_fp"),
  df.stats.r9.24 = combine.stats.fwdrev(all = statlist.r9.all[[grep(pattern = "df.stats.r9.all.24",names(statlist.r9.all))]], 
                                         fwd = statlist.r9.all[[grep(pattern = "df.stats.r9.fwd.24",names(statlist.r9.all))]], 
                                         rev = statlist.r9.all[[grep(pattern = "df.stats.r9.rev.24",names(statlist.r9.all))]],
                                         coords =NA, selectcol = "mean_error_fp"),
  df.stats.r9.25 = combine.stats.fwdrev(all = statlist.r9.all[[grep(pattern = "df.stats.r9.all.25",names(statlist.r9.all))]], 
                                         fwd = statlist.r9.all[[grep(pattern = "df.stats.r9.fwd.25",names(statlist.r9.all))]], 
                                         rev = statlist.r9.all[[grep(pattern = "df.stats.r9.rev.25",names(statlist.r9.all))]],
                                         coords =NA, selectcol = "mean_error_fp"),
  df.stats.r9.e12 = combine.stats.fwdrev(all = statlist.r9.all[[grep(pattern = "df.stats.r9.all.e12",names(statlist.r9.all))]], 
                                        fwd = statlist.r9.all[[grep(pattern = "df.stats.r9.fwd.e12",names(statlist.r9.all))]], 
                                        rev = statlist.r9.all[[grep(pattern = "df.stats.r9.rev.e12",names(statlist.r9.all))]],
                                        coords = coords.e12, selectcol = "mean_error_fp"),
  df.stats.r9.e19 = combine.stats.fwdrev(all = statlist.r9.all[[grep(pattern = "df.stats.r9.all.e19",names(statlist.r9.all))]], 
                                        fwd = statlist.r9.all[[grep(pattern = "df.stats.r9.fwd.e19",names(statlist.r9.all))]], 
                                        rev = statlist.r9.all[[grep(pattern = "df.stats.r9.rev.e19",names(statlist.r9.all))]],
                                        coords = coords.e19, selectcol = "mean_error_fp"),
  df.stats.r10.24 = combine.stats.fwdrev(all = statlist.r10.all[[grep(pattern = "df.stats.r10.all.24",names(statlist.r10.all))]], 
                                        fwd = statlist.r10.all[[grep(pattern = "df.stats.r10.fwd.24",names(statlist.r10.all))]], 
                                        rev = statlist.r10.all[[grep(pattern = "df.stats.r10.rev.24",names(statlist.r10.all))]],
                                        coords =NA, selectcol = "mean_error_fp"),
  df.stats.r10.e12 = combine.stats.fwdrev(all = statlist.r10.all[[grep(pattern = "df.stats.r10.all.e12",names(statlist.r10.all))]], 
                                         fwd = statlist.r10.all[[grep(pattern = "df.stats.r10.fwd.e12",names(statlist.r10.all))]], 
                                         rev = statlist.r10.all[[grep(pattern = "df.stats.r10.rev.e12",names(statlist.r10.all))]],
                                         coords = coords.e12, selectcol = "mean_error_fp"),
  df.stats.flongle.25 = combine.stats.fwdrev(all = statlist.flongle.all[[grep(pattern = "df.stats.flongle.all.25",names(statlist.flongle.all))]], 
                                        fwd = statlist.flongle.all[[grep(pattern = "df.stats.flongle.fwd.25",names(statlist.flongle.all))]], 
                                        rev = statlist.flongle.all[[grep(pattern = "df.stats.flongle.rev.25",names(statlist.flongle.all))]],
                                        coords =NA, selectcol = "mean_error_fp"),
  df.stats.flongle.e19 = combine.stats.fwdrev(all = statlist.flongle.all[[grep(pattern = "df.stats.flongle.all.e19",names(statlist.flongle.all))]], 
                                             fwd = statlist.flongle.all[[grep(pattern = "df.stats.flongle.fwd.e19",names(statlist.flongle.all))]], 
                                             rev = statlist.flongle.all[[grep(pattern = "df.stats.flongle.rev.e19",names(statlist.flongle.all))]],
                                             coords = coords.e19, selectcol = "mean_error_fp")
  
))

#3 Calculate change if you take FOR or REV
for(i in 1:length(statlist.all)) {
  df <- statlist.all[[i]]
  df <- calculate.fwd.rev.change(df = df)
  statlist.all[[i]] <- df
  remove(df)
}
remove(i)

#4 Get positions with a change
#4A Define cutoff
cutoff = -1/1000
#4B Create empty list
posdflist.all <- list()
#4C Get positions and add to list
for(i in 1:length(statlist.all)) {
  df <- get.fwdrev.positions(df = statlist.all[[i]],cutoff)
  posdflist.all <- c(posdflist.all,list(df))
  names(posdflist.all) = c(names(posdflist.all)[-i],
                           paste("posdf.",tail(unlist(strsplit(names(statlist.all)[i], split = "stats.")),n=1),sep=""))
  remove(df)
}
remove(i)

#5 Save positions of posdf
pos.df <- data.frame()
for(i in 1:length(posdflist.all)) {
  tempdf <- posdflist.all[[i]]
  tempdf$flowcell <- unlist(strsplit(names(posdflist.all)[i],split = "[.]"))[2]
  tempdf$REF <- as.character(tempdf$REF)
  pos.df <- rbind(pos.df,tempdf)
  remove(tempdf)
}
remove(i)
write.table(pos.df,file = paste(indir,"/forrevpositions.txt",sep=""),row.names = F,sep="\t")

#6 Correct the filelists
myfiles.r9 <- get.forrev.corrected.filelist (filelist.all = myfiles.r9.all[sample.info.r9.all$analysis =="ALL"],
                                              filelist.fwd = myfiles.r9.all[sample.info.r9.all$analysis =="FOR"],
                                              filelist.rev = myfiles.r9.all[sample.info.r9.all$analysis =="REV"],
                                              posdf = pos.df[which(pos.df$flowcell == "r9"),]
)
myfiles.r10 <- get.forrev.corrected.filelist (filelist.all = myfiles.r10.all[sample.info.r10.all$analysis =="ALL"],
                                               filelist.fwd = myfiles.r10.all[sample.info.r10.all$analysis =="FOR"],
                                               filelist.rev = myfiles.r10.all[sample.info.r10.all$analysis =="REV"],
                                               posdf = pos.df[which(pos.df$flowcell == "r10"),]
)
myfiles.flongle <- get.forrev.corrected.filelist (filelist.all = myfiles.flongle.all[sample.info.flongle.all$analysis =="ALL"],
                                              filelist.fwd = myfiles.flongle.all[sample.info.flongle.all$analysis =="FOR"],
                                              filelist.rev = myfiles.flongle.all[sample.info.flongle.all$analysis =="REV"],
                                              posdf = pos.df[which(pos.df$flowcell == "flongle"),]
)

#7 Create a new sampleinfo table for the corrected filelist
sample.info.r9 <- sample.info.r9.all[which(sample.info.r9.all$analysis == "ALL"),]
sample.info.r9$analysis <- "CORRECTED"
sample.info.r10 <- sample.info.r10.all[which(sample.info.r10.all$analysis == "ALL"),]
sample.info.r10$analysis <- "CORRECTED"
sample.info.flongle <- sample.info.flongle.all[which(sample.info.flongle.all$analysis == "ALL"),]
sample.info.flongle$analysis <- "CORRECTED"

#8 Add names filelists
names(myfiles.r9) = paste(sample.info.r9$samplename,sample.info.r9$type,sample.info.r9$position,sample.info.r9$analysis,sep="__")
names(myfiles.r10) = paste(sample.info.r10$samplename,sample.info.r10$type,sample.info.r10$position,sample.info.r10$analysis,sep="__")
names(myfiles.flongle) = paste(sample.info.flongle$samplename,sample.info.flongle$type,sample.info.flongle$position,sample.info.flongle$analysis,sep="__")

#9 Calculate effect of forrev correction in subset of data 
subsetcorrected.stats.24 <- forrev.correct.fp.subsetdecision(filelist.all = myfiles.r9.all[which(sample.info.r9.all$position =="BB24" & sample.info.r9.all$analysis =="ALL")],
                                                             filelist.fwd = myfiles.r9.all[which(sample.info.r9.all$position =="BB24" & sample.info.r9.all$analysis =="FOR")],
                                                             filelist.rev = myfiles.r9.all[which(sample.info.r9.all$position =="BB24" & sample.info.r9.all$analysis =="REV")],
                                                             coords = NA,
                                                             cutoff = -(1/1000),
                                                             sampleinfo = sample.info.r9.all[which(sample.info.r9.all$position =="BB24"),])

#10 Clean up Rdata
remove(pos.df,posdflist.all)
remove(statlist.r9.all,statlist.r10.all,statlist.flongle.all)
remove(myfiles.r9.all,myfiles.r10.all,myfiles.flongle.all)
remove(sample.info.r9.all,sample.info.r10.all,sample.info.flongle.all)
remove(coords.e12,coords.e19)




# ---- PLOT FORREV CORRECTION EFFECT ----

#1 Improvement in BB24
#1A Plot
plot.improvement <- plot.forrev(plotdf = statlist.all[[which(names(statlist.all) == "df.stats.r9.24")]])
#1B Save
ggsave(paste(outdir,"forrev_improvement_BB24_r9.pdf",sep=""),plot = plot.improvement, width = 6.5, height = 2.666667)
#1C Remove
remove(plot.improvement)

#2 Subsetting plots
#2A Difference in fp rate per base
diff1 <- grid.arrange(plotstatsperpos(plotdf = subsetcorrected.stats.24[[grep(pattern = "OLD_stats",names(subsetcorrected.stats.24))]],title = "uncorrected",plotcolumn = "mean_error_fp",miny=0,maxy=0.01, color = T),
                      plotstatsperpos(plotdf = subsetcorrected.stats.24[[grep(pattern = "CORR_stats",names(subsetcorrected.stats.24))]],title = "corrected",plotcolumn = "mean_error_fp",miny=0,maxy=0.01, color = T),
                      nrow = 2)
#2B Calculate difference in basetype
basetypediff.subset.24 <- prepare.basetypediff.forplotting(basetypediff = subsetcorrected.stats.24[[grep(pattern = "DIFF_basetype",names(subsetcorrected.stats.24))]],
                                                           fporfpdel = "fp")
#2C Generate separate table for data points
basetypediff.subset.24.datapoints <- subsetcorrected.stats.24[[grep(pattern = "DIFF_basetype",names(subsetcorrected.stats.24))]]
basetypediff.subset.24.datapoints <- cbind(basetypediff.subset.24.datapoints[,-grep(pattern = "fp",colnames(basetypediff.subset.24.datapoints))],
                                           basetypediff.subset.24.datapoints[,grep(pattern = "perc_fp_",colnames(basetypediff.subset.24.datapoints))])
basetypediff.subset.24.datapoints <- melt(basetypediff.subset.24.datapoints)
colnames(basetypediff.subset.24.datapoints) <- c("sample","REF","type","value")
basetypediff.subset.24.datapoints$type <- gsub("perc_fp_gold","< 0.1%",basetypediff.subset.24.datapoints$type)
basetypediff.subset.24.datapoints$type <- gsub("perc_fp_silver","0.1-1.0%",basetypediff.subset.24.datapoints$type)
basetypediff.subset.24.datapoints$type <- gsub("perc_fp_bronze","> 1.0%",basetypediff.subset.24.datapoints$type)
#2D Plot Difference in basetype
diff2 <- grid.arrange(plotbasetype(basetype = subsetcorrected.stats.24[[grep(pattern = "OLD_basetype",names(subsetcorrected.stats.24))]],title = "uncorrected",plottype = "fp"),
                      plotbasetype(basetype = subsetcorrected.stats.24[[grep(pattern = "CORR_basetype",names(subsetcorrected.stats.24))]],title = "corrected",plottype = "fp"),
                      ggplot(basetypediff.subset.24, aes(x = type, y = mean, fill = type)) +
                        geom_bar(stat="identity", width = 0.9) +
                        geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.3, size=0.2) +
                        ylab("Positions (%)") +
                        ggtitle("Difference") +
                        xlab("") +
                        coord_cartesian(ylim=c(-5, 5)) +
                        scale_y_continuous(expand = c(0, 0), breaks=seq(-5,5, 1))+
                        theme_bw() +
                        theme(panel.grid.minor = element_blank()) +
                        theme(legend.position = "none") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                        #scale_fill_manual(values=c("goldenrod1","#C0C0C0","#cd7f32"))
                        scale_fill_manual(values=c("gray38","gray38","gray38"))+
                                            geom_jitter(data = basetypediff.subset.24.datapoints, aes(x = type, y = value), size = 0.25), 
                      nrow = 3)

#2D Save plots
ggsave(paste(outdir,"subset_fpratediff_BB24_r9.pdf",sep=""),plot = diff1, width = 6.3, height = 5.33334)
ggsave(paste(outdir,"basetype/subset_basetypediff_BB24_r9.pdf",sep=""),plot = diff2, width = 1.5, height = 7.4)
#2E Clean up
remove(diff1,diff2,basetypediff.subset.24,basetypediff.subset.24.datapoints)




# ---- CALCULATE STATS PER POSITION----

#1 BB: r9 files
df.stats.BB22 <- statsperpos(filelist = myfiles.r9[sample.info.r9$position == "BB22"], stat = "mean")
df.stats.BB24 <- statsperpos(filelist = myfiles.r9[sample.info.r9$position == "BB24"], stat = "mean")
df.stats.BB25 <- statsperpos(filelist = myfiles.r9[sample.info.r9$position == "BB25"], stat = "mean")

#2 BB: flongle files
df.stats.F.BB25 <- statsperpos(filelist = myfiles.flongle[sample.info.flongle$position == "BB25"], stat = "mean")

#3 BB: r10 files
df.stats.R10.BB24 <- statsperpos(filelist = myfiles.r10[sample.info.r10$position == "BB24"], stat = "mean")

#4 I: r9 files
df.stats.INS <- statsperpos(filelist = myfiles.r9[grep(pattern = "TP", sample.info.r9$position)], stat = "mean")

#5 I: flongle files
df.stats.F.INS <- statsperpos(filelist = myfiles.flongle[grep(pattern = "TP", sample.info.flongle$position)], stat = "mean")

#6 I: r10 files
df.stats.R10.INS <- statsperpos(filelist = myfiles.r10[grep(pattern = "TP", sample.info.r10$position)], stat = "mean")

#7 r9 Per table
stats.r9 <- getstats.pertable(filelist = myfiles.r9,sampleinfo = sample.info.r9, mincoverage = 100,blackposlist = sample.info.r9$mutation,perinsert = NA)
stats.r10 <- getstats.pertable(filelist = myfiles.r10,sampleinfo = sample.info.r10, mincoverage = 100,blackposlist = sample.info.r10$mutation,perinsert = NA)
stats.flongle <- getstats.pertable(filelist = myfiles.flongle,sampleinfo = sample.info.flongle, mincoverage = 100,blackposlist = sample.info.flongle$mutation,perinsert = NA)




# ---- PLOT STATS PER POSITION ----

#1 Error rate per position: Base subs
#1A R9 
plot.fp.BB22 <- plotstatsperpos(plotdf = df.stats.BB22, title = "", plotcolumn = "mean_error_fp",miny=0,maxy=0.005, color = T)
plot.fp.BB24 <- plotstatsperpos(plotdf = df.stats.BB24, title = "", plotcolumn = "mean_error_fp",miny=0,maxy=0.005, color = T)
plot.fp.BB25 <- plotstatsperpos(plotdf = df.stats.BB25, title = "", plotcolumn = "mean_error_fp",miny=0,maxy=0.005, color = T)
#1B R10 BB25 with base colors
plot.fp.R10.BB24 <- plotstatsperpos(plotdf = df.stats.R10.BB24,title ="BB24 (R10)", plotcolumn = "mean_error_fp",miny=0,maxy=0.009, color = T)
#1C Flongle
plot.fp.F.BB25 <- plotstatsperpos(plotdf = df.stats.F.BB25,title ="BB25 (Flongle)", plotcolumn = "mean_error_fp",miny=0,maxy=0.009, color = T)
#1D R9 combined
plot.fp <- grid.arrange(plot.fp.BB22,plot.fp.BB24,plot.fp.BB25,nrow=3)
#1E Save plots
ggsave(paste(outdir,"extra/fp_BB22_r9.pdf",sep=""),plot = plot.fp.BB22, width = 9, height = 2)
ggsave(paste(outdir,"extra/fp_BB24_r9.pdf",sep=""),plot = plot.fp.BB24, width = 9, height = 2)
ggsave(paste(outdir,"extra/fp_BB25_r9.pdf",sep=""),plot = plot.fp.BB25, width = 9, height = 2)
ggsave(paste(outdir,"fp_BB24_r10.pdf",sep=""),plot = plot.fp.R10.BB24, width = 6.3, height = 2.666667)
ggsave(paste(outdir,"fp_BB25_flongle.pdf",sep=""),plot = plot.fp.F.BB25, width = 6.3, height = 2.666667)
ggsave(paste(outdir,"fp_BBall_r9.pdf",sep=""),plot = plot.fp, width = 6.3, height = 8)

#2 Error rate per position: base subs + deletions
#2A R9 
plot.error.BB22 <- plotstatsperpos(plotdf = df.stats.BB22, title = "", plotcolumn = "mean_error_fpdel",miny=0,maxy=0.15, color = T)
plot.error.BB24 <- plotstatsperpos(plotdf = df.stats.BB24, title = "", plotcolumn = "mean_error_fpdel",miny=0,maxy=0.15, color = T)
plot.error.BB25 <- plotstatsperpos(plotdf = df.stats.BB25, title = "", plotcolumn = "mean_error_fpdel",miny=0,maxy=0.15, color = T)
#2B R10 
plot.error.R10.BB24 <- plotstatsperpos(plotdf = df.stats.R10.BB24,title ="BB24 (R10)", plotcolumn = "mean_error_fpdel",miny=0,maxy=0.4, color = T)
#2C R10 
plot.error.F.BB25 <- plotstatsperpos(plotdf = df.stats.F.BB25,title ="BB25 (Flongle)", plotcolumn = "mean_error_fpdel",miny=0,maxy=0.2, color = T)
#2D R9 combined
plot.error <- grid.arrange(plot.error.BB22,plot.error.BB24,plot.error.BB25,nrow=3)
#2E Save plots
ggsave(paste(outdir,"extra/error_BB22_r9.pdf",sep=""),plot = plot.error.BB22, width = 9, height = 2)
ggsave(paste(outdir,"extra/error_BB24_r9.pdf",sep=""),plot = plot.error.BB24, width = 9, height = 2)
ggsave(paste(outdir,"extra/error_BB25_r9.pdf",sep=""),plot = plot.error.BB25, width = 9, height = 2)
ggsave(paste(outdir,"error_BB24_r10.pdf",sep=""),plot = plot.error.R10.BB24, width = 6.3, height = 2.666667)
ggsave(paste(outdir,"error_BB25_flongle.pdf",sep=""),plot = plot.error.F.BB25, width = 6.3, height = 2.666667)
ggsave(paste(outdir,"error_BBall_r9.pdf",sep=""),plot = plot.error, width = 6.3, height = 8)

#3 Remove all plots
remove(plot.fp.BB22,plot.fp.BB24,plot.fp.BB25,plot.fp.R10.BB24,plot.fp.F.BB25, plot.fp)
remove(plot.error.BB22,plot.error.BB24,plot.error.BB25,plot.error.R10.BB24,plot.error.F.BB25, plot.error)




# ---- GET GOLD/SILVER/BRONZE BASES ----

#1 Get basetype for BB
#1A BB: r9 files
basetype.BB22 <- get.basetype.pertable(filelist = myfiles.r9[sample.info.r9$position == "BB22"],
                                       blackposlist = sample.info.r9[which(sample.info.r9$position == "BB22"),]$mutation,
                                       mincoverage=100)
basetype.BB24 <- get.basetype.pertable(filelist = myfiles.r9[sample.info.r9$position == "BB24"],
                                       blackposlist = sample.info.r9[which(sample.info.r9$position == "BB24"),]$mutation,
                                       mincoverage=100)
basetype.BB25 <- get.basetype.pertable(filelist = myfiles.r9[sample.info.r9$position == "BB25"],
                                       blackposlist = sample.info.r9[which(sample.info.r9$position == "BB25"),]$mutation,
                                       mincoverage=100)
#1B BB: flongle files
basetype.F.BB25 <- get.basetype.pertable(filelist = myfiles.flongle[sample.info.flongle$position == "BB25"],
                                         blackposlist = sample.info.flongle[which(sample.info.flongle$position == "BB25"),]$mutation,
                                         mincoverage=100)
#1C BB: r10 files
basetype.R10.BB24 <- get.basetype.pertable(filelist = myfiles.r10[sample.info.r10$position == "BB24"],
                                           blackposlist = sample.info.r10[which(sample.info.r10$position == "BB24"),]$mutation,
                                           mincoverage=100)




# ---- PLOT BASETYPE ----
# FP
plot.basetype.fp.BB22 <- plotbasetype(basetype = basetype.BB22,title = "BB22",plottype = "fp")
plot.basetype.fp.BB24 <- plotbasetype(basetype = basetype.BB24,title = "BB24",plottype = "fp")
plot.basetype.fp.BB25 <- plotbasetype(basetype = basetype.BB25,title = "BB25",plottype = "fp")
plot.basetype.fp <- grid.arrange(plot.basetype.fp.BB22,plot.basetype.fp.BB24,plot.basetype.fp.BB25, nrow = 3)
# FPDEL
plot.basetype.error.BB22 <- plotbasetype(basetype = basetype.BB22,title = "BB22",plottype = "fpdel")
plot.basetype.error.BB24 <- plotbasetype(basetype = basetype.BB24,title = "BB24",plottype = "fpdel")
plot.basetype.error.BB25 <- plotbasetype(basetype = basetype.BB25,title = "BB25",plottype = "fpdel")
plot.basetype.error <- grid.arrange(plot.basetype.error.BB22,plot.basetype.error.BB24,plot.basetype.error.BB25, nrow = 3)
# R10
plot.basetype.fp.R10.BB24 <- plotbasetype(basetype = basetype.R10.BB24,title = "BB24",plottype = "fp")
plot.basetype.error.R10.BB24 <- plotbasetype(basetype = basetype.R10.BB24,title = "BB24",plottype = "fpdel")
plot.basetype.R10 <- grid.arrange(plot.basetype.fp.R10.BB24,plot.basetype.error.R10.BB24, nrow = 2)
# Flongle
plot.basetype.fp.F.BB25 <- plotbasetype(basetype = basetype.F.BB25,title = "BB25",plottype = "fp")
plot.basetype.error.F.BB25 <- plotbasetype(basetype = basetype.F.BB25,title = "BB25",plottype = "fpdel")
plot.basetype.F <- grid.arrange(plot.basetype.fp.F.BB25,plot.basetype.error.F.BB25, nrow = 2)
#Save
ggsave(paste(outdir,"basetype/basetype_fp.pdf",sep=""),plot = plot.basetype.fp, width = 1.5, height = 7.4)
ggsave(paste(outdir,"basetype/basetype_error.pdf",sep=""),plot = plot.basetype.error, width = 1.5, height = 7.4)
ggsave(paste(outdir,"basetype/basetype_R10.pdf",sep=""),plot = plot.basetype.R10, width = 1.5, height = 4.9333333)
ggsave(paste(outdir,"basetype/basetype_Flongle.pdf",sep=""),plot = plot.basetype.F, width = 1.5, height = 4.9333333)

remove(plot.basetype.error,plot.basetype.error.BB22,plot.basetype.error.BB24,plot.basetype.error.BB25)
remove(plot.basetype.fp,plot.basetype.fp.BB22,plot.basetype.fp.BB24,plot.basetype.fp.BB25)
remove(plot.basetype.R10,plot.basetype.fp.R10.BB24,plot.basetype.error.R10.BB24)
remove(plot.basetype.F,plot.basetype.fp.F.BB25,plot.basetype.error.F.BB25)




# ---- GET STATS PER TABLE ----

#1 Get error rates for amplicon 12, except the mutated bases in the triple mutant
#1A First get a blacklist and exclude the triple mutant positions from all files (to make it comparable)
black.e12 <- data.frame(position = sample.info.r9$position,
                        mutation = sample.info.r9$mutation)
black.e12[grep(pattern = "12", black.e12$position),]$mutation <- as.character(unique(
  sample.info.r9[which(sample.info.r9$backbone == "PJET" & sample.info.r9$position != "PJET" & sample.info.r9$mutation != "NA"),]$mutation))
#1B Run
stats.r9.e12 <- getstats.pertable(filelist = myfiles.r9,
                                  sampleinfo = sample.info.r9,
                                  mincoverage = 100,
                                  blackposlist = black.e12$mutation,
                                  perinsert = sequences_tp53)
#1C Select amplicon 12
stats.r9.e12 <- stats.r9.e12[grep(pattern = "12", stats.r9.e12$REF),]
#1D Remove blacklist
remove(black.e12)
#1E Add BBtype as analysis
stats.r9.e12[grep(pattern = "PJET", stats.r9.e12$REF),]$ANALYSIS <- c("PCR-free")
stats.r9.e12[grep(pattern = "BB24", stats.r9.e12$REF),]$ANALYSIS <- c("PCR")




# ---- DILUTION EXPERIMENTS ----

#1 Get mutant positions and mutations
#1A Get mutations
positions <- as.character(unique(sample.info.r9[grep(pattern = ";",sample.info.r9$mutation),]$mutation))
positions <- unlist(strsplit(split = " ; ",positions))
mutations <- unlist(strsplit(positions, split = " "))[seq(2,length(positions)*2,2)]
#1B Get ref and mut alleles
refalleles <- unlist(strsplit(mutations,split = ">"))[seq(1,length(mutations)*2,2)]
mutalleles <- unlist(strsplit(mutations,split = ">"))[seq(2,length(mutations)*2,2)]
#1C Get positions
positions <- unlist(strsplit(positions, split = " "))[seq(1,length(positions)*2,2)]
positions <- unlist(strsplit(positions, split = ":"))[seq(2,length(positions)*2,2)]

#2 Get files with amplicon 12 data
myfiles.r9.e12 <- myfiles.r9[which(grepl(sample.info.r9$position, pattern = 12))]

#3 Only focus on these positions
#3A Empty df
e12.mutation.df <- data.frame()
#3B Loop through files
for(i in 1:length(myfiles.r9.e12)) {
  tempdf <- myfiles.r9.e12[[i]]
  tempdf <- tempdf[which(tempdf$POS %in% positions),]
  # don't blacklist for triple mutant
  if(length(grep(pattern = "PJET",tempdf$SAMPLE)) ==3) {
    tempdf <- getstats.perline(filelist = list(tempdf),mincoverage = 100,blackposlist = NA)[[1]]
  }
  e12.mutation.df <- rbind(e12.mutation.df,tempdf)
  remove(tempdf)
}
remove(i)
#3C Check
nrow(e12.mutation.df) == 3*length(myfiles.r9.e12)

#4 Generate FP column
e12.mutation.df$FP <- NA
for(i in 1:length(positions)) {
  e12.mutation.df[which(e12.mutation.df$POS == positions[i]),]$FP <- e12.mutation.df[which(e12.mutation.df$POS == positions[i]),paste("fp_",mutalleles[i],sep="")]
}
remove(i)

#5 Add sample type
#5A All healthy
e12.mutation.df$sampletype <- as.character("Healthy control")
#5B Change patients with mutation
patients <- as.character(metadata.r9[grep(as.character(metadata.r9$Mutation),pattern= "17:7577121 G>A"),]$Cyclomics.Sample.ID)
for(i in patients) {
  e12.mutation.df[which(e12.mutation.df$SAMPLE == i & e12.mutation.df$POS == "7577121"),]$sampletype <- "Patient"
}
remove(i)
#5C Add ratio experiment
e12.mutation.df[grep(pattern = "12WT",e12.mutation.df$SAMPLE),]$sampletype <- "0.000"
e12.mutation.df[grep(pattern = "RATI_0001",e12.mutation.df$SAMPLE),]$sampletype <- "0.001"
e12.mutation.df[grep(pattern = "12MU",e12.mutation.df$SAMPLE),]$sampletype <- "100.0"

#6 Get wt and mut calls for samples with 100% wt and samples with 100%mut
max(e12.mutation.df[grep(pattern = "12WT",e12.mutation.df$SAMPLE),]$FP * 100)
e12.mutation.df[grep(pattern = "12MU",e12.mutation.df$SAMPLE),]$FP * 100
e12.mutation.df[grep(pattern = "12MU",e12.mutation.df$SAMPLE),]$TP * 100

#7 Plot 0.001 mut and 100% wt
#7A Create df
plotdf <- rbind(e12.mutation.df[which(grepl(pattern = "RATI",e12.mutation.df$SAMPLE)),],
                e12.mutation.df[which(grepl(pattern = "12WT",e12.mutation.df$SAMPLE)),])
plotdf$POS <- as.character(plotdf$POS)
#7B Plot
dilplot <- ggplot(plotdf, aes(x = POS, y = FP, colour = SAMPLE, fill = sampletype)) +
  geom_bar(stat = "identity",position=position_dodge(), colour = "white") +
  scale_fill_manual(values = c("black", "gray85")) + 
  scale_x_discrete() +
  scale_y_continuous(expand = c(0,0),labels = comma, limits = c(0,0.0006),breaks = seq(0,0.0006,0.0002)) +
  xlab("Position") + 
  ylab("Observed\nmutation rate") +
  labs(fill = "Mixed mutation rate") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
#7C Save plot
ggsave(plot = dilplot, filename = paste(outdir,"insert/mutationrate.pdf",sep=""),width =4.5, height = 2.5)
#7D Clean up
remove(plotdf,dilplot)

#8 Clean up
remove(positions, refalleles,mutalleles,mutations,patients)

#9 Plot MUT and WT
#9A Create df
plotdf <- rbind(e12.mutation.df[which(grepl(pattern = "12MU",e12.mutation.df$SAMPLE)),],
                e12.mutation.df[which(grepl(pattern = "12WT",e12.mutation.df$SAMPLE)),])
plotdf <- melt(plotdf, id.vars = c("POS","SAMPLE"), measure.vars = c("A","C","G","T"))
plotdf$SAMPLE <- as.character(plotdf$SAMPLE)
plotdf[grep(plotdf$SAMPLE, pattern = "12MU"),]$SAMPLE <- "MUT"
plotdf[grep(plotdf$SAMPLE, pattern = "12WT"),]$SAMPLE <- "WT"
plotdf$variable <- as.character(plotdf$variable)
temp <- data.frame()
for(i in 1:nrow(plotdf)) {
  if(plotdf[i,]$value > 0) {
    temp <- rbind(temp,matrix(unlist(rep(plotdf[i,],plotdf[i,]$value)),ncol = ncol(plotdf), byrow=T))
  }
}
remove(i)
colnames(temp) <- colnames(plotdf)
plotdf <- temp
plotdf$POS <- factor(plotdf$POS,levels = c("7577094","7577120","7577121"))
remove(temp)
#9B Plot
plot1 <- ggplot(plotdf[which(plotdf$SAMPLE == "WT"),], aes(x = POS, y = variable, colour = variable)) +
  geom_point(aes(color = variable),size = 0.000001,shape = '.',position = position_jitter()) +
  scale_color_manual(values = c("royalblue","indianred4","forestgreen","goldenrod")) +
  scale_y_discrete(drop = FALSE) +
  xlab("Position") + 
  ylab(" ") +
  guides(col=guide_legend(" ")) +
  theme_bw() +
  ggtitle("WT") +
  theme(legend.position = "bottom") +
  geom_text(mapping = aes(x = 1, y = "G"),label = nrow(plotdf[which(plotdf$POS == "7577094" & plotdf$SAMPLE == "WT" & plotdf$variable == "G"),]), colour = "white") +
  geom_text(mapping = aes(x = 2, y = "C"),label = nrow(plotdf[which(plotdf$POS == "7577120" & plotdf$SAMPLE == "WT" & plotdf$variable == "C"),]), colour = "white") +
  geom_text(mapping = aes(x = 3, y = "G"),label = nrow(plotdf[which(plotdf$POS == "7577121" & plotdf$SAMPLE == "WT" & plotdf$variable == "G"),]), colour = "white") 
plot2 <- ggplot(plotdf[which(plotdf$SAMPLE == "MUT"),], aes(x = POS, y = variable, colour = variable)) +
  geom_point(aes(color = variable),size = 0.000001,shape = '.',position = position_jitter()) +
  scale_color_manual(values = c("royalblue","indianred4","forestgreen","goldenrod")) +
  xlab("Position") + 
  scale_y_discrete(drop = FALSE) +
  ylab(" ") +
  guides(col=guide_legend(" ")) +
  theme_bw() +
  ggtitle(" MUT") +
  theme(legend.position = "bottom") + 
  geom_text(mapping = aes(x = 1, y = "A"),label = nrow(plotdf[which(plotdf$POS == "7577094" & plotdf$SAMPLE == "MUT" & plotdf$variable == "A"),]), colour = "white") +
  geom_text(mapping = aes(x = 2, y = "T"),label = nrow(plotdf[which(plotdf$POS == "7577120" & plotdf$SAMPLE == "MUT" & plotdf$variable == "T"),]), colour = "white") +
  geom_text(mapping = aes(x = 3, y = "A"),label = nrow(plotdf[which(plotdf$POS == "7577121" & plotdf$SAMPLE == "MUT" & plotdf$variable == "A"),]), colour = "white") +
  geom_text(mapping = aes(x = 1, y = "G"),label = nrow(plotdf[which(plotdf$POS == "7577094" & plotdf$SAMPLE == "MUT" & plotdf$variable == "G"),]), colour = "black") +
  geom_text(mapping = aes(x = 2, y = "C"),label = nrow(plotdf[which(plotdf$POS == "7577120" & plotdf$SAMPLE == "MUT" & plotdf$variable == "C"),]), colour = "black") 
mutplot <- grid.arrange(plot1,plot2,nrow = 1,ncol=2)
#9C Save plot
ggsave(plot = mutplot, filename = paste(outdir,"insert/WTMUT.png",sep=""),width =5.5, height = 3)
#9D Clean up
remove(plotdf,mutplot,plot1,plot2)
  
  


# ---- CALCULATIONS FOR MANUSCRIPT ----

#Nr of samples
length(unique(sample.info.r9[grep(sample.info.r9$position, pattern = "BB2"),]$samplename)) # Total nr of r9 samples
length(unique(sample.info.r9[which(sample.info.r9$position == "BB22"),]$samplename)) # nr of BB22 samples
length(unique(sample.info.r9[which(sample.info.r9$position == "BB24"),]$samplename)) # nr of BB24 samples
length(unique(sample.info.r9[which(sample.info.r9$position == "BB25"),]$samplename)) # nr of BB25 samples

# FP RATE: 10 reps, no FOR/REV correction
paste(
round(mean(rbind(stats.r9.all[which(stats.r9.all$REF == "BB22" & stats.r9.all$ANALYSIS == "ALL"),],
                 stats.r9.all[which(stats.r9.all$REF == "BB24" & stats.r9.all$ANALYSIS == "ALL"),],
                 stats.r9.all[which(stats.r9.all$REF == "BB25" & stats.r9.all$ANALYSIS == "ALL"),])$fp),digits = 4)
,
paste("(",paste(round(min(rbind(stats.r9.all[which(stats.r9.all$REF == "BB22"& stats.r9.all$ANALYSIS == "ALL"),],
                                stats.r9.all[which(stats.r9.all$REF == "BB24"& stats.r9.all$ANALYSIS == "ALL"),],
                                stats.r9.all[which(stats.r9.all$REF == "BB25"& stats.r9.all$ANALYSIS == "ALL"),])$fp),digits = 4),
                round(max(rbind(stats.r9.all[which(stats.r9.all$REF == "BB22"& stats.r9.all$ANALYSIS == "ALL"),],
                                stats.r9.all[which(stats.r9.all$REF == "BB24"& stats.r9.all$ANALYSIS == "ALL"),],
                                stats.r9.all[which(stats.r9.all$REF == "BB25"& stats.r9.all$ANALYSIS == "ALL"),])$fp),digits = 4),sep="-"),")",sep="")
, sep = " ")

# Percentage of positions in the bb that have an FP rate < 0.001
round(((nrow(df.stats.BB22[which(df.stats.BB22$mean_error_fp < 0.001),]) + 
    nrow(df.stats.BB24[which(df.stats.BB24$mean_error_fp < 0.001),]) +
    nrow(df.stats.BB25[which(df.stats.BB25$mean_error_fp < 0.001),])) /
  (nrow(df.stats.BB22) + nrow(df.stats.BB24) + nrow(df.stats.BB25))) *100,digits = 1)

# Number of positions in BB24 that need FOR/REV correction
sum(myfiles.r9[which(sample.info.r9$position == "BB24")][[1]]$TYPE != "10+ ALL")

# FP RATE: 10 reps, including FOR/REV correction
paste(
  round(mean(rbind(stats.r9[which(stats.r9$REF == "BB22"),],
                   stats.r9[which(stats.r9$REF == "BB24"),],
                   stats.r9[which(stats.r9$REF == "BB25"),])$fp),digits = 4)
  ,
  paste("(",paste(round(min(rbind(stats.r9[which(stats.r9$REF == "BB22"),],
                                  stats.r9[which(stats.r9$REF == "BB24"),],
                                  stats.r9[which(stats.r9$REF == "BB25"),])$fp),digits = 4),
                  round(max(rbind(stats.r9[which(stats.r9$REF == "BB22"),],
                                  stats.r9[which(stats.r9$REF == "BB24"),],
                                  stats.r9[which(stats.r9$REF == "BB25"),])$fp),digits = 4),sep="-"),")",sep="")
  , sep = " ")

# Percentage of positions with < 0.1% FP
mean(c(basetype.BB25$perc_fp_gold,basetype.BB24$perc_fp_gold,basetype.BB22$perc_fp_gold))

# Percentage of positions with > 1.0% FP
mean(c(basetype.BB22$perc_fp_bronze,basetype.BB24$perc_fp_bronze,basetype.BB25$perc_fp_bronze))

# Percentage of positions with > 1.0% error (FP+DEL)
mean(c(basetype.BB22$perc_fpdel_bronze,basetype.BB24$perc_fpdel_bronze,basetype.BB25$perc_fpdel_bronze))

#R9 vs Flongle
length(unique(sample.info.flongle$samplename)) # Total nr of flongle samples
paste("Flongle: ",
      round(mean(stats.flongle[which(stats.flongle$REF == "BB25"),]$fp)
            /mean(stats.r9[which(stats.r9$REF == "BB25"),]$fp),digits = 2),
      sep="")

#R9 vs R10
length(unique(sample.info.r10$samplename)) # Total nr of r10 samples
paste("R10: ",
      round(mean(stats.r10[which(stats.r10$REF == "BB24"),]$fp)
            /mean(stats.r9[which(stats.r9$REF == "BB24"),]$fp),digits = 2),
      sep="")

# Improvement PCR-free sample prep
paste(round((mean(stats.r9.e12[which(stats.r9.e12$ANALYSIS == "PCR"),]$fp)) / (mean(stats.r9.e12[which(stats.r9.e12$ANALYSIS == "PCR-free"),]$fp)),digits = 2),"x",sep="")
