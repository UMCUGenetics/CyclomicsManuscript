# @Date created: 19 March 2020
# @Author: Myrthe Jager
# @Description: Plot mutation patient A
# @Version: 1 July 2020




# ---- GET STARTED ----

#1 Load required packages
library(ggplot2)
library(gridExtra)
library(cowplot)

#2 Define output directory
dir = "/.../Cyclomics_manuscript/"
indir = paste(dir,"Data/",sep = "")
outdir = paste(dir,"Results/patients/",sep = "")
datadir = paste(indir,"RCA/",sep="")
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

#2 SUbset metadata
#2A Get Pr12 only
metadata.PR12 <- metadata[grep(pattern = 12,metadata$Amplicon.s.),]
#2B Remove other amplicons
metadata.PR12$Amplicon.s. <- c(12)
#2C remove PJET data
metadata.PR12 <- metadata.PR12[grep(metadata.PR12$Study.hospital.center,pattern ="UMC Utrecht"),]
#2D Remove R10
metadata.PR12 <- metadata.PR12[-grep(metadata.PR12$Flowcell.version,pattern ="R10"),]
#2E Remove saliva sample of patient
metadata.PR12 <- metadata.PR12[-which(metadata.PR12$Cyclomics.Sample.ID == "CY_SS_SC_HN_0001_006_000"),]

#3 Remove general metadata
remove(metadata)

#4 Add 'patient' or 'control'
metadata.PR12$sampletype <- c("Control")
metadata.PR12[grep(pattern = "CY_SS_PC_HN_0001", metadata.PR12$Cyclomics.Sample.ID),]$sampletype <- c("Patient")

#5 Remove muts we don't need
metadata.PR12[which(metadata.PR12$sampletype == "Control"),]$Mutation <- NA
metadata.PR12[which(metadata.PR12$sampletype == "Control"),]$Time.point <- NA




# ---- CHECK IF ALL/FOR/REV ----

#1 Get df with positions
pos.df <- read.table(file = paste(indir,"forrevpositions.txt",sep=""), sep = "\t",header = T)

#2 Get mutation
mutation <- as.character(unique(metadata.PR12[which(metadata.PR12$sampletype == "Patient"),]$Mutation))
mutation <- head(unlist(strsplit(split = " ",mutation)),n=1)
mutation <- unlist(strsplit(split=":",mutation))

#3 Check if the position is in the df, if not: ALL, if so: FOR or REV
position.type <- "ALL"
if(mutation[2] %in% pos.df$POS) {
  temp <- pos.df[which(pos.df$POS == mutation[2]),]
  for(i in 1:nrow(temp)) {
    if(temp[i,]$REF == mutation[1]) {
      position.type <- as.character(temp$type)
    }
  }
  remove(temp,i)
}

#4 Clean up
remove(mutation, pos.df)




# ---- LOAD AND PREPARE I FILES ----

#1 Get sample info of I files
#1A Get info of files
sample.info.PR12 <- rbind(getsampleinfo.10reps(metadata = metadata.PR12),getsampleinfo.forrev(metadata = metadata.PR12))
#1B Select I files only
sample.info.PR12.INS <- sample.info.PR12[grep(pattern = "TP",sample.info.PR12$position),]
#1C Select ALL/FOR/REV files only
sample.info.PR12.INS <- sample.info.PR12.INS[which(sample.info.PR12.INS$analysis == position.type),]
#1D Add patient/control
sample.info.PR12.INS$type <- c("Control")
sample.info.PR12.INS[which(!is.na(sample.info.PR12.INS$mutation)),]$type <- c("Patient")
#1E Remove
remove(sample.info.PR12)

#2 Read I files
#2A I: Open 
myfiles.PR12.INS = lapply(as.character(sample.info.PR12.INS$file), read.delim)
#2B I: Add names
names(myfiles.PR12.INS) <- paste(sample.info.PR12.INS$samplename,sample.info.PR12.INS$type,sample.info.PR12.INS$position,sample.info.PR12.INS$analysis,sep="__")

#3 Prepare I files for calculations
myfiles.PR12.INS <- preparefiles(filelist = myfiles.PR12.INS, sampleinfo = sample.info.PR12.INS)




# ---- CALCULATE STATS PER POSITION----

#1 Get error rates and add them to the original I files, per position
myfiles.PR12.INS <- getstats.perline(filelist = myfiles.PR12.INS, mincoverage = 100, fakeerror = 0.1, blackposlist = NA)

#2 Take the lines with the mutated base and combine into a single table
#2A Get mutation type and base nr
mutbase <- as.character(unique(metadata.PR12[which(metadata.PR12$sampletype == "Patient"),]$Mutation))
muttype <- tail(unlist(strsplit(split = ">",tail(unlist(strsplit(split = " ", mutbase)), n = 1))),n=1)
mutloc <- gsub("[ ACTG>]", "", tail(unlist(strsplit(split = ":", mutbase)),n=1))
remove(mutbase)
#2B Create empty df
pr12.df <- data.frame()
#2C Add mutation lines per file
filelist <- myfiles.PR12.INS
for(file in filelist) {
  file <- file[which(file$POS == mutloc),]
  file <- file[c("SAMPLE","TYPE","REF","POS","REFALLELE","COV","A","C","G","T","DEL",paste("fp_",muttype,sep=""))]
  pr12.df <- rbind(pr12.df,file)
}
#2D Clean
remove(file,mutloc)

#3 Add mutation
pr12.df$mutation <- sample.info.PR12.INS$mutation

#4 Add analysis type
pr12.df$sampletype <- sample.info.PR12.INS$type

#5 Add timepoint
pr12.df$timepoint <- NA
for(i in 1:nrow(pr12.df)) {
  if(pr12.df[i,]$sampletype == "Control") {next}
  pr12.df[i,]$timepoint <- metadata.PR12[c(grep(pattern = pr12.df[i,]$SAMPLE,metadata.PR12$Cyclomics.Sample.ID)),]$Time.point
}
remove(i)

#6 Add ALL,FOR,REV
pr12.df$analysis <- as.character(sample.info.PR12.INS$analysis)




# ---- PLOT ----

plot.patA <- plotpatient(plotdf =pr12.df,miny = 0,maxy = 0.02,mutcolumn = c(paste("fp_",muttype,sep = "")))

ggsave(paste(outdir,"PatientA.pdf",sep=""),plot = plot.patA, width = 3, height = 3)




# ---- ddPCR ----

#1 Table ddPCR
df.dd <- read.delim2(paste(indir,"ddPCR/patientA.txt",sep=""))
df.dd$measurement <- as.numeric(as.character(df.dd$measurement))
df.dd$value <- as.numeric(as.character(df.dd$value))

#2 Add results Cyclomicsseq
for (i in unique(df.dd$time)) {
  tempdf <- data.frame()
  if(length(pr12.df[which(pr12.df$timepoint == i),]$fp) > 0) {
    tempdf <- rbind(tempdf,data.frame(
      time = i,
      Assay = c("CyclomicsSeq"),
      measurement = pr12.df[which(pr12.df$timepoint == i),]$fp,
      value = pr12.df[which(pr12.df$timepoint == i),]$fp
    ))
  }
  df.dd <- rbind(df.dd,tempdf)
  remove(tempdf)
}
remove(i)

#3 Plot
#3A ddPCR
patientplot.ddpcr <- ggplot(data = df.dd[which(df.dd$Assay == "ddPCR"),]) +
  geom_line(aes(x = time, y = value), size=2, colour = "sienna2") +
  geom_point(aes(x = time, y = measurement), size = 2, colour = "sienna4") + 
  scale_x_continuous(limits = c(0,4), expand = c(0, 0),minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,0.02), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  ggtitle("Patient - timeline") +
  theme(legend.position = "none")+
  ylab("") 

#3B Combined
patientplot.combined <- ggplot(data = df.dd) +
  geom_line(aes(x = time, y = value, color = Assay), size=2) +
  geom_point(aes(x = time, y = measurement, fill = Assay), pch=21, cex = 2,color = "black") + 
  scale_x_continuous(limits = c(0,4), expand = c(0, 0),minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,0.02), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  scale_color_manual(values = c("ddPCR" = "sienna2", "CyclomicsSeq" = "skyblue3")) +
  scale_fill_manual(values = c("ddPCR" = "sienna4", "CyclomicsSeq" = "skyblue4")) +
  ggtitle("Patient - timeline") +
  xlab("Time") +
  ylab("Fraction mutation") 

#3c Save plots
ggsave(paste(outdir,"patienta/ddPCR_results.pdf",sep=""), plot = patientplot.ddpcr, width = 2.5, height = 2)
ggsave(paste(outdir,"patienta/ddPCR_comparison.pdf",sep=""), plot = patientplot.combined, width = 4, height = 3)
