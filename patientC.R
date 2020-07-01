# @Date created: 19 March 2020
# @Author: Myrthe Jager
# @Description: Plot mutation patient C
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
#2A Get Pr19 only
metadata.PR19 <- metadata[grep(pattern = 19,metadata$Amplicon.s.),]
#2B Remove other amplicons
metadata.PR19$Amplicon.s. <- c(19)

#3 Remove general metadata
remove(metadata)

#4 Add 'patient' or 'control'
metadata.PR19$sampletype <- c("Control")
metadata.PR19[grep(pattern = "CY_SS_PC_HN_000", metadata.PR19$Cyclomics.Sample.ID),]$sampletype <- c("Patient")
metadata.PR19[which(metadata.PR19$Time.point == 19),]$Time.point <- 5




# ---- DECIDE IF ALL/FOR/REV ----

#1 Get df with positions
pos.df <- read.table(file = paste(indir,"forrevpositions.txt",sep=""), sep = "\t",header = T)

#2 Get mutation
mutation <- as.character(unique(metadata.PR19[which(metadata.PR19$sampletype == "Patient"),]$Mutation))
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
remove(pos.df,mutation)




# ---- LOAD AND PREPARE I FILES ----

#1 Get sample info of I files
#1A Get info of files
sample.info.PR19 <- rbind(getsampleinfo.10reps(metadata = metadata.PR19),getsampleinfo.forrev(metadata = metadata.PR19))
#1B Select I files only
sample.info.PR19.INS <- sample.info.PR19[grep(pattern = "TP",sample.info.PR19$position),]
#1C Select ALL/FOR/REV files only
sample.info.PR19.INS <- sample.info.PR19.INS[which(sample.info.PR19.INS$analysis == position.type),]
#1D Add patient/control
sample.info.PR19.INS$type <- metadata.PR19$sampletype
#1E Remove
remove(sample.info.PR19)

#3 Read I files
#3A I: Open 
myfiles.PR19.INS = lapply(as.character(sample.info.PR19.INS$file), read.delim)
#3B I: Add names
names(myfiles.PR19.INS) <- paste(sample.info.PR19.INS$samplename,sample.info.PR19.INS$type,sample.info.PR19.INS$position,sample.info.PR19.INS$analysis,sep="__")

#4 Prepare I files for calculations
myfiles.PR19.INS <- preparefiles(filelist = myfiles.PR19.INS, sampleinfo = sample.info.PR19.INS)




# ---- CALCULATE STATS PER POSITION----

#1 Get error rates and add them to the original I files, per position
myfiles.PR19.INS <- getstats.perline(filelist = myfiles.PR19.INS, mincoverage = 100, fakeerror = 0.1, blackposlist = NA)

#2 Take the lines with the mutated base and combine into a single table
#2A Get mutation type and base nr
mutbase <- as.character(unique(metadata.PR19[which(metadata.PR19$sampletype == "Patient"),]$Mutation))
muttype <- tail(unlist(strsplit(split = ">",tail(unlist(strsplit(split = " ", mutbase)), n = 1))),n=1)
mutloc <- gsub("[ ACTG>]", "", tail(unlist(strsplit(split = ":", mutbase)),n=1))
remove(mutbase)
#2B Create empty df
pr19.df <- data.frame()
#2C Add mutation lines per file
for(file in myfiles.PR19.INS) {
  file <- file[which(file$POS == mutloc),]
  file <- file[c("SAMPLE","TYPE","REF","POS","REFALLELE","COV","A","C","G","T","DEL",paste("fp_",muttype,sep=""))]
  pr19.df <- rbind(pr19.df,file)
}
#2D Clean
remove(file,mutloc)

#3 Add mutation
pr19.df$mutation <- metadata.PR19$Mutation

#4 Add analysis type
pr19.df$sampletype <- metadata.PR19$sampletype

#5 Add timepoint
pr19.df$timepoint <- metadata.PR19$Time.point

#6 Remove samples with low coverage
if(sum(is.na(pr19.df$fp_T)) > 0) {
  pr19.df <- pr19.df[!is.na(pr19.df$fp_T),]
}




# ---- PLOT ----

plot.patC <- plotpatient(plotdf =pr19.df,miny = 0,maxy = 0.08,mutcolumn = c(paste("fp_",muttype,sep = "")))

ggsave(paste(outdir,"PatientC.pdf",sep=""),plot = plot.patC, width = 3, height = 3)




# ---- ddPCR ----

#1 Table ddPCR
df.dd <-  read.delim2(paste(indir,"ddPCR/patientC.txt",sep=""))
df.dd$measurement <- as.numeric(as.character(df.dd$measurement))
df.dd$value <- as.numeric(as.character(df.dd$value))

#2 Add results Cyclomicsseq
for (i in unique(df.dd$time)) {
  tempdf <- data.frame()
  if(length(pr19.df[which(pr19.df$timepoint == i),]$fp) > 0) {
    tempdf <- rbind(tempdf,data.frame(
      time = i,
      Assay = c("CyclomicsSeq"),
      measurement = pr19.df[which(pr19.df$timepoint == i),]$fp,
      value = pr19.df[which(pr19.df$timepoint == i),]$fp
    ))
  }
  df.dd <- rbind(df.dd,tempdf)
  remove(tempdf)
}
remove(i)
df.dd <- df.dd[!is.na(df.dd$measurement),]

#3 Plot
#3A ddPCR
patientplot.ddpcr <- ggplot(data = df.dd[which(df.dd$Assay == "ddPCR"),]) +
  geom_line(aes(x = time, y = value), size=2, colour = "sienna2") +
  geom_point(aes(x = time, y = measurement), size = 2, colour = "sienna4") + 
  scale_x_continuous(limits = c(0,5), expand = c(0, 0),minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,0.2), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  ggtitle("Patient - timeline") +
  theme(legend.position = "none")+
  ylab("") 

#3B Combined
patientplot.combined <- ggplot(data = df.dd) +
  geom_line(aes(x = time, y = value, color = Assay), size=2) +
  geom_point(aes(x = time, y = measurement, fill = Assay), pch=21, cex = 2,color = "black") + 
  scale_x_continuous(limits = c(0,5), expand = c(0, 0),minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,0.2), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  scale_color_manual(values = c("ddPCR" = "sienna2", "CyclomicsSeq" = "skyblue3")) +
  scale_fill_manual(values = c("ddPCR" = "sienna4", "CyclomicsSeq" = "skyblue4")) +
  ggtitle("Patient - timeline") +
  xlab("Time") +
  ylab("Fraction mutation") 
#6 Save plot
ggsave(paste(outdir,"patientc/ddPCR_results.pdf",sep=""), plot = patientplot.ddpcr, width = 2.5, height = 2)
ggsave(paste(outdir,"patientc/ddPCR_comparison.pdf",sep=""), plot = patientplot.combined, width = 4, height = 3)