# @Date created: 14 May 2019
# @Author: Myrthe Jager
# @Description: Plot ratio/percentage reads supporting a deletion observed in Patient B
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




# ---- Get data ----

#1 Generate df
df <- data.frame(SAMPLE = c("CY_SS_PC_HC_0001_001_000","CY_SM_PC_HC_0002_001_000",
                            "CY_SM_PC_HN_0002_001_000","CY_SM_PC_HN_0002_003_000"),
                 TYPE = rep("10+",4),
                 COV = c(324000,70000,225000,497000),
                 mut = c(0,0,45,0),
                 fp = NA,
                 mutation = c(NA,NA,
                              "17:7577095-7577123 Deletion", "17:7577095-7577123 Deletion"),
                 sampletype = c("Control","Control","Patient","Patient"),
                 timepoint = c(-1,-1,0,1),
                 samplename = c("C1","C2","PB-1b","PB-2b")
)

#2 Calculate fraction mutation/fp
df$fp = df$mut/df$COV




# ---- Plot ----

#1 PLot
plot.patB <- plotpatient(plotdf =df,miny = 0,maxy = 0.001, mutcolumn = "fp")

#2 Save plot
ggsave(paste(outdir,"PatientB.pdf",sep=""), plot = plot.patB, width = 3.4, height = 3)




# ---- ddPCR ----

#1 Table ddPCR
df.dd <- read.delim2(paste(indir,"ddPCR/patientB.txt",sep=""))
df.dd$measurement <- as.numeric(as.character(df.dd$measurement))
df.dd$value <- as.numeric(as.character(df.dd$value))

#2 Add results Cyclomicsseq
for (i in unique(df.dd$time)) {
  tempdf <- data.frame()
  if(length(df[which(df$timepoint == i),]$fp) > 0) {
    tempdf <- rbind(tempdf,data.frame(
      time = i,
      Assay = c("CyclomicsSeq"),
      measurement = df[which(df$timepoint == i),]$fp,
      value = df[which(df$timepoint == i),]$fp
    ))
  }
  df.dd <- rbind(df.dd,tempdf)
  remove(tempdf)
}
remove(i)

#3 Plot
#3A ddPCR
patientplot.ddpcr.alltime <- ggplot(data = df.dd[which(df.dd$Assay == "ddPCR"),]) +
  geom_line(aes(x = time, y = value), size=2, colour = "sienna2") +
  geom_point(aes(x = time, y = measurement), size = 2, colour = "sienna4") + 
  scale_x_continuous(limits = c(0,4), expand = c(0, 0),minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,0.001), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  ggtitle("Patient - timeline") +
  theme(legend.position = "none")+
  ylab("") 
patientplot.ddpcr <- ggplot(data = df.dd[which(df.dd$Assay == "ddPCR"),]) +
  geom_line(aes(x = time, y = value), size=2, colour = "sienna2") +
  geom_point(aes(x = time, y = measurement), size = 2, colour = "sienna4") + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0),minor_breaks = NULL, breaks = c(0,1,1)) +
  scale_y_continuous(limits = c(0,0.001), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  ggtitle("Patient - timeline") +
  theme(legend.position = "none")+
  ylab("") 

#3B Combined
patientplot.combined.alltime <- ggplot(data = df.dd) +
  geom_line(aes(x = time, y = value, color = Assay), size=2) +
  geom_point(aes(x = time, y = measurement, fill = Assay), pch=21, cex = 2,color = "black") + 
  scale_x_continuous(limits = c(0,4), expand = c(0, 0),minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,0.001), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  scale_color_manual(values = c("ddPCR" = "sienna2", "CyclomicsSeq" = "skyblue3")) +
  scale_fill_manual(values = c("ddPCR" = "sienna4", "CyclomicsSeq" = "skyblue4")) +
  ggtitle("Patient - timeline") +
  xlab("Time") +
  ylab("Fraction mutation") 

patientplot.combined.seqtime <- ggplot(data = df.dd) +
  geom_line(aes(x = time, y = value, color = Assay), size=2) +
  geom_point(aes(x = time, y = measurement, fill = Assay), pch=21, cex = 2,color = "black") + 
  scale_x_continuous(limits = c(min(df.dd[which(df.dd$Assay == "CyclomicsSeq"),]$time),max(df.dd[which(df.dd$Assay == "CyclomicsSeq"),]$time)), 
                     expand = c(0, 0),minor_breaks = NULL,breaks = df.dd[which(df.dd$Assay == "CyclomicsSeq"),]$time) +
  scale_y_continuous(limits = c(0,0.001), expand = c(0, 0), minor_breaks = NULL) +
  theme_bw() +
  scale_color_manual(values = c("ddPCR" = "sienna2", "CyclomicsSeq" = "skyblue3")) +
  scale_fill_manual(values = c("ddPCR" = "sienna4", "CyclomicsSeq" = "skyblue4")) +
  ggtitle("Patient - timeline") +
  theme(legend.position = "bottom")+
  xlab("Time") +
  ylab("Fraction mutation") 

#6 Save plot
ggsave(paste(outdir,"patientb/ddPCR_results.pdf",sep=""), plot = patientplot.ddpcr, width = 1.3, height = 2)
ggsave(paste(outdir,"patientb/ddPCR_comparison.pdf",sep=""), plot = patientplot.combined.alltime, width = 4, height = 3)
