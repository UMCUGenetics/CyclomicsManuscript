# @Date created: 15 December 2019
# @Author: Myrthe Jager
# @Description: Functions for statistical analyses and plotting CyclomicsSeq manuscript
# @Version: 1 July 2020




# ---- 1 Create a sample information table ---- 

#1A Get sample information for 0-40 repeats comparison
getsampleinfo.percount <- function (metadata) {
  #Empty table for output
  sampleinfo <- data.frame()
  #Loop through samples
  for (i in 1:nrow(metadata)) {
    tmpmetadata <- metadata[i,]
    # Get Sample dir
    tmpdir <- paste(datadir,metadata$Processed.folder[i],sep="")
    # Get Sample name
    tmpsample <- as.character(metadata$Cyclomics.Sample.ID[i])
    # Get BB type
    tmpbb <- as.character(metadata$Backbone.version[i])
    # Get I type, gene and amplicon
    tmp_igene <- as.character(metadata$Gene.name[i])
    tmp_iamp <- as.character(metadata$Amplicon.s.[i])
    tmpi <- paste(tmp_igene,tmp_iamp,sep=" amplicon ")
    # Get Processed folder name 
    tmpprocessed <- as.character(metadata$Processed.folder[i])
    
    # Get BB files
    bb.split.files <- list.files(paste(tmpdir,"SPLIT_backbone_maxcov_slurm/bin_consensus",sep="/"))
    bb.split.files.full <-list.files(paste(tmpdir,"SPLIT_backbone_maxcov_slurm/bin_consensus",sep="/"), full.names = T)
    bb.nrs <- c()
    # If PJET: get numbers of PJET files
    if(length(grep(pattern = "PJET", x = tmpbb)) > 0) {
      bb.nrs  <- grep(pattern = "pJet", x = bb.split.files)
    }
    # If BB: get numbers of BB24 files
    if(length(grep(pattern = "BB24", x = tmpbb)) > 0) {
      if(length(grep(pattern = "BB200_4n", x = bb.split.files)) > 0 & 
         nrow(tmpmetadata[which(tmpmetadata$Reference.genome.version == "Version 9" | tmpmetadata$Reference.genome.version == "version9_decoy_patient119"),]) == 1) {
        bb.nrs <- grep(pattern = "BB200_4n", x = bb.split.files)  
      }
      if(length(grep(pattern = "BB24", x = bb.split.files)) > 0 & 
         nrow(tmpmetadata[which(tmpmetadata$Reference.genome.version == "Version 9" | tmpmetadata$Reference.genome.version == "version9_decoy_patient119"),]) == 0) {
        bb.nrs <- grep(pattern = "BB24", x = bb.split.files)  
      }
    }
    # If BB25: get numbers of BB25 files
    if(length(grep(pattern = "BB25", x = tmpbb)) > 0) {
      bb.nrs <- grep(pattern = "BB25", x = bb.split.files)
    }
    # If BBCR: get numbers of BBCR files
    if(length(grep(pattern = "BBCR", x = tmpbb)) > 0){
      bb.nrs <- grep(pattern = "BBCR", x = bb.split.files)
    }
    # If BB22: get numbers of BB22 files
    if(length(grep(pattern = "BB22", x = tmpbb)) > 0){
      bb.nrs <- grep(pattern = "BB22", x = bb.split.files)
    }
    # Get BB files
    bb.split.files <- bb.split.files[bb.nrs]
    bb.split.files.full <- bb.split.files.full[bb.nrs]
    # Get rep numbers BB
    bb.temp <- matrix(unlist(strsplit(bb.split.files,"s_")), nrow=length(bb.split.files), byrow=TRUE)
    bb.temp <- matrix(unlist(strsplit(bb.temp[,2],"_f")), nrow=length(bb.split.files), byrow=TRUE)
    bb.reps <- bb.temp[,1]
    
    # Get I files
    i.split.files <- list.files(paste(tmpdir,"SPLIT_insert_maxcov_slurm/bin_consensus",sep="/"))
    i.split.files.full <- list.files(paste(tmpdir,"SPLIT_insert_maxcov_slurm/bin_consensus",sep="/"), full.names = T)
    # Get numbers of I files
    if(length(grep(pattern = "xon",x = i.split.files)) > 0) {
      i.nrs <- grep(pattern = "xon",x = i.split.files)
    }
    if(length(grep(pattern = "TP53",x = i.split.files)) > 0) {
      i.nrs <- grep(pattern = "TP53",x = i.split.files)
    }
    # Get I Files
    i.split.files <- i.split.files[i.nrs]
    i.split.files.full <- i.split.files.full[i.nrs]
    # Get rep numbers I
    i.temp <- matrix(unlist(strsplit(i.split.files,"s_")), nrow=length(i.split.files), byrow=TRUE)
    i.temp <- matrix(unlist(strsplit(i.temp[,2],"_f")), nrow=length(i.split.files), byrow=TRUE)
    i.reps <- i.temp[,1]
    
    # Intermediate cleanup
    remove(i.nrs,bb.nrs,bb.temp,i.temp)
    
    # Get nr of reps
    reps = as.numeric(length(c(bb.split.files.full,i.split.files.full)))
    
    # Should be the same
    if(length(i.reps[i.reps != bb.reps]) != 0) {next}
    # Make sample.info df
    tmpdf <- data.frame(file = c(bb.split.files.full,i.split.files.full),
                        foldername = rep(tmpprocessed,reps),
                        filename = c(bb.split.files,i.split.files), 
                        samplename = rep(tmpsample,reps),
                        type = c(bb.reps,i.reps),
                        analysis = rep("ALL",reps), 
                        position = c(rep(tmpbb,reps/2),rep(tmpi,reps/2)),
                        mutation = c(rep(NA,reps/2),rep(as.character(metadata$Mutation[i]),reps/2)),
                        backbone = rep(tmpbb,reps),
                        insertgene = rep(tmp_igene,reps),
                        insertamplicon = rep(tmp_iamp,reps),
                        flowcell = rep(as.character(metadata$Flowcell.version[i]),reps)
                        )
    
    sampleinfo <- rbind(sampleinfo, tmpdf)
    
    remove(i.split.files,bb.split.files,tmpdir,tmpsample,tmpbb,tmpi,tmp_igene,tmp_iamp,
           bb.split.files.full, i.split.files.full,reps,bb.reps,i.reps,tmpdf,tmpprocessed)
  }
  remove(i)
  
  levels(sampleinfo$type) <- c(levels(sampleinfo$type),"40")
  sampleinfo[which(sampleinfo$type == "40+"),]$type <- c("40")
  
  return(sampleinfo)
}
#1B Get sample information for files with at least 10 repeats
getsampleinfo.10reps <- function (metadata) {
  #Empty table for output
  sampleinfo <- data.frame()
  #Loop through samples
  for (i in 1:nrow(metadata)) {
    tmpmetadata <- metadata[i,]
    # Get Sample dir
    tmpdir <- paste(datadir,metadata$Processed.folder[i],sep="")
    # Get Sample name
    tmpsample <- as.character(metadata$Cyclomics.Sample.ID[i])
    # Get BB type
    tmpbb <- as.character(metadata$Backbone.version[i])
    # Get I type, gene and amplicon
    tmp_igene <- as.character(metadata$Gene.name[i])
    tmp_iamp <- as.character(metadata$Amplicon.s.[i])
    tmpi <- paste(tmp_igene,tmp_iamp,sep=" amplicon ")
    # Get Processed folder name 
    tmpprocessed <- as.character(metadata$Processed.folder[i])
    
    # Get all files
    all.files <- list.files(tmpdir)
    all.files.full <-list.files(tmpdir, full.names = T)
    
    # Get BB files
    bb.nrs <- c()
    # If PJET: get PJET file
    if(length(grep(pattern = "PJET", x = tmpbb)) > 0) {
      bb.nrs  <- grep(pattern = "t_pJet", x = all.files)
    }
    # If BB: get numbers of BB24 files
    if(length(grep(pattern = "BB24", x = tmpbb)) > 0) {
      if(length(grep(pattern = "BB200_4n", x = all.files)) > 0 & 
         nrow(tmpmetadata[which(tmpmetadata$Reference.genome.version == "Version 9" | tmpmetadata$Reference.genome.version == "version9_decoy_patient119"),]) == 1) {
        bb.nrs <- grep(pattern = "t_BB200_4n", all.files)
      }
      if(length(grep(pattern = "BB24", x = all.files)) > 0 & 
         nrow(tmpmetadata[which(tmpmetadata$Reference.genome.version == "Version 9" | tmpmetadata$Reference.genome.version == "version9_decoy_patient119"),]) == 0) {
        bb.nrs <- grep(pattern = "t_BB24", x = all.files)  
      }
    }
    # If BB25: get numbers of BB25 files
    if(length(grep(pattern = "BB25", x = tmpbb)) > 0) {
      bb.nrs <- grep(pattern = "t_BB25.txt", x = all.files)
    }
    # If BBCR: get numbers of BBCR files
    if(length(grep(pattern = "BBCR", x = tmpbb)) > 0){
      bb.nrs <- grep(pattern = "output_BBCR", x = all.files)
    }
    # If BB22: get numbers of BB22 files
    if(length(grep(pattern = "BB22", x = tmpbb)) > 0){
      bb.nrs <- grep(pattern = "t_BB22", x = all.files)
    }
    
    # Get I files
    i.nrs <- c()
    
    # Get numbers of I files
    if(length(grep(pattern = "t_exon",x = all.files)) > 0) {
      i.nrs <- grep(pattern = "t_exon",x = all.files)
    }
    if(length(grep(pattern = "t_Exon",x = all.files)) > 0) {
      i.nrs <- grep(pattern = "t_Exon",x = all.files)
    }
    if(length(grep(pattern = "t_TP53",x = all.files)) > 0) {
      i.nrs <- grep(pattern = "t_TP53",x = all.files)
    }
    # Get I Files
    all.files <- all.files[c(bb.nrs,i.nrs)]
    all.files.full <- all.files.full[c(bb.nrs,i.nrs)]
    
    # Get nr of reps
    reps = as.numeric(length(all.files))
    
    # Get analysis type (F or R)
    analysis <- rep("ALL",reps)
    
    # Get position
    # If both I and BB
    if(length(i.nrs) == length(bb.nrs)) {
      position = c(rep(tmpbb,reps/2),rep(tmpi,reps/2))
      mutation = c(rep(NA,reps/2),rep(as.character(metadata$Mutation[i]),reps/2))
    }
    # If only one of these
    if(length(i.nrs) != length(bb.nrs)) {
      if(length(i.nrs) > 0) {
        position = rep(tmpi,reps)
        mutation = rep(as.character(metadata$Mutation[i]),reps)
      }
      if(length(bb.nrs) > 0) {
        position = rep(tmpbb,reps)
        mutation = rep(NA,reps)
      }
    }

    # Make sample.info df
    tmpdf <- data.frame(file = all.files.full,
                        foldername = rep(tmpprocessed,reps),
                        filename = all.files, 
                        samplename = rep(tmpsample,reps),
                        type = rep("10+",reps),
                        analysis = analysis, 
                        position = position,
                        mutation = mutation,
                        backbone = rep(tmpbb,reps),
                        insertgene = rep(tmp_igene,reps),
                        insertamplicon = rep(tmp_iamp,reps),
                        flowcell = rep(as.character(metadata$Flowcell.version[i]),reps)
    )
    
    sampleinfo <- rbind(sampleinfo, tmpdf)
    
    remove(all.files,all.files.full,tmpdir,tmpsample,tmpbb,tmpi,tmp_igene,tmp_iamp,
           reps,tmpdf,tmpprocessed,analysis,tmpmetadata,i.nrs,bb.nrs,position,mutation)
  }
  remove(i)
  
  levels(sampleinfo$type) <- c(levels(sampleinfo$type),"40")
  
  return(sampleinfo)
} 
#1C For split forward and reverse files with at least 10 repeats
getsampleinfo.forrev <- function (metadata) {
  #Empty table for output
  sampleinfo <- data.frame()
  #Loop through samples
  for (i in 1:nrow(metadata)) {
    tmpmetadata <- metadata[i,]
    # Get Sample dir
    tmpdir <- paste(datadir,metadata$Processed.folder[i],sep="")
    # Get Sample name
    tmpsample <- as.character(metadata$Cyclomics.Sample.ID[i])
    # Get BB type
    tmpbb <- as.character(metadata$Backbone.version[i])
    # Get I type, gene and amplicon
    tmp_igene <- as.character(metadata$Gene.name[i])
    tmp_iamp <- as.character(metadata$Amplicon.s.[i])
    tmpi <- paste(tmp_igene,tmp_iamp,sep=" amplicon ")
    # Get Processed folder name 
    tmpprocessed <- as.character(metadata$Processed.folder[i])
    
    # Get BB files
    forrev.files <- list.files(paste(tmpdir,"for_rev_split",sep="/"))
    forrev.files.full <-list.files(paste(tmpdir,"for_rev_split",sep="/"), full.names = T)
    
    if(length(list.files(paste(tmpdir,"for_rev_split_bb",sep="/")))>0) {
      forrev.files <- list.files(paste(tmpdir,"for_rev_split_bb",sep="/"))
      forrev.files.full <-list.files(paste(tmpdir,"for_rev_split_bb",sep="/"), full.names = T)
    }
    bb.nrs <- c()
    # If PJET: get PJET file
    if(length(grep(pattern = "PJET", x = tmpbb)) > 0) {
      bb.nrs  <- grep(pattern = "t_pJet", x = forrev.files)
    }
    # If BB: get numbers of BB24 files
    if(length(grep(pattern = "BB24", x = tmpbb)) > 0) {
      if(length(grep(pattern = "BB200_4n", x = forrev.files)) > 0 & 
         nrow(tmpmetadata[which(tmpmetadata$Reference.genome.version == "Version 9" | tmpmetadata$Reference.genome.version == "version9_decoy_patient119"),]) == 1) {
        bb.nrs <- grep(pattern = "t_BB200_4n", forrev.files)
      }
      if(length(grep(pattern = "BB24", x = forrev.files)) > 0 & 
         nrow(tmpmetadata[which(tmpmetadata$Reference.genome.version == "Version 9" | tmpmetadata$Reference.genome.version == "version9_decoy_patient119"),]) == 0) {
        bb.nrs <- grep(pattern = "t_BB24", x = forrev.files)  
      }
    }
    # If BB25: get numbers of BB25 files
    if(length(grep(pattern = "BB25", x = tmpbb)) > 0) {
      bb.nrs <- grep(pattern = "t_BB25.txt", x = forrev.files)
    }
    # If BBCR: get numbers of BBCR files
    if(length(grep(pattern = "BBCR", x = tmpbb)) > 0){
      bb.nrs <- grep(pattern = "t_BBCR", x = forrev.files)
    }
    # If BB22: get numbers of BB22 files
    if(length(grep(pattern = "BB22", x = tmpbb)) > 0){
      bb.nrs <- grep(pattern = "t_BB22", x = forrev.files)
    }
    
    forrev.files.bb <- forrev.files[c(bb.nrs)]
    forrev.files.full.bb <- forrev.files.full[c(bb.nrs)]
    
    # Get I files
    forrev.files <- list.files(paste(tmpdir,"for_rev_split",sep="/"))
    forrev.files.full <-list.files(paste(tmpdir,"for_rev_split",sep="/"), full.names = T)
    i.nrs <- c()
    
    # Get numbers of I files
    if(length(grep(pattern = "t_exon",x = forrev.files)) > 0) {
      i.nrs <- grep(pattern = "t_exon",x = forrev.files)
    }
    if(length(grep(pattern = "t_Exon",x = forrev.files)) > 0) {
      i.nrs <- grep(pattern = "t_Exon",x = forrev.files)
    }
    if(length(grep(pattern = "t_TP53",x = forrev.files)) > 0) {
      i.nrs <- grep(pattern = "t_TP53",x = forrev.files)
    }
    # Get I Files
    forrev.files <- c(forrev.files.bb,forrev.files[c(i.nrs)])
    forrev.files.full <- c(forrev.files.full.bb,forrev.files.full[c(i.nrs)])
    remove(forrev.files.bb,forrev.files.full.bb)
    
    # Get nr of reps
    reps = as.numeric(length(forrev.files))
    
    # Get analysis type (F or R)
    analysis <- rep("FOR",reps)
    analysis[grep(pattern = "reverse",forrev.files)] <- c("REV")
    
    # Get position
    # If both I and BB
    if(length(i.nrs) == length(bb.nrs)) {
      position = c(rep(tmpbb,reps/2),rep(tmpi,reps/2))
      mutation = c(rep(NA,reps/2),rep(as.character(metadata$Mutation[i]),reps/2))
    }
    # If only one of these
    if(length(i.nrs) != length(bb.nrs)) {
      if(length(i.nrs) > 0) {
        position = rep(tmpi,reps)
        mutation = rep(as.character(metadata$Mutation[i]),reps)
      }
      if(length(bb.nrs) > 0) {
        position = rep(tmpbb,reps)
        mutation = rep(NA,reps)
      }
    }
    
    # Make sample.info df
    tmpdf <- data.frame(file = forrev.files.full,
                        foldername = rep(tmpprocessed,reps),
                        filename = forrev.files, 
                        samplename = rep(tmpsample,reps),
                        type = rep("10+",reps),
                        analysis = analysis, 
                        position = position,
                        mutation = mutation,
                        backbone = rep(tmpbb,reps),
                        insertgene = rep(tmp_igene,reps),
                        insertamplicon = rep(tmp_iamp,reps),
                        flowcell = rep(as.character(metadata$Flowcell.version[i]),reps)
    )
    
    sampleinfo <- rbind(sampleinfo, tmpdf)
    
    remove(forrev.files,forrev.files.full,tmpdir,tmpsample,tmpbb,tmpi,tmp_igene,tmp_iamp,
           reps,tmpdf,tmpprocessed,analysis,tmpmetadata,position,mutation)
  }
  remove(i)
  
  levels(sampleinfo$type) <- c(levels(sampleinfo$type),"40")
  
  return(sampleinfo)
} 




# ---- 2 Prepare files ---- 
preparefiles <- function(filelist,sampleinfo) {
  for (i in 1:length(filelist)) {
    tempdf <- filelist[[i]]
    #A Change positions (input = 0-based)
    tempdf$POS <- tempdf$POS+1
    #B Change samplename
    tempdf$SAMPLE <- rep(sampleinfo$samplename[i], nrow(tempdf))
    #C Add type
    tempdf$TYPE <- rep(sampleinfo$type[i], nrow(tempdf))
    #D BB ONLY: Change location
    if((length(grep("BB",tempdf$REF)) + length(grep("PJET",toupper(tempdf$REF)))) > 0) {
      tempdf$REF <- rep(sampleinfo$position[i], nrow(tempdf))
    }
    #E Add REF allele + INSERT ONLY: Select regions corresponding to the amplicon simultaneously
    tempdf$REFALLELE <- rep(NA,nrow(tempdf))
    #E1 For BB add REF ALLELE
    if((length(grep("BB",tempdf$REF)) + length(grep("PJET",toupper(tempdf$REF)))) > 0) {
      for (j in 1:nrow(tempdf)) {
        tempdf[j,]$REFALLELE <- toupper(unlist(strsplit(as.character(sequences_backbones[which(sequences_backbones$name == as.character(tempdf$REF[1])),]$sequence),split="",fixed=TRUE))[j])
      }
      remove(j)
    }
    #E2 FOR I 
    if((length(grep("BB",tempdf$REF)) + length(grep("PJET",toupper(tempdf$REF)))) == 0) {
      #E2A Amplicon names
      tmpamplicon <- paste("amplicon",unlist(strsplit(split = ",",as.character(sampleinfo[i,]$insertamplicon))),sep="")
      #E2B Get chr and coords for these amplicons
      tmpcc <- unlist(strsplit(split =":",as.character(sequences_tp53[which(sequences_tp53$name %in% tmpamplicon),]$genomic.coordinates)))
      tmpchr <- tmpcc[seq(1,length(tmpcc),2)]
      tmpcoords <- unlist(strsplit(split ="-",tmpcc[seq(2,length(tmpcc),2)]))
      #E2C Put these in a table, and also add REF allele
      tmpcc <- data.frame()
      for(h in 1:(length(tmpcoords)/2)) {
        tmp <- seq(as.numeric(tmpcoords[1*h+(h-1)]),as.numeric(tmpcoords[(1*h+(h-1))+1]),1)
        tmpseq <- toupper(unlist(strsplit(as.character(sequences_tp53[which(sequences_tp53$name == tmpamplicon[h]),]$sequence),split="",fixed=TRUE)))
        tmpcc <- rbind(tmpcc,
                       data.frame(REF = rep(tmpchr[h], length(tmp)),
                                  POS = tmp,
                                  REFALLELE = tmpseq))
        remove(tmp,tmpseq)
      }
      remove(tmpcoords,tmpamplicon,tmpchr,h)
      #Remove duplicate rows
      if(sum(duplicated(tmpcc)) > 0) {
        tmpcc <- tmpcc[-which(duplicated(tmpcc$POS)),]
      }
      #E2D Loop through new table and only take the right rows + add REF Seq to these
      tmp <- data.frame()
      for(j in 1:nrow(tmpcc)) {
        if(nrow(tempdf[which(tempdf$REF == tmpcc[j,]$REF & tempdf$POS == tmpcc[j,]$POS),]) > 0) {
          tmp <- rbind(tmp,tempdf[which(tempdf$REF == tmpcc[j,]$REF & tempdf$POS == tmpcc[j,]$POS),])
          tmp[nrow(tmp),]$REFALLELE <- as.character(tmpcc[j,]$REFALLELE)
        }
      }
      tempdf <- tmp
      remove(j,tmp,tmpcc)
    }
    #F Change column order
    tempdf <- data.frame(tempdf[,c("SAMPLE","TYPE","REF","POS","REFALLELE","COV","A","C","G","T","DEL")])
    #G Overwrite df in filelist
    filelist[[i]] <- tempdf
    #H close loop
    remove(tempdf)
  }
  remove(i)
  #I Return files
  return(filelist)
} 




# ---- 3 Calculate FP rate & Qscore ---- 

#3A Per line
getstats.perline = function(filelist,fakeerror,mincoverage,blackposlist) {
  for(i in 1:length(filelist)) {
    #A Load file and blacklist position 
    tempdf <- filelist[[i]]
    blackpos <- blackposlist[i]
    #B Add TP, error and qscore columns
    tempdf$TP <- NA
    tempdf$error_fp <- NA
    tempdf$qscore_fp <- NA
    tempdf$error_del <- NA
    tempdf$qscore_del<- NA
    tempdf$error_fpdel <- NA
    tempdf$qscore_fpdel <- NA
    #C Skip empty tables
    if(is.na(tempdf$SAMPLE[1])) {
      tempdf <- tempdf[-which(is.na(tempdf$SAMPLE)),]
      filelist[[i]] <- tempdf
      remove(tempdf)
      next
    }
    #D Skip lines with low cov
    for(j in 1:nrow(tempdf)) {
      if(tempdf[j,]$COV < mincoverage) {next}
    #E fill empty TP and error columns
      #E1 for non-N positions
      if(tempdf$REFALLELE[j] != "N") {
        tempdf[j,]$TP <- tempdf[j,tempdf$REFALLELE[j]]/(tempdf[j,]$COV)
        tempdf[j,]$error_fp <- (tempdf[j,]$COV-tempdf[j,tempdf$REFALLELE[j]]-tempdf[j,]$DEL)/tempdf[j,]$COV
        tempdf[j,]$error_del <- tempdf[j,]$DEL/tempdf[j,]$COV
        tempdf[j,]$error_fpdel <- (tempdf[j,]$COV-tempdf[j,tempdf$REFALLELE[j]])/tempdf[j,]$COV
      }
      #E2 for N-positions
      if(tempdf$REFALLELE[j] == "N") {
        tempdf[j,]$TP <- NA
        tempdf[j,]$error_fp <- NA
        tempdf[j,]$error_del <- tempdf[j,]$DEL/tempdf[j,]$COV
        tempdf[j,]$error_fpdel <- NA
      }
    }
    remove(j)
    #F Calculate the Qscores
    tempdf$qscore_fp <- -10*log10(tempdf$error_fp)
    tempdf$qscore_del <- -10*log10(tempdf$error_del)
    tempdf$qscore_fpdel <- -10*log10(tempdf$error_fpdel)
    #G Substitute Inf values for "likely" qscore
    if(length(tempdf[which(tempdf$qscore_fp == "Inf"),]$qscore_fp) > 0){
      tempdf[which(tempdf$qscore_fp == "Inf"),]$qscore_fp <- -10*log10(fakeerror/tempdf[which(tempdf$qscore_fp == "Inf"),]$COV)}
    if(length(tempdf[which(tempdf$qscore_del == "Inf"),]$qscore_del) > 0){
      tempdf[which(tempdf$qscore_del == "Inf"),]$qscore_del <- -10*log10(fakeerror/tempdf[which(tempdf$qscore_del == "Inf"),]$COV)}
    if(length(tempdf[which(tempdf$qscore_fpdel == "Inf"),]$qscore_fpdel) > 0){
      tempdf[which(tempdf$qscore_fpdel == "Inf"),]$qscore_fpdel <- -10*log10(fakeerror/tempdf[which(tempdf$qscore_fpdel == "Inf"),]$COV)}
    #H Add FP rates of base subs separately
    #H1 Add call rates
    tempdf$fp_A <- tempdf$A/tempdf$COV
    tempdf$fp_C <- tempdf$C/tempdf$COV
    tempdf$fp_G <- tempdf$G/tempdf$COV
    tempdf$fp_T <- tempdf$T/tempdf$COV
    for(j in 1:nrow(tempdf)) {
      #H2 skip lines with low coverage
      if(tempdf[j,]$COV < mincoverage) {
        tempdf[j,]$fp_A <- NA
        tempdf[j,]$fp_C <- NA
        tempdf[j,]$fp_G <- NA
        tempdf[j,]$fp_T <- NA
        next
      }
      #H3 remove TP calls for non-N positions
      if(tempdf$REFALLELE[j] != "N") {
        tempdf[j,paste("fp_",tempdf$REFALLELE[j],sep="")] <- NA
      }
      #H4 remove TP calls for N positions
      if(tempdf$REFALLELE[j] == "N") {
        tempdf[j,]$fp_A <- NA
        tempdf[j,]$fp_C <- NA
        tempdf[j,]$fp_G <- NA
        tempdf[j,]$fp_T <- NA
      }
    }
    remove(j)
    #I Remove data position on blacklist
    if(!is.na(blackpos)) {
      blackpos <- as.character(blackpos)
      # Split if multiple blackpos
      blackpos <- unlist(strsplit(blackpos, split = ";"))
      # Remove del and dup
      if(length(grep(blackpos,pattern ="-")) > 0) {
        blackpos <- blackpos[!grepl(blackpos,pattern ="-")]
      }
      # Only if there is another mutation, change blacklist positions
      if(length(blackpos) > 0) {
        # Remove everything that is not location
        blackpos <- gsub("[ ACTG>]", "", blackpos)
        # Create empty variables for chromosome and location
        blackchr <- c()
        blackloc <- c()
        # Now separate chr and loc
        for(d in 1:length(blackpos)) {
          tmp <- unlist(strsplit(blackpos[d],split = ":"))
          blackchr <- c(blackchr,tmp[1])
          blackloc <- c(blackloc,tmp[2])
          remove(tmp)
        }
        remove(d,blackpos)
        # And remove the data at this position
        for(h in 1:length(blackloc)){
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$TP <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$error_fp <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$qscore_fp <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$error_del <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$qscore_del <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$error_fpdel <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$qscore_fpdel <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$fp_A <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$fp_C <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$fp_G <- NA
          tempdf[which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]$fp_T <- NA
        }
        remove(h,blackchr,blackloc)
      }
    }
    #J Rewrite myfiles
    filelist[[i]] <- tempdf
    #K close loop
    remove(tempdf)
  }
  #L Close
  remove(i)
  return(filelist)
}

#3B Per table
getstats.pertable = function(filelist,sampleinfo, mincoverage,blackposlist,perinsert) {
  rownames(sampleinfo) <- c(1:nrow(sampleinfo))
  # Are there rows that should be duplicated, to get the stats for each insert separately? And if so: which ones?
  temprows <- nrow(sampleinfo)+1
  if(sum(is.na(nchar(perinsert))) == 0 & length(grep(pattern = ",",sampleinfo$position))>0) {
    temprows <- grep(pattern = ",",sampleinfo$position)
  }
  # For all the lines that should not be duplicated, put them in the empty qtable
  qtable <-  data.frame(SAMPLE = paste(sampleinfo[-temprows,]$position,sampleinfo[-temprows,]$samplename,sep="_"),
                        REF = sampleinfo[-temprows,]$position,
                        TYPE = sampleinfo[-temprows,]$type,
                        ANALYSIS = sampleinfo[-temprows,]$analysis,
                        FILENR = row.names(sampleinfo)[-temprows])
  # Convert to character for later on
  qtable$SAMPLE <- as.character(qtable$SAMPLE)
  qtable$REF <- as.character(qtable$REF)
  qtable$ANALYSIS <- as.character(qtable$ANALYSIS)
  qtable$FILENR <- as.character(qtable$FILENR)
  # For all the lines that SHOULD be duplicated, add them
  if(nrow(sampleinfo) != nrow(qtable)) {
    for(row in temprows) {
      # Get the amplicons
      tempampl <- unlist(strsplit(as.character(sampleinfo[row,]$position), split = ","))
      tempampl <- unlist(strsplit(tempampl, split = " "))
      temp <- c()
      for(a in 1:(length(tempampl)-2)) {
        temp <- c(temp,paste(tempampl[1],tempampl[2],tempampl[a+2],sep=" "))
      }
      tempampl <- temp
      remove(a,temp)
      tempnr <- length(tempampl)
      # Create df
      tmpdf <- data.frame(SAMPLE = paste(tempampl,sampleinfo[row,]$samplename,sep="_"),
                          REF = tempampl,
                          TYPE = rep(sampleinfo[row,]$type,tempnr),
                          ANALYSIS = sampleinfo[row,]$analysis,
                          FILENR = rep(row,tempnr))
      # Combine the df with the qtable
      qtable <- rbind(qtable,tmpdf)
      remove(tmpdf,tempnr,tempampl)
    }
    remove(row)
  }
  # Convert to character for later on
  qtable$SAMPLE <- as.character(qtable$SAMPLE)
  qtable$REF <- as.character(qtable$REF)
  qtable$ANALYSIS <- as.character(qtable$ANALYSIS)
  qtable$FILENR <- as.character(qtable$FILENR)
  # For insert: Change the name of the sample (No amplicon, just gene name), and the REF for Insert (to "I in BBx")  
  for(h in 1:nrow(qtable)) {
    if(length(grep(pattern = " ",qtable$SAMPLE[h])) > 0 ) {
      sample <- tail(unlist(strsplit(qtable$SAMPLE[h], split = " ")),n=1)
      gene <- head(unlist(strsplit(qtable$SAMPLE[h], split = " ")),n=1)
      qtable[h,]$SAMPLE <- paste(gene,sample,sep="_")
      remove(sample,gene)
    }
    if((length(grep("BB",qtable$REF[h])) + length(grep("PJET",qtable$REF[h]))) == 0) {
      
      qtable[h,]$REF <- paste(qtable[h,]$REF,sampleinfo[as.numeric(qtable[h,]$FILENR),]$backbone,sep= " in ")
    }
  }
  remove(h)
  # Add columns for coverage, fp and fpdel
  qtable$COV <- NA
  qtable$qscore_fp <- NA
  qtable$qscore_fpdel <- NA
  qtable$fp <- NA
  qtable$error <- NA
  # Change blackpos if NO to NA
  if(length(grep(pattern = "no", tolower(blackposlist))) == 1) {
    blackposlist <- NA
  }
  # Loop through files
  for (i in 1:length(filelist)) {
    tempdf.full <- filelist[[i]]
    # Get the right miniqtable
    miniq <- qtable[which(qtable$FILENR == i),]
    # Loop through inserts
    for(m in 1:nrow(miniq)) {
      k <- as.numeric(row.names(miniq)[m]) # i = filenr, k = linenumber in qtable (and m in miniqtable)
      # Get the table, for each insert separately if necessary
      tempdf <- tempdf.full
      if(nrow(miniq) > 1) {
        #Get amplicon name
        ampl <- unlist(strsplit(miniq[m,]$REF,split = " "))
        ampl <- paste(ampl[2],ampl[3],sep="")
        #Get corresponding coordinates
        tmpcoords <- unlist(strsplit(as.character(perinsert[which(perinsert$name == ampl),]$genomic.coordinates), split = ":"))[-1]
        tmpcoords <- unlist(strsplit(tmpcoords,split = "-"))
        tmpcoords <- c(tmpcoords[1]:tmpcoords[2])
        #Now subset the full tempdf
        tempdf <- data.frame()
        for(pos in tmpcoords) {
          tempdf <- rbind(tempdf,tempdf.full[which(tempdf.full$POS == pos),])
        }
        #Close this part
        remove(pos,ampl,tmpcoords)
      }
      # Fill qtable with NA if table is empty
      if(is.na(tempdf$SAMPLE[1])) {
        qtable[k,]$COV <- NA
        qtable[k,]$qscore_fp <- NA
        qtable[k,]$qscore_fpdel <- NA
        qtable[k,]$fp <- NA
        qtable[k,]$error <- NA
        remove(tempdf)
        next
      }
      # If Coverage is too low, don't calculate scores
      if(max(tempdf$COV) < mincoverage) {
        qtable[k,]$COV <- max(tempdf$COV)
        qtable[k,]$qscore_fp <- NA
        qtable[k,]$qscore_fpdel <- NA
        qtable[k,]$fp <- NA
        qtable[k,]$error <- NA
        remove(tempdf)
        next
      }
      # If coverage is high enough, calculate scores
      if(max(tempdf$COV) >= mincoverage) {
        # remove N positions from backbone
        tempdf <- tempdf[which(tempdf$REFALLELE != "N"),]
        # remove mutated positions from insert
        if(!is.na(blackposlist[i])) {
          #Get blackpos from list
          blackpos <- as.character(blackposlist[i])
          # Split if multiple blackpos
          blackpos <- unlist(strsplit(blackpos, split = ";"))
          # Remove del and dup
          if(length(grep(blackpos,pattern ="-")) > 0) {
            blackpos <- blackpos[!grepl(blackpos,pattern ="-")]
          }
          # And only process the data if there is a mutation left
          if(length(blackpos) > 0) {
            # Remove everything that is not location
            blackpos <- gsub("[ ACTG>]", "", blackpos)
            # Create empty variables for chromosome and location
            blackchr <- c()
            blackloc <- c()
            # Now separate chr and loc
            for(d in 1:length(blackpos)) {
              tmp <- unlist(strsplit(blackpos[d],split = ":"))
              blackchr <- c(blackchr,tmp[1])
              blackloc <- c(blackloc,tmp[2])
              remove(tmp)
            }
            remove(d)
            #remove them, for multiple per row
            for(j in 1:length(blackloc)) {
              if(length(grep(pattern = blackloc[j], tempdf$POS)) > 0) {
                tempdf <- tempdf[-which(tempdf$REF == blackchr[j] & tempdf$POS == blackloc[j]),]
              }
            }
            remove(j, blackchr,blackloc)
          }
          remove(blackpos)
        }
        # Add count TP observations
        tempdf$TP <- NA
        for (j in 1:nrow(tempdf)) {
          tempdf[j,]$TP <- tempdf[j,tempdf$REFALLELE[j]]
        }
        remove(j)
        # Add count FP observations
        tempdf$FP = tempdf$COV-tempdf$TP-tempdf$DEL
        # Calculate qscores
        qtable[k,]$qscore_fp <- -10*log10((sum(tempdf$FP))/(sum(tempdf$COV)))
        qtable[k,]$qscore_fpdel <- -10*log10((sum(tempdf$FP)+sum(tempdf$DEL))/(sum(tempdf$COV)))
        qtable[k,]$fp <- (sum(tempdf$FP))/(sum(tempdf$COV))
        qtable[k,]$error <- (sum(tempdf$FP)+sum(tempdf$DEL))/(sum(tempdf$COV))
        # Add coverage
        qtable[k,]$COV <- max(tempdf$COV)
      }   
      remove(tempdf,k)  
    }
    remove(m,tempdf.full,miniq)
  }
  # Remove filenrs column
  qtable <- qtable[,-grep(colnames(qtable), pattern = "FILENR")]
  # Close
  return(qtable)
  remove(i,temprows)
}

#3C Per position 
statsperpos = function(filelist, stat) {
  if(stat == "average") {stat <- "mean"}
  if(stat != "median" & stat != "mean") {stat <- "median"}
  #A Create empty dfs for error and Qscore (with/without DEL) and for each FP base
  #A1 First loop through files to get all possible positions and ref alleles
  posdf <- data.frame()
  for(h in 1:length(filelist)) {
    posdf <- rbind(posdf, data.frame(REF = filelist[[h]]$REF,
                                     POS = filelist[[h]]$POS,
                                     REFALLELE = filelist[[h]]$REFALLELE))
  }
  remove(h)
  posdf <- posdf[!duplicated(posdf$POS),]
  #A2 Now use this for the dfs
  edf <- posdf
  qdf <- posdf
  edeldf <-posdf
  qdeldf <- posdf
  efpdeldf <- posdf
  qfpdeldf <- posdf
  adf <- posdf
  cdf <- posdf
  gdf <- posdf
  tdf <- posdf
  #B Add counts to countdf
  for(j in 1:length(filelist)) {
    #B1 get columnname
    sample <- names(filelist)[j]
    #B2 Loop through positions
    for(h in 1:nrow(posdf)) {
      tempdf <- filelist[[j]][which(filelist[[j]]$POS == posdf[h,]$POS),]
      #B2 Add data if the position is in the file
      if(nrow(tempdf) == 1) {
        edf[h,sample] <-  tempdf$error_fp
        edeldf[h,sample] <-  tempdf$error_del
        efpdeldf[h,sample] <-  tempdf$error_fpdel
        
        qdf[h,sample] <-  tempdf$qscore_fp
        qdeldf[h,sample] <-  tempdf$qscore_del
        qfpdeldf[h,sample] <-  tempdf$qscore_fpdel
        
        adf[h,sample] <- tempdf$fp_A
        cdf[h,sample] <- tempdf$fp_C
        gdf[h,sample] <- tempdf$fp_G
        tdf[h,sample] <- tempdf$fp_T
      }
      #B3 If the position is not in the file, add NA
      if(nrow(tempdf) == 0) {
        edf[h,sample] <- NA
        edeldf[h,sample] <- NA
        efpdeldf[h,sample] <- NA
        
        qdf[h,sample] <- NA
        qdeldf[h,sample] <- NA
        qfpdeldf[h,sample] <- NA
        
        adf[h,sample] <- NA
        cdf[h,sample] <- NA
        gdf[h,sample] <- NA
        tdf[h,sample] <- NA
      }
      remove(tempdf)
    }
    remove(h,sample)
  }
  remove(j)
  #C Create empty df for error calculations
  errordf <- data.frame(matrix(ncol = 25,nrow=nrow(edf)))
  colnames(errordf) <- c("REF","POS","REFALLELE",
                         paste(stat,"_error_fp",sep =""),"sd_error_fp","se_error_fp",
                         paste(stat,"_qscore_fp",sep=""),"sd_qscore_fp","se_qscore_fp",
                         paste(stat,"_error_del",sep =""),"sd_error_del","se_error_del",
                         paste(stat,"_qscore_del",sep=""),"sd_qscore_del","se_qscore_del",
                         paste(stat,"_error_fpdel",sep =""),"sd_error_fpdel","se_error_fpdel",
                         paste(stat,"_qscore_fpdel",sep=""),"sd_qscore_fpdel","se_qscore_fpdel",
                         paste(stat,"_fp_A",sep=""),paste(stat,"_fp_C",sep=""),paste(stat,"_fp_G",sep=""),paste(stat,"_fp_T",sep=""))
  # And fill the first three columns
  errordf$REF = edf$REF
  errordf$POS = edf$POS
  errordf$REFALLELE = edf$REFALLELE
  #D Calculate statistics
  for(i in 1:nrow(errordf)) {
    #Standard deviation
    errordf[i,]$sd_error_fp = sd(unlist(edf[i,-(1:3)]),na.rm=TRUE)
    errordf[i,]$sd_error_del = sd(unlist(edeldf[i,-(1:3)]),na.rm=TRUE)
    errordf[i,]$sd_error_fpdel = sd(unlist(efpdeldf[i,-(1:3)]),na.rm=TRUE)
    errordf[i,]$sd_qscore_fp = sd(unlist(qdf[i,-(1:3)]),na.rm=TRUE)
    errordf[i,]$sd_qscore_del = sd(unlist(qdeldf[i,-(1:3)]),na.rm=TRUE)
    errordf[i,]$sd_qscore_fpdel = sd(unlist(qfpdeldf[i,-(1:3)]),na.rm=TRUE)
    #Standard error
    errordf[i,]$se_error_fp = errordf[i,]$sd_error_fp/sqrt(sum(!is.na(unlist(edf[i,-(1:3)]))))
    errordf[i,]$se_error_del = errordf[i,]$sd_error_del/sqrt(sum(!is.na(unlist(edf[i,-(1:3)]))))
    errordf[i,]$se_error_fpdel = errordf[i,]$sd_error_fpdel/sqrt(sum(!is.na(unlist(edf[i,-(1:3)]))))
    errordf[i,]$se_qscore_fp = errordf[i,]$sd_qscore_fp/sqrt(sum(!is.na(unlist(edf[i,-(1:3)]))))
    errordf[i,]$se_qscore_del = errordf[i,]$sd_qscore_del/sqrt(sum(!is.na(unlist(edf[i,-(1:3)]))))
    errordf[i,]$se_qscore_fpdel = errordf[i,]$sd_qscore_fpdel/sqrt(sum(!is.na(unlist(edf[i,-(1:3)]))))
    #Median (if stat = median)
    if(stat == "median") {
      errordf[i,]$median_error_fp = median(unlist(edf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_error_del = median(unlist(edeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_error_fpdel = median(unlist(efpdeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_qscore_fp = median(unlist(qdf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_qscore_del = median(unlist(qdeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_qscore_fpdel = median(unlist(qfpdeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_fp_A = median(unlist(adf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_fp_C = median(unlist(cdf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_fp_G = median(unlist(gdf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$median_fp_T = median(unlist(tdf[i,-(1:3)]),na.rm=TRUE)
    }
    #Mean (if stat = mean)
    if(stat == "mean") {
      errordf[i,]$mean_error_fp = mean(unlist(edf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_error_del = mean(unlist(edeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_error_fpdel = mean(unlist(efpdeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_qscore_fp = mean(unlist(qdf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_qscore_del = mean(unlist(qdeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_qscore_fpdel = mean(unlist(qfpdeldf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_fp_A = mean(unlist(adf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_fp_C = mean(unlist(cdf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_fp_G = mean(unlist(gdf[i,-(1:3)]),na.rm=TRUE)
      errordf[i,]$mean_fp_T = mean(unlist(tdf[i,-(1:3)]),na.rm=TRUE)
    }
  }
  remove(i, edf, edeldf,efpdeldf,qdf,qdeldf,qfpdeldf,adf,cdf,gdf,tdf,posdf)
  return(errordf)
}




# ---- 4 FOR/REV functions ---- 

#4A: Join ALL, FWD and REV tables
combine.stats.fwdrev <- function(all,fwd,rev,coords,selectcol) {
  if(sum(!is.na(coords)) > 0) {
    rows <- seq(rownames(all[which(all$POS == coords[1]),]),
                    rownames(all[which(all$POS == coords[2]),]))
    all <- all[rows,]
    fwd <- fwd[rows,]
    rev <- rev[rows,]
    remove(rows)
  }
  df <- data.frame(REF = all$REF,
                   POS = all$POS,
                   REFALLELE = all$REFALLELE,
                   ALL = all[,which(colnames(all) == selectcol)],
                   FWD = fwd[,which(colnames(fwd) == selectcol)],
                   REV = rev[,which(colnames(rev) == selectcol)])
  colnames(df) <- c(colnames(df[1:3]),
                    paste(colnames(df[4]),tail(unlist(strsplit(selectcol,split = "_")),n=1),sep="_"),
                    paste(colnames(df[5]),tail(unlist(strsplit(selectcol,split = "_")),n=1),sep="_"),
                    paste(colnames(df[6]),tail(unlist(strsplit(selectcol,split = "_")),n=1),sep="_"))
  return(df)
}

#4B: Calculate change if you take FOR or REV
calculate.fwd.rev.change <- function(df) {
  df$change_fwd <- df[,grep(pattern = "FWD_",colnames(df))] - df[,grep(pattern = "ALL_",colnames(df))]
  df$change_rev <- df[,grep(pattern = "REV_",colnames(df))] - df[,grep(pattern = "ALL_",colnames(df))]
  
  fwdlow <- df[,grep(pattern = "FWD_",colnames(df))] < df[,grep(pattern = "REV_",colnames(df))]
  revlow <- df[,grep(pattern = "FWD_",colnames(df))] > df[,grep(pattern = "REV_",colnames(df))]
  bothlow <- df[,grep(pattern = "FWD_",colnames(df))] == df[,grep(pattern = "REV_",colnames(df))]
  
  df$type <- NA
  df$change <- NA
  df$newvalue <- NA
  if(sum(fwdlow,na.rm=T) > 0) {
    df[which(fwdlow),]$type <- c("FOR")
    df[which(fwdlow),]$change <- df[which(fwdlow),]$change_fwd
    df[which(fwdlow),]$newvalue <- df[which(fwdlow),grep(pattern = "FWD_",colnames(df))]
  }
  if(sum(revlow,na.rm=T) > 0) {
    df[which(revlow),]$type <- c("REV")
    df[which(revlow),]$change <- df[which(revlow),]$change_rev
    df[which(revlow),]$newvalue <- df[which(revlow),grep(pattern = "REV_",colnames(df))]
  }
  if(sum(bothlow,na.rm=T) > 0) {
    df[which(bothlow),]$type <- c("BOTH")
    df[which(bothlow),]$change <- df[which(bothlow),]$change_rev
    df[which(bothlow),]$newvalue <- df[which(bothlow),grep(pattern = "REV_",colnames(df))]
  }
  remove(fwdlow,revlow,bothlow)
  return(df)
}

#4C: Select positions which really change
get.fwdrev.positions <- function(df,cutoff) {
  positions <- unique(c(df[which(df$change_fwd <= cutoff),]$POS,
                        df[which(df$change_rev <= cutoff),]$POS))
  temp_df <- data.frame()
  for(pos in positions) {
    temp_df <- rbind(temp_df,df[which(df$POS == pos),])
  }
  df <- temp_df
  remove(temp_df)
  return(df)
}

#4D: Get new filelist with corrected forrev fp rates
get.forrev.corrected.filelist <- function(filelist.all, filelist.fwd,filelist.rev,posdf) {
  # Create new filelist
  corr.filelist <- list()
  for(i in 1:length(filelist.all)) {
    # For each sample, load the filelists all/for/rev
    tempall <- filelist.all[[i]]
    tempfor <- filelist.fwd[[i]]
    temprev <- filelist.rev[[i]]
    # Create a df with position and change
    positions <- data.frame(POS = posdf[which(posdf$REF == unique(tempall$REF)),]$POS,
                            TYPE = posdf[which(posdf$REF == unique(tempall$REF)),]$type)
    # Now generate an empty data frame for the output --> the corrected & combined df
    newdf <- data.frame()
    # Loop through positions and change, if it exists in posdf
    for(j in tempall$POS) {
      # If it does not occur in the posdf, just take the 'ALL' with changed type
      if(sum(positions$POS == j) == 0) {
        templine <- tempall[which(tempall$POS == j),]
        templine$TYPE <- paste(templine$TYPE,"ALL",sep=" ")
        newdf <- rbind(newdf,templine)
        remove(templine)
      }
      if(sum(positions$POS == j) > 0) {
        if(positions[which(positions$POS == j),]$TYPE == "FOR") {
          templine <- tempfor[which(tempfor$POS == j),]
          templine$TYPE <- paste(templine$TYPE,"FOR",sep=" ")
          newdf <- rbind(newdf,templine)
          remove(templine)
        }
        if(positions[which(positions$POS == j),]$TYPE == "REV") {
          templine <- temprev[which(temprev$POS == j),]
          templine$TYPE <- paste(templine$TYPE,"REV",sep=" ")
          newdf <- rbind(newdf,templine)
          remove(templine)
        }
      }
    }
    remove(j)
    remove(tempall,tempfor,temprev)
    corr.filelist <- c(corr.filelist,list(newdf))
    remove(newdf)
  }
  remove(i, positions)
  return(corr.filelist)
}

#4E: Calculate difference in basetype before & after correcting FORREV positions
get.basetype.diff <- function(basetype.old, basetype.new) {
  basetype.diff <- basetype.old
  
  basetype.diff$fp_gold <- as.numeric(as.character(basetype.new$fp_gold)) -basetype.old$fp_gold 
  basetype.diff$fp_silver <- as.numeric(as.character(basetype.new$fp_silver)) - basetype.old$fp_silver
  basetype.diff$fp_bronze <- as.numeric(as.character(basetype.new$fp_bronze)) - basetype.old$fp_bronze
  
  basetype.diff$fpdel_gold <- as.numeric(as.character(basetype.new$fpdel_gold)) -basetype.old$fpdel_gold 
  basetype.diff$fpdel_silver <- as.numeric(as.character(basetype.new$fpdel_silver)) - basetype.old$fpdel_silver
  basetype.diff$fpdel_bronze <- as.numeric(as.character(basetype.new$fpdel_bronze)) - basetype.old$fpdel_bronze
  
  basetype.diff$perc_fp_bronze = basetype.diff$fp_bronze/basetype.diff$fp_total*100
  basetype.diff$perc_fp_silver = basetype.diff$fp_silver/basetype.diff$fp_total*100
  basetype.diff$perc_fp_gold = basetype.diff$fp_gold/basetype.diff$fp_total*100
  
  basetype.diff$perc_fpdel_bronze = basetype.diff$fpdel_bronze/basetype.diff$fpdel_total*100
  basetype.diff$perc_fpdel_silver = basetype.diff$fpdel_silver/basetype.diff$fpdel_total*100
  basetype.diff$perc_fpdel_gold = basetype.diff$fpdel_gold/basetype.diff$fpdel_total*100
  return(basetype.diff)
}

#4F: Subset dataset and define FORREV positions on this subset, test in the remaining samples how this affects the fprate
forrev.correct.fp.subsetdecision <- function(filelist.all,filelist.fwd,filelist.rev,coords,cutoff,sampleinfo) {
  # Calculate how many files for for/rev analysis and how many for new fp rate
  nr.forrev <- floor(length(filelist.all)/2)
  nr.fp <- length(filelist.all) - nr.forrev
  # Select random files for forrev
  files.forrev <- sample(1:length(filelist.all), nr.forrev)
  # New filelist with only forrev files
  forrev.all <- filelist.all[files.forrev]
  forrev.fwd <- filelist.fwd[files.forrev]
  forrev.rev <- filelist.rev[files.forrev]
  # New filelist with only fp files
  fp.all <- filelist.all[-files.forrev]
  fp.fwd <- filelist.fwd[-files.forrev]
  fp.rev <- filelist.rev[-files.forrev]
  # Calculate df stats for forrevfiles
  forrev.df.stats.all <- statsperpos(filelist = forrev.all, stat = "mean")
  forrev.df.stats.fwd <- statsperpos(filelist = forrev.fwd, stat = "mean")
  forrev.df.stats.rev <- statsperpos(filelist = forrev.rev, stat = "mean")
  # Get for/rev positions
  forrev.df.stats <- combine.stats.fwdrev(all = forrev.df.stats.all, fwd = forrev.df.stats.fwd, rev = forrev.df.stats.rev, 
                                          coords = coords, selectcol = "mean_error_fp")
  forrev.df.stats <- calculate.fwd.rev.change(df = forrev.df.stats)
  pos.df <- get.fwdrev.positions(df = forrev.df.stats,cutoff = cutoff)
  # Intermediate cleanup
  remove(nr.forrev,nr.fp)
  remove(forrev.all,forrev.fwd,forrev.rev)
  remove(forrev.df.stats.fwd,forrev.df.stats.rev,forrev.df.stats.all)
  # Now make new files for fp with the all/for/rev positions
  fp.filelist <- get.forrev.corrected.filelist(filelist.all = fp.all,
                                               filelist.fwd = fp.fwd,
                                               filelist.rev = fp.rev,
                                               posdf = pos.df)
  names(fp.filelist) <- gsub("ALL","CORR",names(fp.all))
  # Get fp stats for corrected files
  fp.df.stats <- statsperpos(filelist = fp.filelist, stat = "mean")
  fp.df.stats.uncorrected <- statsperpos(filelist = fp.all,stat = "mean")
  # Calculate basetype
  basetype.old <- get.basetype.pertable(filelist = fp.all,
                                        blackposlist = sampleinfo[-files.forrev,]$mutation,
                                        mincoverage=100)
  basetype.new <- get.basetype.pertable(filelist = fp.filelist,
                                        blackposlist = sampleinfo[-files.forrev,]$mutation,
                                        mincoverage=100)
  # Calculate the difference in basetype
  basetype.diff <- get.basetype.diff(basetype.old = basetype.old, basetype.new = basetype.new)
  # Create info table
  df.info <- data.frame(files.usedfor.forrev.calculations = paste(names(filelist.all)[files.forrev],collapse = ","),
                        files.usedfor.correctedfp.calculations = paste(names(filelist.all)[-files.forrev],collapse = ","),
                        file.order = paste(c("info","stats_files_used_for_forrev_decision","forrev_positions",
                                             "OLD_stats_otherfiles","CORR_stats_otherfiles","OLD_basetype_otherfiles","CORR_basetype_otherfiles","DIFF_basetype_otherfiles",rep("OLD_otherfiles",length(fp.all)),rep("CORR_otherfiles",length(fp.filelist))),collapse=","))
  # Combine
  stats <- c(list(df.info,
                  forrev.df.stats,
                  pos.df,
                  fp.df.stats.uncorrected,
                  fp.df.stats,
                  basetype.old,
                  basetype.new,
                  basetype.diff),
             fp.all,
             fp.filelist)
  names(stats)[1:8] <- c("info","stats_files_used_for_forrev_decision","forrev_positions",
                         "OLD_stats_otherfiles","CORR_stats_otherfiles","OLD_basetype_otherfiles","CORR_basetype_otherfiles","DIFF_basetype_otherfiles")
  
  # Final cleanup
  remove(df.info,forrev.df.stats,pos.df,fp.df.stats.uncorrected,fp.df.stats,basetype.old,basetype.new,fp.filelist)
  remove(files.forrev)
  remove(fp.all,fp.fwd,fp.rev)
  # return stat df
  return(stats)
}

#4G Prepare basetype difference table for plotting
prepare.basetypediff.forplotting <- function(basetypediff,fporfpdel) {
  basetypediff <- cbind(basetypediff[,-grep(pattern = fporfpdel,colnames(basetypediff))],
                        basetypediff[,grep(pattern = paste(c("perc",fporfpdel,""),collapse="_"),colnames(basetypediff))])
  # Calculate mean and se values
  plotdf <- data.frame(type = factor(c("< 0.1%","0.1-1.0%","> 1.0%"), levels = c("< 0.1%","0.1-1.0%","> 1.0%")),
                       mean = c(mean(basetypediff[,3]),
                                mean(basetypediff[,4]),
                                mean(basetypediff[,5])),
                       se = c((sd(basetypediff[,3]))/(sqrt(nrow(basetypediff))),
                              (sd(basetypediff[,4]))/(sqrt(nrow(basetypediff))),
                              (sd(basetypediff[,5]))/(sqrt(nrow(basetypediff))))
  )
  return(plotdf)
}




# ---- 5 GET BASETYPE ---- 
get.basetype.pertable = function(filelist,blackposlist,mincoverage) { 
  # Generate a df to write the data to
  basedf <- data.frame()
  # Define cutoff
  goldcutoff <- 1/1000
  silvercutoff <- 1/100
  # For each file
  for(i in 1:length(filelist)) {
    tempdf <- filelist[[i]]
    # Remove lines that contain mutations
    blackposlist <- as.character(blackposlist)
    if(!is.na(blackposlist[i]) & length(grep(pattern = ">",blackposlist[i])) > 0) {
      blackpos <- blackposlist[i]
      # Split for multiple mutations
      if(length(grep(pattern = ";",blackposlist[i])) > 0) {
        blackpos <- unlist(strsplit(blackposlist[i], split = ";"))
        blackpos <- blackpos[grep(pattern = ">",blackpos)]
      }
      blackpos <- gsub("[ ACTG>]", "", blackpos)
      blackchr <- c()
      blackloc <- c()
      for(d in 1:length(blackpos)) {
        tmp <- unlist(strsplit(blackpos[d],split = ":"))
        blackchr <- c(blackchr,tmp[1])
        blackloc <- c(blackloc,tmp[2])
        remove(tmp)
      }
      remove(d,blackpos)
      # now remove the lines
      for(h in 1:length(blackloc)){
        tempdf <- tempdf[-which(tempdf$REF == blackchr[h] & tempdf$POS == blackloc[h]),]
      }
      remove(h)
    }
    # Add empty fptype columns
    tempdf$fptype <- NA
    tempdf$fpdeltype <- NA
    # Remove lines with low coverage
    tempdf <- tempdf[which(tempdf$COV > mincoverage),]
    if(nrow(tempdf) > 0) {
      # Add type: bronze
      tempdf[which(!is.na(tempdf$error_fp)),]$fptype <- c("Bronze")
      tempdf[which(!is.na(tempdf$error_fpdel)),]$fpdeltype <- c("Bronze")
      # Add type: silver
      tempdf[which(tempdf$error_fp <= silvercutoff),]$fptype <- c("Silver")
      tempdf[which(tempdf$error_fpdel <= silvercutoff),]$fpdeltype <- c("Silver")
      # Add type: gold
      tempdf[which(tempdf$error_fp < goldcutoff),]$fptype <- c("Gold")
      tempdf[which(tempdf$error_fpdel < goldcutoff),]$fpdeltype <- c("Gold")
    }
    # Add counts to basedf table
    if(nrow(tempdf) > 0) {
      basedf <- rbind(basedf,data.frame(
        SAMPLE = unique(tempdf$SAMPLE),
        REF = unique(tempdf$REF),
        fp_gold = length(grep(pattern = "Gold",tempdf$fptype)),
        fp_silver = length(grep(pattern = "Silver",tempdf$fptype)),
        fp_bronze = length(grep(pattern = "Bronze",tempdf$fptype)),
        fp_total = nrow(tempdf[!is.na(tempdf$fptype),]),
        fpdel_gold = length(grep(pattern = "Gold",tempdf$fpdeltype)),
        fpdel_silver = length(grep(pattern = "Silver",tempdf$fpdeltype)),
        fpdel_bronze = length(grep(pattern = "Bronze",tempdf$fpdeltype)),
        fpdel_total = nrow(tempdf[!is.na(tempdf$fpdeltype),])
      ))
    }
    if(nrow(tempdf) == 0) {
      temp <- unlist(strsplit(names(filelist)[i], split = "__"))
      basedf <- rbind(basedf,data.frame(
        SAMPLE = temp[1],
        REF = temp[3],
        fp_gold = 0,
        fp_silver = 0,
        fp_bronze = 0,
        fp_total = 0,
        fpdel_gold = 0,
        fpdel_silver = 0,
        fpdel_bronze = 0,
        fpdel_total = 0
      ))
      remove(temp)
    }
    remove(tempdf)
  }
  basedf$perc_fp_gold <- (basedf$fp_gold/basedf$fp_total)*100
  basedf$perc_fp_silver <- (basedf$fp_silver/basedf$fp_total)*100
  basedf$perc_fp_bronze <- (basedf$fp_bronze/basedf$fp_total)*100
  basedf$perc_fpdel_gold <- (basedf$fpdel_gold/basedf$fpdel_total)*100
  basedf$perc_fpdel_silver <- (basedf$fpdel_silver/basedf$fpdel_total)*100
  basedf$perc_fpdel_bronze <- (basedf$fpdel_bronze/basedf$fpdel_total)*100
  
  # Cleanup
  remove(i, goldcutoff,silvercutoff)
  # Close
  return(basedf)
}




# ---- 6 PLOT ----

#6A FP/Qscore Across 0-40 reps
plotstats <- function (plotdf,groupon,title,plotcolumn,miny,maxy,col) {
  #If 1-fp, split the name
  do.minus <- c("NO")
  if(length(grep(pattern = "1-",plotcolumn))>0) {
    plotcolumn <- tail(unlist(strsplit(plotcolumn,split = "-")),n=1)
    do.minus <- c("YES")
  }
  # Define error bars
  errorbar <- paste("se",tail(unlist(strsplit(plotcolumn, split = "an_")),n=1),sep="_")
  if(length(grep(pattern = errorbar, colnames(plotdf))) == 0) {
    errorbar <- NA
  }
  # Get the values for grouping
  groupon.value <- plotdf[,colnames(plotdf) == groupon]
  # Exclude columns that are not necessary
  plotdf <- cbind(plotdf[,-grep(colnames(plotdf),pattern = "fp|error")],plotdf[,which(colnames(plotdf) == plotcolumn)],plotdf[,which(colnames(plotdf) == errorbar)])
  # Change colnames, table WITHOUT errorbars
  if(ncol(plotdf) == 6) {
    colnames(plotdf)[ncol(plotdf)] <- c("stat")
    #Define plottype (y-axis label)
    plottype <- strsplit(plotcolumn,split = "_")[[1]][1]
  }
  # Change colnames, table WITH errorbars
  if(ncol(plotdf) == 7) {
    colnames(plotdf)[(ncol(plotdf)-1):ncol(plotdf)] <- c("stat","def")
    #Define plottype (y-axis label)
    plottype <- paste(strsplit(plotcolumn,split = "_")[[1]][1],strsplit(plotcolumn,split = "_")[[1]][2],sep=" ")
  }
  # Longer title
  if(length(grep(pattern = "fpdel",plotcolumn)) == 1) {title <- paste(title,"(including DEL)",sep=" ")}
  # Define breaks (per 5 for qscores, per 0.01 for fp/error)
  perbreak <- 5
  if(length(grep(pattern = "qscore",plotcolumn)) == 0 & miny >= 0.90) {
    perbreak <- 0.01
  }
  if(length(grep(pattern = "qscore",plotcolumn)) == 0 & miny >= 0.965) {
    perbreak <- 0.005
  }
  if(length(grep(pattern = "qscore",plotcolumn)) == 0 & miny >= 0.995) {
    perbreak <- 0.0005
  }
  if(length(grep(pattern = "qscore",plotcolumn)) == 0 & miny == 0.925) {
    perbreak <- 0.015
  }
  # Now change the data if 1-fp
  if(do.minus == "YES") {
    plotdf$stat <- 1-plotdf$stat
    plottype <- paste("1-sn",toupper(plottype),sep="")
  }
  # Plot stats (unbinned files, without colors)
  if(length(grep(plotdf$TYPE,pattern =","))==0 & is.na(errorbar) & is.na(col)) {
    plotdf$TYPE <- as.numeric(as.character(plotdf$TYPE))
    plot <- ggplot(plotdf,aes(x=TYPE,y=stat)) +
      geom_line(aes(color= groupon.value, group = SAMPLE), size = 0.3) +
      ylab(plottype) +
      xlab("Number of repeats") +
      coord_cartesian(ylim=c(miny,maxy)) +
      scale_y_continuous(expand = c(0, 0), breaks=seq(miny,maxy, perbreak),minor_breaks = NULL)+
      ggtitle(title) +
      scale_x_continuous(expand = c(0,0))+
      theme_bw() +
      scale_colour_grey()+
      theme(legend.title=element_blank()) +
      geom_vline(xintercept=10, size = 0.3, colour = "Black", linetype = "dotted")
  }
  # Plot stats (unbinned files, with colors)
  if(length(grep(plotdf$TYPE,pattern =","))==0 & is.na(errorbar) & !is.na(col) & is.numeric(col)) {
    plotdf$TYPE <- as.numeric(as.character(plotdf$TYPE))
    cols <- c("mediumpurple3", "indianred2","seagreen3","Black","orange3")
    cols <- cols[1:col]
    plot <- ggplot(plotdf,aes(x=TYPE,y=stat)) +
      geom_line(aes(color= groupon.value, group = SAMPLE), size = 0.3) +
      ylab(plottype) +
      xlab("Number of repeats") +
      coord_cartesian(ylim=c(miny,maxy)) +
      scale_y_continuous(expand = c(0, 0), breaks=seq(miny,maxy, perbreak),minor_breaks = NULL)+
      ggtitle(title) +
      scale_x_continuous(expand = c(0,0))+
      theme_bw() +
      scale_color_manual(values = cols) +
      theme(legend.title=element_blank()) +
      geom_vline(xintercept=10, size = 0.3, colour = "Black", linetype = "dotted")
  }
  if(length(grep(plotdf$TYPE,pattern =","))==0 & is.na(errorbar) & !is.na(col) & col == T) {
    plotdf$TYPE <- as.numeric(as.character(plotdf$TYPE))
    plot <- ggplot(plotdf,aes(x=TYPE,y=stat)) +
      geom_line(aes(color= groupon.value, group = SAMPLE), size = 0.3) +
      ylab(plottype) +
      xlab("Number of repeats") +
      coord_cartesian(ylim=c(miny,maxy)) +
      scale_y_continuous(expand = c(0, 0), breaks=seq(miny,maxy, perbreak),minor_breaks = NULL)+
      ggtitle(title) +
      scale_x_continuous(expand = c(0,0))+
      theme_bw() +
     # scale_color_manual(values = cols) +
      theme(legend.title=element_blank()) +
      geom_vline(xintercept=10, size = 0.3, colour = "Black", linetype = "dotted")
  }
  if(length(grep(plotdf$TYPE,pattern =","))==0 & is.na(errorbar) & !is.na(col) & col == "boxplot") {
    plotdf$TYPE <- as.numeric(as.character(plotdf$TYPE))
    plotdf$GROUP <- paste(plotdf$TYPE,groupon.value,sep="_")
    plot <- ggplot(plotdf,aes(x=TYPE,y=stat,group = GROUP)) +
      geom_boxplot(size = 0.3,aes(fill = groupon.value)) +
      ylab(plottype) +
      xlab("Number of repeats") +
      coord_cartesian(ylim=c(miny,maxy)) +
      scale_y_continuous(expand = c(0, 0), breaks=seq(miny,maxy, perbreak),minor_breaks = NULL)+
      ggtitle(title) +
      scale_x_continuous(expand = c(0,0))+
      theme_bw() +
      scale_fill_grey()+
      theme(legend.title=element_blank()) 
  }
  return(plot)
}

#6B Plot stats per position
plotstatsperpos <- function(plotdf,title,plotcolumn,miny,maxy,color) {
  #A Define error bars
  errorbar <- paste("se",tail(unlist(strsplit(plotcolumn, split = "an_")),n=1),sep="_")
  #B Define plottype (y-axis label)
  plottype <- paste(strsplit(plotcolumn,split = "_")[[1]][1]," sn",toupper(strsplit(plotcolumn,split = "_")[[1]][3]),sep="")
  if(length(grep(pattern = "fpdel", plotcolumn)) > 0) {
    plottype <- paste(strsplit(plotcolumn,split = "_")[[1]][1],strsplit(plotcolumn,split = "_")[[1]][2],sep=" ")
  }
  #C Define max x (length bb/i)
  minx = head(plotdf$POS,n=1)
  maxx = tail(plotdf$POS,n=1)
  if(maxx != minx+maxx) {maxx <- plotdf$POS[nrow(plotdf)]}
  #D Check whether color should be included in the plot
  if(toupper(color) == "Y" | toupper(color) == "YES" | toupper(color) == "T" | toupper(color) == "TRUE") {color <- "T"}
  if(toupper(color) != "T") {color <- "F"}
  #E Exclude columns that are not necessary
  #E1 for plots without color
  if(toupper(color) == "F") {
    plotdf <- cbind(plotdf[,-grep(colnames(plotdf),pattern = "fp|se|sd|_del")],plotdf[,which(colnames(plotdf) == plotcolumn)],plotdf[,which(colnames(plotdf) == errorbar)])
    colnames(plotdf)[grep("plotcolumn", colnames(plotdf))] <- "stat"
    colnames(plotdf)[grep("errorbar", colnames(plotdf))] <- "def"
  }
  #E2 For plots with color, fp
  if(toupper(color) == "T" & length(grep("fpdel", plotcolumn)) == 0) {
    # Define columns to look for ATCG
    cols <- paste(strsplit(plotcolumn, split = "_")[[1]][-2],collapse="_")
    # generate table with right columns
    plotdf <- cbind(plotdf[,-grep(colnames(plotdf),pattern = "fp|se|sd|_del")],plotdf[,which(colnames(plotdf) == plotcolumn)],plotdf[,which(colnames(plotdf) == errorbar)],
                    plotdf[,grep(cols,colnames(plotdf))])
    # change colnames
    colnames(plotdf)[grep("plotcolumn", colnames(plotdf))] <- "stat"
    colnames(plotdf)[grep("errorbar", colnames(plotdf))] <- "def"
    # Melt for plotting
    plotdf <- melt(plotdf,id.vars = colnames(plotdf)[-grep("an_fp",colnames(plotdf))], 
                   measure.vars = colnames(plotdf)[grep("an_fp",colnames(plotdf))])
    #Remove "median"or "mean from ATCG names
    if(length(grep("median", plotcolumn)) == 1) {plotdf$variable <- data.frame(lapply(plotdf, function(x) gsub("median_fp_",'', x)))$variable}
    if(length(grep("mean", plotcolumn)) == 1) {plotdf$variable <- data.frame(lapply(plotdf, function(x) gsub("mean_fp_",'', x)))$variable}
    remove(cols)
  }
  #E3 For plots with color, fpdel
  if(toupper(color) == "T" & length(grep("fpdel", plotcolumn)) == 1) {
    # Define columns to look for ATCG AND DEL
    cols <- paste(strsplit(strsplit(plotcolumn, split = "de")[[1]][1], split = "_")[[1]][-2],collapse="_")
    delcol <- paste(strsplit(plotcolumn,"_fp")[[1]],collapse = "_")
    # generate table with right columns
    plotdf <- cbind(plotdf[,-grep(colnames(plotdf),pattern = "fp|se|sd|_del")],plotdf[,which(colnames(plotdf) == plotcolumn)],plotdf[,which(colnames(plotdf) == errorbar)],
                    plotdf[,grep(cols,colnames(plotdf))], plotdf[,grep(delcol,colnames(plotdf))])
    # change colnames
    colnames(plotdf)[grep("plotcolumn", colnames(plotdf))] <- "stat"
    colnames(plotdf)[grep("errorbar", colnames(plotdf))] <- "def"
    colnames(plotdf)[grep("delcol", colnames(plotdf))] <- "DEL"
    # Melt for plotting
    plotdf <- melt(plotdf,id.vars = colnames(plotdf)[-grep("an_fp|DEL",colnames(plotdf))], 
                   measure.vars = colnames(plotdf)[grep("an_fp|DEL",colnames(plotdf))])
    #Remove "median"or "mean from ATCG names
    if(length(grep("median", plotcolumn)) == 1) {plotdf$variable <- data.frame(lapply(plotdf, function(x) gsub("median_fp_",'', x)))$variable}
    if(length(grep("mean", plotcolumn)) == 1) {plotdf$variable <- data.frame(lapply(plotdf, function(x) gsub("mean_fp_",'', x)))$variable}
    remove(cols, delcol)
  }
  #F Get name for x axis
  xname <- as.character(plotdf[1,"REF"])
  if(xname == "17") {xname <- "Chr 17"}
  #title <- paste(title,plotcolumn, sep = " ")
  perbreak <- 0.025
  if(maxy <= 0.2) {
    perbreak <- 0.05
  }
  if(maxy <= 0.06) {
    perbreak <- 0.02
  }
  if(maxy <= 0.015) {
    perbreak <- 0.002
  }
  if(maxy == 0.009){
    perbreak <- 0.003
  }
  if(maxy <= 0.005) {
    perbreak <- 0.001
  }
  if(maxy > 0.2) {
    perbreak <- 0.1
  }
  #G For 'N'positions: do not plot anything
  if(length(grep("N", plotdf$REFALLELE)) > 0 & length(grep("DEL", plotdf$variable)) > 0 ) {
    plotdf[which(plotdf$REFALLELE == "N"),]$value <- NA
  }
  #H PLot
  #H1 Qscore
  if(length(grep(plottype,pattern = "qscore")) ==1) {
    plot <- ggplot(plotdf,aes(x = POS, y = stat)) +
      geom_point(stat="identity", fill="black", size = 0.1 ) +
      geom_errorbar(aes(ymin=stat-def, ymax=stat+def),width=.7, size=0.2) +
      ylab(plottype) +
      ggtitle(title) +
      xlab(paste("Position in ",xname,sep="")) +
      coord_cartesian(ylim=c(miny, maxy),xlim=c(minx,maxx)) +
      scale_y_continuous(expand = c(0, 0))+
      scale_x_continuous(expand = c(0, 0))+
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(panel.grid.major = element_blank())
  }
  #H2 fp/fpdel without colors
  if (length(grep(plottype,pattern = "error|snFP")) ==1 & length(grep("variable",colnames(plotdf))) == 0) {
    plot <- ggplot(plotdf,aes(x = POS, y = stat)) +
      geom_bar(stat="identity", width = 1, fill="black" ) +
      geom_errorbar(aes(ymin=stat-def, ymax=stat+def),width=.7, size=0.2) +
      ylab(plottype) +
      ggtitle(title) +
      xlab(paste("Position in ",xname,sep="")) +
      coord_cartesian(ylim=c(miny, maxy),xlim=c(minx,maxx)) +
      scale_y_continuous(expand = c(0, 0), breaks=seq(miny,maxy, perbreak))+
      scale_x_continuous(expand = c(0, 0))+
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(panel.grid.minor = element_blank())
  }
  #H3 fp/fpdel with colors
  if(length(grep(plottype,pattern = "error|snFP")) ==1 & length(grep("variable",colnames(plotdf))) == 1) {
    #Define colors for either fp or fpdel
    if(length(grep("DEL", plotdf$variable)) == 0) {colors <- c("royalblue","indianred4","forestgreen","whitesmoke","goldenrod")}
    if(length(grep("DEL", plotdf$variable)) > 0)  {colors <- c("royalblue","indianred4","grey","forestgreen","whitesmoke","goldenrod")}
    #For sequences without N: define colors
    if(length(grep("DEL", plotdf$variable)) == 0 & length(grep(pattern = "N", plotdf$REFALLELE)) == 0) {colors <- c("royalblue","indianred4","forestgreen","goldenrod")}
    if(length(grep("DEL", plotdf$variable)) > 0 & length(grep(pattern = "N", plotdf$REFALLELE)) == 0)  {colors <- c("royalblue","indianred4","grey","forestgreen","goldenrod")}
    #plot 
    plot <- ggplot(plotdf,aes(x = POS, y = value, fill = variable)) +
      geom_bar(stat="identity", width = 1) +
      geom_col(aes(y = c(-1*(perbreak/5)), fill = REFALLELE), width = 1) +
      geom_errorbar(aes(ymin=stat-def, ymax=stat+def),width=.7, size=0.2) +
      ylab(plottype) +
      ggtitle(title) +
      xlab(paste("Position in ",xname,sep="")) +
      scale_fill_manual(values=colors) +
      scale_y_continuous(expand = c(0,0), breaks=seq(miny,maxy, perbreak))+
      coord_cartesian(ylim=c(c(-1*(perbreak/7)),maxy)) + 
      scale_x_continuous(expand = c(0, 0))+
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "bottom")+
      geom_hline(yintercept=0, size = 0.3, colour = "Black", linetype = "solid")
    
  }
  #I Close
  remove(errorbar,plottype,plotdf,maxx,xname)
  return(plot)
} 

#6C Plot basetype
plotbasetype <- function(basetype, title, plottype){
  # Take either fp or fpdel
  basetype <- cbind(basetype[,-grep(pattern = "fp",colnames(basetype))],
                    basetype[,grep(pattern = paste(c("perc",plottype,""),collapse="_"),colnames(basetype))])
  # Calculate mean and se values
  plotdf <- data.frame(type = factor(c("< 0.1%","0.1-1.0%","> 1.0%"), levels = c("< 0.1%","0.1-1.0%","> 1.0%")),
                       mean = c(mean(basetype[,3]),
                                mean(basetype[,4]),
                                mean(basetype[,5])),
                       se = c((sd(basetype[,3]))/(sqrt(nrow(basetype))),
                              (sd(basetype[,4]))/(sqrt(nrow(basetype))),
                              (sd(basetype[,5]))/(sqrt(nrow(basetype))))
  )
  # prepare basetype table
  basetype <- melt(basetype)
  colnames(basetype) <- c("sample","REF","type","value")
  basetype$type <- gsub("perc_fp_gold","< 0.1%",basetype$type)
  basetype$type <- gsub("perc_fp_silver","0.1-1.0%",basetype$type)
  basetype$type <- gsub("perc_fp_bronze","> 1.0%",basetype$type)
  basetype$type <- gsub("perc_fpdel_gold","< 0.1%",basetype$type)
  basetype$type <- gsub("perc_fpdel_silver","0.1-1.0%",basetype$type)
  basetype$type <- gsub("perc_fpdel_bronze","> 1.0%",basetype$type)
  # Plot
  plot <- ggplot(plotdf, aes(x = type, y = mean, fill = type)) +
    geom_bar(stat="identity", width = 0.9) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.3, size=0.2) +
    ylab("Positions (%)") +
    ggtitle(title) +
    xlab("") +
    coord_cartesian(ylim=c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks=seq(0,100, 25))+
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #scale_fill_manual(values=c("goldenrod1","#C0C0C0","#cd7f32"))
    scale_fill_manual(values=c("gray38","gray38","gray38"))
  plot <- plot + geom_jitter(data = basetype, aes(x = type, y = value), size = 0.25)
  # Clean
  remove(plotdf)
  return(plot)
} 

#6D Plot patient
plotpatient <- function(plotdf,miny,maxy,mutcolumn) {
  #Define min and max values in controls
  max <- max(plotdf[which(plotdf$sampletype == "Control"),c(mutcolumn)])
  median <- median(plotdf[which(plotdf$sampletype == "Control"),c(mutcolumn)])
  #Change NA to -1
  if(length(plotdf[which(is.na(plotdf$timepoint)),]$timepoint) > 0) {
    plotdf[which(is.na(plotdf$timepoint)),]$timepoint <- c(-1)
  }
  # Change colname of mutbase
  colnames(plotdf)[grep(mutcolumn,colnames(plotdf))] <- c("fp")
  #Controlplot
  controlplot <- ggplot(data = plotdf[which(plotdf$sampletype == "Control"),], aes(x = timepoint, y = fp, fill = sampletype)) +
    geom_boxplot(colour = "darkgoldenrod", fill = "goldenrod") +
    geom_point(size = 2, colour = "darkgoldenrod4",position = position_jitter(w = 0.05, h = 0)) +
    scale_x_continuous(limits = c(-1.5,-0.5), breaks = c(-1)) +
    scale_y_continuous(limits = c(-maxy/500,maxy), expand = c(0, 0), minor_breaks = NULL) +
    theme_bw()+
    xlab("") +
    ylab("Fraction mutation") +
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(colour="white"),axis.ticks = element_blank())
  #Patientplot
  patientplot <- ggplot(data = plotdf[which(plotdf$sampletype == "Patient"),], aes(x = timepoint, y = fp, fill = sampletype)) +
    geom_line(colour = "skyblue3", size=2) +
    geom_point(colour = "skyblue4", size=2) + 
    scale_x_continuous(limits = c(min(plotdf[which(plotdf$sampletype == "Patient"),]$timepoint),max(plotdf[which(plotdf$sampletype == "Patient"),]$timepoint)),
                       breaks = seq(min(plotdf[which(plotdf$sampletype == "Patient"),]$timepoint),
                                    max(plotdf[which(plotdf$sampletype == "Patient"),]$timepoint)), minor_breaks = NULL,expand = c(0, 0)) +
    scale_y_continuous(limits = c(-maxy/500,maxy), expand = c(0, 0), minor_breaks = NULL) +
    theme_bw() +
    geom_hline(yintercept = median, size = 1, colour = "darkgoldenrod", linetype = "dashed") +
    geom_rect(aes(ymin = 0, ymax = max, xmin = -Inf, xmax = Inf),
              fill = "darkgoldenrod", alpha = 0.05) +
    ylab("") +
    xlab("Time") +
    guides(fill=FALSE) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  #Combined plot
  plot <- plot_grid(controlplot,patientplot,nrow=1, rel_widths = c(1/3,2/3))
  return(plot)
}

#6E Plot coverage
plotcov <- function(plotdf,maxy,pertype,breaks) {
  # Change type for plotting
  plotdf$TYPE <- as.numeric(as.character(plotdf$TYPE))
  # calculate where the major breaks will be
  temp <- c()
  for(i in 1:breaks) {
    temp <- c(temp,maxy/breaks*i)
  }
  breaks <- temp
  remove(i,temp)
  #PLOT
  if(is.na(pertype)) {
    plot <- ggplot(plotdf, aes(x = TYPE, y = COV, group = TYPE)) + 
      geom_boxplot(fill = "gray69",outlier.shape = NA,color = "gray69") +
      coord_cartesian(ylim=c(0,maxy)) +
      scale_y_continuous(position = "right",expand = c(0, 0), breaks = breaks, minor_breaks = NULL)+
      scale_x_continuous(expand = c(0,0))+
      ylab("Coverage") +
      xlab("Number of repeats") +
      ggtitle("") +
      theme_bw() +
      theme(legend.title=element_blank())
  }
  if(!is.na(pertype)) {
    plotdf$PT <- paste(plotdf$TYPE,plotdf$REF,sep= " ")
    if(length(grep(pattern = "PCR",plotdf$ANALYSIS)) > 0) {
      plotdf$PT <- paste(plotdf$ANALYSIS,plotdf$TYPE,sep= " ")
      plotdf$REF <- plotdf$ANALYSIS
    }
    plot <- ggplot(plotdf, aes(x = TYPE, y = COV, group = PT)) + 
      geom_boxplot(aes(fill = REF,color = REF),outlier.shape = NA) +
      coord_cartesian(ylim=c(0,maxy)) +
      scale_y_continuous(position = "right",expand = c(0, 0), breaks = breaks, minor_breaks = NULL)+
      scale_x_continuous(expand = c(0,0))+
      ylab("Coverage") +
      xlab("Number of repeats") +
      ggtitle("") +
      scale_fill_grey()+
      scale_color_grey()+
      theme_bw() 
  }
  return(plot)
}

#6F Plot forrev
plot.forrev <- function(plotdf) {
  # Remove columns that we don't need
  plotdf <- plotdf[,-grep(colnames(plotdf),pattern ="ALL_|FWD_|REV_|type|newvalue")]
  plotdf <- plotdf[,-which(colnames(plotdf) == "change")]
  # Melt for ggplot
  plotdf <- melt(plotdf, id.vars=c("REF","POS","REFALLELE"))
  # Add type: For or Rev
  plotdf$type <- NA
  plotdf[grep(pattern = "fwd",plotdf$variable),]$type <- "Forward"
  plotdf[grep(pattern = "rev",plotdf$variable),]$type <- "Reverse"
  # Select either forward or reverse
  for(i in unique(plotdf$POS)) {
    tempdf <- plotdf[which(plotdf$POS==i),]
    rows <- rownames(tempdf)
    tempdf <- tempdf[which(tempdf$value == min(tempdf$value)),]
    if(nrow(tempdf) == 2) {
      tempdf <- tempdf[1,]
    }
    if(nrow(tempdf) > 0) {
      rows <- rows[rows != rownames(tempdf)]
    }
    if(nrow(tempdf)==0) {
      rows <- rows[1]
    }
    plotdf <- plotdf[-which(rownames(plotdf) == rows),]
    remove(tempdf,rows)
  }
  plotdf <- plotdf[order(plotdf$POS),]
  # Now reverse the change to improvement
  plotdf$value = plotdf$value*-1
  # Calculate ylims
  ylims <- round(max(abs(plotdf$value),na.rm=T),3)
  if(max(abs(plotdf$value),na.rm=T) > ylims) {
    ylims <- ylims + 0.001
  }
  # Plot
  plot <- ggplot(plotdf,aes (x = POS, y = value, fill = type))+
    geom_bar(stat="identity", position = "dodge") +
    ylab("Improvement in snFP") + 
    scale_y_continuous(expand = c(0,0), limits = c(0,ylims))+
    scale_x_continuous(expand = c(0,0))+
    xlab("Position in backbone") +
    theme_bw()+
    geom_hline(yintercept = 1/1000, size = 0.3, colour = "Black", linetype = "dashed")+
    theme(legend.title = element_blank(),
            legend.position = "bottom",
          panel.grid.minor = element_blank())
  # Close
  remove(ylims)
  return(plot)
} 
