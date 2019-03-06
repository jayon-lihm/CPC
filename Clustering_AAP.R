## All chr1 positions, #(AAP>0 samples)>55
## 4,969,820

args <- commandArgs(TRUE)
index <- args[1]
tmpDir <- args[2]
chrom <- args[3]

last_letter <- substr(tmpDir, nchar(tmpDir), nchar(tmpDir))
if ( last_letter != "/" ){
  tmpDir <- paste(tmpDir, "/", sep="")
}

## index <- 781
## chrom <- 2
## tmpDir <- "./tmp"

library("cluster")
samplefile <- "./familyinfo_parents.txt"
segdupfile <- "./UCSC_genomicSuperDups_short.sorted_merged.sorted.bed.merged_noCHR"

## Text file with position and # of samples with AAP>0 where # of samples are greater than 55
file <- paste("./chr", chrom, "/chr", chrom, "_d30AAP_sampleCounts.sti", index, ".bedg.gz", sep="")

outDir <- paste("./chr", chrom, "/", sep="")

aap_mat_file <- paste("chr", chrom, "_aap_mat_sti", index, "_test.RData", sep="") ## AAP matrix R object
outfile <- paste("chr", chrom, "_d30AAP_clusterSummary.sti", index, ".txt", sep="")

sample_info <- read.table(samplefile, header=F, as.is=T, sep="\t")
colnames(sample_info) <- c("famID", "sampleID", "parents")

base_data <- read.table(file, header=F, as.is=T, sep="\t")
colnames(base_data) <- c("chrom", "start", "end", "aap_count")

####### 1. COLLECT AAP from all samples
AAPdata <- base_data[,c(1,2,3)]
cat("1. Collect AAP from all samples and make AAP matrix\n")
cat(" - number of positions: ", dim(base_data)[1], "\n")
startTime <- Sys.time()

for (i in 1:dim(sample_info)[1]){
  AAPfile <- paste("./", sample_info$famID[i], "/", sample_info$sampleID[i], "/chr", chrom, "_q20_d30.aap.bed.gz", sep="") ## Specify
  REFfile <- paste("./", sample_info$famID[i], "/", sample_info$sampleID[i], "/chr", chrom, "_q20_d30.ref.bed.gz", sep="") ## Specify

  ## collect i-th sample's AAP values overlapping with SGE_TASK_ID's "bedg" file, 
  test_cmd <- paste("zcat ", AAPfile, REFfile, " | bedtools intersect -a ", file, " -b stdin -wb | awk 'BEGIN{OFS=\"\t\"}{print $5, $6, $7, $8}' | wc -l", sep=" ")
  a <- system(test_cmd, intern=TRUE)
  a <- as.numeric(a)
  cat(i, ",", a, "\n")
  if (a>0){
    cmd <- paste("zcat ", AAPfile, REFfile, " | bedtools intersect -a ", file, " -b stdin -wb | awk 'BEGIN{OFS=\"\t\"}{print $5, $6, $7, $8}'", sep=" ")
    tmp <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
    colnames(tmp) <- c("chrom", "start", "end", sample_info$sampleID[i])

    tmp <- tmp[order(tmp$chrom, tmp$start),]
    AAPdata <- merge(AAPdata, tmp, all=TRUE, by=c("chrom", "start", "end"))
  } else if (a==0){
    AAPdata[,sample_info$sampleID[i]] <- rep(NA, dim(AAPdata)[1])
  }
}
cat("\n")
endTime <- Sys.time()
print(paste("Start: ", startTime, sep=""))
print(paste("End: ", endTime, sep=""))
cat("Done\n")
save(AAPdata, file=paste(outDir, aap_mat_file, sep=""))

######### 2. CLUSTERING
base_data$count2 <- apply(as.matrix(AAPdata[,4:dim(AAPdata)[2]]), 1, function(x){ length(which(x>0 & x<=1)) })
base_data$d30_count <- apply(as.matrix(AAPdata[,4:dim(AAPdata)[2]]), 1, function(x){ length(which(x>=0 & x<=1)) })
base_data$index <- 1:dim(base_data)[1] ## row number
base_data$job_id <- index ## SGE_TASK_ID

base_data$sil3.Inc0_NTr <- -9
base_data$sil3.Exc0_NTr <- -9
base_data$sil3.Inc0_Tr <- -9
base_data$sil3.Exc0_Tr <- -9

base_data$best.case <- -9
base_data$best.sil <- -9
base_data$medoid1 <- -9
base_data$medoid2 <- -9
base_data$medoid3 <- -9

base_data$lower_bound <- -9
base_data$upper_bound <- -9

base_data$count_c1 <- -9
base_data$count_c2 <- -9
base_data$count_c3 <- -9

cat("2. Record Clustering Results\n")
startTime <- Sys.time()

for (i in 1:dim(AAPdata)[1]){
  cat(i, ",")
  chr <- AAPdata$chrom[i]
  pos <- AAPdata$end[i]
  tmp <- AAPdata[i, 4:dim(AAPdata)[2]]
  aap_count <- base_data$count[i]
  d30_count <- base_data$d30_count[i]
  
  x_ge0 <- as.numeric(t(tmp[which(!is.na(tmp) & tmp>=0)]))
  x_gt0 <- as.numeric(t(tmp[which(!is.na(tmp) & tmp>0)]))
  x_ge0_tr <- x_ge0[which( x_ge0>quantile(x_ge0, 0.05) & x_ge0<quantile(x_ge0, 0.95))]
  x_gt0_tr <- x_gt0[which( x_gt0>quantile(x_gt0, 0.05) & x_gt0<quantile(x_gt0, 0.95))]

  x_list <- list(x_ge0, x_gt0, x_ge0_tr, x_gt0_tr)
  x_list.name <- c("Inc0,NoTrim", "Exc0,NoTrim", "Inc0,Trim", "Exc0,Trim")

  sil.vec <- NULL
  for (j in 1:length(x_list)){
    x2 <- x_list[[j]]
    if ( length(x2)>3 ){
      clus3 <- pam(as.vector(x2), k=3)
      sil.vec <- c(sil.vec, clus3$silinfo$avg.width)
    } else {
      sil.vec <- c(sil.vec, -999)
    }
  }

  col_idx <- grep("^sil3", colnames(base_data))
  base_data[i, col_idx] <- round(sil.vec, 4)
  base_data$best.case[i] <- which.max(base_data[i, col_idx])
  base_data$best.sil[i] <- max(base_data[i, col_idx])
  
  final_x <- x_list[[base_data$best.case[i]]]
  final.clus <- pam(final_x, k=3)
  clus.order <- order(final.clus$medoids)

  min.vec <- apply(as.matrix(1:3), 1, function(k, clus=final.clus){ min(final_x[which(clus$clustering==k)]) })
  max.vec <- apply(as.matrix(1:3), 1, function(k, clus=final.clus){ max(final_x[which(clus$clustering==k)]) })
  min.vec <- min.vec[clus.order]
  max.vec <- max.vec[clus.order]

  med.vec <- final.clus$medoids[clus.order]

  base_data$medoid1[i] <- round(med.vec[1], 4)
  base_data$medoid2[i] <- round(med.vec[2], 4)
  base_data$medoid3[i] <- round(med.vec[3], 4)
  
  base_data$lower_bound[i] <- round((max.vec[1] + min.vec[2])/2, 4)
  base_data$upper_bound[i] <- round((max.vec[2] + min.vec[3])/2, 4)

  base_data$count_c1[i] <- length(which(!is.na(x_ge0) & x_ge0<=base_data$lower_bound[i]))
  base_data$count_c2[i] <- length(which(!is.na(x_ge0) & x_ge0>base_data$lower_bound[i] & x_ge0<=base_data$upper_bound[i]))
  base_data$count_c3[i] <- length(which(!is.na(x_ge0) & x_ge0>base_data$upper_bound[i]))
}
cat("\n")
endTime <- Sys.time()
print(paste("Start: ", startTime, sep=""))
print(paste("End: ", endTime, sep=""))
cat("Done\n")

############ 3. MARK SEGDUP
cat("3. Mark Segdup\n")
startTime <- Sys.time()

base_data$segdup <- FALSE
tmpfile <- paste(tmpDir, "tmp.chr", chrom, "_sti", index, ".txt", sep="")
write.table(base_data[,1:3], tmpfile, col.names=F, row.names=F, quote=F, sep="\t") ## write tmp file for segdup check

cmd <- paste("bedtools intersect -a ", tmpfile, " -b ", segdupfile, " -wao ", sep=" ")
segdupOverlaps <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
colnames(segdupOverlaps) <- c("chrom", "start", "end", "segdup.chrom", "segdup.start", "segdup.end", "overlaps")
tmp <- segdupOverlaps[which(segdupOverlaps$overlaps>0), 1:3]

if (dim(tmp)[1]>0){
  segdupPos <- paste(tmp$chrom, tmp$end, sep=":")
  base_data$segdup[which(paste(base_data$chrom, base_data$end, sep=":") %in% segdupPos)] <- TRUE
}

system(paste("rm -f ", tmpfile, sep=" "))

write.table(base_data, paste(tmpDir, outfile, sep=""), col.names=T, row.names=F, quote=F, sep="\t")
system(paste("gzip -f ", paste(tmpDir, outfile, sep=""), sep=" "))
system(paste("cp -f ", paste(tmpDir, outfile, ".gz", sep=""), paste(outDir, outfile, ".gz", sep=""), sep=" "))

endTime <- Sys.time()
print(paste("Start: ", startTime, sep=""))
print(paste("End: ", endTime, sep=""))
cat("Done\n")

