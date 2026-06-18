##### SCRIPT : shallowHRD_hg19 #####


args = commandArgs(trailingOnly=TRUE)

path_to_ratio_file = args[1]
inputPath = normalizePath(dirname(path_to_ratio_file))

NAMEEE_intermediary1 = gsub("\\", "/", normalizePath(path_to_ratio_file), fixed = TRUE) # for window path
NAMEEE_intermediary2 = sub("*/*.bam_ratio.txt", "", NAMEEE_intermediary1)

NAMEEE = sub(".*/", "", NAMEEE_intermediary2)

output_relative_Path = args[2]
outputPath = normalizePath(output_relative_Path)


path_to_cyto = args[3]
cytoFile = normalizePath(path_to_cyto)

# Optional custom gene list path (supports both arg4 and arg5 for compatibility)
custom_gene_file <- NULL
if (length(args) >= 5 && nzchar(args[5])) {
  custom_gene_file <- args[5]
} else if (length(args) >= 4 && nzchar(args[4])) {
  custom_gene_file <- args[4]
}

if (!is.null(custom_gene_file)) {
  if (file.exists(custom_gene_file)) {
    cat("Custom gene file detected:", custom_gene_file, "\n")
  } else {
    warning("Custom gene file does not exist and will be ignored: ", custom_gene_file)
    custom_gene_file <- NULL
  }
} else {
  cat("No custom gene file provided; using built-in hg19 gene list.\n")
}

continue_on_error <- function() { 
  paste() 
}
options(error=continue_on_error) 

cat("======================================================== \n")
cat("Sample : ", NAMEEE, "\n")
cat("======================================================== \n")


##### Cytoband_hg19 #####

cyt_Annot<-read.csv(cytoFile,header=T,fill=T)

cyt_Annot<-cyt_Annot[,1:6]
for(k in 2:dim(cyt_Annot)[[1]]){
  str1<-as.numeric(cyt_Annot[k-1,2:3])
  str2<-as.numeric(cyt_Annot[k,2:3])
  if(sum(abs(str1-str2))==0){cyt_Annot[k,1]<-0}
}

cyt_Annot[1:(dim(cyt_Annot)[[1]]-1),1]<-cyt_Annot[1:(dim(cyt_Annot)[[1]]-1),1]+cyt_Annot[2:dim(cyt_Annot)[[1]],1]
cyt_Annot[dim(cyt_Annot)[[1]],1]<-1
cyt_Annot<-cyt_Annot[-which(cyt_Annot[,1]==0),]
cyt_Annot[,1]<-cyt_Annot[,6]-cyt_Annot[,5]
cyt_Annot[,2]<-cyt_Annot[,2]+(cyt_Annot[,3]-1)/2

for(k in 2:dim(cyt_Annot)[[1]]){
  if(cyt_Annot[k,2]==cyt_Annot[k-1,2]){cyt_Annot[k,1]<-cyt_Annot[k,6]-cyt_Annot[k-1,5];cyt_Annot[k-1,1]<-0}
}

cyt_Annot<-cyt_Annot[-which(cyt_Annot[,1]==0),]
arm_len<-cyt_Annot[,1:2]
chr_len<-arm_len
chr_len[,2]<-trunc(chr_len[,2])
for(k in 2:dim(chr_len)[[1]]){
  if(chr_len[k,2]==chr_len[k-1,2]){chr_len[k,1]==chr_len[k-1,1]+chr_len[k,1]
    chr_len[k-1,1]<-0
  }
}
chr_len<-chr_len[-which(chr_len[,1]==0),]
chr_len<-chr_len[-dim(chr_len)[[1]],]
GG<-sum(chr_len[,1])
chr_p<-round(chr_len[,1]/GG,3)


##### GENE LIST LOADER #####

load_gene_list <- function(custom_path = NULL) {
  # Built-in gene list coordinates are hg19 and match shallowHRD assumptions
  built_in_genes <- data.frame(
    gene = c("DOCK7", "HCN3", "KLHL12", "RBBP5", "TSC22D2", "PIK3CA", "ANKRD17", "CDKN2AIP",
             "RTN4IP1", "AIM1", "INTS10", "PPP2R2A", "BRF2", "ZNF703", "MYC", "CD274_PDL1",
             "MTAP", "SEPHS1", "ZMIZ1", "WAPAL", "PTEN", "CSTF3", "BBS1", "CTTN", "CCND1",
             "ATG16L2", "INTS4", "CCDC77", "FOXM1", "YEATS4", "MDM2", "BRCA2", "C13orf23",
             "TRIM13", "SDCCAG1", "SNAP23", "IGF1R", "CYB5B", "P53", "ELAC2", "MAP2K4",
             "ERAL1", "NF1", "ERBB2", "BRCA1", "PHB", "SUPT4H1", "RAD51C", "GALK1",
             "AKAP8", "BRD4", "PIK3R2", "CCNE1", "NOSIP", "C20orf111", "ZNF217",
             "TSHZ2", "SAMD10", "PCNT"),
    chr = c(1, 1, 1, 1, 3, 3, 4, 4, 6, 6, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10,
            11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 15, 15, 16, 17, 17,
            17, 17, 17, 17, 17, 17, 17, 17, 19, 19, 19, 19, 19, 20, 20, 20, 20, 21),
    start = c(63037227, 155253447, 202878988, 205073188, 150155147, 178912013, 74031804,
              184368003, 107048506, 106914242, 19692253, 26189610, 37704083, 37555419,
              128751439, 5460548, 21834858, 13374861, 80828723, 88238281, 89677535,
              33144578, 66289588, 70263658, 69462583, 72533068, 77647740, 525161, 2976593,
              69769087, 69223209, 32932025, 39598114, 50581891, 50284154, 42806544,
              99349764, 69479345, 7581274, 12908137, 11985670, 27185057, 29563320,
              37864629, 41236847, 47486829, 56426497, 56791453, 73754412, 15477397,
              15394840, 18272658, 30309059, 50071269, 42831995, 52196994, 51850383,
              62608232, 47804876),
    stringsAsFactors = FALSE
  )

  if (is.null(custom_path) || !nzchar(custom_path)) {
    cat("Using built-in gene list (", nrow(built_in_genes), " genes).\n", sep = "")
    return(built_in_genes)
  }

  if (!file.exists(custom_path)) {
    warning("Custom gene list not found: ", custom_path, ". Falling back to built-in genes.")
    return(built_in_genes)
  }

  cat("Loading custom gene list from: ", custom_path, "\n", sep = "")

  custom_genes <- tryCatch(read.csv(custom_path, stringsAsFactors = FALSE),
                           error = function(e) NULL)

  if (is.null(custom_genes)) {
    warning("Unable to read custom gene list; falling back to built-in genes.")
    return(built_in_genes)
  }

  required_cols <- c("gene", "chr", "start")
  missing_cols <- setdiff(required_cols, colnames(custom_genes))
  if (length(missing_cols) > 0) {
    warning("Missing required custom gene columns: ", paste(missing_cols, collapse = ", "),
            ". Falling back to built-in genes.")
    return(built_in_genes)
  }

  custom_genes <- custom_genes[, required_cols]
  custom_genes$chr <- as.numeric(sub("^chr", "", as.character(custom_genes$chr), ignore.case = TRUE))
  custom_genes$start <- as.numeric(custom_genes$start)
  custom_genes <- custom_genes[!is.na(custom_genes$chr) & !is.na(custom_genes$start), ]
  # Exclude chrX/chrY/chrM; this no-chrX workflow uses autosomes only (chr 1-22)
  custom_genes <- custom_genes[custom_genes$chr >= 1 & custom_genes$chr <= 22, ]

  if (nrow(custom_genes) == 0) {
    warning("No valid custom genes after filtering chr/start; falling back to built-in genes.")
    return(built_in_genes)
  }

  cat("Using custom gene list with ", nrow(custom_genes), " genes.\n", sep = "")
  custom_genes
}

lookup_gene_ratios <- function(gene_list, B, C, THR) {
  out <- data.frame(
    gene = gene_list$gene,
    chr = gene_list$chr,
    start = gene_list$start,
    ratio_point_initial = NA_real_,
    CN_point_initial = NA_real_,
    ratio_segment_initial = NA_real_,
    CN_segment_initial = NA_real_,
    ratio_segment_final = NA_real_,
    CN_segment_final = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(gene_list))) {
    chr_i <- gene_list$chr[i]
    start_i <- gene_list$start[i]

    B_chr <- B[B$chr == chr_i, ]
    if (nrow(B_chr) > 0) {
      idx_initial <- which.min(abs(B_chr$start - start_i))
      initial_row <- B_chr[idx_initial, ]
      out$ratio_point_initial[i] <- initial_row$ratio
      out$CN_point_initial[i] <- round(initial_row$ratio / THR, 3)
      out$ratio_segment_initial[i] <- initial_row$ratio_median
      out$CN_segment_initial[i] <- round(initial_row$ratio_median / THR, 3)
    }

    C_chr <- C[C$chr == chr_i, ]
    if (nrow(C_chr) > 0) {
      idx_final <- which.min(abs(C_chr$start - start_i))
      final_row <- C_chr[idx_final, ]
      out$ratio_segment_final[i] <- final_row$ratio_median
      out$CN_segment_final[i] <- round(final_row$ratio_median / THR, 3)
    }
  }

  out
}

##### LIBRARIES #####

suppressMessages(library("ggpubr"))
suppressMessages(library("gridExtra"))
suppressMessages(library("DescTools"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("ks"))
suppressMessages(library("ggrepel"))

##### FUNCTIONS #####
## Function Local minima ##

localMinima <- function(x) { 
  y <- diff(c(.Machine$integer.max, x)) < 0L 
  rle(y)$lengths                             
  y <- cumsum(rle(y)$lengths)                 
  y <- y[seq.int(1L, length(y), 2L)]         
  if (x[[1]] == x[[2]]) {                    
    y <- y[-1]
  }
  y                                          
}


## Function Local maxima ##

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}


## Function readSegmFile ##

readSegmFile<-function(segFileName){
  aa<-scan(segFileName,nlines=1,what="character",quiet = T)
  k<-1
  
  while((length(grep("#",aa))!=0)&&(grep("#",aa)==1)){
    aa<-scan(segFileName,skip=k,nlines=1,what="character",quiet = T)
    k<-k+1
  }
  
  tmp<-read.delim(segFileName,skip=k-1, header=T,sep="\t")
  tmp<-data.matrix(tmp)
  if(k==2){
    results<-getInfoSegm(scan(segFileName,nlines=1,what="character",quiet = T))
    homoConst<-results$homoConst
    sdH<-results$sdH	
    p_BAF<-results$p_BAF
    q_LRR<-results$q_LRR
    Delta<-results$Delta
  }else{homoConst<-1;sdH<-NA;p_BAF<-NA;q_LRR<-NA;Delta<-c(2,NA)}
  
  out<-list(tmp=tmp,homoConst=homoConst,sdH=sdH,p_BAF=p_BAF,q_LRR=q_LRR,Delta=Delta)
}


## getSegmentID : gather in the same name segment of the length  ##

getSegmentID<-function(THR,tmp,c_chr,c_cn,c_conf){
  
  tmp[,c_conf]<-0;tmp[1,c_conf]<-1
  
  for(k in 2:dim(tmp)[[1]]){
    
    if(tmp[k,c_chr]==tmp[k-1,c_chr] & abs(tmp[k,c_cn]-tmp[k-1,c_cn])<THR){
      
      tmp[k,c_conf]<-tmp[k-1,c_conf]
      
    }else{tmp[k,c_conf]<-tmp[k-1,c_conf]+1}
  } 
  
  out<-tmp
  
}


## shrinkReprTMP

shrinkReprTMP<-function(tmp,c_posS,c_posE,c_cn,c_conf){
  
  B <- data.frame(ratio_file_tsv)
  B=B[,-1] 
  colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")
  
  GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                       seqnames.field = "chr",
                                                       start.field = "start", end.field = "end")
  
  for(k in 1:(dim(tmp)[[1]]-1)){
    
    if(tmp[k,c_conf]==tmp[k+1,c_conf]){
      
      tmp[k+1,c_posS]<-tmp[k,c_posS]
      
      # w<-c(tmp[k+1,c_posE]-tmp[k+1,c_posS],tmp[k,c_posE]-tmp[k,c_posS])
      
      # tmp[k+1,c_cn]<-weighted.mean(c(tmp[k+1,c_cn],tmp[k,c_cn]),w)
      
      gr=GRanges(seqnames=c(tmp[k+1,c_chr]),
                 ranges=IRanges(start=c(tmp[k+1,c_posS]),end=c(tmp[k+1,c_posE])),
                 strand=c("*"))
      subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
      
      tmp[k+1,c_cn]<-median(subsetGRobject$ratio)
      
      tmp[k,1]<-0
    }
  }
  
  tt<-which(tmp[,1]==0);if(length(tt)!=0){tmp<-tmp[-tt,]}
  out<-tmp
  
}


## Function getInfoSegm

getInfoSegm<-function(infoString){
  
  tt<-unlist(strsplit(infoString,";"))
  hh<-grep("homoConst=",tt)
  if(length(hh)==1){homoConst<-as.numeric(sub("homoConst=","",tt[hh]))}else{homoConst<-1}
  hh<-grep("sdH=",tt)
  if(length(hh)==1){sdH<-as.numeric(sub("sdH=","",tt[hh]))}else{sdH<-NA}
  
  hh<-grep("p=",tt)
  if(length(hh)==1){p_BAF<-as.numeric(sub("p=","",tt[hh]))}else{p_BAF<-NA}
  
  hh<-grep("q=",tt)
  if(length(hh)==1){q_LRR<-as.numeric(sub("q=","",tt[hh]))}else{q_LRR<-NA}
  
  hh<-grep("2LR=",tt)
  if(length(hh)==1){Delta<-as.numeric(sub("2LR=","",tt[hh]))}else{Delta<-NA}
  Delta<-c(2,Delta)
  
  
  if((is.na(homoConst))||(!is.numeric(homoConst))){homoConst<-1
  }else{if((homoConst<0.8)||(homoConst>1)){homoConst<-1}
  }
  if((is.na(sdH))||(!is.numeric(sdH))){sdH<-NA
  }else{if((sdH<0.05)||(sdH>1)){sdH<-NA}
  }
  
  out<-list(homoConst=homoConst,sdH=sdH,p_BAF=p_BAF,q_LRR=q_LRR,Delta=Delta)
  
}


## breakSmoothToLGA

breakSmoothToLGA<-function(THR,tmp,c_ind,c_chr,c_posS,c_posE,c_cn,c_CN,c_conf){
  
  B <- data.frame(ratio_file_tsv)
  B=B[,-1] 
  colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")
  
  GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                       seqnames.field = "chr",
                                                       start.field = "start", end.field = "end")
  
  tmp<-getSegmentID(THR=THR,tmp=tmp,c_chr=c_chr+1,c_cn=c_cn,c_conf=c_conf)
  tmp<-shrinkReprTMP(tmp=tmp,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  
  tmp[,c_ind]<-seq(1,dim(tmp)[[1]])
  
  BreaksSmall<-0
  
  # print("ARRIVE_HERE")
  
  pass=1
  
  FL<-T
  while(FL){
    
    tmp[,c_ind]<-seq(1,dim(tmp)[[1]])
    
    kk<-which(round((tmp[,c_posE]-tmp[,c_posS])/10^6,1)<3) # segments smaller than 3Mb
    
    if(length(kk)>0 & pass < 10000){
      
      kk<-kk[order((tmp[kk,c_posE]-tmp[kk,c_posS])/10^6)] 
      
      pass=pass+1
      
      for(k in 1:length(kk)){
        
        if(kk[k]==1){                                        
          if(tmp[kk[k]+1,c_chr+1]==tmp[kk[k],c_chr+1]){      
            if(tmp[kk[k]+1,c_ind]!=0){                       
              tmp[kk[k],c_ind]<-0                               
            }
          }                          
        }
        
        else{
          if(kk[k]==dim(tmp)[[1]]){                                                
            if(tmp[kk[k]-1,c_chr+1]==tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_ind]!=0){    
              tmp[kk[k],c_ind]<-0}                                                 
          }
          else{
            
            if(tmp[kk[k]-1,c_chr+1]==tmp[kk[k]+1,c_chr+1] & tmp[kk[k]-1,c_ind]!=0 & tmp[kk[k]+1,c_ind]!=0){tmp[kk[k],1]<-0}
            
            if(tmp[kk[k]-1,c_chr+1]==tmp[kk[k],c_chr+1] & tmp[kk[k]+1,c_chr+1]!=tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_ind]!=0){tmp[kk[k],c_ind]<-0}
            
            if(tmp[kk[k]+1,c_chr+1]==tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_chr+1]!=tmp[kk[k],c_chr+1] & tmp[kk[k]+1,c_ind]!=0){tmp[kk[k],c_ind]<-0}
            
            if(tmp[kk[k]+1,c_chr+1]!=tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_chr+1]!=tmp[kk[k],c_chr+1]){tmp[kk[k],c_ind]<-0}
            
          }
        }
      } 
      
      # print("ARRIVE_HERE 222222 \n")
      
      tt<-which(tmp[,c_ind]==0) 
      
      if(length(tt)>0){tmp <- tmp[-tt,];BreaksSmall<-BreaksSmall+length(tt)} 
      
      for(k in 1:(dim(tmp)[[1]]-1)){ 
        
        # print(k)
        
        if(round((tmp[k+1,c_posS]-tmp[k,c_posE])/10^6,1)<3){    
          
          if(tmp[k,c_chr+1]==tmp[k+1,c_chr+1] & abs(tmp[k,c_cn]-tmp[k+1,c_cn])<THR){ 
            
            tmp[k+1,c_posS]<-tmp[k,c_posS]        
            
            # w<-c(tmp[k+1,c_posE]-tmp[k+1,c_posS],tmp[k,c_posE]-tmp[k,c_posS]) 
            # 
            # tmp[k+1,c_cn]<-weighted.mean(c(tmp[k+1,c_cn],tmp[k,c_cn]),w) 
            
            gr=GRanges(seqnames=c(tmp[k+1,c_chr]),
                       ranges=IRanges(start=c(tmp[k+1,c_posS]),end=c(tmp[k+1,c_posE])),
                       strand=c("*"))
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp[k+1,c_cn]<-median(subsetGRobject$ratio)
            
            tmp[k,c_ind]<-0                                              
            
          }
          
          
        }
      }
      
      # print("kk")
      # print(length(kk)) # can't remove all small segments?
      # 
      # print("pass")
      # print(pass)
      
      tt<-which(tmp[,c_ind]==0);if(length(tt)>0){tmp<-tmp[-tt,]} 
      
      
    }else{FL<-F}  
  }
  
  
  
  tmp<-getSegmentID(THR=THR,tmp=tmp,c_chr=c_chr+1,c_cn=c_cn,c_conf=c_conf)
  tmp<-shrinkReprTMP(tmp=tmp,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  
  tmp[1,c_ind]<-BreaksSmall
  
  out<-tmp
  
  
}


## LGA_control

LGA_control<-function(THR,lenBIN,lenMB,tmp,c_ind,c_chr,c_posS,c_posE,c_cn,c_CN,c_conf){
  
  tmp[,c_ind]<-seq(1,dim(tmp)[[1]])
  
  
  WC<-tmp[,c_CN]*0
  
  j<-5
  
  for(k in 2:dim(tmp)[[1]]){
    
    if(tmp[k,c_chr+1]==tmp[k-1,c_chr+1] & round((tmp[k,c_posS]-tmp[k-1,c_posE])/10^6,1)<3){
      
      if(round((tmp[k,c_posE]-tmp[k,c_posS]+1)/10^6,0)>=lenMB & round((tmp[k-1,c_posE]-tmp[k-1,c_posS]+1)/10^6,0)>=lenMB){
        
        
        if(abs(tmp[k,c_cn]-tmp[k-1,c_cn])>THR*coefficient){WC[k]<-1}
        
      }
      
    }
  }
  
  
  out<-WC
  
  
}


##### Fast gathering 1 #####

cat("Fast gathering segments 1... \n")

dataTable <- read.table(paste0(inputPath,"/",NAMEEE,".bam_ratio.txt"), header = TRUE)
dataTable = dataTable[,1:4]

dataTable = dataTable[which(!dataTable[,1] == "X"),] 
dataTable = dataTable[which(!dataTable[,1] == "Y"),] 

dataTable[,1] = as.numeric(as.character(dataTable[,1]))

dataTable = dataTable[order(dataTable[,1]),]

L = dim(dataTable)[1]

size_window = dataTable[2,2] - dataTable[1,2]

dataTable = cbind(1:L, dataTable[,1], dataTable[,2], dataTable[,2]+size_window-1, 
                  dataTable[,3], dataTable[,4]) 

dataTable = as.data.frame(dataTable)
colnames(dataTable) = c("feature", "chromosome", "start", "end", "ratio", "ratio_median")

# log2 transformation
# dataTable[,5] = log2(dataTable[,5])
# dataTable[,6] = log2(dataTable[,6])

ratio <- data.frame(dataTable)

short_size_window = size_window/1000

X=ratio

ratio_file_tsv = ratio # save in the moment

X=X[,-1] # feature
X=X[,-4] 


## Remove spurious regions based on telomere and centromere from UCSC

X <- X[which(!(X[,1]==1 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==1 & X[,2]>=121535434 & X[,3]<=124535434)),]
X <- X[which(!(X[,1]==1 & X[,2]>=249240621 & X[,3]<=249250621)),]
X <- X[which(!(X[,1]==2 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==2 & X[,2]>=92326171 & X[,3]<=95326171)),]
X <- X[which(!(X[,1]==2 & X[,2]>=243189373 & X[,3]<=243199373)),]
X <- X[which(!(X[,1]==3 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==3 & X[,2]>=90504854 & X[,3]<=93504854)),]
X <- X[which(!(X[,1]==3 & X[,2]>=198012430 & X[,3]<=198022430)),]
X <- X[which(!(X[,1]==4 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==4 & X[,2]>=49660117 & X[,3]<=52660117)),]
X <- X[which(!(X[,1]==4 & X[,2]>=191144276 & X[,3]<=191154276)),]
X <- X[which(!(X[,1]==5 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==5 & X[,2]>=46405641 & X[,3]<=9405641)),]
X <- X[which(!(X[,1]==5 & X[,2]>=180905260 & X[,3]<=180915260)),]
X <- X[which(!(X[,1]==6 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==6 & X[,2]>=58830166 & X[,3]<=61830166)),]
X <- X[which(!(X[,1]==6 & X[,2]>=171105067 & X[,3]<=171115067)),]
X <- X[which(!(X[,1]==7 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==7 & X[,2]>=58054331 & X[,3]<=61054331)),]
X <- X[which(!(X[,1]==7 & X[,2]>=159128663 & X[,3]<=159138663)),]
X <- X[which(!(X[,1]==8 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==8 & X[,2]>=43838887 & X[,3]<=46838887)),]
X <- X[which(!(X[,1]==8 & X[,2]>=146354022 & X[,3]<=146364022)),]
X <- X[which(!(X[,1]==9 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==9 & X[,2]>=47367679 & X[,3]<=50367679)),]
X <- X[which(!(X[,1]==9 & X[,2]>=141203431 & X[,3]<=141213431)),]
X <- X[which(!(X[,1]==10 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==10 & X[,2]>=39254935 & X[,3]<=42254935)),]
X <- X[which(!(X[,1]==10 & X[,2]>=135524747 & X[,3]<=135534747)),]
X <- X[which(!(X[,1]==11 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==11 & X[,2]>=51644205 & X[,3]<=54644205)),]
X <- X[which(!(X[,1]==11 & X[,2]>=134996516 & X[,3]<=135006516)),]
X <- X[which(!(X[,1]==12 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==12 & X[,2]>=34856694 & X[,3]<=37856694)),]
X <- X[which(!(X[,1]==12 & X[,2]>=133841895 & X[,3]<=133851895)),]
X <- X[which(!(X[,1]==13 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=115159878 & X[,3]<=115169878)),]
X <- X[which(!(X[,1]==14 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=107339540 & X[,3]<=107349540)),]
X <- X[which(!(X[,1]==15 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=17000000 & X[,3]<=20000000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=102521392 & X[,3]<=102531392)),]
X <- X[which(!(X[,1]==16 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==16 & X[,2]>=35335801 & X[,3]<=38335801)),]
X <- X[which(!(X[,1]==16 & X[,2]>=90344753 & X[,3]<=90354753)),]
X <- X[which(!(X[,1]==17 & X[,2]>=22263006 & X[,3]<=25263006)),]
X <- X[which(!(X[,1]==18 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==18 & X[,2]>=15460898 & X[,3]<=18460898)),]
X <- X[which(!(X[,1]==18 & X[,2]>=78067248 & X[,3]<=78077248)),]
X <- X[which(!(X[,1]==19 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==19 & X[,2]>=24681782 & X[,3]<=27681782)),]
X <- X[which(!(X[,1]==19 & X[,2]>=59118983 & X[,3]<=59128983)),]
X <- X[which(!(X[,1]==20 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==20 & X[,2]>=26369569 & X[,3]<=29369569)),]
X <- X[which(!(X[,1]==20 & X[,2]>=63015520 & X[,3]<=63025520)),]
X <- X[which(!(X[,1]==21 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==21 & X[,2]>=11288129 & X[,3]<=14288129)),]
X <- X[which(!(X[,1]==21 & X[,2]>=48119895 & X[,3]<=48129895)),]
X <- X[which(!(X[,1]==22 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=13000000 & X[,3]<=16000000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=51294566 & X[,3]<=51304566)),]


## Add chromosome arm 

options(show.error.messages = FALSE)

X$chr_arm <- rep(0,nrow(X))
colnames(X) <- c("chr", "start", "end", "ratio_median", "chr_arm")
X=X[c("chr", "chr_arm", "start", "end", "ratio_median")]

X[,3]=as.numeric(as.character(X[,3])) 
X[,4]=as.numeric(as.character(X[,4])) 

tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] < 125000000
X[tt,2] = 1
tt <- X[,1] == 1 & X[,3] > 125000000 & X[,4] > 125000000
X[tt,2] = 2
tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] > 125000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 125000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 125000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])  
rownames(X) <- NULL 

tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] < 93300000
X[tt,2] = 1
tt <- X[,1] == 2 & X[,3] > 93300000 & X[,4] > 93300000
X[tt,2] = 2
tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] > 93300000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 93300000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 93300001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] < 91000000
X[tt,2] = 1
tt <- X[,1] == 3 & X[,3] > 91000000 & X[,4] > 91000000
X[tt,2] = 2
tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] > 91000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 91000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 91000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] < 50400000
X[tt,2] = 1
tt <- X[,1] == 4 & X[,3] > 50400000 & X[,4] > 50400000
X[tt,2] = 2
tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] > 50400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 50400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 50400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] < 48400000
X[tt,2] = 1
tt <- X[,1] == 5 & X[,3] > 48400000 & X[,4] > 48400000
X[tt,2] = 2
tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] > 48400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 48400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 48400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] < 61000000
X[tt,2] = 1
tt <- X[,1] == 6 & X[,3] > 61000000 & X[,4] > 61000000
X[tt,2] = 2
tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] > 61000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 61000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 61000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] < 59900000
X[tt,2] = 1
tt <- X[,1] == 7 & X[,3] > 59900000 & X[,4] > 59900000
X[tt,2] = 2
tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] > 59900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 59900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 59900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] < 45600000
X[tt,2] = 1
tt <- X[,1] == 8 & X[,3] > 45600000 & X[,4] > 45600000
X[tt,2] = 2
tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] > 45600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 45600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 45600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] < 49000000
X[tt,2] = 1
tt <- X[,1] == 9 & X[,4] > 49000000 & X[,4] > 49000000
X[tt,2] = 2
tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] > 49000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 49000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 49000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] < 40200000
X[tt,2] = 1
tt <- X[,1] == 10 & X[,3] > 40200000 & X[,4] > 40200000
X[tt,2] = 2
tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] > 40200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 40200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 40200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] < 53700000
X[tt,2] = 1
tt <- X[,1] == 11 & X[,3] > 53700000 & X[,4] > 53700000
X[tt,2] = 2
tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] > 53700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 53700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 53700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] < 35800000
X[tt,2] = 1
tt <- X[,1] == 12 & X[,3] > 35800000 & X[,4] > 35800000
X[tt,2] = 2
tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] > 35800000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 35800000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 35800001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] < 17900000
X[tt,2] = 1
tt <- X[,1] == 13 & X[,3] > 17900000 & X[,4] > 17900000
X[tt,2] = 2
tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] > 17900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] < 17600000
X[tt,2] = 1
tt <- X[,1] == 14 & X[,3] > 17600000 & X[,4] > 17600000
X[tt,2] = 2
tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] > 17600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] < 19000000
X[tt,2] = 1
tt <- X[,1] == 15 & X[,3] > 19000000 & X[,4] > 19000000
X[tt,2] = 2
tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] > 19000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 19000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 19000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] < 36600000
X[tt,2] = 1
tt <- X[,1] == 16 & X[,3] > 36600000 & X[,4] > 36600000
X[tt,2] = 2
tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] > 36600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 36600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 36600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] < 24000000
X[tt,2] = 1
tt <- X[,1] == 17 & X[,3] > 24000000 & X[,4] > 24000000
X[tt,2] = 2
tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] > 24000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 24000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 24000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] < 17200000
X[tt,2] = 1
tt <- X[,1] == 18 & X[,3] > 17200000 & X[,4] > 17200000
X[tt,2] = 2
tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] > 17200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] < 26500000
X[tt,2] = 1
tt <- X[,1] == 19 & X[,3] > 26500000 & X[,4] > 26500000
X[tt,2] = 2
tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] > 26500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 26500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 26500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] < 27500000
X[tt,2] = 1
tt <- X[,1] == 20 & X[,3] > 27500000 & X[,4] > 27500000
X[tt,2] = 2
tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] > 27500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 27500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 27500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] < 13200000
X[tt,2] = 1
tt <- X[,1] == 21 & X[,3] > 13200000 & X[,4] > 13200000
X[tt,2] = 2
tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] > 13200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 13200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 13200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] < 14700000
X[tt,2] = 1
tt <- X[,1] == 22 & X[,3] > 14700000 & X[,4] > 14700000
X[tt,2] = 2
tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] > 14700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 14700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 14700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 


X[,2]<-X[,1]+(X[,2]-1)/2      

tt=which(X$ratio_median==-Inf)

if (length(tt) >= 1){
  X=X[-tt,]  
}


X=as.matrix(X)

L = dim(X)[1]
A = matrix(0, ncol=6, nrow=L)

i=1    
c=1    

while (i < L){
  if(X[i,2]==X[i+1,2]){             
    if (X[i,5]==X[i+1,5]){          
      n=1
      while (X[i,2]==X[i+n,2] & X[i,5]==X[i+n,5] & i+n < L){
        n=n+1
      }
      if(i+n == L){
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n,4], X[i,5], X[i+n,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
      else{
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n-1,4], X[i,5], X[i+n-1,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
    }
    else{
      A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
      i=i+1
      c=c+1
    }
  }
  else{
    A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
    i=i+1                                        
    c=c+1                         
  }  
}

A=subset(A, A[,1] != 0)
rownames(A) <- NULL 
colnames(A) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")


write.table(A, file = paste0(outputPath,"/",NAMEEE,"_ratio_median_gathered.txt"), sep = "\t", row.names = FALSE)

options(show.error.messages = TRUE)


##### Find Threshold 1 #####

cat("Find Threshold 1... \n")

X = read.table(paste0(outputPath,"/",NAMEEE,"_ratio_median_gathered.txt"), sep = "\t", header = TRUE)

X = X[which(X$chr != 23),]


# true value

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")


GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")
i=1
L=dim(X)[1]

while (i < L){
  gr=GRanges(seqnames=c(X[i,1]),
             ranges=IRanges(start=c(X[i,3]),end=c(X[i,4])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  X[i,5] = median(subsetGRobject$ratio)
  i=i+1
}


X = X[which(X[,6] > 2999999),]

L=dim(X)[1]
X=data.matrix(X)

test=c()
i=1
for (i in 1:L){
  v=i+1
  for (y in v:L-1){
    test = c(test, abs(X[i,5] - X[y,5]))}
}  


list_THR = c()
c = 0

nb_simulations = 100000


while (c < nb_simulations){
  n = length(test)
  n
  
  rand_n = sample(2:n,1)
  rand_n
  
  sampled_test = sample(test,rand_n)
  
  ## First maxima
  
  maxx = localMaxima(density(sampled_test)$y)[1]
  First_local_maxima = density(sampled_test)$x[maxx]
  
  
  ## First minima 
  
  minx = localMinima(density(sampled_test)$y)
  First_local_minima = density(sampled_test)$x[minx][which(density(sampled_test)$x[minx] > First_local_maxima)][1]
  
  
  ## Second maxima
  
  maxx = localMaxima(density(sampled_test)$y)[2]
  Second_local_maxima = density(sampled_test)$x[maxx]
  
  if (is.na(Second_local_maxima) == FALSE){
    if (abs((First_local_minima-First_local_maxima)-(Second_local_maxima-First_local_minima)) < 0.10){
      list_THR = c(list_THR,First_local_minima)
      test_optimized = sampled_test
    }
  }
  c = c + 1
}



if (is.null(list_THR) == TRUE){
  
  ## First maxima
  
  maxx = localMaxima(density(test)$y)[1]
  First_maxima = density(test)$x[maxx]
  
  ## First minima
  
  minx = localMinima(density(test)$y)
  First_local_minima = density(test)$x[minx][which(density(test)$x[minx] > First_maxima)][1]
  
  h <- hns(test, deriv.order = 2)
  den <- kdde(test, h=h, deriv.order = 2)
  
  flections<-c()
  for(i in 2:length(den$estimate)){
    if(sign(den$estimate[i])!=sign(den$estimate[i-1])){
      flections<-c(flections,i)
    }
  }
  
  All_inflexion_point = den$x[flections][which(den$x[flections] > First_maxima)]
  Threshold_change_sign_second_derivative_2nd_value = sort(All_inflexion_point)[1]
  
  Threshold = min(Threshold_change_sign_second_derivative_2nd_value, First_local_minima)
  
}


if (is.na(mean(list_THR)) == FALSE){
  number_positive = length(list_THR)
  
  Threshold = mean(list_THR)
  
  Threshold_first_round_simulation = Threshold
}



# THR_intermediary

test = as.data.frame(test)

ploty <- ggplot(test, aes(x = test)) + geom_density() +
  geom_vline(aes(xintercept = Threshold), color = "blue") + 
  ggtitle(paste("Density pairwise difference large segment")) + xlab("difference") +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=13),
        axis.text.y = element_text(size=15),
        axis.title.y = element_blank(),
        panel.background = element_blank())
ggsave(paste0(outputPath,"/",NAMEEE,"_THR_intermediary",".jpeg"), plot = ploty, device = "jpeg", width = 20, height = 10)


##### Reading and initialization 1 #####

cat("Reading and initialisation 1... \n")

segFiles <- list.files(outputPath, pattern = paste0(NAMEEE, "_ratio_median_gathered.txt"),full.names=T)

nSample<-1

coefficient = 1.0

ColNamesTMP<-c("index", "chr","chr_arm","posStart", "posEnd", "ratio")

c_ind<-1; c_chr<-2;c_posS<-4; c_posE<-5; c_cn<-6; c_conf<-8

LLBB<-c(3,4,5,6,7,8,9,10,11)          # size LGA tested

TAB<-matrix(0,nSample,20)              # matrix to stock results 
rownames(TAB)<-seq(1:nSample)          
colnames(TAB)<-seq(1,20)               

colnames(TAB)[1:6]<-c("p_BAF","q_LRR","2copyLRR","DNA_ind","All_breaks","less3Mb_breaks")

ii<-1

rownames(TAB)[ii]<-"result"

results<-readSegmFile(segFileName=segFiles[ii])  
tmp<-results$tmp                                  

tmp <- cbind(seq(1,dim(tmp)[[1]]),tmp) 
colnames(tmp)[1]<-c("index")

tmp <- cbind(tmp,rep(0,nrow(tmp))) 
colnames(tmp)[8]<-c("level")

graphe_I_tab = tmp



# recalculate inital

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")
i=1
L=dim(tmp)[1]

while (i < L){
  gr=GRanges(seqnames=c(tmp[i,2]),
             ranges=IRanges(start=c(tmp[i,4]),end=c(tmp[i,5])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  tmp[i,6] = median(subsetGRobject$ratio)
  i=i+1
}




tt<-which(tmp[,c_chr+1]>23.6)                            # exclusion chr24 and more   
if(length(tt)>0){tmp<-tmp[-tt,]}

tt<-which(tmp[,c_chr+1]==21)                             # exclusion short arm chr21
if(length(tt)>0){tmp<-tmp[-tt,]}

tt<-which(tmp[,c_chr+1]==22)                             # exclusion short arm chr22
if(length(tt)>0){tmp<-tmp[-tt,]}

graphe_II_tab = tmp

tmp[,7] = tmp[,5] - tmp[,4] + 1

tmp_3mb = tmp[which(tmp[,7] > 2999999),]   
cop_tmp = tmp_3mb                    
THR=Threshold

if (THR > 0.45){
  THR = 0.45}

if (THR < 0.025){
  THR = 0.025
}

level=1


##### Gather big segment (> 3Mb) 1 ####

cat("Gather big segments 1... \n")

options(show.error.messages = FALSE)

while (dim(cop_tmp)[1] > 1){  
  
  line_largest_index = which.max(cop_tmp[,7])             # biggest segment
  largest_index = cop_tmp[line_largest_index,1]           # associated index
  ratio_largest = cop_tmp[line_largest_index,6]           # median_ratio segment
  ratio_largest = ratio_largest + 0.00001                 # avoid picking same value
  
  close_index = cop_tmp[which.min(abs(cop_tmp[,6] - ratio_largest)),1]
  closest_index = c(largest_index, close_index)
  
  ## Find all value that matches the two closest index fit with THR 
  
  L=dim(cop_tmp)[1]
  
  for (i in 1:L){
    validation = 0
    n = length(closest_index)
    for (g in 1:n){
      if (abs(cop_tmp[i,6] - cop_tmp[which((cop_tmp[,1]==closest_index[g])),6]) < THR){   
        validation = validation + 1}
    }
    if (validation == n){
      closest_index = c(closest_index, cop_tmp[i,1])}}
  
  closest_index = unique(closest_index)
  
  ##  Assign level based on the index of close segment in genome    
  
  n = length(closest_index)
  for (i in 1:n){
    tmp_3mb[which((tmp_3mb[,1] == closest_index[i])),8] = level
    cop_tmp <- cop_tmp[which(!(cop_tmp[,1]==closest_index[i])),]
  }
  dim(tmp_3mb)
  dim(cop_tmp)
  level = level + 1
} 


B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")

##  Gather with level

options(show.error.messages = TRUE)

L1 = dim(tmp_3mb)[1]

A = matrix(0, ncol=8, nrow=L1)

i=1
c=1   
while (i<L1+1){
  if (i==L1){
    A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])    # segment pas sur deux lignes, a stocker seul 
    i=i+1
    c=c+1}
  else{
    if(tmp_3mb[i,3] == tmp_3mb[i+1,3]){                                           
      if (tmp_3mb[i,8] == tmp_3mb[i+1,8]){                                        
        n=1
        vector_chr = c(tmp_3mb[i,2])
        vector_segment_start = c(tmp_3mb[i,4])
        vector_segment_end = c(tmp_3mb[i,5])
        vector_strand = c("*")
        
        while (tmp_3mb[i,8] == tmp_3mb[i+n,8] && tmp_3mb[i,3] == tmp_3mb[i+n,3]){ # while same chr_arm & level 
          vector_chr = c(vector_chr, tmp_3mb[i+n,2])
          vector_segment_start = c(vector_segment_start, tmp_3mb[i+n,4])
          vector_segment_end = c(vector_segment_end, tmp_3mb[i+n,5])
          vector_strand = c(vector_strand, "*")
          n = n + 1
          if(i+n == L1+1){
            break}}
        
        gr=GRanges(seqnames=vector_chr,
                   ranges=IRanges(start=vector_segment_start,end=vector_segment_end),
                   strand=vector_strand) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr) # overlapping
        
        A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i+n-1,5], 
                median(subsetGRobject$ratio), tmp_3mb[i+n-1,5]- tmp_3mb[i,4]+1, tmp_3mb[i,8]) # recalculating new_median
        i=i+n
        c=c+1}
      else{                                            
        A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])    # segment pas sur deux lignes, a stocker seul 
        i=i+1
        c=c+1}}
    else{                                              
      A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])
      i=i+1                                                              
      c=c+1}
  }  
}

A=subset(A, A[,1] != 0)

rownames(A) <- NULL 
colnames(A) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

tmp_3mb=A

graphe_III_tab = tmp_3mb 


##### Reput small segment & Smoothing to final segmentation 1 #####

cat("Reput small segments 1... \n")

options(show.error.messages = FALSE)

tmp_0.1_3mb=tmp[which(tmp[,7] > 99999),]
tmp_0.1_3mb=tmp_0.1_3mb[which(tmp_0.1_3mb[,7] < 2999999),]  

L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

values = unique(tmp_0.1_3mb[,3][!tmp_0.1_3mb[,3] %in% tmp_3mb[,3]])  
length=length(values)

for (missing_chr_arm in values){
  c=1
  L_3mb=dim(tmp_3mb)[1]
  while(c < L_3mb+1){
    if (missing_chr_arm > max(tmp_3mb[,3])){
      tmp_3mb=rbind(tmp_3mb, tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),])}
    else if (missing_chr_arm < tmp_3mb[c,3]){
      if (c==1){
        tmp_3mb=rbind(tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),] , tmp_3mb[1:L_3mb,])
        c=L_3mb+1}
      else{
        tmp_3mb=rbind(tmp_3mb[1:(c-1),], tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),] , tmp_3mb[c:L_3mb,])}
      c=L_3mb+1}
    else{
      c=c+1
    }
  }
}

for (missing_chr_arm in values){
  tmp_0.1_3mb <- tmp_0.1_3mb[which(!(tmp_0.1_3mb[,3]==missing_chr_arm)),]}


B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")

## initialisation

halt = 0 

i=1
c=1

L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

while (halt < 1){  
  if (tmp_0.1_3mb[i,3] > tmp_3mb[c,3] | (tmp_0.1_3mb[i,3] == tmp_3mb[c,3] && tmp_0.1_3mb[i,4] > tmp_3mb[c,5])){
    halt = halt + 1
  }
  
  else{                                            
    if (tmp_0.1_3mb[i,5] < tmp_3mb[c,4]){                       ########## segment avant
      if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR){          # No gathering 
        tmp_3mb=rbind(c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_0.1_3mb[i,5], 
                        tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]) , tmp_3mb[c:L_3mb,])
        i=i+1
      }
      
      else{                                                    # gathering
        
        gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                   ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                   strand=c("*")) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
        
        
        tmp_3mb=rbind(c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                        median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
        i=i+1}
    }
    else if (tmp_0.1_3mb[i,4] < tmp_3mb[c,4] && tmp_0.1_3mb[i,5] >= tmp_3mb[c,4]  && tmp_0.1_3mb[i,5] <= tmp_3mb[c,5]){  ########## segment overlapping mais debut avant
      if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR) {
        gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                   ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,4])),
                   strand=c("*")) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
        tmp_3mb=rbind(c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_3mb[c,4], 
                        median(subsetGRobject$ratio),tmp_3mb[c,4]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]), tmp_3mb[c:L_3mb,])
        i=i+1}
      else{
        gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                   ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                   strand=c("*")) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
        
        tmp_3mb=rbind(c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                        median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
        i=i+1}  
    }
    else{
      i=i+1
    }
  }
}

tmp_0.1_3mb = tmp_0.1_3mb[i:l_0.1_3mb,]  



L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

i=1
c=1

while (i < l_0.1_3mb + 1){          
  c=1                                                                            # reset 
  L_3mb=dim(tmp_3mb)[1]
  while (c < L_3mb+1){                                                            
    if (tmp_0.1_3mb[i,3] == tmp_3mb[c,3]){                                            
      if (tmp_0.1_3mb[i,5] < tmp_3mb[c,4]){                         ############### 1 - Before all segments chr_arm
        
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR){                     # not gathering with THR
          tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_0.1_3mb[i,5], 
                                             tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]) , tmp_3mb[c:L_3mb,])
          i=i+1
          c=L_3mb+1}
        else{                                                                # gathering
          gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                     ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                     strand=c("*")) # to overlap
          subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
          
          tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                                             median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
          i=i+1
          c=L_3mb+1}
      }
      
      else if (tmp_0.1_3mb[i,4] < tmp_3mb[c,4] && tmp_0.1_3mb[i,5] >= tmp_3mb[c,4]  && tmp_0.1_3mb[i,5] <= tmp_3mb[c,5] && tmp_0.1_3mb[i,4] >= tmp_3mb[c-1,5]){    ########## 2 - Inside segments beginning before but after previous segment
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR) {
          if (c+1 > L_3mb) {
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,4])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_3mb[c,4], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,4]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]), tmp_3mb[c:L_3mb,])
            i=i+1
            c=L_3mb+1}
          else{
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,4])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_3mb[c,4], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,4]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]))
            i=i+1
            c=L_3mb+1}
        }
        
        else{
          if (c+1 > L_3mb) {
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]))
            i=i+1
            c=L_3mb+1
          }
          else{
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}
        }
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,4] && tmp_0.1_3mb[i,5] <= tmp_3mb[c,5]){    ########## 3 - Inside segments pile
        i=i+1
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,4] && tmp_0.1_3mb[i,4] <= tmp_3mb[c,5] && tmp_0.1_3mb[i,5] > tmp_3mb[c,5] && tmp_0.1_3mb[i,5] <= tmp_3mb[c+1,5]){    ########## 4 - Inside segments end after but before next segment
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR) {
          if (c+1 > L_3mb) {
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_0.1_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c),], c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_3mb[c,5], tmp_0.1_3mb[c,5], 
                                             median(subsetGRobject$ratio),tmp_0.1_3mb[c,5]-tmp_3mb[c,5]+1,tmp_0.1_3mb[i,8]))
            i=i+1
            c=L_3mb+1 
          }
          else{
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_0.1_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c),], c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_3mb[c,5], tmp_0.1_3mb[c,5], 
                                             median(subsetGRobject$ratio),tmp_0.1_3mb[c,5]-tmp_3mb[c,5]+1,tmp_0.1_3mb[i,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1
          }
        }
        
        else{
          if (c+1 > L_3mb){
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], 
                                               median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]))
            i=i+1
            c=L_3mb+1
          }
          else{
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], 
                                               median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1
          }
        }
      }
      
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,5] && tmp_0.1_3mb[i,5] <= tmp_3mb[c+1,4] && tmp_3mb[c,3] == tmp_3mb[c+1,3]){     ########## 5 - between two segments same chr_arm
        
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR && abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6]) <= THR){       # two in range THR
          if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6])){     # c closer
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], median(subsetGRobject$ratio),
                                                tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]) , tmp_3mb[(c+1):L_3mb,])  
            i=i+1
            c=L_3mb+1}
          
          else {                      # c+1 closer
            if (c+2 > L_3mb){         # anticipate end of file
              gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8])) 
              i=i+1
              c=L_3mb+1}
            else{
              gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]) , tmp_3mb[(c+2):L_3mb,])  
              i=i+1
              c=L_3mb+1}
          }
        }
        
        else if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR){                   # only c
          gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                     ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                     strand=c("*")) # to overlap
          subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
          
          tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], median(subsetGRobject$ratio),
                                              tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]) , tmp_3mb[(c+1):L_3mb,])  
          i=i+1
          c=L_3mb+1}
        
        else if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6]) <= THR){                  # only c+1
          if (c+2 > L_3mb){                                                       # anticipate end of file
            gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                            tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]))  
            i=i+1
            c=L_3mb+1}
          
          else{
            if (tmp_0.1_3mb[i+1,5] < tmp_3mb[c+1,4] && tmp_0.1_3mb[i+1,3] == tmp_3mb[c+1,3]){  # small segments afterward
              gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_0.1_3mb[i,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_0.1_3mb[i,5], median(subsetGRobject$ratio),
                                              tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1,tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])  
              i=i+1
              c=L_3mb+1}
            else{
              gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]) , tmp_3mb[(c+2):L_3mb,])  
              i=i+1
              c=L_3mb+1}
          }
        }
        
        else{ # aucun                                      
          tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                          tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])
          i=i+1
          c=L_3mb+1}
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,5] && (tmp_3mb[c,3] != tmp_3mb[c+1,3] | is.null(tmp_3mb[c+1,3]))){     ############## 6 - after last segment chr_arm
        
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR){   # in range THR
          if (c+1 > L_3mb){                                 # anticipate end of file
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3],tmp_3mb[c,4],
                                                tmp_0.1_3mb[i,5],median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8])) 
            i=i+1
            c=L_3mb+1}
          else{                                             # not end of file
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3],tmp_3mb[c,4],
                                                tmp_0.1_3mb[i,5],median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}}
        
        else{                                       # not in range THR
          if (c+1 > L_3mb){                         # anticipate end of file
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                            tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]))
            i=i+1
            c=L_3mb+1}
          else{                                     # not end of file
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                            tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}}}
      
      else {             # none of the cases listed up there
        c=c+1
      }
    }    
    else{                                            
      c=c+1
    }
  }
}

rownames(tmp_3mb) <- NULL



## add leftovers

L_3mb=dim(tmp_3mb)[1]
c = 1

while (c < L_3mb){
  if (tmp_3mb[c,3] == tmp_3mb[c+1,3]){
    if (tmp_3mb[c+1,4] - tmp_3mb[c,5] > 1){
      
      gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                 ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_3mb[c+1,4])),
                 strand=c("*")) # to overlap
      subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
      
      if (abs(tmp_3mb[c,6] - median(subsetGRobject$ratio)) < abs(tmp_3mb[c+1,6] - median(subsetGRobject$ratio))){
        tmp_3mb[c,5] = tmp_3mb[c+1,4] - 1
      }
      else{
        tmp_3mb[c+1,4] = tmp_3mb[c,5] + 1
      }
    }
  }
  c=c+1
}

rownames(tmp_3mb) <- NULL

options(show.error.messages = TRUE)

colnames(tmp_3mb) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

tmp_3mb[,7] = tmp_3mb[,5] - tmp_3mb[,4] + 1

graphe_IV_tab = tmp_3mb

write.table(graphe_IV_tab, file = paste0(outputPath,"/",NAMEEE,"_IV.txt"), sep = "\t", row.names = FALSE)



## V

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")


tmp_3mb<-breakSmoothToLGA(THR,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)


# add leftovers

rownames(tmp_3mb) <- NULL

L_3mb=dim(tmp_3mb)[1]
c = 1

while (c < L_3mb){
  if (tmp_3mb[c,3] == tmp_3mb[c+1,3]){
    if (tmp_3mb[c+1,4] - tmp_3mb[c,5] > 1){
      
      gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                 ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_3mb[c+1,4])),
                 strand=c("*")) # to overlap
      subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
      
      if (abs(tmp_3mb[c,6] - median(subsetGRobject$ratio)) < abs(tmp_3mb[c+1,6] - median(subsetGRobject$ratio))){
        tmp_3mb[c,5] = tmp_3mb[c+1,4] - 1
      }
      else{
        tmp_3mb[c+1,4] = tmp_3mb[c,5] + 1
      }
    }
  }
  c=c+1
}

rownames(tmp_3mb) <- NULL

colnames(tmp_3mb) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

tmp_3mb[,7] = tmp_3mb[,5] - tmp_3mb[,4] + 1



# go to begin and end of chromosome

B <- data.frame(ratio_file_tsv)
B=B[,-1]
B=B[,-5]
colnames(B) <- c("chr", "start", "end", "readcount")
attach(B)

df=merge(aggregate(start ~ chr, B, min), aggregate(end ~ chr, B, max))[seq(from=1, to=22, by = 1),]
df=as.data.frame(df)
detach(B)

# first
tmp_3mb[1,4] = df[1,2]

# inside
L_3mb=dim(tmp_3mb)[1]
c = 1
n = 1

while (c < L_3mb){
  if (tmp_3mb[c,2] != tmp_3mb[c+1,2]){ # m?me chr
    tmp_3mb[c,5] = df[n,3]
    n=n+1
    tmp_3mb[c+1,4] = df[n,2]     
  }
  c=c+1
}

# final
tmp_3mb[L_3mb,5] = df[22,3]

tmp_3mb[,7] = tmp_3mb[,5] - tmp_3mb[,4] + 1


# recalculate complete value

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")

c=1
L_3mb=dim(tmp_3mb)[1]

while (c < L_3mb){
  gr=GRanges(seqnames=c(tmp_3mb[c,2]),
             ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_3mb[c,5])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  tmp_3mb[c,6] = median(subsetGRobject$ratio)
  c=c+1
}


graphe_V_tab = tmp_3mb

## graphe representative of the final segmentation 

test_data_frame = as.data.frame(tmp_3mb)
colnames(test_data_frame) <- make.unique(names(test_data_frame))

test_ordered = test_data_frame[order(test_data_frame$ratio_median),]
test_ordered$num_line <- seq.int(nrow(test_ordered))

test_ploty_CN_level <- ggplot(test_ordered, aes(x = ratio_median, y = num_line)) + 
  geom_point(shape = 1, size = 1, color = "#0072B2", na.rm = TRUE) + geom_line() + ggtitle("Final segmentation diagnostic") + 
  xlab("ratio median") + xlim(-2.5,2.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())
suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_final_segmentation_visual",".jpeg"), plot = test_ploty_CN_level, device = "jpeg", width = 20, height = 10))


##### Call LGAs 1 #####

cat("Call Large Genomic Alterations 1... \n")

LGAs_data_frame <- as.data.frame(matrix(0, ncol = 2, nrow = 9))
LGAs_data_frame[,1] = c(3:11)
colnames(LGAs_data_frame) <- c("Size_LGA", "Number_LGA")


tmp_3mb = data.frame(tmp_3mb)

tmp_3mb = tmp_3mb[which(tmp_3mb$chr != 23),]

tmp_3mb = as.matrix(tmp_3mb)


for (i in (3:11)){
  WC<-LGA_control(THR,lenBIN=500,lenMB=i,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  LGAs_data_frame[i-2,2] = sum(WC[,1])}
write.table(LGAs_data_frame, file = paste0(outputPath,"/",NAMEEE,"_number_LGAs.txt"), sep = "\t", row.names = FALSE)


# For graphe 10Mb 

WC<-LGA_control(THR,lenBIN=500,lenMB=10,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)

test=cbind(tmp_3mb,WC[,1])
colnames(test) <- c("index", "chr","chr_arm", "start", "end", "ratio_median", "size", "level", "WC")


L=dim(test)[1]

i=1
while(i + 1 < L){
  L=dim(test)[1]
  if (test[i,9] == 0 & test[i+1,9] != 1){
    test=test[-i,]}
  else{
    i=i+1}}

l=dim(test)[1]

if (is.null(l) == FALSE){ 
  if (test[l,9] == 0){ 
    while (test[l,9] == 0){
      test=test[-l,]
      l=l-1}}} 

graphe_VI_tab = test



##### Graphe different steps LGAs calling procedure intermediary #####

table_graphe_I_tab = as.data.frame(graphe_I_tab)

table_graphe_I_tab = table_graphe_I_tab[which(table_graphe_I_tab$chr != 23),]


lower_limit_graphe = -2
higher_limit_graphe = 2

if (max(abs(table_graphe_I_tab$ratio_median)) > 2){
  lower_limit_graphe = -max(abs(table_graphe_I_tab$ratio_median))
  higher_limit_graphe = max(abs(table_graphe_I_tab$ratio_median))
}


## preliminary_preparation

B <- data.frame(ratio_file_tsv)
B=B[,-1]
B=B[,-5]
colnames(B) <- c("chr", "start", "end", "readcount")
attach(B)

df=merge(aggregate(start ~ chr, B, min), aggregate(end ~ chr, B, max))[seq(from=1, to=22, by = 1),]
df=as.data.frame(df)
detach(B)



## add chr_arm information

adding_centromere = c(125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000, 49000000,
                      40200000, 53700000, 35800000, 36600000, 24000000, 17200000, 26500000,
                      27500000)

adding_centromere_chr = c(1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,20)


table_adding_centromere = data.frame(adding_centromere,adding_centromere_chr)
colnames(table_adding_centromere) = c("start_centromere", "chr")



## true part V

graphe_V_tab = data.frame(graphe_V_tab)
graphe_V_tab$size = graphe_V_tab$end - graphe_V_tab$start + 1

copy_graphe_V_tab = graphe_V_tab
copy_graphe_V_tab = copy_graphe_V_tab[,-8]
copy_graphe_V_tab = copy_graphe_V_tab[,-7]
copy_graphe_V_tab = copy_graphe_V_tab[,-1]

write.table(copy_graphe_V_tab, file = paste0(outputPath,"/",NAMEEE,"_final_segmentation.txt"), sep = "\t", row.names = FALSE)

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr!=23),]



C <- data.frame(graphe_V_tab)

C = C[which(C$chr!=23),] 



closest_higlight = Closest(B$start, 30302805)[1]

higlight_CCNE1 = B[B$chr == 19 & B$start == closest_higlight,]

data.segm = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.4, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)
data.text = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.5, chr = 19, label = "CCNE1")

if (higlight_CCNE1[1,5] >= 0.5){
  data.segm = data.frame(x=higlight_CCNE1[1,2], y=1.85, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)
  data.text = data.frame(x=higlight_CCNE1[1,2], y=1.95, chr = 19, label = "CCNE1")
}  else if(higlight_CCNE1[1,5] >= 1.7) {
  data.segm = data.frame(x=higlight_CCNE1[1,2], y=1.1, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] - 0.05, chr = 19)
  data.text = data.frame(x=higlight_CCNE1[1,2], y=1, chr = 19, label = "CCNE1")
}

V <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60", na.rm=TRUE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000", na.rm = TRUE) +
  geom_segment(data=data.segm, mapping=aes(x=x, y=y, xend=xend, yend=yend), 
               arrow=arrow(unit(0.30,"cm"), angle = 20), size=0.6, color="black", inherit.aes = FALSE, na.rm=TRUE) +
  ylim(lower_limit_graphe,higher_limit_graphe) +   scale_x_continuous(expand = c(0, 0)) + ggtitle(NAMEEE) +
  geom_vline(data=table_adding_centromere, mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(data = data.text, mapping = aes(x = x, y = y, label = label), size = 3.3, inherit.aes = TRUE, na.rm=TRUE) +
  geom_point(data = higlight_CCNE1, aes(x=start, y=ratio_median), color = "orange", size = 2, na.rm=TRUE) +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")


V <- V + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=15),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               legend.position = "none",
               panel.background = element_blank())

V = ggplotGrob(x = V)
V$layout$clip = "off"

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_final_segmentation_intermediary",".jpeg"), plot = V, device = "jpeg", width = 23, height = 13))


##  True Part VI : LGAs if called

if (sum(WC[,1]) != 0){
  
  B <- data.frame(ratio_file_tsv)
  B=B[,-1] 
  colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")
  
  B = B[which(B$chr!=23),]
  
  C <- data.frame(graphe_VI_tab)
  
  C = C[which(C$chr != 23),]
  
  
  graphe_VI_tab = data.frame(graphe_VI_tab)
  graphe_VI_tab$size = graphe_VI_tab$end - graphe_VI_tab$start + 1
  
  copy_graphe_VI_tab = graphe_VI_tab
  copy_graphe_VI_tab = copy_graphe_VI_tab[,-9]
  copy_graphe_VI_tab = copy_graphe_VI_tab[,-8]
  copy_graphe_VI_tab = copy_graphe_VI_tab[,-1]
  
  write.table(copy_graphe_VI_tab, file = paste0(outputPath,"/",NAMEEE,"_LGAs.txt"), sep = "\t", row.names = FALSE) 
  
  closest_higlight = Closest(B$start, 30302805)[1]
  
  higlight_CCNE1 = B[B$chr == 19 & B$start == closest_higlight,]
  
  data.segm = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.4, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)
  data.text = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.5, chr = 19, label = "CCNE1")
  
  if (higlight_CCNE1[1,5] >= 0.5){
    data.segm = data.frame(x=higlight_CCNE1[1,2], y=1.85, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)
    data.text = data.frame(x=higlight_CCNE1[1,2], y=1.95, chr = 19, label = "CCNE1")
  }  else if(higlight_CCNE1[1,5] >= 1.7) {
    data.segm = data.frame(x=higlight_CCNE1[1,2], y=1.1, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] - 0.05, chr = 19)
    data.text = data.frame(x=higlight_CCNE1[1,2], y=1, chr = 19, label = "CCNE1")
  }
  
  VI <- ggplot() +
    geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
    geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60", na.rm=TRUE) +
    geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#006600", na.rm=TRUE) +
    geom_segment(data=data.segm, mapping=aes(x=x, y=y, xend=xend, yend=yend), 
                 arrow=arrow(unit(0.30,"cm"), angle = 20), size=0.6, color="black", inherit.aes = FALSE, na.rm=TRUE) +
    ylim(lower_limit_graphe,higher_limit_graphe) +   scale_x_continuous(expand = c(0, 0)) + ggtitle(NAMEEE) +
    geom_vline(data=table_adding_centromere, mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
    geom_point(data = higlight_CCNE1, aes(x=start, y=ratio_median), color = "orange", size = 2, na.rm=TRUE) +
    geom_text(data = data.text, mapping = aes(x = x, y = y, label = label), size = 4, inherit.aes = TRUE, na.rm=TRUE) +
    geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
    scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
    facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")
  
  VI <- VI + theme(plot.title = element_text(hjust = 0.5),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_text(size=15),
                   panel.spacing = unit(0, "lines"),
                   strip.text.x = element_blank(),
                   line = element_blank(),
                   legend.position = "none",
                   panel.background = element_blank())
  
  VI = ggplotGrob(x = VI)
  VI$layout$clip = "off"
  
  suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_LGAs_intermediary",".jpeg"), plot = VI, device = "jpeg", width = 23, height = 13))
}

##### NEW Find Threshold 2 #####

cat("Find new Threshold 2... \n")


graphe_V_tab = data.frame(graphe_V_tab)

graphe_V_tab$size = graphe_V_tab$end - graphe_V_tab$start + 1

copy_graphe_V_tab = graphe_V_tab
copy_graphe_V_tab = copy_graphe_V_tab[,-8]
copy_graphe_V_tab = copy_graphe_V_tab[,-7]
copy_graphe_V_tab = copy_graphe_V_tab[,-1]

copy_graphe_V_tab$size <- copy_graphe_V_tab$end - copy_graphe_V_tab$start + 1

X = copy_graphe_V_tab

Y = X

Y = Y[which(Y$chr != 23),]


L=dim(Y)[1]
X=data.matrix(Y)

test=c()
i=1
for (i in 1:L){
  v=i+1
  for (y in v:L-1){
    test = c(test, abs(Y[i,5] - Y[y,5]))}
}  



### multiple test THR

list_THR = c()
list_max2 = c()
c = 0


while (c < nb_simulations){
  n = length(test)
  n
  
  rand_n = sample(2:n,1)
  rand_n
  
  sampled_test = sample(test,rand_n)
  
  ## First maxima
  
  maxx = localMaxima(density(sampled_test)$y)[1]
  First_local_maxima = density(sampled_test)$x[maxx]
  
  
  ## First minima 
  
  minx = localMinima(density(sampled_test)$y)
  First_local_minima = density(sampled_test)$x[minx][which(density(sampled_test)$x[minx] > First_local_maxima)][1]
  
  
  ## Second maxima
  
  maxx = localMaxima(density(sampled_test)$y)[2]
  Second_local_maxima = density(sampled_test)$x[maxx]
  
  if (is.na(Second_local_maxima) == FALSE){
    if (abs((First_local_minima-First_local_maxima)-(Second_local_maxima-First_local_minima)) < 0.05){
      list_THR = c(list_THR,First_local_minima)
      list_max2 = c(list_max2, Second_local_maxima)
      test_optimized = sampled_test
    }
  }
  c = c + 1
}


### Faire second inflexion point if no THR found


if (is.null(list_THR) == TRUE){
  
  ## First maxima
  
  maxx = localMaxima(density(test)$y)[1]
  First_maxima = density(test)$x[maxx]
  
  
  cat("First_minima... \n")
  
  ## First minima
  
  minx = localMinima(density(test)$y)
  First_local_minima = density(test)$x[minx][which(density(test)$x[minx] > First_maxima)][1]
  
  h <- hns(test, deriv.order = 2)
  den <- kdde(test, h=h, deriv.order = 2)
  
  flections<-c()
  for(i in 2:length(den$estimate)){
    if(sign(den$estimate[i])!=sign(den$estimate[i-1])){
      flections<-c(flections,i)
    }
  }
  
  All_inflexion_point = den$x[flections][which(den$x[flections] > First_maxima)]
  Threshold_change_sign_second_derivative_2nd_value = sort(All_inflexion_point)[1]
  
  Threshold = min(Threshold_change_sign_second_derivative_2nd_value, First_local_minima)
  
}


if (is.na(mean(list_THR)) == FALSE){
  number_positive = length(list_THR)
  
  Threshold = mean(list_THR)
  
  Threshold_second_round_simulation = Threshold
  
  MAX2 = mean(list_max2)
}



###  Classical THR + graphe

X = copy_graphe_V_tab

X = X[which(X[,6] > 2999999),]

X = X[which(X[,6] > ((quantile(X[,6])[[4]] - quantile(X[,6])[[2]])/2)),]  


X = X[which(X$chr != 23),]

L=dim(X)[1]
X=data.matrix(X)

test=c()
i=1
for (i in 1:L){
  v=i+1
  for (y in v:L-1){
    test = c(test, abs(X[i,5] - X[y,5]))}
}  

## First maxima

maxx = localMaxima(density(test)$y)[1]
First_local_maxima = density(test)$x[maxx]


## First minima 

minx = localMinima(density(test)$y)
First_local_minima = density(test)$x[minx][which(density(test)$x[minx] > First_local_maxima)][1]


## Second maxima

maxx = localMaxima(density(test)$y)[2]
Second_local_maxima = density(test)$x[maxx]


if (is.na(Second_local_maxima) == FALSE){
  if (abs((First_local_minima-First_local_maxima)-(Second_local_maxima-First_local_minima)) < 0.10){
    if (First_local_minima < Threshold){
      Threshold = First_local_minima
      Threshold_global_view = First_local_minima
      MAX2 = Second_local_maxima
    }
  }
}



test = as.data.frame(test)

ploty <- ggplot(test, aes(x = test)) + geom_density() +
  geom_vline(aes(xintercept = Threshold), color = "blue") + 
  ggtitle(paste("Density pairwise difference large segment")) + xlab("difference") +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=13),
        axis.text.y = element_text(size=15),
        axis.title.y = element_blank(),
        panel.background = element_blank())
ggsave(paste0(outputPath,"/",NAMEEE,"_THR",".jpeg"), plot = ploty, device = "jpeg", width = 20, height = 10)



##### NEW Fast gathering 2 #####

cat("Find fast gathering 2... \n")

dataTable <- read.table(paste0(inputPath,"/",NAMEEE,".bam_ratio.txt"), header = TRUE)
dataTable = dataTable[,1:4]

dataTable = dataTable[which(!dataTable[,1] == "X"),] 
dataTable = dataTable[which(!dataTable[,1] == "Y"),] 

dataTable[,1] = as.numeric(as.character(dataTable[,1]))

dataTable = dataTable[order(dataTable[,1]),]

L = dim(dataTable)[1]

size_window = dataTable[2,2] - dataTable[1,2]

dataTable = cbind(1:L, dataTable[,1], dataTable[,2], dataTable[,2]+size_window-1, 
                  dataTable[,3], dataTable[,4]) 

dataTable = as.data.frame(dataTable)
colnames(dataTable) = c("feature", "chromosome", "start", "end", "ratio", "ratio_median")

# log2 transformation
# dataTable[,5] = log2(dataTable[,5])
# dataTable[,6] = log2(dataTable[,6])

ratio <- data.frame(dataTable)

short_size_window = size_window/1000

X=ratio

ratio_file_tsv = ratio # save for later

X=X[,-1] # feature
X=X[,-4] 


## Remove spurious regions based on telomere and centromere from UCSC

X <- X[which(!(X[,1]==1 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==1 & X[,2]>=121535434 & X[,3]<=124535434)),]
X <- X[which(!(X[,1]==1 & X[,2]>=249240621 & X[,3]<=249250621)),]
X <- X[which(!(X[,1]==2 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==2 & X[,2]>=92326171 & X[,3]<=95326171)),]
X <- X[which(!(X[,1]==2 & X[,2]>=243189373 & X[,3]<=243199373)),]
X <- X[which(!(X[,1]==3 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==3 & X[,2]>=90504854 & X[,3]<=93504854)),]
X <- X[which(!(X[,1]==3 & X[,2]>=198012430 & X[,3]<=198022430)),]
X <- X[which(!(X[,1]==4 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==4 & X[,2]>=49660117 & X[,3]<=52660117)),]
X <- X[which(!(X[,1]==4 & X[,2]>=191144276 & X[,3]<=191154276)),]
X <- X[which(!(X[,1]==5 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==5 & X[,2]>=46405641 & X[,3]<=9405641)),]
X <- X[which(!(X[,1]==5 & X[,2]>=180905260 & X[,3]<=180915260)),]
X <- X[which(!(X[,1]==6 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==6 & X[,2]>=58830166 & X[,3]<=61830166)),]
X <- X[which(!(X[,1]==6 & X[,2]>=171105067 & X[,3]<=171115067)),]
X <- X[which(!(X[,1]==7 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==7 & X[,2]>=58054331 & X[,3]<=61054331)),]
X <- X[which(!(X[,1]==7 & X[,2]>=159128663 & X[,3]<=159138663)),]
X <- X[which(!(X[,1]==8 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==8 & X[,2]>=43838887 & X[,3]<=46838887)),]
X <- X[which(!(X[,1]==8 & X[,2]>=146354022 & X[,3]<=146364022)),]
X <- X[which(!(X[,1]==9 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==9 & X[,2]>=47367679 & X[,3]<=50367679)),]
X <- X[which(!(X[,1]==9 & X[,2]>=141203431 & X[,3]<=141213431)),]
X <- X[which(!(X[,1]==10 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==10 & X[,2]>=39254935 & X[,3]<=42254935)),]
X <- X[which(!(X[,1]==10 & X[,2]>=135524747 & X[,3]<=135534747)),]
X <- X[which(!(X[,1]==11 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==11 & X[,2]>=51644205 & X[,3]<=54644205)),]
X <- X[which(!(X[,1]==11 & X[,2]>=134996516 & X[,3]<=135006516)),]
X <- X[which(!(X[,1]==12 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==12 & X[,2]>=34856694 & X[,3]<=37856694)),]
X <- X[which(!(X[,1]==12 & X[,2]>=133841895 & X[,3]<=133851895)),]
X <- X[which(!(X[,1]==13 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=115159878 & X[,3]<=115169878)),]
X <- X[which(!(X[,1]==14 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=107339540 & X[,3]<=107349540)),]
X <- X[which(!(X[,1]==15 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=17000000 & X[,3]<=20000000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=102521392 & X[,3]<=102531392)),]
X <- X[which(!(X[,1]==16 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==16 & X[,2]>=35335801 & X[,3]<=38335801)),]
X <- X[which(!(X[,1]==16 & X[,2]>=90344753 & X[,3]<=90354753)),]
X <- X[which(!(X[,1]==17 & X[,2]>=22263006 & X[,3]<=25263006)),]
X <- X[which(!(X[,1]==18 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==18 & X[,2]>=15460898 & X[,3]<=18460898)),]
X <- X[which(!(X[,1]==18 & X[,2]>=78067248 & X[,3]<=78077248)),]
X <- X[which(!(X[,1]==19 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==19 & X[,2]>=24681782 & X[,3]<=27681782)),]
X <- X[which(!(X[,1]==19 & X[,2]>=59118983 & X[,3]<=59128983)),]
X <- X[which(!(X[,1]==20 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==20 & X[,2]>=26369569 & X[,3]<=29369569)),]
X <- X[which(!(X[,1]==20 & X[,2]>=63015520 & X[,3]<=63025520)),]
X <- X[which(!(X[,1]==21 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==21 & X[,2]>=11288129 & X[,3]<=14288129)),]
X <- X[which(!(X[,1]==21 & X[,2]>=48119895 & X[,3]<=48129895)),]
X <- X[which(!(X[,1]==22 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=13000000 & X[,3]<=16000000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=51294566 & X[,3]<=51304566)),]


## Add chromosome arm 

options(show.error.messages = FALSE)

X$chr_arm <- rep(0,nrow(X))
colnames(X) <- c("chr", "start", "end", "ratio_median", "chr_arm")
X=X[c("chr", "chr_arm", "start", "end", "ratio_median")]

X[,3]=as.numeric(as.character(X[,3])) 
X[,4]=as.numeric(as.character(X[,4])) 

tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] < 125000000
X[tt,2] = 1
tt <- X[,1] == 1 & X[,3] > 125000000 & X[,4] > 125000000
X[tt,2] = 2
tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] > 125000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 125000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 125000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])  
rownames(X) <- NULL 

tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] < 93300000
X[tt,2] = 1
tt <- X[,1] == 2 & X[,3] > 93300000 & X[,4] > 93300000
X[tt,2] = 2
tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] > 93300000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 93300000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 93300001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] < 91000000
X[tt,2] = 1
tt <- X[,1] == 3 & X[,3] > 91000000 & X[,4] > 91000000
X[tt,2] = 2
tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] > 91000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 91000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 91000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] < 50400000
X[tt,2] = 1
tt <- X[,1] == 4 & X[,3] > 50400000 & X[,4] > 50400000
X[tt,2] = 2
tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] > 50400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 50400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 50400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] < 48400000
X[tt,2] = 1
tt <- X[,1] == 5 & X[,3] > 48400000 & X[,4] > 48400000
X[tt,2] = 2
tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] > 48400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 48400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 48400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] < 61000000
X[tt,2] = 1
tt <- X[,1] == 6 & X[,3] > 61000000 & X[,4] > 61000000
X[tt,2] = 2
tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] > 61000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 61000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 61000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] < 59900000
X[tt,2] = 1
tt <- X[,1] == 7 & X[,3] > 59900000 & X[,4] > 59900000
X[tt,2] = 2
tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] > 59900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 59900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 59900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] < 45600000
X[tt,2] = 1
tt <- X[,1] == 8 & X[,3] > 45600000 & X[,4] > 45600000
X[tt,2] = 2
tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] > 45600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 45600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 45600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] < 49000000
X[tt,2] = 1
tt <- X[,1] == 9 & X[,4] > 49000000 & X[,4] > 49000000
X[tt,2] = 2
tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] > 49000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 49000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 49000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] < 40200000
X[tt,2] = 1
tt <- X[,1] == 10 & X[,3] > 40200000 & X[,4] > 40200000
X[tt,2] = 2
tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] > 40200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 40200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 40200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] < 53700000
X[tt,2] = 1
tt <- X[,1] == 11 & X[,3] > 53700000 & X[,4] > 53700000
X[tt,2] = 2
tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] > 53700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 53700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 53700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] < 35800000
X[tt,2] = 1
tt <- X[,1] == 12 & X[,3] > 35800000 & X[,4] > 35800000
X[tt,2] = 2
tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] > 35800000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 35800000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 35800001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] < 17900000
X[tt,2] = 1
tt <- X[,1] == 13 & X[,3] > 17900000 & X[,4] > 17900000
X[tt,2] = 2
tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] > 17900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] < 17600000
X[tt,2] = 1
tt <- X[,1] == 14 & X[,3] > 17600000 & X[,4] > 17600000
X[tt,2] = 2
tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] > 17600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] < 19000000
X[tt,2] = 1
tt <- X[,1] == 15 & X[,3] > 19000000 & X[,4] > 19000000
X[tt,2] = 2
tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] > 19000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 19000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 19000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] < 36600000
X[tt,2] = 1
tt <- X[,1] == 16 & X[,3] > 36600000 & X[,4] > 36600000
X[tt,2] = 2
tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] > 36600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 36600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 36600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] < 24000000
X[tt,2] = 1
tt <- X[,1] == 17 & X[,3] > 24000000 & X[,4] > 24000000
X[tt,2] = 2
tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] > 24000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 24000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 24000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] < 17200000
X[tt,2] = 1
tt <- X[,1] == 18 & X[,3] > 17200000 & X[,4] > 17200000
X[tt,2] = 2
tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] > 17200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] < 26500000
X[tt,2] = 1
tt <- X[,1] == 19 & X[,3] > 26500000 & X[,4] > 26500000
X[tt,2] = 2
tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] > 26500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 26500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 26500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] < 27500000
X[tt,2] = 1
tt <- X[,1] == 20 & X[,3] > 27500000 & X[,4] > 27500000
X[tt,2] = 2
tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] > 27500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 27500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 27500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] < 13200000
X[tt,2] = 1
tt <- X[,1] == 21 & X[,3] > 13200000 & X[,4] > 13200000
X[tt,2] = 2
tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] > 13200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 13200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 13200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] < 14700000
X[tt,2] = 1
tt <- X[,1] == 22 & X[,3] > 14700000 & X[,4] > 14700000
X[tt,2] = 2
tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] > 14700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 14700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 14700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 


X[,2]<-X[,1]+(X[,2]-1)/2      

tt=which(X$ratio_median==-Inf)

if (length(tt) >= 1){
  X=X[-tt,]  
}

X=as.matrix(X)

L = dim(X)[1]
A = matrix(0, ncol=6, nrow=L)

i=1    
c=1    

while (i < L){
  if(X[i,2]==X[i+1,2]){             
    if (X[i,5]==X[i+1,5]){          
      n=1
      while (X[i,2]==X[i+n,2] & X[i,5]==X[i+n,5] & i+n < L){
        n=n+1
      }
      if(i+n == L){
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n,4], X[i,5], X[i+n,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
      else{
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n-1,4], X[i,5], X[i+n-1,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
    }
    else{
      A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
      i=i+1
      c=c+1
    }
  }
  else{
    A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
    i=i+1                                        
    c=c+1                         
  }  
}

A=subset(A, A[,1] != 0)
rownames(A) <- NULL 
colnames(A) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")


write.table(A, file = paste0(outputPath,"/",NAMEEE,"_ratio_median_gathered.txt"), sep = "\t", row.names = FALSE)

options(show.error.messages = TRUE)


##### NEW Reading and initialization 2 #####

cat("Reading and initialisation 2... \n")

segFiles <- list.files(outputPath, pattern = paste0(NAMEEE, "_ratio_median_gathered.txt"),full.names=T)

nSample<-1

coefficient = 1.0

ColNamesTMP<-c("index", "chr","chr_arm","posStart", "posEnd", "ratio")

c_ind<-1; c_chr<-2;c_posS<-4; c_posE<-5; c_cn<-6; c_conf<-8

LLBB<-c(3,4,5,6,7,8,9,10,11)          # size LGA tested

TAB<-matrix(0,nSample,20)              # matrix to stock results 
rownames(TAB)<-seq(1:nSample)          
colnames(TAB)<-seq(1,20)               

colnames(TAB)[1:6]<-c("p_BAF","q_LRR","2copyLRR","DNA_ind","All_breaks","less3Mb_breaks")

ii<-1

rownames(TAB)[ii]<-"result"

results<-readSegmFile(segFileName=segFiles[ii])  
tmp<-results$tmp                                  

tmp <- cbind(seq(1,dim(tmp)[[1]]),tmp) 
colnames(tmp)[1]<-c("index")

tmp <- cbind(tmp,rep(0,nrow(tmp))) 
colnames(tmp)[8]<-c("level")

graphe_I_tab = tmp



# recalculate inital

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")
i=1
L=dim(tmp)[1]

while (i < L){
  gr=GRanges(seqnames=c(tmp[i,2]),
             ranges=IRanges(start=c(tmp[i,4]),end=c(tmp[i,5])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  tmp[i,6] = median(subsetGRobject$ratio)
  i=i+1
}


tt<-which(tmp[,c_chr+1]>23.6)                            # exclusion chr24 and more   
if(length(tt)>0){tmp<-tmp[-tt,]}

tt<-which(tmp[,c_chr+1]==21)                             # exclusion short arm chr21
if(length(tt)>0){tmp<-tmp[-tt,]}

tt<-which(tmp[,c_chr+1]==22)                             # exclusion short arm chr22
if(length(tt)>0){tmp<-tmp[-tt,]}

graphe_II_tab = tmp

tmp[,7] = tmp[,5] - tmp[,4] + 1

tmp_3mb = tmp[which(tmp[,7] > 2999999),]   
cop_tmp = tmp_3mb                    
THR=Threshold

if (THR > 0.45){
  THR = 0.45}

if (THR < 0.025){
  THR = 0.025
}

level=1


##### NEW Gather big segment (> 3Mb) #####

cat("Gather big segments 2... \n")

options(show.error.messages = FALSE)

while (dim(cop_tmp)[1] > 1){  
  
  line_largest_index = which.max(cop_tmp[,7])             # biggest segment
  largest_index = cop_tmp[line_largest_index,1]           # associated index
  ratio_largest = cop_tmp[line_largest_index,6]           # median_ratio segment
  ratio_largest = ratio_largest + 0.00001                 # avoid picking same value
  
  close_index = cop_tmp[which.min(abs(cop_tmp[,6] - ratio_largest)),1]
  closest_index = c(largest_index, close_index)
  
  ## Find all value that matches the two closest index fit with THR 
  
  L=dim(cop_tmp)[1]
  
  for (i in 1:L){
    validation = 0
    n = length(closest_index)
    for (g in 1:n){
      if (abs(cop_tmp[i,6] - cop_tmp[which((cop_tmp[,1]==closest_index[g])),6]) < THR){   
        validation = validation + 1}
    }
    if (validation == n){
      closest_index = c(closest_index, cop_tmp[i,1])}}
  
  closest_index = unique(closest_index)
  
  ##  Assign level based on the index of close segment in genome    
  
  n = length(closest_index)
  for (i in 1:n){
    tmp_3mb[which((tmp_3mb[,1] == closest_index[i])),8] = level
    cop_tmp <- cop_tmp[which(!(cop_tmp[,1]==closest_index[i])),]
  }
  dim(tmp_3mb)
  dim(cop_tmp)
  level = level + 1
} 


B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")

##  Gather with level

options(show.error.messages = TRUE)

L1 = dim(tmp_3mb)[1]

A = matrix(0, ncol=8, nrow=L1)

i=1
c=1   
while (i<L1+1){
  if (i==L1){
    A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])    # segment pas sur deux lignes, a stocker seul 
    i=i+1
    c=c+1}
  else{
    if(tmp_3mb[i,3] == tmp_3mb[i+1,3]){                                           
      if (tmp_3mb[i,8] == tmp_3mb[i+1,8]){                                        
        n=1
        vector_chr = c(tmp_3mb[i,2])
        vector_segment_start = c(tmp_3mb[i,4])
        vector_segment_end = c(tmp_3mb[i,5])
        vector_strand = c("*")
        
        while (tmp_3mb[i,8] == tmp_3mb[i+n,8] && tmp_3mb[i,3] == tmp_3mb[i+n,3]){ # while same chr_arm & level 
          vector_chr = c(vector_chr, tmp_3mb[i+n,2])
          vector_segment_start = c(vector_segment_start, tmp_3mb[i+n,4])
          vector_segment_end = c(vector_segment_end, tmp_3mb[i+n,5])
          vector_strand = c(vector_strand, "*")
          n = n + 1
          if(i+n == L1+1){
            break}}
        
        gr=GRanges(seqnames=vector_chr,
                   ranges=IRanges(start=vector_segment_start,end=vector_segment_end),
                   strand=vector_strand) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr) # overlapping
        
        A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i+n-1,5], 
                median(subsetGRobject$ratio), tmp_3mb[i+n-1,5]- tmp_3mb[i,4]+1, tmp_3mb[i,8]) # recalculating new_median
        i=i+n
        c=c+1}
      else{                                            
        A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])    # segment pas sur deux lignes, a stocker seul 
        i=i+1
        c=c+1}}
    else{                                              
      A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])
      i=i+1                                                              
      c=c+1}
  }  
}

A=subset(A, A[,1] != 0)

rownames(A) <- NULL 
colnames(A) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

tmp_3mb=A

graphe_III_tab = tmp_3mb 


##### NEW Reput small segment & Smoothing to final segmentation 2 #####

cat("Reput small segments 2... \n")

options(show.error.messages = FALSE)

tmp_0.1_3mb=tmp[which(tmp[,7] > 99999),]
tmp_0.1_3mb=tmp_0.1_3mb[which(tmp_0.1_3mb[,7] < 2999999),]  

L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

values = unique(tmp_0.1_3mb[,3][!tmp_0.1_3mb[,3] %in% tmp_3mb[,3]])

length=length(values)


for (missing_chr_arm in values){
  c=1
  L_3mb=dim(tmp_3mb)[1]
  while(c < L_3mb+1){
    if (missing_chr_arm > max(tmp_3mb[,3])){
      tmp_3mb=rbind(tmp_3mb, tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),])}
    else if (missing_chr_arm < tmp_3mb[c,3]){
      if (c==1){
        tmp_3mb=rbind(tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),] , tmp_3mb[1:L_3mb,])
        c=L_3mb+1}
      else{
        tmp_3mb=rbind(tmp_3mb[1:(c-1),], tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),] , tmp_3mb[c:L_3mb,])}
      c=L_3mb+1}
    else{
      c=c+1
    }
  }
}


for (missing_chr_arm in values){
  tmp_0.1_3mb <- tmp_0.1_3mb[which(!(tmp_0.1_3mb[,3]==missing_chr_arm)),]}


B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")

## initialisation

halt = 0 

i=1
c=1

L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

while (halt < 1){  
  if (tmp_0.1_3mb[i,3] > tmp_3mb[c,3] | (tmp_0.1_3mb[i,3] == tmp_3mb[c,3] && tmp_0.1_3mb[i,4] > tmp_3mb[c,5])){
    halt = halt + 1
  }
  
  else{                                            
    if (tmp_0.1_3mb[i,5] < tmp_3mb[c,4]){                       ########## segment avant
      if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR){          # No gathering 
        tmp_3mb=rbind(c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_0.1_3mb[i,5], 
                        tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]) , tmp_3mb[c:L_3mb,])
        i=i+1
      }
      
      else{                                                    # gathering
        
        gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                   ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                   strand=c("*")) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
        
        
        tmp_3mb=rbind(c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                        median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
        i=i+1}
    }
    else if (tmp_0.1_3mb[i,4] < tmp_3mb[c,4] && tmp_0.1_3mb[i,5] >= tmp_3mb[c,4]  && tmp_0.1_3mb[i,5] <= tmp_3mb[c,5]){  ########## segment overlapping mais debut avant
      if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR) {
        gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                   ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,4])),
                   strand=c("*")) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
        tmp_3mb=rbind(c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_3mb[c,4], 
                        median(subsetGRobject$ratio),tmp_3mb[c,4]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]), tmp_3mb[c:L_3mb,])
        i=i+1}
      else{
        gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                   ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                   strand=c("*")) # to overlap
        subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
        
        tmp_3mb=rbind(c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                        median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
        i=i+1}  
    }
    else{
      i=i+1
    }
  }
}

tmp_0.1_3mb = tmp_0.1_3mb[i:l_0.1_3mb,]  



L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

i=1
c=1

while (i < l_0.1_3mb + 1){          
  c=1                                                                            # reset 
  L_3mb=dim(tmp_3mb)[1]
  while (c < L_3mb+1){                                                            
    if (tmp_0.1_3mb[i,3] == tmp_3mb[c,3]){                                            
      if (tmp_0.1_3mb[i,5] < tmp_3mb[c,4]){                         ############### 1 - Before all segments chr_arm
        
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR){                     # not gathering with THR
          tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_0.1_3mb[i,5], 
                                             tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]) , tmp_3mb[c:L_3mb,])
          i=i+1
          c=L_3mb+1}
        else{                                                                # gathering
          gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                     ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                     strand=c("*")) # to overlap
          subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
          
          tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                                             median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
          i=i+1
          c=L_3mb+1}
      }
      
      else if (tmp_0.1_3mb[i,4] < tmp_3mb[c,4] && tmp_0.1_3mb[i,5] >= tmp_3mb[c,4]  && tmp_0.1_3mb[i,5] <= tmp_3mb[c,5] && tmp_0.1_3mb[i,4] >= tmp_3mb[c-1,5]){    ########## 2 - Inside segments beginning before but after previous segment
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR) {
          if (c+1 > L_3mb) {
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,4])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_3mb[c,4], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,4]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]), tmp_3mb[c:L_3mb,])
            i=i+1
            c=L_3mb+1}
          else{
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,4])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_3mb[c,4], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,4]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]))
            i=i+1
            c=L_3mb+1}
        }
        
        else{
          if (c+1 > L_3mb) {
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]))
            i=i+1
            c=L_3mb+1
          }
          else{
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                                               median(subsetGRobject$ratio),tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}
        }
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,4] && tmp_0.1_3mb[i,5] <= tmp_3mb[c,5]){    ########## 3 - Inside segments pile
        i=i+1
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,4] && tmp_0.1_3mb[i,4] <= tmp_3mb[c,5] && tmp_0.1_3mb[i,5] > tmp_3mb[c,5] && tmp_0.1_3mb[i,5] <= tmp_3mb[c+1,5]){    ########## 4 - Inside segments end after but before next segment
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR) {
          if (c+1 > L_3mb) {
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_0.1_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c),], c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_3mb[c,5], tmp_0.1_3mb[c,5], 
                                             median(subsetGRobject$ratio),tmp_0.1_3mb[c,5]-tmp_3mb[c,5]+1,tmp_0.1_3mb[i,8]))
            i=i+1
            c=L_3mb+1 
          }
          else{
            gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_0.1_3mb[c,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c),], c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_3mb[c,5], tmp_0.1_3mb[c,5], 
                                             median(subsetGRobject$ratio),tmp_0.1_3mb[c,5]-tmp_3mb[c,5]+1,tmp_0.1_3mb[i,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1
          }
        }
        
        else{
          if (c+1 > L_3mb){
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], 
                                               median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]))
            i=i+1
            c=L_3mb+1
          }
          else{
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], 
                                               median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1
          }
        }
      }
      
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,5] && tmp_0.1_3mb[i,5] <= tmp_3mb[c+1,4] && tmp_3mb[c,3] == tmp_3mb[c+1,3]){     ########## 5 - between two segments same chr_arm
        
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR && abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6]) <= THR){       # two in range THR
          if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6])){     # c closer
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], median(subsetGRobject$ratio),
                                                tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]) , tmp_3mb[(c+1):L_3mb,])  
            i=i+1
            c=L_3mb+1}
          
          else {                      # c+1 closer
            if (c+2 > L_3mb){         # anticipate end of file
              gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8])) 
              i=i+1
              c=L_3mb+1}
            else{
              gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]) , tmp_3mb[(c+2):L_3mb,])  
              i=i+1
              c=L_3mb+1}
          }
        }
        
        else if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR){                   # only c
          gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                     ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                     strand=c("*")) # to overlap
          subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
          
          tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], median(subsetGRobject$ratio),
                                              tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]) , tmp_3mb[(c+1):L_3mb,])  
          i=i+1
          c=L_3mb+1}
        
        else if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6]) <= THR){                  # only c+1
          if (c+2 > L_3mb){                                                       # anticipate end of file
            gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
                       ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                            tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]))  
            i=i+1
            c=L_3mb+1}
          
          else{
            if (tmp_0.1_3mb[i+1,5] < tmp_3mb[c+1,4] && tmp_0.1_3mb[i+1,3] == tmp_3mb[c+1,3]){  # small segments afterward
              gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_0.1_3mb[i,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_0.1_3mb[i,5], median(subsetGRobject$ratio),
                                              tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1,tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])  
              i=i+1
              c=L_3mb+1}
            else{
              gr=GRanges(seqnames=c(tmp_0.1_3mb[i,2]),
                         ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
                         strand=c("*")) # to overlap
              subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
              
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]) , tmp_3mb[(c+2):L_3mb,])  
              i=i+1
              c=L_3mb+1}
          }
        }
        
        else{ # aucun                                      
          tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                          tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])
          i=i+1
          c=L_3mb+1}
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,5] && (tmp_3mb[c,3] != tmp_3mb[c+1,3] | is.null(tmp_3mb[c+1,3]))){     ############## 6 - after last segment chr_arm
        
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR){   # in range THR
          if (c+1 > L_3mb){                                 # anticipate end of file
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3],tmp_3mb[c,4],
                                                tmp_0.1_3mb[i,5],median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8])) 
            i=i+1
            c=L_3mb+1}
          else{                                             # not end of file
            gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                       ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_0.1_3mb[i,5])),
                       strand=c("*")) # to overlap
            subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
            
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3],tmp_3mb[c,4],
                                                tmp_0.1_3mb[i,5],median(subsetGRobject$ratio),tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}}
        
        else{                                       # not in range THR
          if (c+1 > L_3mb){                         # anticipate end of file
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                            tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]))
            i=i+1
            c=L_3mb+1}
          else{                                     # not end of file
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                            tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}}}
      
      else {             # none of the cases listed up there
        c=c+1
      }
    }    
    else{                                            
      c=c+1
    }
  }
}

rownames(tmp_3mb) <- NULL



## add leftovers

L_3mb=dim(tmp_3mb)[1]
c = 1

while (c < L_3mb){
  if (tmp_3mb[c,3] == tmp_3mb[c+1,3]){
    if (tmp_3mb[c+1,4] - tmp_3mb[c,5] > 1){
      
      gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                 ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_3mb[c+1,4])),
                 strand=c("*")) # to overlap
      subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
      
      if (abs(tmp_3mb[c,6] - median(subsetGRobject$ratio)) < abs(tmp_3mb[c+1,6] - median(subsetGRobject$ratio))){
        tmp_3mb[c,5] = tmp_3mb[c+1,4] - 1
      }
      else{
        tmp_3mb[c+1,4] = tmp_3mb[c,5] + 1
      }
    }
  }
  c=c+1
}

rownames(tmp_3mb) <- NULL

options(show.error.messages = TRUE)

colnames(tmp_3mb) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

tmp_3mb[,7] = tmp_3mb[,5] - tmp_3mb[,4] + 1

graphe_IV_tab = tmp_3mb

write.table(graphe_IV_tab, file = paste0(outputPath,"/",NAMEEE,"_IV.txt"), sep = "\t", row.names = FALSE)




## V

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")

tmp_3mb<-breakSmoothToLGA(THR,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)



## missing chr_arm

L_3mb=dim(tmp_3mb)[1]
L_graphe_IV_tab=dim(graphe_IV_tab)[1]

values = unique(graphe_IV_tab[,3][!graphe_IV_tab[,3] %in% tmp_3mb[,3]])



for (missing_chr_arm in values){
  c=1
  L_3mb=dim(tmp_3mb)[1]
  while(c < L_3mb+1){
    if (missing_chr_arm > max(tmp_3mb[,3])){
      tmp_3mb=rbind(tmp_3mb, graphe_IV_tab[which(graphe_IV_tab[,3] == missing_chr_arm),])}
    else if (missing_chr_arm < tmp_3mb[c,3]){
      if (c==1){
        tmp_3mb=rbind(graphe_IV_tab[which(graphe_IV_tab[,3] == missing_chr_arm),] , tmp_3mb[1:L_3mb,])
        c=L_3mb+1}
      else{
        tmp_3mb=rbind(tmp_3mb[1:(c-1),], graphe_IV_tab[which(graphe_IV_tab[,3] == missing_chr_arm),] , tmp_3mb[c:L_3mb,])}
      c=L_3mb+1}
    else{
      c=c+1
    }
  }
}



## add leftovers

rownames(tmp_3mb) <- NULL

L_3mb=dim(tmp_3mb)[1]
c = 1

while (c < L_3mb){
  if (tmp_3mb[c,3] == tmp_3mb[c+1,3]){
    if (tmp_3mb[c+1,4] - tmp_3mb[c,5] > 1){
      
      gr=GRanges(seqnames=c(tmp_3mb[c,2]),
                 ranges=IRanges(start=c(tmp_3mb[c,5]),end=c(tmp_3mb[c+1,4])),
                 strand=c("*")) # to overlap
      subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
      
      if (abs(tmp_3mb[c,6] - median(subsetGRobject$ratio)) < abs(tmp_3mb[c+1,6] - median(subsetGRobject$ratio))){
        tmp_3mb[c,5] = tmp_3mb[c+1,4] - 1
      }
      else{
        tmp_3mb[c+1,4] = tmp_3mb[c,5] + 1
      }
    }
  }
  c=c+1
}

rownames(tmp_3mb) <- NULL

colnames(tmp_3mb) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

tmp_3mb[,7] = tmp_3mb[,5] - tmp_3mb[,4] + 1




## go to begin and end of chromosome

B <- read.table(paste0(outputPath,"/",NAMEEE, "_ratio_median_gathered.txt"), header = TRUE)

B <- data.frame(B)

colnames(B) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")
attach(B)

df=merge(aggregate(start ~ chr_arm, B, min), aggregate(end ~ chr_arm, B, max))[seq(from=1, to=39, by = 1),]
df=as.data.frame(df)
detach(B)



# first
tmp_3mb[1,4] = df[1,2]

# inside
L_3mb=dim(tmp_3mb)[1]
c = 1
n = 1

while (c < L_3mb){
  if (tmp_3mb[c,3] != tmp_3mb[c+1,3]){ # pas même chr_arm
    tmp_3mb[c,5] = df[n,3]
    n=n+1
    tmp_3mb[c+1,4] = df[n,2]
  }
  c=c+1
}

# final
tmp_3mb[L_3mb,5] = df[39,3]

tmp_3mb[,7] = tmp_3mb[,5] - tmp_3mb[,4] + 1




# recalculate complete value

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")

c=1
L_3mb=dim(tmp_3mb)[1]

while (c < L_3mb){
  gr=GRanges(seqnames=c(tmp_3mb[c,2]),
             ranges=IRanges(start=c(tmp_3mb[c,4]),end=c(tmp_3mb[c,5])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  tmp_3mb[c,6] = median(subsetGRobject$ratio)
  c=c+1
}

graphe_V_tab = tmp_3mb


## graphe representative of the final segmentation 

test_data_frame = as.data.frame(tmp_3mb)

test_data_frame = test_data_frame[which(test_data_frame$chr != 23),]

colnames(test_data_frame) <- make.unique(names(test_data_frame))

test_ordered = test_data_frame[order(test_data_frame$ratio_median),]
test_ordered$num_line <- seq.int(nrow(test_ordered))

test_ploty_CN_level <- ggplot(test_ordered, aes(x = ratio_median, y = num_line)) + 
  geom_point(shape = 1, size = 1, color = "#0072B2", na.rm = TRUE) + geom_line() + ggtitle("Final segmentation diagnostic") + 
  xlab("ratio median") + xlim(-2.5,2.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())
suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_final_segmentation_visual",".jpeg"), plot = test_ploty_CN_level, device = "jpeg", width = 20, height = 10))



##### NEW Call LGAs 2 #####

cat("Call Large Genomic Alterations 2... \n")

LGAs_data_frame <- as.data.frame(matrix(0, ncol = 2, nrow = 9))
LGAs_data_frame[,1] = c(3:11)
colnames(LGAs_data_frame) <- c("Size_LGA", "Number_LGA")

tmp_3mb = data.frame(tmp_3mb)

tmp_3mb = tmp_3mb[which(tmp_3mb$chr != 23),]

tmp_3mb = as.matrix(tmp_3mb)


for (i in (3:11)){
  WC<-LGA_control(THR,lenBIN=500,lenMB=i,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  LGAs_data_frame[i-2,2] = sum(WC[,1])}
write.table(LGAs_data_frame, file = paste0(outputPath,"/",NAMEEE,"_number_LGAs.txt"), sep = "\t", row.names = FALSE)


# For graphe 10Mb 

WC<-LGA_control(THR,lenBIN=500,lenMB=10,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)

test=cbind(tmp_3mb,WC[,1])
colnames(test) <- c("index", "chr","chr_arm", "start", "end", "ratio_median", "size", "level", "WC")


L=dim(test)[1]

i=1
while(i + 1 < L){
  L=dim(test)[1]
  if (test[i,9] == 0 & test[i+1,9] != 1){
    test=test[-i,]}
  else{
    i=i+1}}

l=dim(test)[1]

if (is.null(l) == FALSE){ 
  if (test[l,9] == 0){ 
    while (test[l,9] == 0){
      test=test[-l,]
      l=l-1}}} 

graphe_VI_tab = test



##### Graphe different steps LGAs calling procedure #####
## limit ##

cat("Creating graphes... \n")

table_graphe_I_tab = as.data.frame(graphe_I_tab)

table_graphe_I_tab = table_graphe_I_tab[which(table_graphe_I_tab$chr != 23),]

lower_limit_graphe = -2
higher_limit_graphe = 2

if (max(abs(table_graphe_I_tab$ratio_median)) > 2){
  lower_limit_graphe = -max(abs(table_graphe_I_tab$ratio_median))
  higher_limit_graphe = max(abs(table_graphe_I_tab$ratio_median))
}


## preliminary_preparation


B <- data.frame(ratio_file_tsv)
B=B[,-1]
B=B[,-5]
colnames(B) <- c("chr", "start", "end", "readcount")
attach(B)

B = B[which(B$chr != 23),]

df=merge(aggregate(start ~ chr, B, min), aggregate(end ~ chr, B, max))[seq(from=1, to=22, by = 1),]
df=as.data.frame(df)
detach(B)


## add chr_arm information

adding_centromere = c(125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000, 49000000,
                      40200000, 53700000, 35800000, 36600000, 24000000, 17200000, 26500000,
                      27500000)

adding_centromere_chr = c(1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,20)


table_adding_centromere = data.frame(adding_centromere,adding_centromere_chr)
colnames(table_adding_centromere) = c("start_centromere", "chr")


## Normalized read count

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]


Z <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "black") + ggtitle(NAMEEE) + 
  ylim(lower_limit_graphe,higher_limit_graphe) + scale_x_continuous(expand = c(0, 0)) + geom_vline(data=table_adding_centromere, 
                                                                                                   mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x") 

  
Z <- Z + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=20),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               legend.position = "none",
               panel.background = element_blank())

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_normalised_read_count",".jpeg"), plot = Z, device = "jpeg", width = 23, height = 13))


## Part I

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]

C <- data.frame(graphe_I_tab)

C = C[which(C$chr != 23),]

I <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") + ggtitle(NAMEEE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(lower_limit_graphe,higher_limit_graphe) + scale_x_continuous(expand = c(0, 0)) + geom_vline(data=table_adding_centromere, 
                                                                  mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x") 

I <- I + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=20),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               legend.position = "none",
               panel.background = element_blank())

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_beginning_segmentation",".jpeg"), plot = I, device = "jpeg", width = 23, height = 13))



##### Part II to IV graphes(debug) #####
## Part II ##


B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]

C <- data.frame(graphe_II_tab)

C = C[which(C$chr != 23),]

II <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") + ggtitle(NAMEEE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(lower_limit_graphe,higher_limit_graphe) + scale_x_continuous(expand = c(0, 0)) + geom_vline(data=table_adding_centromere, 
                                                                  mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x") 

II <- II + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=20),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               legend.position = "none",
               panel.background = element_blank())

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_II",".jpeg"), plot = II, device = "jpeg", width = 23, height = 13))


## Part III ##

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]

C <- data.frame(graphe_III_tab)

C = C[which(C$chr != 23),]

III <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") + ggtitle(NAMEEE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(lower_limit_graphe,higher_limit_graphe) + scale_x_continuous(expand = c(0, 0)) + geom_vline(data=table_adding_centromere, 
                                                                  mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x") 

III <- III + theme(plot.title = element_text(hjust = 0.5),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_text(size=20),
                 panel.spacing = unit(0, "lines"),
                 strip.text.x = element_blank(),
                 line = element_blank(),
                 legend.position = "none",
                 panel.background = element_blank())

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_III",".jpeg"), plot = III, device = "jpeg", width = 23, height = 13))


## Part IV ##

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]

C <- data.frame(graphe_IV_tab)

C = C[which(C$chr != 23),]

IV <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") + ggtitle(NAMEEE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(lower_limit_graphe,higher_limit_graphe) + scale_x_continuous(expand = c(0, 0)) + geom_vline(data=table_adding_centromere, 
                                                                  mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x") 

IV <- IV + theme(plot.title = element_text(hjust = 0.5),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_text(size=20),
                   panel.spacing = unit(0, "lines"),
                   strip.text.x = element_blank(),
                   line = element_blank(),
                   legend.position = "none",
                   panel.background = element_blank())

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_IV",".jpeg"), plot = IV, device = "jpeg", width = 23, height = 13))



##### Graphe different steps LGAs calling procedure (second part) #####
## true part V

graphe_V_tab = data.frame(graphe_V_tab)
graphe_V_tab$size = graphe_V_tab$end - graphe_V_tab$start + 1

copy_graphe_V_tab = graphe_V_tab
copy_graphe_V_tab = copy_graphe_V_tab[,-8]
copy_graphe_V_tab = copy_graphe_V_tab[,-7]
copy_graphe_V_tab = copy_graphe_V_tab[,-1]

write.table(copy_graphe_V_tab, file = paste0(outputPath,"/",NAMEEE,"_final_segmentation.txt"), sep = "\t", row.names = FALSE)

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]

C <- data.frame(graphe_V_tab)

C = C[which(C$chr != 23),]

V <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60", na.rm=TRUE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000", na.rm = TRUE) +
  ylim(lower_limit_graphe,higher_limit_graphe) +   scale_x_continuous(expand = c(0, 0)) + ggtitle(NAMEEE) +
  geom_vline(data=table_adding_centromere, mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")


V <- V + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=15),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               legend.position = "none",
               panel.background = element_blank())

V = ggplotGrob(x = V)
V$layout$clip = "off"

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_final_segmentation",".jpeg"), plot = V, device = "jpeg", width = 23, height = 13))



##  True Part VI : LGAs if called

if (sum(WC[,1]) != 0){    
  B <- data.frame(ratio_file_tsv)
  B=B[,-1] 
  colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")
  
  B = B[which(B$chr != 23),]
  
  C <- data.frame(graphe_VI_tab)
  
  C = C[which(C$chr != 23),]
  
  graphe_VI_tab = data.frame(graphe_VI_tab)
  graphe_VI_tab$size = graphe_VI_tab$end - graphe_VI_tab$start + 1
  
  copy_graphe_VI_tab = graphe_VI_tab
  copy_graphe_VI_tab = copy_graphe_VI_tab[,-9]
  copy_graphe_VI_tab = copy_graphe_VI_tab[,-8]
  copy_graphe_VI_tab = copy_graphe_VI_tab[,-1]
  
  write.table(copy_graphe_VI_tab, file = paste0(outputPath,"/",NAMEEE,"_LGAs.txt"), sep = "\t", row.names = FALSE) 
  
  VI <- ggplot() +
    geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
    geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60", na.rm=TRUE) +
    geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#006600", na.rm=TRUE) +
    ylim(lower_limit_graphe,higher_limit_graphe) +   scale_x_continuous(expand = c(0, 0)) + ggtitle(NAMEEE) +
    geom_vline(data=table_adding_centromere, mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
    geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
    scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
    facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")
  
  VI <- VI + theme(plot.title = element_text(hjust = 0.5),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_text(size=15),
                   panel.spacing = unit(0, "lines"),
                   strip.text.x = element_blank(),
                   line = element_blank(),
                   legend.position = "none",
                   panel.background = element_blank())
  
  VI = ggplotGrob(x = VI)
  VI$layout$clip = "off"
  
  suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_LGAs",".jpeg"), plot = VI, device = "jpeg", width = 23, height = 13))
}





### ZOOMED

B <- data.frame(ratio_file_tsv)
B=B[,-1]
B=B[,-5]
colnames(B) <- c("chr", "start", "end", "readcount")
attach(B)

df=merge(aggregate(start ~ chr, B, min), aggregate(end ~ chr, B, max))[seq(from=1, to=22, by = 1),]
df=as.data.frame(df)
detach(B)



## add chr_arm information

adding_centromere = c(125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000, 49000000,
                      40200000, 53700000, 35800000, 36600000, 24000000, 17200000, 26500000,
                      27500000)

adding_centromere_chr = c(1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,20)

table_adding_centromere = data.frame(adding_centromere,adding_centromere_chr)
colnames(table_adding_centromere) = c("start_centromere", "chr")



## true part V ZOOMED

graphe_V_tab = data.frame(graphe_V_tab)
graphe_V_tab$size = graphe_V_tab$end - graphe_V_tab$start + 1

copy_graphe_V_tab = graphe_V_tab
copy_graphe_V_tab = copy_graphe_V_tab[,-8]
copy_graphe_V_tab = copy_graphe_V_tab[,-7]
copy_graphe_V_tab = copy_graphe_V_tab[,-1]

B <- data.frame(ratio_file_tsv)
B=B[,-1] 
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

C <- data.frame(graphe_V_tab)

V <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60", na.rm=TRUE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000", na.rm = TRUE) +
  ylim(-2,2) +   scale_x_continuous(expand = c(0, 0)) + ggtitle(NAMEEE) +
  geom_vline(data=table_adding_centromere, mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_text(aes(x=start,y=2, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) + 
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")


V <- V + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=15),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               legend.position = "none",
               panel.background = element_blank())

V = ggplotGrob(x = V)
V$layout$clip = "off"

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_final_segmentation_zoomed",".jpeg"), plot = V, device = "jpeg", width = 23, height = 13))



##### Quality assessment #####

X_Ratio <- ratio_file_tsv
X_I = read.table(paste0(outputPath,"/",NAMEEE, "_ratio_median_gathered.txt"), header = TRUE)

X_I = X_I[which(X_I$chr != 23),]


X_Ratio["ratio_new"] <- NA

l = dim(X_Ratio)[1]
j=1 # X_Ratio

L = dim(X_I)[1]
i=1 # X_I

colnames(X_Ratio) <- c("feature", "chr", "start", "end", "ratio", "ratio_median", "ratio_new")

X_Ratio = X_Ratio[which(X_Ratio$chr != 23),]

options(show.error.messages = FALSE)

cat("Quality assessment... \n")

while (j < l + 1){
  while (i < L + 1){
    
    if (X_Ratio$chr[j] == X_I$chr[i]){
      if (X_Ratio$start[j] >= X_I$start[i] && X_Ratio$end[j] <= X_I$end[i]){
        X_Ratio$ratio_new[j] = X_Ratio$ratio[j] - X_I$ratio_median[i]
        j = j + 1
      }
      else if (X_Ratio$start[j] >= X_I$end[i]) {
        i = i + 1
      }
      else {
        j = j + 1
      }
    }
    
    else if (X_Ratio$chr[j] > X_I$chr[i]) {
      i = i + 1
    }
    
    else{
      j = j + 1
    }
  }
}

X_Ratio$ratio_new = abs(X_Ratio$ratio_new)

vector_ratio_new = c(na.omit(X_Ratio$ratio_new))

QC_MAD_point_corrected = mad(vector_ratio_new, constant = 1)


# quality

quality = "Good"

if (QC_MAD_point_corrected > 0.50){
  quality = "Bad"
} else if (QC_MAD_point_corrected > 0.14){
  if (Threshold > 0.45){
    quality = "Bad"
  }
  else{
    quality = "Average" 
  }
} else{
  if (Threshold > 0.45){
    quality = "Average"
  }
}

if (Threshold < 0.025){
  quality = "Normal or low cellularity"
}

##### amplifications & deletions #####
### loading

options(show.error.messages = TRUE)

cat("Amplifications & deletions... \n")

## Initial file

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]

## Final segmentation

C <- data.frame(graphe_V_tab)

C=C[,-8]
C=C[,-7]
C=C[,-1]

colnames(C) <- c("chr", "chr_arm", "start", "end", "ratio_median")

C = C[which(C$chr != 23),]


gene_list <- load_gene_list(custom_gene_file)
lookup_gene_table <- lookup_gene_ratios(gene_list, B, C, THR)

amplification_deletion_table = data.frame(
  gene = lookup_gene_table$gene,
  chr = lookup_gene_table$chr,
  start = lookup_gene_table$start,
  ratio_point_initial = lookup_gene_table$ratio_point_initial,
  CN_to_baseline_point_initial = lookup_gene_table$CN_point_initial,
  ratio_segment_initial = lookup_gene_table$ratio_segment_initial,
  CN_to_baseline_segment_initial = lookup_gene_table$CN_segment_initial,
  ratio_segment_final = lookup_gene_table$ratio_segment_final,
  CN_to_baseline_segment_final = lookup_gene_table$CN_segment_final
)

write.table(amplification_deletion_table, file = paste0(outputPath,"/",NAMEEE,"_amplification_deletion_table.txt"), sep = "\t", row.names = FALSE)





### Figure summary deletions and amplifications

# Amp en haut
# Deletions en bas

# limit

table_graphe_I_tab = as.data.frame(graphe_I_tab)

table_graphe_I_tab = table_graphe_I_tab[which(table_graphe_I_tab$chr != 23),]

lower_limit_graphe = -2
higher_limit_graphe = 2

if (max(abs(table_graphe_I_tab$ratio_median)) > 2){
  lower_limit_graphe = -max(abs(table_graphe_I_tab$ratio_median))
  higher_limit_graphe = max(abs(table_graphe_I_tab$ratio_median))
}



# preliminary_preparation

B <- data.frame(ratio_file_tsv)
B=B[,-1]
B=B[,-5]
colnames(B) <- c("chr", "start", "end", "readcount")

B = B[which(B$chr != 23),]

df=merge(aggregate(start ~ chr, B, min), aggregate(end ~ chr, B, max))[seq(from=1, to=22, by = 1),]
df=as.data.frame(df)


# add chr_arm information

adding_centromere = c(125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000, 49000000,
                      40200000, 53700000, 35800000, 36600000, 24000000, 17200000, 26500000,
                      27500000)

adding_centromere_chr = c(1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,20)

table_adding_centromere = data.frame(adding_centromere,adding_centromere_chr)
colnames(table_adding_centromere) = c("start_centromere", "chr")



# Normalised read count

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

B = B[which(B$chr != 23),]

# Plot-only highlighting for key HRD genes (not written to output table)
hrd_genes <- c("BRCA1", "BRCA2", "RAD51C", "CDKN2AIP", "PIK3CA", "MYC", "CD274_PDL1", "PTEN", "CCND1", "CCNE1", "NF1", "ERBB2")
amplification_deletion_table$label_color <- ifelse(amplification_deletion_table$gene %in% hrd_genes, "red", "black")

Z <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = chr %% 2 == 0)) +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") + ggtitle(NAMEEE) +
  ylim(lower_limit_graphe,higher_limit_graphe) + scale_x_continuous(expand = c(0, 0)) + geom_vline(data=table_adding_centromere,
                                                                                                   mapping=aes(xintercept=start_centromere), color="black", linetype="dotted") +
  geom_point(data=amplification_deletion_table, aes(x=start, y=ratio_segment_initial), color = "orange", size = 2) +
  coord_cartesian(clip = "off") +
  geom_text_repel(data = amplification_deletion_table,
                  mapping = aes(x = start, y = ratio_segment_initial, label = gene, size = 7, fontface = 'bold', color = label_color),
                  min.segment.length = 0, box.padding = 2, xlim = c(-Inf, NA),
                  ylim = c(-Inf, Inf), max.overlaps = Inf, max.iter = Inf, max.time = 60, force = 30) +
  scale_color_identity() +
  geom_text(aes(x=start,y=higher_limit_graphe, hjust = "left", vjust = "top", label = chr), data = df, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")


Z <- Z + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=20),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               legend.position = "none",
               panel.background = element_blank())


suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_amplifications_deletions",".jpeg"), plot = Z, device = "jpeg", width = 23, height = 13))
options(show.error.messages = FALSE)


##### Summary figure #####
## graphe

cat("Creating summary_plot... \n")

graphe = V
if (sum(WC[,1]) != 0){
  graphe = VI
}

## Threshold

Threshold = round(Threshold, 3)

## THR (corrected)

THR = round(THR, 3)

## number LGAs

number_LGAs = sum(WC[,1])

## BRCAness

if (number_LGAs >= 20){
  HRD = "Yes (>= 20)"
} else if (number_LGAs >= 15){
  HRD = "Borderline [15;19]"
} else{
  HRD = "No (< 15)"
}

## MAX2_diag

MAX2_diag = "Correct"

if (MAX2 <0.16){
  MAX2_diag = "Low cellularity"
}

## QC_MAD_point

QC_MAD_point_corrected = round(QC_MAD_point_corrected,3)


## quality simulations

quality = "good"

if (number_positive/nb_simulations < 0.1){
  quality = "Low quality"
}


## cMAD inclusion 

if (cMAD > 0.05){
  quality = "Discard sample"
}


## Table

a = c("CNA cut-off value", "Number simulations for cut-off", "Number successfull simulations", "MAX2", "cMAD",  "Number LGAs 10Mb", "HRD")
b = c(THR, format(nb_simulations,scientific = F), number_positive, round(MAX2,3), QC_MAD_point_corrected,   number_LGAs, HRD)


A = data.frame(a,b)

colnames(A) <- NULL
row.names(A) <- NULL

testy <- qplot(1:10, 1:10, geom = "blank") + 
  theme(line = element_blank(), text = element_blank(), 
        panel.background = element_blank()) +
  annotation_custom(grob = tableGrob(A, rows = NULL, theme = ttheme_default(base_size = 14, base_colour = "black", base_family = "",
                                                                            parse = FALSE)))


## Summary figure

summary_plot <- suppressWarnings(ggarrange(graphe,
                                           ggarrange(ploty, test_ploty_CN_level, testy, labels = c("B", "C", "D"), ncol = 3, font.label = list(size = 14, face = "bold")),
                                           nrow = 2, labels = "A", font.label = list(size = 14, face = "bold")))

suppressWarnings(ggsave(paste0(outputPath,"/",NAMEEE,"_summary_plot",".jpeg"), plot = summary_plot, device = "jpeg", width = 18, height = 10.17, dpi = 300))
