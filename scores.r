rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)




#-------------------------------------read matrices ------------------------------------------


pfm_ARF<- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF <- round((t(as.matrix(pfm_ARF)))*nRegion)+1 ;pfm_ARF
pfm_ARF <- cbind(rep(151,4),rep(151,4),pfm_ARF)
maxi_ARF <- apply(pfm_ARF,FUN=max, 2)
maxi_ARF <- matrix(nrow=4, rep(maxi_ARF,4),byrow=TRUE)
pwm_ARF <- log(pfm_ARF/maxi_ARF)
pwm_ARF_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 
pwm_ARF <-  reverseComplement(pwm_ARF_rev) ; pwm_ARF

#-------------------------------------read fasta files-----------------------------------------

ARF_pos <- readDNAStringSet('ARF5.fas')
width_pos <- width(ARF_pos)
seq_pos <- as.character(ARF_pos)
seq_rev_pos <- as.character(reverseComplement(ARF_pos))

#-------------------------------------read bed files-----------------------------------------

bed <- read.csv('ARF5.bed',sep="\t", header=FALSE)
DAP_score <- rep(bed[,5],2)


#-------------------------------------Compute Scores-----------------------------------------
th <- -6
#
scores_ARF_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
scores_ARF_rev_pos<- mapply(seq_rev_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
pos_max <- apply(FUN=which.max,scores_ARF_pos,2)
pos_max_rev<- apply(FUN=which.max,scores_ARF_rev_pos,2)
#
middle <- median(c(pos_max,pos_max_rev))
#
seq_plus <- seq_pos[pos_max==middle]
seq_rev <- seq_rev_pos[pos_max_rev==middle]
DAP_score <- DAP_score[c(pos_max,pos_max_rev)==middle]
#
scores_plus <- scores_ARF_pos[,pos_max==middle]
scores_rev <- scores_ARF_rev_pos[,pos_max_rev==middle]
#
seq <- c(seq_plus,seq_rev)
scores <- cbind(scores_plus,scores_rev)
#

seq_ind <- apply(FUN=max,scores,2) > th & apply(FUN=max,scores,2) < (th + 1)
seq_sel <- seq[seq_ind]
DAP_score <- DAP_score[seq_ind]

########################### Compute Dinucleotide ###############################

pwm <- matrix(c(1,5,9,13,0,1,2,3),byrow=FALSE,nrow=4)
rownames(pwm) <- c("A","C","G","T")
seq_set <- DNAStringSet(seq_sel)
reg_size <- mean(width(seq_set))
#----------------------- read table ------------------------------------

dinuc <- read.table("table.txt",header=TRUE,sep="\t")
names_dinuc <- dinuc[,2]
dinuc <- data.matrix(dinuc[,-(1:2)])
rownames(dinuc) <- names_dinuc
dinuc <- dinuc[!sapply(names_dinuc,FUN=str_detect,pattern="RNA"),]

#---------------------- Compute Dinucleotides -----------------------


codes <- t(sapply(seq_sel,FUN=PWMscoreStartingAt,pwm=pwm,starting.at=(1:(reg_size-1))))
dim_codes <- dim(codes)

pearson <- numeric()
for (j in (1:dim(dinuc)[1]))
{
    M <- numeric()
    for (i in 1:(reg_size-1))
    {
        M <- cbind(M,dinuc[j,codes[,i]])
    }
    corel <- apply(FUN=cor,M,y=DAP_score,2)
    pearson <- rbind(pearson,corel)
    print(j)
}
rownames(pearson) <- rownames(dinuc)


library(gplots)
library(RColorBrewer)


my_palette <- colorRampPalette(c("blue","white","red"))(n = 400)


heatmap.2(pearson[6:15,(50:150)],
  ## cellnote = bound,  # same data set for cell labels
  ## main = "Correlation", # heat map title
  ## notecol="black",      # change font color of cell labels to black
  density.info="none",
  trace="none",                                      # turns off density plot inside color legend
  ## trace="none",
  Rowv=FALSE,                      # turns off trace lines inside the heat map
  symkey=FALSE,
  symbreaks=FALSE,
  margins =c(5,15),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  ## #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  Colv="NA",
  keysize=1
  )            # turn off column clustering


