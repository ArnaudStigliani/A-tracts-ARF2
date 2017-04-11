rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)




#-------------------------------------read matrices ------------------------------------------


pfm_ARF<- read.table("m_ARF2.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF <- round((t(as.matrix(pfm_ARF)))*nRegion)+1 ;pfm_ARF
maxi_ARF <- apply(pfm_ARF,FUN=max, 2)
maxi_ARF <- matrix(nrow=4, rep(maxi_ARF,4),byrow=TRUE)
pwm_ARF <- log(pfm_ARF/maxi_ARF)
pwm_ARF_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 
pwm_ARF <-  reverseComplement(pwm_ARF_rev) ; pwm_ARF

#-------------------------------------read fasta-----------------------------------------------

# pos
ARF_pos <- readDNAStringSet('ARF2.fas')
width_pos <- width(ARF_pos)
seq_pos <- as.character(ARF_pos)
seq_rev_pos <- as.character(reverseComplement(ARF_pos))

# neg
ARF_neg <- readDNAStringSet('ARF2_neg.fas')
width_neg <- width(ARF_neg)
seq_neg <- as.character(ARF_neg)
seq_rev_neg <- as.character(reverseComplement(ARF_neg))




#-------------------------------------Compute Scores-----------------------------------------

# pos
#
scores_ARF_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
scores_ARF_rev_pos<- mapply(seq_rev_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
maxi <- apply(FUN=max,scores_ARF_pos,2)
#
seq_pos_final <- ifelse(apply(FUN=max,scores_ARF_pos,2) > apply(FUN=max,scores_ARF_rev_pos,2) , seq_pos , seq_rev_pos)
#
scores_pos<- mapply(seq_pos_final,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
seq_pos <- seq_pos[apply(FUN=which.max,scores_pos,2) < 150]
#
scores_pos <- scores_pos[,apply(FUN=which.max,scores_pos,2) < 150]
#
pos_pos <- apply(FUN=which.max,scores_pos,2)
# neg
#
scores_ARF_neg <- mapply(seq_neg,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
scores_ARF_rev_neg <- mapply(seq_rev_neg,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
maxi <- apply(FUN=max,scores_ARF_neg,2)
#
seq_neg_final <- ifelse(apply(FUN=max,scores_ARF_neg,2) > apply(FUN=max,scores_ARF_rev_neg,2) , seq_neg , seq_rev_neg)
#
scores_neg <- mapply(seq_neg_final,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
seq_neg <- seq_neg[apply(FUN=which.max,scores_neg,2) < 150]
#
scores_neg <- scores_neg[,apply(FUN=which.max,scores_neg,2) < 150]
#
pos_neg <- apply(FUN=which.max,scores_neg,2)

#-------------------------------------retrieve seq after scores---------------------------------------

after_pos_site <-  mapply(str_sub,seq_pos,pos_pos+10,pos_pos+30)
after_neg_site <-  mapply(str_sub,seq_neg,pos_neg+10,pos_neg+30)

list_motifs <- list('A','AA','AAA','AAAA','T','TT','TTT','TTTT','AT','TA','AAT','TTA','ATT','TAA','AAAT','TTTA','AATT','TTAA','ATTT','TAAA')
mult_pos <- numeric()
mult_neg <- numeric()
for (elt in list_motifs)
{
    Atracts_pos <- unlist(str_locate_all(after_pos_site,elt))
    Atracts_neg <- unlist(str_locate_all(after_neg_site,elt))
    list_pos <- Atracts_pos[!is.na(Atracts_pos)]  
    list_neg <- Atracts_neg[!is.na(Atracts_neg)] 
    mult_pos <- c(mult_pos,length(list_pos)/length(seq_pos))
    mult_neg <- c(mult_neg,length(list_neg)/length(seq_neg))
}
posneg <- round(cbind(mult_pos,mult_neg),4)
rownames(posneg) <- list_motifs
colnames(posneg) <- c('pos','neg')
write.table(posneg,"AT_content_10_to_30_bases_after_site.csv",quote=FALSE,sep="\t")
