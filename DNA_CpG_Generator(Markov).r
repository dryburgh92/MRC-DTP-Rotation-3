## CpG Island and Non-CpG island Sequence Generator 
###################################################################################################################
#Libraries 
library(compiler)
library(rtracklayer)
library(GenomicFeatures)


###################################################################################################################
#    Example of a Markov Chain of first order (the current nucleotide only
#    depends on the previous nucleotide).

#    The multinomial model per base is:
#    p(A|A)+p(C|A)+p(G|A)+p(T|A) = 1.0
#    p(A|C)+p(C|C)+p(G|C)+p(T|C) = 1.0
#    p(A|G)+p(C|G)+p(G|G)+p(T|G) = 1.0
#    p(A|T)+p(C|T)+p(G|T)+p(T|T) = 1.0  


alp <- c("A", "C", "G","T")
#Order = A, C, G, T (row1-4) & (col1-4)
#### Create the matrix that will store the probability distribution given 
# a certain nucleotide:
# Add the probability distribution per base:
pos_mod <- matrix(NA,nr=4,nc=4)
colnames(pos_mod) <- rownames(pos_mod) <- alp
pos_mod[1,] <- c(0.18,0.27,0.43,0.12)# A probability transitions
pos_mod[2,] <- c(0.17,0.37,0.27,0.19)# C probability transitions
pos_mod[3,] <- c(0.16,0.34,0.37,0.13)# G probability transitions
pos_mod[4,] <- c(0.08,0.36,0.38,0.18)# T probability transitions

#### Matrix for Non_CpG island Probability~Negative Model
neg_mod <- matrix(NA,nr=4,nc=4)
colnames(neg_mod) <- rownames(neg_mod) <- alp
neg_mod[1,] <- c(0.30,0.21,0.28,0.21)# A probability transitions
neg_mod[2,] <- c(0.32,0.30,0.08,0.30)# C probability transitions
neg_mod[3,] <- c(0.25,0.25,0.30,0.20)# G probability transitions
neg_mod[4,] <- c(0.18,0.24,0.29,0.29)# T probability transitions

#### Log_likelihood of finding a CpG island 
likelihood.values <- round(log2(pos_mod/neg_mod), 3)
#logx >0 = more likely to generate a CpG island

#(Pos mod = >prob of CG vs neg mod = <GC)=1
#Probabilities for bases starting sequence
inProb <- c(0.4,0.1,0.1,0.4); names(inProb) <- alp
 

generateFirstOrderSeq <- function(lengthSeq,
                                  alphabet,  
                                  initialProb,
                                  firstOrderMatrix){
  #    lengthSeq = length of the sequence
  #    alphabet = alphabet that compounds the sequence 
  #    initialProb   = initial probability distribution
  #    firstOrderMatrix = matrix that stores the probability distribution of a 
  #                       first order Markov Chain
  
  # Construct the object that stores the sequence
  outputSeq <- rep(NA,lengthSeq)
  # Which is the first base:
  outputSeq[1]  <- sample(alphabet,1,prob=initialProb) 
  # Let the computer decide:
  for(i in 2:length(outputSeq)){
    prevNuc <- outputSeq[i-1]    
    currentProb <- firstOrderMatrix[prevNuc,]
    outputSeq[i] <- sample(alp,1,prob=currentProb)
  } 
  cat("** DONE: Sequence computation is complete **\n")
  return(outputSeq)
}

# Use the generateFirstOrderSeq function to generate a sequence of 1000 bases for CpG island and non-CpG sequences 
Pos_mod_seq <- generateFirstOrderSeq(1000,alp,inProb,pos_mod)
Neg_mod_seq <- generateFirstOrderSeq(1000,alp,inProb,neg_mod)













 
 





