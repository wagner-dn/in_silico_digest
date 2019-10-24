#A script to automate the digestion of multiple protein sequences 
#using multiple restriction enzymes.
#Written by Dominique N. Wagner & Nathan R. Vaughan, May 2014
#Hope it works :)

#Here you set the directory to read your data file from
path="/Users/dom/scripts/R"
setwd(paste(path,"",sep="/"))

library(seqRFLP)

#All you need to input are your reference protein list and a list of enzyme cut sequences
Data<-read.fasta("RNA_est_cgbAssembly100_original.txt")

ENZYME_SEQUENCES <- c("GCGC")

#These are vectors to summarize the digested sequence output data
SEQUENCE<-vector(length=(length(Data)/2))
SEQUENCE.LENGTH<-vector(length=(length(Data)/2))
N.FRAGS<-vector(length=(length(Data)/2))
TAIL.FRAG<-vector(length=(length(Data)/2))
TAIL.LENGTH<-vector(length=(length(Data)/2))
MEAN.LENGTH<-vector(length=(length(Data)/2))
MEDIAN.LENGTH<-vector(length=(length(Data)/2))
MIN.LENGTH<-vector(length=(length(Data)/2))
MAX.LENGTH<-vector(length=(length(Data)/2))
ALL.FRAGS<-vector(length=(length(Data)/2))
ALL.FRAGS.LENGTH<-vector(length=(length(Data)/2))


MULTI.DIGEST<-function(SEQUENCE_IN, ENZYMES_IN)
{
  ### code to digest a RNA sequence with multiple restriction enzymes
  ### Written by Dominique N. Wagner & Nathan R. Vaughan, May 2014
  ### Use code at your own risk. No guarantees.
  #5' to 3' cut site
  
  for(i in 1:(length(ENZYMES_IN)))
  {
    if(i==1)
    {
      DIGEST.SEQUENCE.WHOLE<-SEQUENCE_IN
    }else{
      DIGEST.SEQUENCE.WHOLE<-FRAGS.DIGESTED
    }
    ENZYME.SEQUENCE<-ENZYMES_IN[i]
    Enzyme.length=nchar(ENZYME.SEQUENCE, type = "chars", allowNA = FALSE)
    CUT_SITE_3prime <- substring(ENZYME.SEQUENCE, first=1, last=(Enzyme.length/2))
    CUT_SITE_5prime <- substring(ENZYME.SEQUENCE, first=((Enzyme.length/2)+1), last=Enzyme.length) 
    N<-length(DIGEST.SEQUENCE.WHOLE)
    for(j in 1:N)
    {
      ADD.3PRIME=FALSE
      ADD.5PRIME=FALSE
      DIGEST.SEQUENCE<-DIGEST.SEQUENCE.WHOLE[i]  
      Sequence.length=nchar(DIGEST.SEQUENCE, type = "chars", allowNA = FALSE)
      end.piece<-substring(DIGEST.SEQUENCE, first=(Sequence.length+1-Enzyme.length), last=Sequence.length)
      start.piece<-substring(DIGEST.SEQUENCE, first=1, last=Enzyme.length)
        
      if(grepl(ENZYME.SEQUENCE,end.piece,ignore.case=TRUE)){ADD.3PRIME=TRUE}
      if(grepl(ENZYME.SEQUENCE,start.piece,ignore.case=TRUE)){ADD.5PRIME=TRUE}
        
      FRAG.LIST <- strsplit(DIGEST.SEQUENCE, split=ENZYME.SEQUENCE,  fixed = FALSE, perl = FALSE)
        
      RAW.FRAGS <- FRAG.LIST[[1]]
      N.TEMP <- length(RAW.FRAGS)
      if(N.TEMP > 1)
      {
        FRAGS.3PRIME <- paste(RAW.FRAGS[1:(N.TEMP-1)], CUT_SITE_3prime, sep='')#Add 3' to end of cut segments
        FRAGS <- c(FRAGS.3PRIME, RAW.FRAGS[N.TEMP])#insert the tail cut
        FRAGS.5PRIME <- paste(CUT_SITE_5prime, FRAGS[2:N.TEMP], sep='')#Add 5' to start of cut segments
        FRAGS.TEMP <- c(FRAGS[1], FRAGS.5PRIME)
      }else{
        FRAGS.TEMP <- RAW.FRAGS
      }
          
      if(ADD.3PRIME==TRUE)
      {
        FRAGS.TEMP[length(FRAGS.TEMP)] <- paste(FRAGS.TEMP[length(FRAGS.TEMP)], CUT_SITE_3prime, sep='')
        FRAGS.TEMP <- c(FRAGS.TEMP, CUT_SITE_5prime)
      }
      if(ADD.5PRIME==TRUE)
      {
        FRAGS.TEMP[1] <- paste(CUT_SITE_5prime, FRAGS.TEMP[1], sep='')
        FRAGS.TEMP <- c(CUT_SITE_3prime, FRAGS.TEMP)         
      }
          
      #if this is the first fragment we start the digested sequence fresh
      if(j==1)
      {
        FRAGS.DIGESTED <- c(FRAGS.TEMP) #insert the digested frags
      }else{ #Otherwise the fragments are appended to the sequence
        FRAGS.DIGESTED <- c(FRAGS.DIGESTED, FRAGS.TEMP) #append the digested frags
      }
    } 
  }
  #Return the fully digested sequence
  return(DIGEST.SEQUENCE.WHOLE) 
}

#This section digests each protein individually and compliles the results into the summary vectors
lastFrag=0
for(i in 1:(length(Data)/2))
{
  Data[[(2*i)]]<-gsub("a","A",Data[[(2*i)]])
  Data[[(2*i)]]<-gsub("g","G",Data[[(2*i)]])
  Data[[(2*i)]]<-gsub("c","C",Data[[(2*i)]])
  Data[[(2*i)]]<-gsub("t","T",Data[[(2*i)]])
  SEQUENCE[i]<-Data[[(2*i)]]
  FRAGMENTS <- MULTI.DIGEST(SEQUENCE[1,1], ENZYME_SEQUENCES)
  N.FRAGS[i] <- length(FRAGMENTS)
  TAIL.FRAG[i] <- FRAGMENTS[N.FRAGS[i]]
  TAIL.LENGTH[i] <- nchar(TAIL.FRAG[i], type = "chars", allowNA = FALSE)
  MEAN.LENGTH[i] <- mean(nchar(FRAGMENTS, type = "chars", allowNA = FALSE)) 
  MEDIAN.LENGTH[i] <- median(nchar(FRAGMENTS, type = "chars", allowNA = FALSE)) 
  MIN.LENGTH[i] <- min(nchar(FRAGMENTS, type = "chars", allowNA = FALSE)) 
  MAX.LENGTH[i] <- max(nchar(FRAGMENTS, type = "chars", allowNA = FALSE))
  SEQUENCE.LENGTH[i] <- nchar(SEQUENCE[i], type = "chars", allowNA = FALSE)
  for(j in 1:N.FRAGS[i])
  {
    ALL.FRAGS[(lastFrag+j)]=FRAGMENTS[j] 
    ALL.FRAGS.LENGTH[(lastFrag+j)]=nchar(FRAGMENTS[j], type = "chars", allowNA = FALSE)
  }
  lastFrag=lastFrag+N.FRAGS[i]
}

N.MIN <- min(N.FRAGS)
N.MAX <- max(N.FRAGS)
N.MEAN <- mean(N.FRAGS)
N.MEDIAN <- median(N.FRAGS)

TAIL.MIN <- min(TAIL.LENGTH)
TAIL.MAX <- max(TAIL.LENGTH)
TAIL.MEAN <- mean(TAIL.LENGTH)
TAIL.MEDIAN <- median(TAIL.LENGTH)

GLOBAL.MIN <- min(ALL.FRAGS.LENGTH)
GLOBAL.MAX <- max(ALL.FRAGS.LENGTH)
GLOBAL.MEAN <- mean(ALL.FRAGS.LENGTH)
GLOBAL.MEDIAN <- median(ALL.FRAGS.LENGTH)

SEQUENCE.MIN <- min(SEQUENCE.LENGTH)
SEQUENCE.MAX <- max(SEQUENCE.LENGTH)
SEQUENCE.MEAN <- mean(SEQUENCE.LENGTH)
SEQUENCE.MEDIAN <- median(SEQUENCE.LENGTH)


#---You can work on this to get the plots and output Matrix you want--------------------------

TAIL.LENGTH.DIST <- hist(TAIL.LENGTH, br=c(0,100, 180, 260, 340, 420, 500, 580, 660, 740, 820, 900, 980,TAIL.MAX+1), plot =FALSE)
Number_260to340 <- gg$counts[5]
NumberA <- gg$counts
jj <- hist(FRAGMENT_SIZE, br=c(0,1000, 1050, 1100, 1180, 1260, 1340, 1420, 1580, 1660, 1740, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,MAX_Size+1), plot =TRUE)
NumberB <- jj$counts

plot(x=TAIL.LENGTH.DIST$mids, y=TAIL.LENGTH.DIST$density)

RESULTS_MATRIX[i, 1] <- CHROMOSOME_LIST[i]
RESULTS_MATRIX[i, 2] <- CHROMOSOME_SIZE
RESULTS_MATRIX[i, 3] <- NUMBER_FRAGMENTS
RESULTS_MATRIX[i, 4] <- MEAN_Size
RESULTS_MATRIX[i, 5] <- MEDIAN_Size
RESULTS_MATRIX[i, 6] <- MIN_Size
RESULTS_MATRIX[i, 7] <- MAX_Size
RESULTS_MATRIX[i, 8] <- NA
RESULTS_MATRIX[i, 9:21] <- NumberA
RESULTS_MATRIX[i, 22]   <- NA
RESULTS_MATRIX[i, 23:42] <- NumberB
