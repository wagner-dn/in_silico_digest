###The following function performs an in silico double digest of a DNA sequence. 
###It has been tested on a single thread of a desktop up to the size of a large chromosome ###from the Zebra Finch.
#SEQUENCE = the DNA sequence in fasta format
#RECOGNITION_CODE = the recognition site for restriction enzyme 1
#CUT_SITE_3prime = the 3' cut site for restriction enzyme 1
#CUT_SITE_5prime = the 5' cut site
#RECOGNITION_CODE2 = enzyme 2 etc.

library(seqRFLP)

DOUBLE.DIGEST <- function(SEQUENCE, RECOGNITION_CODE, CUT_SITE_3prime, CUT_SITE_5prime, RECOGNITION_CODE2, CUT_SITE_3prime2, CUT_SITE_5prime2){
### code to digest a DNA sequence with one restriction enzyme
### Written by Jason T. Weir, 6 July 2011
### Modified by Dominique N. Wagner, May 2014
### Use code at your own risk. No guarantees.
#5' to 3' cut site
FRAGMENTS <- strsplit(SEQUENCE, split=RECOGNITION_CODE,  fixed = FALSE, perl = FALSE)
FRAGMENTS1 <- FRAGMENTS[[1]]
N <- length(FRAGMENTS1) #number of fragments after cutting
FRAGMENTS2 <- paste(FRAGMENTS1[1:(N-1)], CUT_SITE_3prime, sep='')
FRAGMENTS3 <- c(FRAGMENTS2, FRAGMENTS1[N])
FRAGMENTS4 <- paste(CUT_SITE_5prime, FRAGMENTS3[2:N], sep='')
FRAGMENTS5 <- c(FRAGMENTS3[1], FRAGMENTS4)

nchar(FRAGMENTS5, type = "chars", allowNA = FALSE)

	### code to now digest with a second restriction enzyme
FRAGMENTS6 <- strsplit(FRAGMENTS5, split=RECOGNITION_CODE2,  fixed = FALSE, perl = FALSE)
for(i in 1:N){
	Ni <- length(FRAGMENTS6[[i]])
	if(Ni > 1){
		temp <- paste(FRAGMENTS6[[i]][1:(Ni-1)], CUT_SITE_3prime2, sep='')
		FRAGMENTS6[[i]] <- c(temp, FRAGMENTS6[[i]][Ni])
		temp2 <- paste(CUT_SITE_5prime2, FRAGMENTS6[[i]][2:Ni], sep='')
		FRAGMENTS6[[i]] <- c(FRAGMENTS6[[i]][1], temp2)
	}
}
FRAGMENTS7 <- unlist(FRAGMENTS6)
FRAGMENTS7

return(FRAGMENTS7)

}
#Restriction Enzyme 1
RECOGNITION_CODE <- "CCGC" #AciI
CUT_SITE_3prime <- "CC"
CUT_SITE_5prime <- "GC"

#Restriction Enzyme 2
RECOGNITION_CODE2 <- "TCGA" #Taq_aI
CUT_SITE_3prime2 <- "TC"
CUT_SITE_5prime2 <- "GA"

CHROMOSOME_LIST <- c(
"RNA_est_cgbAssembly100_2.fa"
)

RESULTS_MATRIX <- matrix(data = NA, nrow = length(CHROMOSOME_LIST), ncol = 45, byrow = FALSE, dimnames = NULL)

for(i in 1:length(CHROMOSOME_LIST)){
	SEQUENCE <- read.fasta(file = CHROMOSOME_LIST[i]) #this is a fasta format DNA sequence for the whole chromosome
	SEQUENCE <- SEQUENCE[[2]] #The sequence is in the second element of the list after reading in the FAST format file

	OUTPUT <- DOUBLE.DIGEST(SEQUENCE, RECOGNITION_CODE, CUT_SITE_3prime, CUT_SITE_5prime, RECOGNITION_CODE2, CUT_SITE_3prime2, CUT_SITE_5prime2)

	FRAGMENT_SIZE <- nchar(OUTPUT)
	FRAGMENT_SIZE <- sort(FRAGMENT_SIZE)
	MEAN_Size <- mean(FRAGMENT_SIZE) 
	MEDIAN_Size <- median(FRAGMENT_SIZE) 
	MIN_Size <- min(FRAGMENT_SIZE) 
	MAX_Size <- max(FRAGMENT_SIZE) 
	NUMBER_FRAGMENTS <- length(OUTPUT)
	#hist(FRAGMENT_SIZE, breaks = 3000)
	CHROMOSOME_SIZE <- nchar(SEQUENCE)

	gg <- hist(FRAGMENT_SIZE, br=c(0,100, 180, 260, 340, 420, 500, 580, 660, 740, 820, 900, 980,MAX_Size+1), plot =FALSE)
	Number_260to340 <- gg$counts[5]
	NumberA <- gg$counts
	jj <- hist(FRAGMENT_SIZE, br=c(0,1000, 1050, 1100, 1180, 1260, 1340, 1420, 1580, 1660, 1740, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,MAX_Size+1), plot =TRUE)
	NumberB <- jj$counts

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
}
RESULTS_MATRIX
