linearModel <- function (filename)
{
	mydata <- read.csv (file = filename, h=T)

	mydata$gene = as.factor(mydata$gene)
	mydata$pas = as.factor(mydata$pas)
	mydata$tissue = as.factor(mydata$tissue)	

	#fit original model
	result <-gls(log(count) ~1 + tissue + gene + (tissue * gene) ,  na.action="na.omit" ,data=mydata, method="ML")        
        return (result)
}
###############################################
library(nlme)
args <- commandArgs(TRUE)
permute_start <- as.numeric(args[1])
permute_end <- as.double(args[2])

 print (permute_start)
 print (permute_end)
#print (args)

#fit permuted files depending on the range (jobs are divided)
#permute_start <- 501
#permute_end <- 550
result <- list()

infile = "permutedFiles/PAS"
outfile = "resultFiles/PAS"
for (i in permute_start: permute_end)
{
        file_permute <- paste(infile,i, sep=".")

        permute.fit <- linearModel(file_permute)

	residual_permute <- permute.fit$residuals
       
	outfilename <- paste(outfile,i,"residual",sep=".")
	write.csv(residual_permute ,file = outfilename, row.names = FALSE)

}


