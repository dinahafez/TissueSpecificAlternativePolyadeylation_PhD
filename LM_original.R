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
#args <- commandArgs(TRUE)
#permute_start <- as.numeric(args[1])
#permute_end <- as.double(args[2])

 #print (permute_start)
 #print (permute_end)



result <- list()

infile = "PAS"
outfile = "PAS_out"



permute.fit <- linearModel(infile)

	residual_permute <- permute.fit$residuals
       

	write.csv(residual_permute ,file = outfile, row.names = FALSE)




