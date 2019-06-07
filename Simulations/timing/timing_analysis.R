library(gdata)

d <- read.table("timingresults_sd.txt", sep = ',')
columns <- c("num.shapes", "num.cones", "dir.per.cone", "ec.curve", "time", "sd")
names(d) <- columns

times <- data.frame(matrix(vector(), 0, 6,
                       dimnames=list(c(), columns)),
                stringsAsFactors=F)
times

for(i in 1:(dim(d)[1])){
  for(j in 1:(dim(d)[2])){
    times[i,j] <- as.numeric(strsplit(as.character(d[i,j]), ':')[[1]][2])
  }
}


write.csv(times, file = "timings_sd.csv")
