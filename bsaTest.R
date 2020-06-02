# BSA test on counts
# ============================================================

args <- commandArgs(TRUE)

# read count data
cov <- read.table(args[1], header = FALSE, as.is = TRUE)

# function to perform the test
bsa.test <- function(n1, n2, n3, n4, nchr1, nchr2) {
  p1 <- n2/(n1 + n2)
  p2 <- n4/(n3 + n4)
  p0 <- (p1 + p2)/2
  chisq.score <- ((p1 - p2)^2) / ( p0*(1-p0)*(1/(n1 + n2) + 1/(n3 + n4) + (1/nchr1 + 1/nchr2)) )
  pval <- pchisq(chisq.score, 1, lower.tail = FALSE)
  return(pval)
}

nchr1 <- as.numeric(args[2])
nchr2 <- as.numeric(args[3])

pval <- numeric(nrow(cov))
f1 <- numeric(nrow(cov))
f2 <- numeric(nrow(cov))

for (i in 1:nrow(cov)) {
  sample1.count <- as.numeric(unlist(strsplit(cov[i, 7], split = ",")))
  sample2.count <- as.numeric(unlist(strsplit(cov[i, 9], split = ",")))
  if (cov[i, 6] == "PASS" | cov[i, 8] == "PASS") {
    pval[i] <- bsa.test(sample1.count[1], sample1.count[2], sample2.count[1], sample2.count[2], nchr1, nchr2)
  } else {
    pval[i] <- NA
  }
  f1[i] <- sample1.count[2]/sum(sample1.count)
  f2[i] <- sample2.count[2]/sum(sample2.count)
}

write.table(cbind(cov[, 1:7], f1, cov[, 8:9], f2, pval), file = args[4], sep = " ", col.names = F, row.names = F, quote = F)
