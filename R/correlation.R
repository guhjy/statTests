#' Compute a p-value for a pearson correlation using permutation.
pearson.pvalue <- function(a, b, R=10000) {
    pear.obs <- cor(a,b) # Pearson Correlation
    pear <- c()
    for (i in 1:R){
        new.a  <- permute(a)
        pear[i] <- cor(new.a,b)
    }
    return(mean(abs(pear)>=abs(pear.obs)))
}