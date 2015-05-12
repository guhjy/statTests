# Log-like transformation that works with positive and negative values.
# b = a coefficient that controls the linearity of HL(y) with respect to y. Larger
# values create a more linear relationship.
# d = dynamic ranges for y values (-1,1)
# r = dynamic range for y values (-Inf,-1] and [1,Inf)
# See Bagwell, 2005 http://www.cores.utah.edu/labs/flowcytometry/HyperLog.pdf
# Adapted from PyFCM http://code.google.com/p/py-fcm/source/browse/src/core/transforms.py
hyperlog <- function(y, b, d, r, intervals=1000.0) {
    ub = log(max(y) + 1 - min(y))
    xx = exp(seq(0, ub, ub / intervals)) - 1 + min(y)
    yy = sapply(xx, function(x) uniroot(.EH, c(-10^6, 10^6), x, b, d, r)$root)
    spline(xx, yy, xout=y)
}

.EH <- function(x, y, b, d, r) {
    e = as.numeric(d) / r
    sgn = sign(x)
    expval <- 10 ^ (sgn * e * x)
    if (is.infinite(expval)) {
        sgn * .Machine$double.xmax
    }
    else {
        (sgn * expval) + (b * e * x) - sgn - y
    }
}

discretize <- function(m, levels) {
    for (i in 2:length(levels)) {
        m[m > levels[i-1] & m <= levels[i]] <- levels[i]
    }
    m[m == levels[1]] <- levels[2]
    m
}