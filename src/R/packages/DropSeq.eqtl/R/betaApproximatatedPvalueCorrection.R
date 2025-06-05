# This uses the fastEQTL implementation in tensorQTL
# Because of how likelihoods are computed pvalues must be > 0 and < 1.


#' Compute a pvalue correction using a Beta approximation (fastQTL)
#'
#' Fit permuted p-value distribution to beta distribution.
#' This characterizes the extreme tail of the null without directly sampling it.
#' See the paper: https://academic.oup.com/bioinformatics/article/32/10/1479/1742545
#' See the code: https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/core.py
#'
#' @param empiricP The empiric pvalue to correct
#' @param permutedPDist A distribution of permuted pvalues.  This is used to estimate the shape of the beta distribution.
#' 1000-10000 permutations should be sufficient for a well calibrated result.
#'
#' @return A data frame containing the direct permutation result and the beta approximation
#' @import fitdistrplus
#' @export
betaApproximationPvalue<-function (empiricP, permutedPDist) {

    mean=mean(permutedPDist)
    var=var(permutedPDist)
    beta_shape1=mean*(mean*(1-mean) /var-1)
    beta_shape2=beta_shape1*(1/mean-1)

    #beta parameter estimation from fastQTL implementation

    #This is validated against the python code.
    beta_likelihood<-function (x, shape1, shape2) {
        logbeta = log(gamma(shape1)) + log(gamma(shape2)) - log(gamma(shape1+shape2))
        r= (1.0-shape1)*sum(log(x)) + (1.0-shape2)*sum(log(1.0-x)) + length(x)*logbeta
        return (r)
    }

    #wrapper around the function to bundle parameters together in a single vector for optim.
    beta_likelihood_wrapper<-function (shapes, x) {
        z=beta_likelihood(x, shapes[1], shapes[2])
        return (z)
    }

    #note: p-values must be greater than 0 and less than 1.
    maxValue=0.999
    minValue=1e-200
    permutedPDistFixed=permutedPDist

    idx=which(permutedPDistFixed==0)
    if (length(idx)>0) permutedPDistFixed[idx]<-minValue

    idx=which(permutedPDistFixed==1)
    if (length(idx)>0) permutedPDistFixed[idx]<-0.999

    #optim must take a vector of parameters.
    #note: the optimization should be bounded between 0 and infinity for each parameter.
    # This and fitdistrplus seem to converge on the same answer once input data is properly bounded.
    # r=optim(par=c(beta_shape1, beta_shape2), fn=beta_likelihood_wrapper, x=permutedPDistFixed, method="L-BFGS-B", lower=c(0,0), upper=c(Inf,Inf))
    # beta_shape1=r$par[1]
    # beta_shape2=r$par[2]

    #fit the p-values to the beta distribution with the initial shape estimates
    f=fitdistrplus::fitdist(data = permutedPDistFixed, "beta", start = list(shape1 = beta_shape1, shape2 = beta_shape2))
    beta_shape1=f$estimate[1]
    beta_shape2=f$estimate[2]

    pval_perm=(length(which(permutedPDist<=empiricP))+1)/(length(permutedPDist)+1)
    #draw the new p-value from the beta distribution
    pval_beta=pbeta(empiricP, beta_shape1, beta_shape2, lower.tail = TRUE)
    data.frame(pval_perm=pval_perm, pval_beta=pval_beta)

}
