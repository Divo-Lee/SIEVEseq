#' Fit the skew-normal distribution to CLR-transformed RNA-Seq data
#'
#' @description Estimate the mean, standard deviation, and skewness parameters of
#'              the skew-normal distribution using centered log-ratio (CLR)
#'              transformed RNA-Seq data.
#'
#' @param data A table of CLR-transformed count data, with genes/transcripts on the rows and
#'             samples on columns.
#'
#' @return
#' \item{mu}{The maximum likelihood estimate of the mean parameter.}
#' \item{se.mu}{The standard error of the maximum likelihood estimate of \code{mu}.}
#' \item{z.mu}{The Wald statistic for \code{mu}.}
#' \item{p.mu}{The p-value of the Wald statistic for \code{mu}.}
#' \item{sigma}{The maximum likelihood estimate of the standard deviation parameter.}
#' \item{se.sigma}{The standard error of the maximum likelihood estimate of \code{sigma}.}
#' \item{z.sigma}{The Wald statistic for \code{sigma}.}
#' \item{p.sigma}{The p-value of the Wald statistic for \code{sigma}.}
#' \item{gamma}{The maximum likelihood estimate of the skewness parameter.}
#' \item{se.gamma}{The standard error of the maximum likelihood estimate of \code{gamma}.}
#' \item{z.gamma}{The Wald statistic for \code{gamma}.}
#' \item{p.gamma}{The p-value of the Wald statistic for \code{gamma}.}
#'
#'
#'
#'
#' @references
#'  Azzalini, A. (1985). \emph{A class of distributions which
#'  includes the normal ones}. \emph{Scandinavian Journal of Statistics} \bold{12}(2),
#'  171--178, JSTOR
#'
#'  Azzalini, A. and Capitanio, A. (2014).
#'  \emph{The Skew-Normal and Related Families}.
#'  Cambridge University Press, IMS Monographs series.
#'
#'  Azzalini, A. and Arellano-Valle, R. B. (2013).
#'  Maximum penalized likelihood estimation for skew-normal and skew-\emph{t}
#'  distributions. \emph{Journal of Statistical Planning and Inference} \bold{143}, 419--433.
#'
#'  Azzalini, A. (2022). \emph{The R package \bold{sn}: The skew-normal
#'  and related distribution such as the skew-t and the SUN (version 2.0.2)}.
#'  Universit\`a degli Studi di Padova, Italia.
#'  Home page: \url{https://cran.r-project.org/package=sn}.
#'
#'  Aitchison, J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#'  Chapman & Hall, London.
#'
#'
#' @import sn
#'
#' @examples
#'  library(SIEVEseq)
#'  data(clrCounts1)
#'  clr.SN.fit(clrCounts1[1:2, ])
#'  clr.SN.fit(clrCounts1[1, ])
#'
#'
#'@export
clr.SN.fit <- function(data){
  # data is only for one group
  # data must be clr-transformed counts

  if (is.matrix(data) == FALSE) {
    # only one particular gene
    esti._mat <- c(rep(NA, 12))
    names(esti._mat) <- c("mu", "se.mu", "z.mu", "p.mu",
                         "sigma", "se.sigma", "z.sigma", "p.sigma",
                         "gamma", "se.gamma.", "z.gamma", "p.gamma")
    fit_sn <- selm(data ~ 1,
                   family="SN") # regular MLE without penalty, centered parameters (CP) skew-normal (SN)
    if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
      esti._mat[1:4] <- summary(fit_sn)@param.table[1,]
      esti._mat[5:8] <- summary(fit_sn)@param.table[2,]
      esti._mat[9:12] <- summary(fit_sn)@param.table[3,]
    } else {
      fit_sn <- selm(data ~ 1,
                     method = "MPLE", # "Qpenalty"
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        esti._mat[1:4] <- summary(fit_sn)@param.table[1,]
        esti._mat[5:8] <- summary(fit_sn)@param.table[2,]
        esti._mat[9:12] <- summary(fit_sn)@param.table[3,]
      } else {
        fit_sn <- selm(data ~ 1,
                       method = "MPLE",
                       penalty = "MPpenalty",
                       family="SN")
        esti._mat[1:4] <- summary(fit_sn)@param.table[1,]
        esti._mat[5:8] <- summary(fit_sn)@param.table[2,]
        esti._mat[9:12] <- summary(fit_sn)@param.table[3,]
      }
    }
    return(esti._mat)
  } else {
    d1 <- dim(data)[1]
    # number of genes >= 2
    esti._mat <- matrix(nrow = d1, ncol = 12)
    colnames(esti._mat) <-c("mu", "se.mu", "z.mu", "p.mu",
                            "sigma", "se.sigma", "z.sigma", "p.sigma",
                            "gamma", "se.gamma.", "z.gamma", "p.gamma")
    for (i in 1:d1) {
      fit_sn <- selm(data[i,] ~ 1,
                     family="SN")
      if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
        esti._mat[i, 1:4] <- summary(fit_sn)@param.table[1,]
        esti._mat[i, 5:8] <- summary(fit_sn)@param.table[2,]
        esti._mat[i, 9:12] <- summary(fit_sn)@param.table[3,]
      } else {
        fit_sn <- selm(data[i,] ~ 1,
                       method = "MPLE",
                       family="SN")
        if(sum(is.na(summary(fit_sn)@param.table[, 4])) == 0){
          esti._mat[i, 1:4] <- summary(fit_sn)@param.table[1,]
          esti._mat[i, 5:8] <- summary(fit_sn)@param.table[2,]
          esti._mat[i, 9:12] <- summary(fit_sn)@param.table[3,]
        } else {
          fit_sn <- selm(data[i,] ~ 1,
                         method = "MPLE",
                         penalty = "MPpenalty",
                         family="SN")
          esti._mat[i, 1:4] <- summary(fit_sn)@param.table[1,]
          esti._mat[i, 5:8] <- summary(fit_sn)@param.table[2,]
          esti._mat[i, 9:12] <- summary(fit_sn)@param.table[3,]
        }
      }
    }
    row.names(esti._mat) <- row.names(data)
    return(esti._mat)
  }
}
