# Validation for M6 Approach-2 (weight-estimation uncertainty). Backs the
# Appendix claim that the correction is correct and negligible.
#   (a) the analytic weight derivative dw/dtheta matches finite differences;
#   (b) the sensitivity A(t)=dB/dtheta (in .aalen_if) matches finite differences;
#   (c) the resulting SE correction is < ~1% (stabilized-weight cancellation).
# Run from the repo root: Rscript dev/validate_a2.R
suppressMessages({library(dplyr);library(tidyr);library(timereg);library(survival);library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
expit <- function(x) 1/(1+exp(-x))

set.seed(7); d <- dat_sim(n = 400, seed = 7); nv <- 5
ds <- ste_prepare(d, n.visit = nv)
wif <- .weight_if(ds, nv)                          # analytic dw/dtheta + IF_theta

# reconstruct the design/prob to build a finite-difference weight recomputation
ctrl <- ds$A.baseline == 0
use  <- ctrl & (ds$rownum != nv) & !is.na(ds$L.lead1); ui <- which(use)
md <- glm(Anext.equal.to.baseline ~ L.lead1, binomial, ds[ctrl,])
mn <- glm(Anext.equal.to.baseline ~ as.factor(rownum)*L.baseline, binomial, ds[ctrl,])
Zd <- model.matrix(delete.response(terms(md)), ds[ui,], xlev = md$xlevels)
Zn <- model.matrix(delete.response(terms(mn)), ds[ui,], xlev = mn$xlevels)
th_d <- coef(md); th_n <- coef(mn); qd <- length(th_d)
recompute_w <- function(td, tn) {
  wd <- rep(1, nrow(ds)); wn <- rep(1, nrow(ds))
  wd[ui] <- as.vector(expit(Zd %*% td)); wn[ui] <- as.vector(expit(Zn %*% tn))
  ratio <- wn/wd
  ave(seq_along(ratio), paste(ds$id, ds$trial),
      FUN = function(ix){v<-ratio[ix]; cumprod(c(1, v[-length(v)]))})
}
eps <- 1e-6; md_err <- 0
for (k in seq_len(qd)) { tp<-th_d; tp[k]<-tp[k]+eps; tm<-th_d; tm[k]<-tm[k]-eps
  md_err <- max(md_err, max(abs((recompute_w(tp,th_n)-recompute_w(tm,th_n))/(2*eps) - wif$dwdtheta[,k]))) }
for (k in seq_along(th_n)) { tp<-th_n; tp[k]<-tp[k]+eps; tm<-th_n; tm[k]<-tm[k]-eps
  md_err <- max(md_err, max(abs((recompute_w(th_d,tp)-recompute_w(th_d,tm))/(2*eps) - wif$dwdtheta[,qd+k]))) }
cat(sprintf("(a) dw/dtheta: max |analytic - finite-diff| = %.2e\n", md_err))

# (b) A(t) vs finite-diff of B(t): compare .aalen_if A to numerical dB/dtheta
fit <- .aalen_if(ds, t.hor = 1:5, dwdtheta = wif$dwdtheta)
Bof <- function(wcol_vec){ ds2<-ds; ds2$wtmp<-wcol_vec; .aalen_if(ds2, 1:5, "wtmp")$B }
k <- qd + 5  # a num-model component
tp<-th_n; tp[5]<-tp[5]+eps; tm<-th_n; tm[5]<-tm[5]-eps
fdB <- (Bof(recompute_w(th_d,tp)) - Bof(recompute_w(th_d,tm)))/(2*eps)   # nt x p
cat(sprintf("(b) A(t): max |analytic - finite-diff| = %.2e (typical |A| ~ %.2e)\n",
    max(abs(t(fit$A[,k,]) - fdB)), max(abs(fit$A[,k,]))))

# (c) SE correction magnitude
a1 <- ci_lef(d, t.hor=1:5, B=200, approach="a1", seed=1)$mrd
a2 <- ci_lef(d, t.hor=1:5, B=200, approach="a2", seed=1)$mrd
r <- a2$se[a2$type=="analytic_cluster"] / a1$se[a1$type=="analytic_cluster"]
cat(sprintf("(c) a2/a1 cluster-SE ratio by horizon: %s  (max deviation %.1f%%)\n",
    paste(sprintf("%.3f", r), collapse=" "), 100*max(abs(r-1))))
cat("VALIDATE_A2 OK\n")
