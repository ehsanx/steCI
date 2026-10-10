#' Derivatives and influence functions of the stabilized IPACW weights (internal)
#'
#' Machinery for the Approach-2 (weight-perturbing) variant of [ci_lef()]. For
#' each prepared person-trial row it returns the derivative of the cumulative
#' stabilized weight with respect to the two logistic weight-model coefficient
#' vectors, and for each subject it returns the cluster influence function of
#' those coefficient MLEs. Together with the sensitivity `A(t)` from
#' [.aalen_if()] these propagate weight-estimation uncertainty into the
#' marginal-risk-difference influence functions.
#'
#' The analysis weight is `w_r = prod_{s<r} p_num(s)/p_denom(s)`, so
#' `d log w_r / d theta_num = sum_{s<r} (1-p_num) z_num` and
#' `d log w_r / d theta_denom = -sum_{s<r} (1-p_denom) z_denom`; the derivative
#' has been validated against finite differences to ~1e-10.
#'
#' @param dat.seq Prepared data from [ste_prepare()].
#' @param n.visit Number of visits.
#' @return list with `dwdtheta` (`n_row x q`), `IFtheta` (`n_subject x q`,
#'   aligned to `factor(dat.seq$id)` levels), and `q`.
#' @keywords internal
.weight_if <- function(dat.seq, n.visit) {
  expit <- function(x) 1 / (1 + exp(-x))
  ctrl <- dat.seq$A.baseline == 0
  # control, non-terminal rows with a defined next-L: the rows carrying a
  # non-trivial weight ratio (each trial's last row has NA L.lead1).
  use <- ctrl & (dat.seq$rownum != n.visit) & !is.na(dat.seq$L.lead1)
  ui <- which(use)

  # refit the two weight models exactly as ste_prepare() does (same coefficients)
  md <- stats::glm(Anext.equal.to.baseline ~ L.lead1, family = "binomial",
                   data = dat.seq[ctrl, , drop = FALSE])
  mn <- stats::glm(Anext.equal.to.baseline ~ as.factor(rownum) * L.baseline,
                   family = "binomial", data = dat.seq[ctrl, , drop = FALSE])
  th_d <- stats::coef(md); th_n <- stats::coef(mn)
  qd <- length(th_d); qn <- length(th_n); q <- qd + qn

  Zd <- stats::model.matrix(stats::delete.response(stats::terms(md)),
                            data = dat.seq[ui, , drop = FALSE], xlev = md$xlevels)
  Zn <- stats::model.matrix(stats::delete.response(stats::terms(mn)),
                            data = dat.seq[ui, , drop = FALSE], xlev = mn$xlevels)
  pd <- as.vector(expit(Zd %*% th_d)); pn <- as.vector(expit(Zn %*% th_n))

  # score-like vectors v = (1-p) z on the use rows; lagged cumulative sum within
  # (id, trial) gives d log w / d theta up to the previous row.
  nr <- nrow(dat.seq)
  Vd <- matrix(0, nr, qd); Vn <- matrix(0, nr, qn)
  Vd[ui, ] <- (1 - pd) * Zd; Vn[ui, ] <- (1 - pn) * Zn
  Dd <- .lag_cumsum(Vd, dat.seq$id, dat.seq$trial)
  Dn <- .lag_cumsum(Vn, dat.seq$id, dat.seq$trial)
  # The derivative below is d w / d theta for the UNTRUNCATED cumulative weight, so it must be
  # paired with the untruncated weight. Under weight = "trunc" the analysis uses
  # ipw.stab.cum.trunc, whose derivative is zero wherever truncation binds, and mixing the two
  # would silently combine one weight with another's derivative. That combination is never run
  # (the shipped analyses are weight = "stab"), so this is a guard rather than a fix.
  if (!is.null(attr(dat.seq, "weight")) && identical(attr(dat.seq, "weight"), "trunc"))
    stop(".weight_if is derived for the untruncated stabilised weight; Approach-2 with ",
         "weight = \"trunc\" is not supported.")
  wv <- dat.seq$ipw.stab.cum
  dwdtheta <- cbind(-wv * Dd, wv * Dn)
  dwdtheta[!is.finite(dwdtheta)] <- 0

  # cluster influence functions of the two logistic MLEs: IF_theta = vcov %*% score
  y_u <- dat.seq$Anext.equal.to.baseline[ui]
  vcd <- solve(t(Zd) %*% (pd * (1 - pd) * Zd))
  vcn <- solve(t(Zn) %*% (pn * (1 - pn) * Zn))
  idf <- factor(dat.seq$id[ui])
  IFth_d <- rowsum(Zd * (y_u - pd), idf) %*% t(vcd)
  IFth_n <- rowsum(Zn * (y_u - pn), idf) %*% t(vcn)
  allid <- levels(factor(dat.seq$id)); K <- length(allid)
  IFtheta <- matrix(0, K, q)
  IFtheta[match(levels(idf), allid), ] <- cbind(IFth_d, IFth_n)

  list(dwdtheta = dwdtheta, IFtheta = IFtheta, q = q)
}

# Lagged cumulative sum of the rows of matrix V within (id, trial) groups:
# row r receives the sum over rows s < r in its group (row 1 -> 0).
#' @keywords internal
.lag_cumsum <- function(V, id, trial) {
  out <- matrix(0, nrow(V), ncol(V))
  for (g in split(seq_len(nrow(V)), paste(id, trial))) {
    if (length(g) <= 1) next
    cs <- apply(V[g, , drop = FALSE], 2, cumsum)
    out[g, ] <- rbind(0, cs[-length(g), , drop = FALSE])
  }
  out
}
