# Robustness study for the appendix: repeat the E1 (constant hazard
# difference) coverage comparison under a DGM with TWO time-varying confounders,
# to show the row-level SE anticonservatism (and the cluster-robust fix) is not an
# artefact of the single-confounder mechanism. Self-contained (does not modify the
# single-confounder steCI package). Reduced grid; run in background.
#   Rscript dev/robustness_2confounder.R
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
expit <- function(x) 1/(1+exp(-x))
K <- 6L                         # visits
G0 <- -1; GL <- 0.4             # treatment model intercept / confounder slope (per L)
AA <- -0.04; AL1 <- 0.015; AL2 <- 0.010; AU <- 0.015   # hazard coefficients

# --- 2-confounder DGM (structurally different: two confounders, distinct dynamics) ---
dat_sim2 <- function(n, alpha0, seed) {
  set.seed(seed)
  U <- rnorm(n, 0, 0.1)
  L1 <- L2 <- A <- matrix(nrow = n, ncol = K)
  L1[,1] <- rnorm(n, U, 1); L2[,1] <- rnorm(n, 0.5*U, 1)
  A[,1] <- rbinom(n, 1, expit(G0 + GL*L1[,1] + GL*L2[,1]))
  for (k in 2:K) {
    L1[,k] <- rnorm(n, 0.8*L1[,k-1] - A[,k-1] + 0.1*(k-1) + U, 1)
    L2[,k] <- rnorm(n, 0.6*L2[,k-1] - 0.5*A[,k-1] + U, 1)
    A[,k] <- ifelse(A[,k-1]==1, 1, rbinom(n, 1, expit(G0 + GL*L1[,k] + GL*L2[,k])))
  }
  Tt <- rep(NA_real_, n)
  for (k in 1:K) { u <- runif(n); h <- alpha0 + AA*A[,k] + AL1*L1[,k] + AL2*L2[,k] + AU*U
    nt <- -log(u)/h; Tt <- ifelse(is.na(Tt) & nt<1 & h>0, k-1+nt, Tt) }
  D <- ifelse(is.na(Tt), 0, 1); Tt <- ifelse(is.na(Tt), K, Tt)
  cn <- function(p) paste0(p, ".", 0:(K-1))
  colnames(A)<-cn("A"); colnames(L1)<-cn("L1"); colnames(L2)<-cn("L2")
  d <- data.frame(id=1:n, T.obs=Tt, D.obs=D, A, L1, L2)
  dl <- reshape(d, varying=c(cn("A"),cn("L1"),cn("L2")), direction="long", idvar="id")
  dl <- dl[order(dl$id, dl$time),]
  dl$time.stop <- dl$time+1; dl <- dl[dl$time < dl$T.obs,]
  dl$time.stop <- pmin(dl$time.stop, dl$T.obs)
  dl$event <- as.integer(dl$time.stop==dl$T.obs & dl$D.obs==1)
  dl <- dl %>% group_by(id) %>%
    mutate(Alag1 = replace_na(lag(A),0), visit=row_number(),
           L1.baseline=first(L1), L2.baseline=first(L2)) %>% ungroup()
  as.data.frame(dl)
}

# --- sequential trials + IPACW with two confounders; fit const-aalen -> gamma + SEs ---
fit2 <- function(dl) {
  trials <- lapply(0:(K-1), function(g) {
    x <- dl[dl$time>=g,]; x <- x %>% mutate(trial=g) %>% group_by(id) %>%
      mutate(rownum=row_number(), time.new=time-g, time.stop.new=time.stop-g,
             Alag1.baseline=first(Alag1))
    if (g>0) x <- filter(x, Alag1.baseline==0)
    x %>% mutate(A.baseline=first(A), L1.baseline=first(L1), L2.baseline=first(L2)) %>% ungroup()
  })
  ds <- bind_rows(trials) %>% group_by(id,trial) %>%
    mutate(A.lead1=lead(A), L1.lead1=lead(L1), L2.lead1=lead(L2),
           Anext=ifelse(A.lead1==A.baseline,1,0)) %>% ungroup()
  ds <- as.data.frame(ds[ds$A==ds$A.baseline,])
  ctrl <- ds$A.baseline==0; pr <- ds$rownum!=K & !is.na(ds$L1.lead1)
  md <- glm(Anext ~ L1.lead1 + L2.lead1, binomial, ds[ctrl,])
  mn <- glm(Anext ~ as.factor(rownum)*(L1.baseline+L2.baseline), binomial, ds[ctrl,])
  ds$wd <- ds$wn <- 1
  ds$wd[pr] <- predict(md, ds[pr,], type="response"); ds$wd <- ifelse(ds$A.baseline==1,1,ds$wd)
  ds$wn[pr] <- predict(mn, ds[pr,], type="response"); ds$wn <- ifelse(ds$A.baseline==1,1,ds$wn)
  ds$rat <- ds$wn/ds$wd
  ds <- ds %>% group_by(id,trial) %>%
    mutate(lg=ifelse(rownum==1,1,lag(rat)), w=cumprod(lg)) %>% ungroup()
  ds <- as.data.frame(ds)
  f0 <- aalen(Surv(time.new,time.stop.new,event)~const(A.baseline)+L1.baseline+L2.baseline,
              data=ds, n.sim=0, weights=ds$w, robust=1)
  f1 <- aalen(Surv(time.new,time.stop.new,event)~const(A.baseline)+L1.baseline+L2.baseline,
              data=ds, n.sim=0, weights=ds$w, robust=1, clusters=ds$id)
  c(gamma=unname(f0$gamma[1,1]), model=sqrt(f0$var.gamma[1,1]),
    robust=sqrt(f0$robvar.gamma[1,1]), cluster=sqrt(f1$robvar.gamma[1,1]))
}

# --- true marginal gamma via a 2-confounder pseudo-RCT ---
true_gamma2 <- function(alpha0, n=6e4, seed=2024) {
  set.seed(seed); U <- rnorm(n,0,0.1); L10 <- rnorm(n,U,1); L20 <- rnorm(n,0.5*U,1)
  arm <- function(av){ L1<-L2<-matrix(nrow=n,ncol=K); A<-matrix(av,n,K); L1[,1]<-L10; L2[,1]<-L20
    for (k in 2:K){ L1[,k]<-rnorm(n,0.8*L1[,k-1]-A[,k-1]+0.1*(k-1)+U,1); L2[,k]<-rnorm(n,0.6*L2[,k-1]-0.5*A[,k-1]+U,1) }
    Tt<-rep(NA,n); for(k in 1:K){u<-runif(n);h<-alpha0+AA*A[,k]+AL1*L1[,k]+AL2*L2[,k]+AU*U;nt<--log(u)/h;Tt<-ifelse(is.na(Tt)&nt<1&h>0,k-1+nt,Tt)}
    D<-ifelse(is.na(Tt),0,1);Tt<-ifelse(is.na(Tt),K,Tt); list(Tt=Tt,D=D) }
  a1<-arm(1); a0<-arm(0)
  dd<-data.frame(Ts=c(a1$Tt,a0$Tt),ev=c(a1$D,a0$D),A=rep(c(1,0),each=n))
  unname(aalen(Surv(Ts,ev)~const(A),data=dd,n.sim=0)$gamma[1,1])
}

grid <- expand.grid(n=c(300L,1000L), alpha0=c(0.15,0.25,0.40), event=NA)
grid$event <- c("low","med","high")[match(grid$alpha0, c(0.15,0.25,0.40))]
REPS <- 300L; z <- qnorm(0.975)
truths <- sapply(c(0.15,0.25,0.40), true_gamma2); names(truths) <- c("0.15","0.25","0.4")
cat("true marginal gamma:", paste(sprintf("%s=%.4f",names(truths),truths),collapse=" "), "\n\n")
cat("n\tevent\tmodel_cov\trobust_cov\tcluster_cov\tmodel_ser\trobust_ser\tcluster_ser\n")
for (gi in seq_len(nrow(grid))) {
  cell <- grid[gi,]; tg <- truths[[as.character(cell$alpha0)]]
  est <- matrix(NA, REPS, 4, dimnames=list(NULL, c("gamma","model","robust","cluster")))
  for (r in seq_len(REPS)) {
    d <- dat_sim2(cell$n, cell$alpha0, seed = 1000*gi + r)
    v <- tryCatch(fit2(d), error=function(e) NULL); if (!is.null(v)) est[r,] <- v
  }
  ok <- is.finite(est[,"gamma"]); e <- est[ok,,drop=FALSE]
  cov <- sapply(c("model","robust","cluster"), function(m) mean(abs(e[,"gamma"]-tg) <= z*e[,m]))
  ser <- sapply(c("model","robust","cluster"), function(m) mean(e[,m])/sd(e[,"gamma"]))
  cat(sprintf("%d\t%s\t%.3f\t\t%.3f\t\t%.3f\t\t%.2f\t\t%.2f\t\t%.2f\n",
      cell$n, cell$event, cov["model"], cov["robust"], cov["cluster"],
      ser["model"], ser["robust"], ser["cluster"]))
}
cat("\nROBUSTNESS_2CONF DONE\n")
