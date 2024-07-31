makelegaldat <- function(scendat,subjdat,clickdat,pthresh=NULL) {
  msubj <- subjdat[!str_detect(uid,"test") & task_complete,uid]
  msubj <- msubj[msubj %in% subjdat[,.N,by=uid][N==1,uid]] #exclude subjects who repeated the task
  msubj <- msubj[msubj %in% clickdat[,uniqueN(scenario)==31 & .N==31,by=uid][V1==T,uid]] #exclude subjects with wrong number of reponses
  scendat <- scendat[uid %in% msubj]
  # clickdat <- clickdat[!is.na(guilty)]
  
  if ("legalguilt" %in% names(clickdat)) names(clickdat)[which(names(clickdat) == "legalguilt")] <- "bardguilt"
  
  dat <- merge(scendat,clickdat,by=c("uid","scenario"))
  
  # were there conditions?
  cond <- str_subset(colnames(subjdat),"cond")
  if (length(cond)>0) dat <- merge(dat,subjdat[,c("uid",cond),with=F],by="uid")
  # relable evidence conditions 
  if ("cond_evidence" %in% cond) {
    dat[cond_evidence == "random", cond_evidence:="balanced"]
    dat[cond_evidence == "inculpatory", cond_evidence:="defenseless"]
  }
  if ("cond_rating" %in% cond) {
    dat$cond_rating <- as.character(dat$cond_rating)
    dat[cond_rating == "0", cond_rating:="without"]
    dat[cond_rating == "1", cond_rating:="with"]
    dat[cond_rating == "without", rating:=NA]
  }
  # are there trial types?
  wcond <- str_subset(colnames(scendat),"cond")
  # if (length(wcond)>0) dat <- merge(dat,scendat[,c("uid","scenario",wcond),with=F],by=c("uid","scenario"))
  if ("cond_capstone" %in% wcond) dat[,cond_capstone:=(cond_capstone==1)]
  
  # orderize the evidence levels 
  evidord <- c("none", "clear_ex", "ambiguous", "clear_in")
  dat[,physical:=factor(physical,evidord)]
  dat[,document:=factor(document,evidord)]
  dat[,witness:=factor(witness,evidord)]
  dat[,character:=factor(character,evidord[-3])]
  # dat[, severity := scenario/max(scenario)]
  
  dat[,start:=as.POSIXct(start)]
  dat[,stop:=as.POSIXct(stop)]
  
  if (!is.null(pthresh)) {
    fail <- dat[!is.na(rating), .(cor.test(as.numeric(bardguilt), rating)$p.value, cor(as.numeric(bardguilt),rating)), by=uid][is.na(V1) | V1>pthresh | V2<0, uid]
    dat <- dat[!(uid %in% fail)]
  }
  return(dat)
}

makebalancedat <- function(dat) {
  return(dat[,.(n_exculp=str_detect(sapply(.(physical,document,witness,character),as.character),"ex") %>% sum(),
                n_inculp=str_detect(sapply(.(physical,document,witness,character),as.character),"in") %>% sum(),
                n_ambig=str_detect(sapply(.(physical,document,witness,character),as.character),"amb") %>% sum()),
             by=.(uid,scenario)])
}

readindat <- function(whichtask, dobalance=T, pthresh=NULL) {
  fpath <- paste0("data/", whichtask)
  
  catchdat <- fread(paste0(fpath, "_catches.csv"))
  clickdat <- fread(paste0(fpath, "_answers.csv"))
  scendat <- fread(paste0(fpath, "_realized_scen.csv"))
  subjdat <- fread(paste0(fpath, "_subjectinfo.csv"))[,-1,with=F]
  subjdat[,start:=as.POSIXct(start,format="%Y-%m-%d %H:%M:%S")]
  # issue: those who completed no catch trials will not show up
  setkey(subjdat,"start")
  
  ldat <- makelegaldat(scendat, subjdat, clickdat, pthresh)
  setkey(ldat,uid,question)
  ldat[,evconf:=paste0(physical,document,witness,character)]
  if (dobalance) {
    ldat <- merge(makebalancedat(ldat),ldat,by=c("uid","scenario"))
    ldat[, nev:= n_inculp + n_exculp + n_ambig]
  }
  return(ldat)
}

post_summary_dt <- function(samps,name = colnames(samps), ci=.95) {
  if (length(dim(samps)) %in% c(0,1))
    samps <- matrix(samps,ncol=1)
  return(data.table(mean=colMeans(samps),lb=apply(samps,2,quantile,prob=(1-ci)/2),
                    ub=apply(samps,2,quantile,prob=0.5+ci/2),level=name))
}

ybarify <- function(bdat, resp="binary") {
  ybar_dat <- bdat[,.(ybar=ifelse(resp=="binary", mean(bardguilt), mean(rating)/100),
                      se=ifelse(resp=="binary", sd(bardguilt)/sqrt(.N), sd(rating/100)/sqrt(.N)),
                      M=.N), by=.(physical,document,witness,character,n_exculp,n_inculp,n_ambig)]
  ybar_dat[,evconf:=paste0(physical,document,witness,character)]
  ybar_dat[,nev:= n_inculp + n_exculp + n_ambig]
  # fit <- lmer(rating/100 ~ 0 + factor(evconf, levels=unique(evconf)) + (1|uid) + (1|scenario), data=bdat)
  # ybar_dat[,`:=`(ybar=fixef(fit), se=vcov(fit) |> diag() |> sqrt())]
  return(ybar_dat)
}

fit_onesided <- function(ybar_dat, model) {
  ybar_onesided <- get_onesided(ybar_dat)
  X <- model.matrix(~ physical + document + witness + character, ybar_onesided)[,-1]
  X_pred <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_onesided$ybar, M=ybar_onesided$M, m0=1, se=ybar_onesided$se,
                  X_pred=X_pred, N_pred=nrow(X_pred))
  return(optimizing(model, standat, as_vector=F, hessian=T))
}

fit_singleton <- function(ybar_dat, model) {
  ybar_single <- ybar_dat[nev < 2]
  X <- model.matrix(~ physical + document + witness + character, ybar_single)[,-1]
  X_pred <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_single$ybar, M=ybar_single$M, se=ybar_single$se, m0=1,
                  X_pred=X_pred, N_pred=nrow(X_pred))
  return(optimizing(model, standat, as_vector=F))
}

get_onesided <- function(df) return(df[((n_exculp>0) + (n_inculp>0) + (n_ambig>0)) %in% c(0,1)])

plot_ypred_poly <- function(ybar_dat, ypred, degree=3) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  plt <- ggplot(ybar_comp[nev>0], aes(y=ybar, x=yhat, color=ordered(nev))) + geom_point(alpha=0.3) + geom_abline() + 
    geom_smooth(method="lm", formula=y ~ poly(x, degree), se=F) + xlab("Sum of evidence weights (compressed)") + ylab("Observed case strength") +
    scale_colour_viridis_d("Evidence \n count", end=0.9)
  return(plt)
}

plot_ypred_poly_onesided <- function(ybar_dat, ypred, degree=3, evidscheme=NULL) {
  if (is.null(evidscheme)) evidscheme <- c("Baseline"="black","Inculpatory"="#ba4040ff","Exculpatory"="#406bbaff", "Ambiguous"="#765884ff")
  ybar_dat <- valence_from_evcounts(ybar_dat)
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  plt <- ggplot(ybar_comp, aes(y=ybar, x=yhat)) + geom_pointrange(aes(ymin=ybar-2*se, ymax=ybar+2*se, color=valence)) + geom_abline(alpha=0.5, linetype=2) +
    geom_smooth(method="lm", formula= y ~ poly(x, 3), se=F, color="black") + scale_color_manual("Evidence valence", breaks=names(evidscheme)[-1], values=evidscheme) +
    ylab("Observed case strength") + xlab("Sum of weights")
  return(plt)
}

lm_est_int <- function(ybar_dat, ypred, degree=1) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  return(ybar_comp[nev>0,lm(ybar-yhat ~ 1 + poly(yhat, degree) + poly(yhat, degree):nev)])
}

lm_est_resid <- function(ybar_dat, ypred, degree=1) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  return(ybar_comp[nev>0,lm(ybar-yhat ~ 1 + poly(yhat, degree))])
}

test_int <- function(ybar_dat, ypred, degree=1) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  return(anova(lm_est_resid(ybar_dat, ypred, degree), lm_est_int(ybar_dat, ypred, degree)))
}

fit_full <- function(ybar_dat, model) {
  X <- X_pred <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_dat$ybar, M=ybar_dat$M, m0=1, se=ybar_dat$se,
                  X_pred=X_pred, N_pred=nrow(X_pred))
  return(optimizing(model, standat, as_vector=F, hessian=T))
}

fit_divglm <- function(ybar_dat, model, initfrom, type="count") {
  X <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  if (type=="count") {
    Z <- model.matrix(~ factor(nev), data=ybar_dat)[, -1]
  } else if (type=="flex") {
    Z <- X
  } else if (type=="linear_count") {
    Z <- model.matrix(~ nev, data=ybar_dat)[, -1]
    dim(Z) <- c(length(Z), 1)
  } else if (type=="type") {
    Z <- model.matrix(~ (physical!="none") + (document!="none") + (character!="none") + (witness!="none"), ybar_dat)[,-1]
  } else if (type=="valence") {
    Z <- model.matrix(~ n_exculp + n_ambig + n_inculp, data=ybar_dat)[, -1]
  } else if (type=="type+valence") {
    Z <- model.matrix(~ (physical!="none") + (document!="none") + (character!="none") + (witness!="none") + 
                        n_exculp + n_ambig + n_inculp, data=ybar_dat)[, -1]
  } else {
    Z <- matrix()
  }
  
  if (!is.vector(initfrom)) initfrom <- fit_full(ybar_dat, initfrom)$par[1:3]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_dat$ybar, M=ybar_dat$M, m0=1, se=ybar_dat$se, Z=Z, Q=ncol(Z))
  return(optimizing(model, standat, init=initfrom, as_vector=F, hessian=T))
}

fit_divglm2 <- function(ybar_dat, model1, model2, byval=F) {
  fit1 <- fit_divglm(ybar_dat, model1, "linear_count")
  init <-fit1$par[1:4]
  if (byval) {
    init$lambda_0 <- rep(init$gamma, 3)
    init$lambda_1 <- rep(0, 3)
  } else {
    init$lambda <- c(init$gamma, 0)
  }
  
  X <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_dat$ybar, M=ybar_dat$M, m0=1, se=ybar_dat$se)
  return(optimizing(model2, standat, as_vector=F, hessian=T, init=init))
}

valence_from_evcounts <- function(dt) {
  dt[n_exculp>0, valence := "Exculpatory"]
  dt[n_inculp>0, valence := "Inculpatory"]
  dt[n_ambig>0, valence := "Ambiguous"]
  dt[is.na(valence), valence := "Baseline"]
  dt[, valence := factor(valence, levels=c("Baseline","Exculpatory","Ambiguous","Inculpatory"))]
  return(dt)
}

compare_models <- function(model1, model2, parind1=4, parind2=3) {
  kdiff <- purrr::list_c(model1$par[1:parind1]) |> length() - purrr:::list_c(model2$par[1:parind2]) |> length()
  lldiff <- model1$value - model2$value
  lrtest <- 1-pchisq(2*lldiff, kdiff)
  lr_pred <- exp(lldiff - kdiff) 
  return(data.table(kdiff=kdiff, lldiff=lldiff, pval=lrtest, lr_pred=lr_pred))
}

pdiv <- function(q, pars, nev) return(pt(q/c(exp(pars$gamma*nev)), df=exp(pars$lognu)))
divreg <- function(x, a, b) return(1- 1/exp(a + abs(x)*b))

maketable <- function(ybar_dat, shrinky_model, divnorm_model, divreg_model, divregbyval_model, response) {
  shrinky_fit <- fit_full(ybar_dat, shrinky_model)
  initfrom <- shrinky_fit$par[1:3]
  div_count_fit <- fit_divglm(ybar_dat, divnorm_model, c(initfrom, list(gamma=array(0, dim=1))), type="linear_count")
  div_reg_fit <- fit_divglm(ybar_dat, divreg_model, c(initfrom, list(l0=0, l1=0)),  type="dontworryaboutit")
  div_regbyval_fit <- fit_divglm(ybar_dat, divregbyval_model, c(initfrom, list(l0=rep(0, 3), l1=rep(0, 2))), type="dontworryaboutit")
  div_flex_fit <- fit_divglm(ybar_dat, divnorm_model, c(initfrom, list(gamma=rep(0, 11))), type="flex")
  
  foo <- rbind(compare_models(shrinky_fit, shrinky_fit, parind1 = 3),
               compare_models(div_count_fit, shrinky_fit),
               compare_models(div_reg_fit, div_count_fit, 5, 4),
               compare_models(div_regbyval_fit, div_reg_fit, 5, 5),
               compare_models(div_flex_fit, div_regbyval_fit, parind2 = 5))
  foo <- cbind(response=response, model=c("SC (no normalization)", "SNC Averaging","SNC W Regression", "SNC WxV Regression", "SNC Flexible"), foo)
  foo[, pval:=formatC(pval, digits=1, format="g") |> str_replace("\\+0?", "\\+") |> str_replace("\\-0?", "\\-")]
  foo[, lr_pred:=formatC(lr_pred, digits=2, format="g") |> str_replace("\\+0?", "\\+") |> str_replace("\\-0?", "\\-")]
  return(foo)
}

model_comps_rating <- function(data) return(
  maketable(ybarify(data, resp="rating"), shrinky, shrinky_div, shrinky_divreg, shrinky_divreg_byval, "Case strength"))
model_comps_binary <- function(data) return(
  maketable(ybarify(data, resp="binary"), shrinky_beta, shrinky_beta_div, shrinky_beta_divreg, shrinky_beta_divreg_byval, "Guilt judgment"))

H2S <- function(hessian, stdev=T) {
  S <- solve(-hessian)
  if (stdev) S <- diag(S) |> sqrt()
  return(S)
}

parse_evidence <- function(X,baseline=T) {
  X[str_detect(evidence,"physical"),type:="Physical"]
  X[str_detect(evidence,"document"),type:="Document"]
  X[str_detect(evidence,"witness"),type:="Witness"]
  X[str_detect(evidence,"character"),type:="Character"]
  X[str_detect(evidence,"clear_ex"),valence:="Exculpatory"]
  X[str_detect(evidence,"ambiguous"),valence:="Ambiguous"]
  X[str_detect(evidence,"clear_in"),valence:="Inculpatory"]
  # X[,c(".chain",".iteration","evidence",".variable"):=NULL]
  if (baseline) X[is.na(type),c("type","valence"):="Baseline"]
  
  X[,type:=factor(type,levels=c("Baseline","Physical","Document","Witness","Character"))]
  X[,valence:=factor(valence,levels=c("Baseline","Exculpatory","Ambiguous","Inculpatory"))]
  return(X)
}

