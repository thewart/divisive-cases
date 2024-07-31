onesided_dat <- get_onesided(ybarify(bdat, "rating"))
singleton_fit <- lmer(rating/100 ~ 1 + physical + document + witness + character + (1|uid) + (1|scenario), data=bdat[nev %in% c(0,1)])
onesided_fit <- lmer(rating/100 ~ 0 + factor(evconf, levels=unique(evconf)) + (1|uid) + (1|scenario), data=get_onesided(bdat))

# panel a
effdt <- data.table(mean=fixef(singleton_fit), se=sqrt(diag(vcov(singleton_fit))), evidence=names(fixef(singleton_fit))) |> parse_evidence()
effplt <- ggplot(effdt,aes(y=mean,x=type, color=valence)) + geom_pointrange(aes(ymin=mean-2*se, ymax=mean+2*se), position = position_dodge(width=0.3)) +
  xlab("Evidence type") + ylab("Evidence weigths (points)") + geom_hline(yintercept = 0) + scale_x_discrete(drop=F) +
  scale_color_manual("Evidence valence",breaks = names(evidscheme)[-1],values=evidscheme) 

# panel b
onesided_lmer <- copy(onesided_dat)
onesided_lmer[, `:=` (ybar=fixef(onesided_fit), se=vcov(onesided_fit) |> diag() |> sqrt())] |> valence_from_evcounts()
saturated <- onesided_lmer[n_exculp==nev | n_inculp==nev][(nev %in% c(0,4)) | (nev==1 & physical %in% c("clear_ex", "clear_in"))]
saturated[, `:=` (valence1=valence %in% c("Baseline","Inculpatory"), valence2=valence %in% c("Baseline","Exculpatory"))]
phalf <-  ggplot(saturated, aes(y=ybar, x=ordered(nev), color=valence)) + 
  geom_line(data=saturated[valence1==T], aes(group=valence1), color=evidscheme[2]) + 
  geom_line(data=saturated[valence2==T], aes(group=valence2), color=evidscheme[3]) + 
  geom_pointrange(aes(ymin=ybar-2*se, ymax=ybar+2*se)) + scale_color_manual("Evidence valence", breaks=names(evidscheme)[-1], values=evidscheme) + 
  scale_x_discrete(NULL,labels=c("No evidence", "Physical only", "All evidence")) + ylab("Case strength rating") + theme(legend.position = "none")

# panel c
p1 <- plot_ypred_poly_onesided(onesided_dat, predict(singleton_fit, newdata=get_onesided(onesided_dat), re.form=NA)) + xlab("Sum of weights \n (singleton fit)")

ypred <- lmer(rating/100 ~ 1 + physical + document + witness + character + (1 | uid), data=get_onesided(bdat)) |>
  predict(newdata=onesided_dat, re.form=NA)
p2 <- plot_ypred_poly_onesided(onesided_dat, ypred) + xlab("Sum of weights \n (single-valence fit)")

shrinky_max_onesided <- fit_onesided(onesided_dat, shrinky)
p3 <- plot_ypred_poly_onesided(onesided_dat, shrinky_max_onesided$par$yhat) + xlab("Compressed sum of weights \n (single-valence fit)")

pb123 <- egg::ggarrange(p1 + theme(legend.position = "none") + scale_y_continuous(limits = c(0,1), breaks=c(0, .25, .5, .75, 1), expand=c(.01,.01)),
                        p2 + theme(legend.position = "none") + scale_y_continuous(NULL, limits = c(0,1), breaks=c(0, .25, .5, .75, 1), expand=c(.01,.01), labels = NULL) + scale_x_continuous(breaks=c(0.33, 0.67, 1)),
                        p3 + theme(legend.position = "none") + scale_y_continuous(NULL, limits = c(0,1), breaks=c(0, .25, .5, .75, 1), labels = NULL, expand=c(.01,.01)), nrow=1)

# panel d
p4 <- ggplot(NULL) + xlim(c(-3, 3)) + geom_function(fun=pnorm, args=list(sd=sqrt(3)), aes(linetype="Probit")) +
  geom_function(fun=pt, args=list(df=3), aes(linetype="t CDF")) + geom_function(fun=plogis, args=list(scale=3/pi), aes(linetype="Logit")) + 
  xlab("Sum of weights \n (uncompressed)") + scale_y_continuous("Case strength rating", breaks=c(0, .5, 1)) +
  scale_linetype_manual(NULL, values=c("t CDF"="solid", "Probit"="dashed", "Logit"="dotted")) + theme(legend.position = "left")

# assembly
pb <- plot_grid(pb123, p4, nrow=1, rel_widths = c(0.65, 0.35), labels = c("c","d"), vjust = c(0.5, 1.5))
ptop <- plot_grid(effplt + theme(legend.position = "top"), phalf, nrow=1, labels="auto", rel_widths = c(1, 0.6), axis="tblr", align = "h")
plot_grid(ptop, pb, ncol=1)
