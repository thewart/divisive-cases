ybar_dat <- ybarify(bdat, resp="rating")
juh <- bdat[,.(rating, nev, evconf)][
  data.table(evconf=ybar_dat$evconf, ybar=ybar_dat$ybar, ypred=shrinky_max$par$ypred), on="evconf"][
    , .(se=sd(rating/100)), by=.(nev, ypred, ybar, evconf)]
juh[, ybar_bin:=cut(ybar, breaks=c(0, .2, .3, .5, .65, 1))]
ggplot(juh[nev>0], aes(y=se, x=nev, color=ordered(nev))) + geom_smooth(aes(color=NULL), method="lm", se=F, color="black") + geom_jitter(width=0.1, height=0) + 
  facet_wrap(vars(ybar_bin), nrow=1, strip.position = "bottom") + ylab("Rating SD") + xlab("") + scale_color_viridis_d("Evidence count", end=0.9) + theme(legend.position = "top")
