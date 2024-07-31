ybar_dat <- ybarify(bdat, resp="rating")

shrinky_max <- fit_full(ybar_dat, shrinky)
alligator <- bdat[(n_exculp==0 | n_inculp==0) & (character=="none"), .(ybar=mean(rating/100), se=sd(rating/100)/sqrt(.N)), by=.(n_inculp, n_exculp, n_ambig)]
alligator <- rbind(alligator[n_inculp==1][, `:=` (valence="Inculpatory", n_clear=n_inculp)],
                   alligator[n_exculp==1][, `:=` (valence="Exculpatory", n_clear=n_exculp)]) |> cbind(source="Observed")

ybar_dat$yhat <- shrinky_max$par$yhat
alligator_shrinky <- ybar_dat[(n_exculp==0 | n_inculp==0) & (character=="none"), .(ybar=mean(yhat), se=0), by=.(n_inculp, n_exculp, n_ambig)]
alligator_shrinky <- rbind(alligator_shrinky[n_exculp==1][, `:=` (valence="Exculpatory", n_clear=n_exculp)],
                           alligator_shrinky[n_inculp==1][, `:=` (valence="Inculpatory", n_clear=n_inculp)]) |> cbind(source="Summation-compression")

divglm <- fit_divglm(ybar_dat, shrinky_div, shrinky, type="flex")
ybar_dat$yhat <- divglm$par$yhat
alligator_div <- ybar_dat[(n_exculp==0 | n_inculp==0) & (character=="none"), .(ybar=mean(yhat), se=0), by=.(n_inculp, n_exculp, n_ambig)]
alligator_div <- rbind(alligator_div[n_exculp==1][, `:=` (valence="Exculpatory", n_clear=n_exculp)],
                       alligator_div[n_inculp==1][, `:=` (valence="Inculpatory", n_clear=n_inculp)]) |> cbind(source="Sum-normalize-compress")
alligator <- rbind(alligator, alligator_shrinky, alligator_div)

ggplot(alligator, aes(y=ybar, x=n_ambig, ymin=ybar-2*se, ymax=ybar+2*se, color=valence)) + 
  facet_grid(factor(valence, levels=unique(valence)) ~ factor(source, levels=unique(source)), scales = "free_y") + 
  geom_smooth(method="lm", se=F, color="black") + geom_pointrange() + theme(legend.position = "none") +
  scale_color_manual("Evidence valence",breaks = names(evidscheme)[-1],values=evidscheme) + scale_x_continuous("Ambiguoius evidence count", breaks=seq(0,2)) + scale_y_continuous("Case strength rating", breaks=scales::breaks_extended(n=4))
