##### top half ####
ybar_dat <- ybarify(bdat, resp="rating")

shrinky_max_onesided <- fit_onesided(ybar_dat, shrinky)
p1 <- plot_ypred_poly(ybar_dat, shrinky_max_onesided$par$ypred, degree = 3)

shrinky_max <- fit_full(ybar_dat, shrinky)
p2 <- plot_ypred_poly(ybar_dat, shrinky_max$par$ypred, 1)

divglm <- fit_divglm(ybar_dat, shrinky_div, shrinky, type="linear_count")
p3 <- plot_ypred_poly(ybar_dat, divglm$par$yhat, 1)

p4 <- ggplot(NULL) + geom_function(fun=pdiv, args=list(pars=divglm$par, nev=1), aes(color="1"), linewidth=0.75) + 
  geom_function(fun=pdiv, args=list(pars=divglm$par, nev=2), aes(color="2"), linewidth=0.75) + 
  geom_function(fun=pdiv, args=list(pars=divglm$par, nev=3), aes(color="3"), linewidth=0.75) + 
  geom_function(fun=pdiv, args=list(pars=divglm$par, nev=4), aes(color="4"), linewidth=0.75) + 
  scale_color_viridis_d("Evidence \n count", end=0.9) + scale_y_continuous("Predicted case strength", breaks=scales::extended_breaks(3)) + 
  scale_x_continuous("Sum of weights", limits = c(-4, 4), labels=c("", "\U2190 Exculp.", "", "Inculp. \U2192", ""))

p123 <- egg::ggarrange(p1 + theme(legend.position = "none") + xlab("Compressed predictions \n (single-valence fit)"), 
                       p2 + scale_y_continuous(NULL, labels=NULL) + xlab("Compressed predictions \n (all cases fit)") + theme(legend.position = "none"), 
                       p3 + scale_y_continuous(NULL, labels=NULL) + xlab("Normalized predictions \n (all cases fit)"), nrow=1)
plt_count <- plot_grid(p123, p4 + theme(legend.position = "none"), rel_widths = c(0.75, 0.25), labels = "auto")

##### this is gonna get ugly for the bottom half ####
library(Matrix)
fitdiv_exp1 <- fit_divglm(ybarify(bdat, resp="rating"), shrinky_div, shrinky, type="flex")
fitdiv_exp2 <- fit_divglm(ybarify(bdat_cond, resp="rating"), shrinky_div, shrinky, type="flex")
fitreg_exp1 <- fit_divglm(ybarify(bdat, resp="rating"), shrinky_divreg_byval, shrinky, type="juh")
fitreg_exp2 <- fit_divglm(ybarify(bdat_cond, resp="rating"), shrinky_divreg_byval, shrinky, type="juh")

lentheta <- nrow(fitdiv_exp1$hessian)
lengamma <- length(fitdiv_exp1$par$gamma)
X <- rbind(diag(0.5, lentheta, lentheta), diag(0.5, lentheta, lentheta))
bigSig <- bdiag(solve(-fitdiv_exp1$hessian), solve(-fitdiv_exp2$hessian))
mu <- t(X) %*% c(purrr::list_c(fitdiv_exp1$par[1:4]), purrr::list_c(fitdiv_exp2$par[1:4]))
Sigma <- t(X) %*% bigSig %*% X

gamma <- mu[(lentheta-lengamma+1):(lentheta)]
sigma <- diag(Sigma) |> sqrt()
sigma_gamma <- sigma[(lentheta-lengamma+1):(lentheta)]
sigma_beta <- sigma[2:(lengamma+1)]

gamma <- fitdiv_exp2$par$gamma |> as.vector()
sigma <- H2S(fitdiv_exp2$hessian)
sigma_gamma <- sigma[names(sigma) |> str_detect("gamma")] |> as.vector()
sigma_beta <- sigma[names(sigma) |> str_detect("beta")] |> as.vector()
divdt <- data.table(gamma=1 - 1/exp(gamma), gamma.lb = 1 - 1/exp(gamma - 2*sigma_gamma), gamma.ub = 1 - 1/exp(gamma + 2*sigma_gamma),
                    beta=fitdiv_exp2$par$beta, beta.lb=fitdiv_exp2$par$beta - 2*sigma_beta, beta.ub=fitdiv_exp2$par$beta + 2*sigma_beta,
                    evidence=evnames, fromwhich="Experiment 2") |> parse_evidence()

gamma <- fitdiv_exp1$par$gamma |> as.vector()
sigma <- H2S(fitdiv_exp1$hessian)
sigma_gamma <- sigma[names(sigma) |> str_detect("gamma")] |> as.vector()
sigma_beta <- sigma[names(sigma) |> str_detect("beta")] |> as.vector()
divdt <- rbind(divdt, data.table(gamma=1 - 1/exp(gamma), gamma.lb = 1 - 1/exp(gamma - 2*sigma_gamma), gamma.ub = 1 - 1/exp(gamma + 2*sigma_gamma),
                                 beta=fitdiv_exp1$par$beta, beta.lb=fitdiv_exp1$par$beta - 2*sigma_beta, beta.ub=fitdiv_exp1$par$beta + 2*sigma_beta,
                                 evidence=evnames, fromwhich="Experiment 1") |> parse_evidence())

x1 <- seq(-1, -.2, .025)
x2 <- seq(.25, 1.5, .025)
x3 <- seq(-.2, .25, .025)
curve1 <- rbind(
  data.table(beta=x1, gamma=with(fitreg_exp2$par, sapply(x1, divreg, a=l0[1], b=l1[1])), valence="Exculpatory"),
  data.table(beta=x2, gamma=with(fitreg_exp2$par, sapply(x2, divreg, a=l0[3], b=l1[2])), valence="Inculpatory"),
  data.table(beta=x3, gamma=with(fitreg_exp2$par, sapply(x3, divreg, a=l0[2], b=0)), valence="Ambiguous"))
curve1[, `:=` (type="Physical", fromwhich="Experiment 2")]
curve2 <- rbind(
  data.table(beta=x1, gamma=with(fitreg_exp1$par, sapply(x1, divreg, a=l0[1], b=l1[1])), valence="Exculpatory"),
  data.table(beta=x2, gamma=with(fitreg_exp1$par, sapply(x2, divreg, a=l0[3], b=l1[2])), valence="Inculpatory"),
  data.table(beta=x3, gamma=with(fitreg_exp1$par, sapply(x3, divreg, a=l0[2], b=0)), valence="Ambiguous"))
curve2[, `:=` (type="Physical", fromwhich="Experiment 1")]
curve <- rbind(curve1, curve2)
curve[, `:=` (gamma.lb = gamma, gamma.ub = gamma)]

plt1 <- ggplot(divdt, aes(y=gamma, x=beta, fill=valence, shape=type)) + geom_line(data=curve, aes(color=valence), linewidth=0.8) + 
  geom_errorbar(aes(ymin=gamma.lb, ymax=gamma.ub), alpha=0.5, width=0) + geom_errorbarh(aes(xmin=beta.lb, xmax=beta.ub), alpha=0.5) + geom_point(size=3.5) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + facet_wrap(vars(fromwhich), nrow = 2, scales = "free_y") +
  xlab("Evidence weight") + scale_y_continuous("Suppression", breaks=scales::breaks_extended(3), labels = scales::percent) + 
  scale_fill_manual("Valence", breaks=names(evidscheme)[-1], values=evidscheme) + scale_shape_manual("Type", values=c(21:24)) + 
  guides(fill=guide_legend(override.aes=list(shape=21))) + scale_color_manual("Valence", breaks=names(evidscheme)[-1], values=evidscheme) + 
  theme(panel.grid=element_blank(), axis.line.x = element_line("gray"), axis.line.y = element_line("gray")) + coord_cartesian(xlim=c(-1.25, 1.75))

coefplt <- rbind(
  with(fitreg_exp1, rbind(
    data.table(mean = par$l0, se=H2S(hessian)[colnames(hessian) |> str_detect("l0")],
               valence=names(evidscheme)[c(3, 4, 2)], type="Intercept", fromwhich="Experiment 1"),
    data.table(mean = par$l1, se=H2S(hessian)[colnames(hessian) |> str_detect("l1")],
               valence=names(evidscheme)[c(3, 2)], type="Slope", fromwhich="Experiment 1"))),
  with(fitreg_exp2, rbind(
    data.table(mean = par$l0, se=H2S(hessian)[colnames(hessian) |> str_detect("l0")],
               valence=names(evidscheme)[c(3, 4, 2)], type="Intercept", fromwhich="Experiment 2"),
    data.table(mean = par$l1, se=H2S(hessian)[colnames(hessian) |> str_detect("l1")],
               valence=names(evidscheme)[c(3, 2)], type="Slope", fromwhich="Experiment 2"))))

plt2 <- ggplot(coefplt[type=="Intercept"], aes(y=mean, x=valence, ymin=mean-2*se, ymax=mean+2*se, color=valence)) + geom_hline(yintercept = 0) + geom_pointrange() + 
  facet_wrap(vars(fromwhich)) + scale_x_discrete(NULL, labels=c("Amb.", "Exc.", "Inc.")) + scale_y_continuous("Intercept", breaks=scales::extended_breaks(3)) +
  scale_color_manual("Valence", breaks=names(evidscheme)[-1], values=evidscheme, guide=NULL) 

plt3 <- ggplot(coefplt[type=="Slope"], aes(y=mean, x=valence, ymin=mean-2*se, ymax=mean+2*se, color=valence)) + geom_hline(yintercept = 0) + geom_pointrange() + 
  facet_wrap(vars(fromwhich)) + scale_x_discrete(NULL, labels=c("Exc.", "Inc.")) + scale_y_continuous("Slope magnitude", breaks=scales::extended_breaks(3)) +
  scale_color_manual("Valence", breaks=names(evidscheme)[-1], values=evidscheme, guide=NULL) 
plt23 <- plot_grid(plt2, plt3, rel_widths = c(1, 0.75), axis = "tbrl", nrow=2, labels = c("d", "e"))

plt_wreg <- plot_grid(plt1, plt23 + theme(plot.margin = margin(l=5)), nrow = 1, rel_widths = c(1, 0.5), labels = c("c", ""))

##### put the two halves together ####
plot_grid(plt_count, plt_wreg, ncol=1, rel_heights = c(0.5, 0.6)) #9x6 is a good aspect ratio for this one

