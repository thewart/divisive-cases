bvr <- stan_glmer(bardguilt ~ 1 + scale(rating) + (1 + scale(rating) | uid),iter=400,data=bdat,family = binomial)
foo <- apply(bdat[,((0:100)-mean(rating))/sd(rating)] %*%
               t(extract(bvr$stanfit,"scale(rating)")[[1]]), 1, function(x) x + extract(bvr$stanfit,"(Intercept)")[[1]])
foo <- post_summary_dt(1/(1+exp(-foo)))
foo$rating <- 0:100

bvr_c <- stan_glmer(bardguilt ~ 0 + cut_interval(rating,7) + (0 + cut_interval(rating,7) | uid), iter=1000, data=bdat, family = binomial)
juh <- 1/(1+exp(-extract(bvr_c$stanfit,"beta")[[1]])) %>% post_summary_dt()
juh$rating <- seq(0, 100, length.out=15)[seq(2, 15, 2)]
rvb_plt <- ggplot(foo, aes(y=mean,ymax=ub,ymin=lb,x=rating)) + geom_line() + geom_ribbon(alpha=0.25) + geom_pointrange(data=juh) +
  ylab("Guilt probability") + scale_x_continuous("Case strength rating (single response)", labels=c(0, 0.25, 0.5, 0.75, 1))

rvb_case_plt <- ggplot(NULL, aes(y=ybarify(bdat, resp = "binary")$ybar, x=ybarify(bdat, resp = "rating")$ybar)) + 
  geom_point(alpha=0.3) + geom_smooth(method="lm", formula=y ~ poly(x,1), color="black", se=F) + 
  xlab("Case strength rating (case average)") + ylab("Guilt probability") + ylim(c(0,1)) + xlim(c(0,1))

juh <- data.table(x=1:4, y=c(summary(lm(rating/100 ~ evconf, data=bdat))$r.squared,
                             summary(lm(rating/100 ~ uid, data=bdat))$r.squared,
                             summary(lm(rating/100 ~ scenario, data=bdat))$r.squared,
                             summary(lm(rating/100 ~ evconf + uid + scenario, data=bdat))$r.squared))
r2comp <- ggplot(juh, aes(y=y, x=x)) + geom_bar(stat = "identity") + 
  ylab(expression("R"^"2")) + scale_x_continuous(NULL, breaks=1:4, labels = c("Evidence set", "Participant", "Crime", "Combined"), minor_breaks = NULL)

plot_grid(rvb_plt, rvb_case_plt, r2comp, nrow=1)
