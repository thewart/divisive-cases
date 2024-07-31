t1 <- rbind(model_comps_rating(bdat), model_comps_binary(bdat))
t2 <- rbind(model_comps_rating(bdat_cond), model_comps_binary(bdat_cond))

good <- bdat_rate[,cor(n_inculp-n_exculp, as.numeric(bardguilt)), by=uid][V1>0 & !is.na(V1), uid]
bdat_rate <- bdat_rate[uid %in% good]
t3_1 <- model_comps_binary(bdat_rate[cond_rating=="with"])
t3_2 <- model_comps_binary(bdat_rate[cond_rating=="without"])

# print table 1, experiment 1
t1[response=="Case strength"]

# print table 1, experiment 2
t2[response=="Case strength"]

# print supplementary table (in segments)
# experiment 1
t1[response=="Guilt judgment"]
# experiment 2
t2[response=="Guilt judgment"]
# experiment 3 (w/ rating)
t3_1
# experiment 3 (w/0 rating)
t3_2
