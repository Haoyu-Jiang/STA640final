# cnames <- c("Sample Size", "Bias", "SE_ACM", "SD", "SE_ACM/SD",
#                        "SE_ACM Based", "95% Lower", "95% Upper", "SE_SD Based",
#                        "95% Lower", "95% Upper", "Percentile Based", "95% Lower", "95% Upper")
# 
# colnames(res1.df) <- colnames(res2.df) <- colnames(res3.df) <- cnames
# 
# res_all <- rbind(res1.df, res2.df, res3.df)
# 
# 
# res_all[, 2:14] <- round(res_all[, 2:14], 3)
# res_all[, 6:14] <- 100 * res_all[, 6:14]
# View(res_all)

save(res_all, file = "res_all.RData")
