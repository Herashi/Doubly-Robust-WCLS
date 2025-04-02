out = omit[["all_data"]]
out_wcls = out[,c("sec")]
out_rwcls = out[,c("se_ML")]
out_drwcls = out[,c("se_DR")]

out_wcls_oracle = out[,c("sec_oracle")]
out_rwcls_oracle = out[,c("se_ML_oracle")]
out_drwcls_oracle = out[,c("se_DR_oracle")]


# omit[["out"]]

# compare WCLS vs R-WCLS & DR-WCLS

df = out_wcls ^2 / (out_rwcls ^2)
sum(df>=1)
summary(df)

(omit[["out"]]$sdc/omit[["out"]]$sd_ML)^2

#### 
df = out_wcls ^2 / (out_drwcls ^2)
sum(df>=1)
summary(df)

(omit[["out"]]$sdc/omit[["out"]]$sd_DR)^2

#####
df = out_wcls_oracle ^2 / (out_rwcls_oracle ^2)
# sum(df>=1)
summary(df)

# (omit[["out"]]$sdc_oracle/omit[["out"]]$sd_ML_oracle)^2

######
df = out_wcls_oracle ^2 / (out_drwcls_oracle ^2)
# sum(df>=1)
summary(df)

# (omit[["out"]]$sdc_oracle/omit[["out"]]$sd_DR_oracle)^2



# df = out_drwcls[,2]^2 / (out_rwcls[,2]^2)
# sum(df>=1)
# summary(df)
# hist(df)
# 
# df = out_rwcls[,2]^2 / (out_drwcls[,2]^2)
# sum(df>=1)
# summary(df)
# hist(df)


# omit[["out"]][["sec_oracle"]]^2/omit[["out"]][["se_ML_oracle"]]^2
# omit[["out"]][["sec_oracle"]]^2/omit[["out"]][["se_DR_oracle"]]^2
# 
# omit[["out"]][["sec"]]^2/omit[["out"]][["se_ML"]]^2
# omit[["out"]][["sec"]]^2/omit[["out"]][["se_DR"]]^2
