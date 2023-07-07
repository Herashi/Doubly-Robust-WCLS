## specify group structure
group_all = list()
group_all[["100"]] = list()


# 625 (euqal group size, different variances)
group_all[["100"]][["n"]] = 100
group_all[["100"]][["group_id"]] = rep(1:5,each=20)
group_all[["100"]][["baseline sigma2"]] = rep(0,5)
group_all[["100"]][["bg sigma2"]] = rep(0,5)


 


