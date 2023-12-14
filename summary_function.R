
library(tidyverse)
library(beepr)

asd_over_list <- function(list){
  
  ## df 1
    df1_weight_sum  <-
    list[[1]][[1]]  %>%
    as.tibble() %>%
      reframe(
      asd_x1 = 
        abs(
          sum(x1*trt)/sum(trt) - sum(x1*(1-trt))/sum(1-trt)
        )/sqrt(
          (((x1*trt) - mean(x1*trt))^2 + ((x1*(1 - trt)) - mean(x1*(1 - trt)))^2)/2
        ),
      asd_x2 = 
        abs(
          sum(x2*trt)/sum(trt) - sum(x2*(1-trt))/sum(1-trt)
        )/sqrt(
          (((x2*trt) - mean(x1*trt))^2 + ((x2*(1 - trt)) - mean(x2*(1 - trt)))^2)/2
        ),
      asd_x3 = 
        abs(
          sum(x3*trt)/sum(trt) - sum(x3*(1-trt))/sum(1-trt)
        )/sqrt(
          (((x3*trt) - mean(x3*trt))^2 + ((x3*(1 - trt)) - mean(x3*(1 - trt)))^2)/2
        ),
      weights = "None"
    )
  
  ## df 2
  
  df2_weight_sum  <-
    list[[2]][[1]]  %>%
    as.tibble() %>%
    reframe(
      asd_x1 = 
        abs(
          sum(x1*trt)/sum(trt) - sum(x1*(1-trt))/sum(1-trt)
        )/sqrt(
          (((x1*trt) - mean(x1*trt))^2 + ((x1*(1 - trt)) - mean(x1*(1 - trt)))^2)/2
        ),
      asd_x2 = 
        abs(
          sum(x2*trt)/sum(trt) - sum(x2*(1-trt))/sum(1-trt)
        )/sqrt(
          (((x2*trt) - mean(x1*trt))^2 + ((x2*(1 - trt)) - mean(x2*(1 - trt)))^2)/2
        ),
      asd_x3 = 
        abs(
          sum(x3*trt)/sum(trt) - sum(x3*(1-trt))/sum(1-trt)
        )/sqrt(
          (((x3*trt) - mean(x3*trt))^2 + ((x3*(1 - trt)) - mean(x3*(1 - trt)))^2)/2
        ),
      weights = "IPTW"
    )
  
  ## df 3
  
  df3_weight_sum  <-
    list[[3]][[1]]  %>%
    as.tibble() %>%
    mutate(
      x1 = x1*cens_weights,
      x2 = x2*cens_weights,
      x3 = x3*cens_weights
    ) %>%
    reframe(
      asd_x1 = 
        abs(
          sum(x1*trt)/sum(trt*cens_weights) - sum(x1*(1-trt))/sum((1-trt)*cens_weights)
        )/sqrt(
          (((x1*trt) - mean(x1*trt))^2 + ((x1*(1 - trt)) - mean(x1*(1 - trt)))^2)/2
        ),
      asd_x2 = 
        abs(
          sum(x2*trt)/sum(trt*cens_weights) - sum(x2*(1-trt))/sum((1-trt)*cens_weights)
        )/sqrt(
          (((x2*trt) - mean(x1*trt))^2 + ((x2*(1 - trt)) - mean(x2*(1 - trt)))^2)/2
        ),
      asd_x3 = 
        abs(
          sum(x3*trt)/sum(trt*cens_weights) - sum(x3*(1-trt))/sum((1-trt)*cens_weights)
        )/sqrt(
          (((x3*trt) - mean(x3*trt))^2 + ((x3*(1 - trt)) - mean(x3*(1 - trt)))^2)/2
        ),
      weights = "IPCW * IPTW"
    )
  
  
  ## df 4
  
  df4_weight_sum  <-
    list[[4]][[1]]  %>%
    as.tibble() %>%
    mutate(
      x1 = x1*cens_weights,
      x2 = x2*cens_weights,
      x3 = x3*cens_weights
    ) %>%
    reframe(
      asd_x1 = 
        abs(
          sum(x1*trt)/sum(trt*cens_weights) - sum(x1*(1-trt))/sum((1-trt)*cens_weights)
        )/sqrt(
          (((x1*trt) - mean(x1*trt))^2 + ((x1*(1 - trt)) - mean(x1*(1 - trt)))^2)/2
        ),
      asd_x2 = 
        abs(
          sum(x2*trt)/sum(trt*cens_weights) - sum(x2*(1-trt))/sum((1-trt)*cens_weights)
        )/sqrt(
          (((x2*trt) - mean(x1*trt))^2 + ((x2*(1 - trt)) - mean(x2*(1 - trt)))^2)/2
        ),
      asd_x3 = 
        abs(
          sum(x3*trt)/sum(trt*cens_weights) - sum(x3*(1-trt))/sum((1-trt)*cens_weights)
        )/sqrt(
          (((x3*trt) - mean(x3*trt))^2 + ((x3*(1 - trt)) - mean(x3*(1 - trt)))^2)/2
        ),
      weights = "IPCW"
    )
  
  ## df 5
  
  df5_weight_sum  <-
    list[[5]][[1]]  %>%
    as.tibble() %>%
    mutate(
      x1 = x1*trt_weight,
      x2 = x2*trt_weight,
      x3 = x3*trt_weight
    ) %>%
    reframe(
      asd_x1 = 
        abs(
          sum(x1*trt)/sum(trt*trt_weight) - sum(x1*(1-trt))/sum((1-trt)*trt_weight)
        )/sqrt(
          (((x1*trt) - mean(x1*trt))^2 + ((x1*(1 - trt)) - mean(x1*(1 - trt)))^2)/2
        ),
      asd_x2 = 
        abs(
          sum(x2*trt)/sum(trt*trt_weight) - sum(x2*(1-trt))/sum((1-trt)*trt_weight)
        )/sqrt(
          (((x2*trt) - mean(x1*trt))^2 + ((x2*(1 - trt)) - mean(x2*(1 - trt)))^2)/2
        ),
      asd_x3 = 
        abs(
          sum(x3*trt)/sum(trt*trt_weight) - sum(x3*(1-trt))/sum((1-trt)*trt_weight)
        )/sqrt(
          (((x3*trt) - mean(x3*trt))^2 + ((x3*(1 - trt)) - mean(x3*(1 - trt)))^2)/2
        ),
      weights = "IPTW * IPCW"
    )
  
  return(bind_rows(df1_weight_sum, df2_weight_sum, df3_weight_sum, df4_weight_sum, df5_weight_sum))
}

library(tidyverse)

#all_asd_avgs <- tibble()
#all_sum_res <- tibble()
sample_size <- 100

load("scen1_sum_100.RData")
tot_obs <- length(scen1_sum[[1]]$data)
asd_scen1 <-
  tibble(df_list = scen1_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen1"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen1)
all_sum_res <- bind_rows(all_sum_res, scen1_sum[[2]])
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
rm(asd_scen1)
rm(scen1_sum)
beep()

load("scen2_sum_100.RData")
asd_scen2 <-
  tibble(df_list = scen2_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen2"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen2)
all_sum_res <- bind_rows(all_sum_res, scen2_sum[[2]])
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
rm(asd_scen2)
rm(scen2_sum)
beep()

load("scen3_sum_100.RData")
asd_scen3 <-
  tibble(df_list = scen3_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen3"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen3)
all_sum_res <- bind_rows(all_sum_res, scen3_sum[[2]])
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
rm(asd_scen3)
rm(scen3_sum)
beep()


load("scen4_sum_100.RData")
asd_scen4 <-
  tibble(df_list = scen4_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen4"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen4)
all_sum_res <- bind_rows(all_sum_res, scen4_sum[[2]])
rm(asd_scen4)
rm(scen4_sum)
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
beep()

load("scen5_sum_100.RData")
asd_scen5 <-
  tibble(df_list = scen5_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen5"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen5)
all_sum_res <- bind_rows(all_sum_res, scen5_sum[[2]])
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
rm(asd_scen5)
rm(scen5_sum)
beep()

load("scen6_sum_100.RData")
asd_scen6 <-
  tibble(df_list = scen6_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen6"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen6)
all_sum_res <- bind_rows(all_sum_res, scen6_sum[[2]])
rm(asd_scen6)
rm(scen6_sum)
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
beep()

load("scen7_sum_100.RData")
asd_scen7 <-
  tibble(df_list = scen7_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen7"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen7)
all_sum_res <- bind_rows(all_sum_res, scen7_sum[[2]])
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
rm(asd_scen7)
rm(scen7_sum)
beep()

load("scen8_sum_100.RData")
asd_scen8 <-
  tibble(df_list = scen8_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen8"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen8)
all_sum_res <- bind_rows(all_sum_res, scen8_sum[[2]])
rm(asd_scen8)
rm(scen8_sum)
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
beep()

load("scen9_sum_100.RData")
asd_scen9 <-
  tibble(df_list = scen9_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen9"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen9)
all_sum_res <- bind_rows(all_sum_res, scen9_sum[[2]])
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
rm(asd_scen9)
rm(scen9_sum)
beep()

load("scen10_sum_100.RData")
asd_scen10 <-
  tibble(df_list = scen10_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen10"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen10)
all_sum_res <- bind_rows(all_sum_res, scen10_sum[[2]])
rm(asd_scen10)
rm(scen10_sum)
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
beep()

load("scen11_sum_100.RData")
asd_scen11 <-
  tibble(df_list = scen11_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen11"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen11)
all_sum_res <- bind_rows(all_sum_res, scen11_sum[[2]])
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
rm(asd_scen11)
rm(scen11_sum)
beep()

load("scen12_sum_100.RData")
asd_scen12 <-
  tibble(df_list = scen12_sum[[1]]$data[sample(1:40000, sample_size, replace=FALSE)]) %>%
  mutate(
    asds = map(df_list, asd_over_list, .progress=TRUE)
  ) %>%
  select(-df_list) %>%
  unnest(asds) %>%
  group_by(weights) %>%
  summarise(
    avg_asd_x1 = mean(asd_x1, na.rm=TRUE),
    avg_asd_x2 = mean(asd_x2, na.rm=TRUE),
    avg_asd_x3 = mean(asd_x3, na.rm=TRUE),
    scenario = "scen12"
  )
all_asd_avgs <- bind_rows(all_asd_avgs, asd_scen12)
all_sum_res <- bind_rows(all_sum_res, scen12_sum[[2]])
rm(asd_scen12)
rm(scen12_sum)
save(all_asd_avgs, file="all_asd_avgs.RData")
save(all_sum_res, file="all_sum_res.RData")
beep()

beep(8)



wider_asd_sum <- all_asd_avgs %>%
  group_by(scenario) %>%
  rename(
    x1 = avg_asd_x1,
    x2 = avg_asd_x2,
    x3 = avg_asd_x3
  ) %>%
  pivot_wider(
      names_from = weights,
      values_from = c(x1, x2, x3)
      ) %>%
  janitor::clean_names() %>%
  relocate(
    scenario,
    x1_none,
    x1_iptw,
    x1_iptw_ipcw,
    x1_ipcw,
    x1_ipcw_iptw,
    x2_none,
    x2_iptw,
    x2_iptw_ipcw,
    x2_ipcw,
    x2_ipcw_iptw,
    x3_none,
    x3_iptw,
    x3_iptw_ipcw,
    x3_ipcw,
    x3_ipcw_iptw
  )

wider_res_sum <- 
  all_sum_res %>%
  group_by(scenario) %>%
  pivot_wider(
    names_from = est_name,
    values_from = c(avg_te, avg_se, relative_bias, coverage_prop)
  ) %>%
  janitor::clean_names() %>%
  relocate(
    scenario,
    avg_te_trt_adj_km,
    avg_se_trt_adj_km,
    relative_bias_trt_adj_km,
    coverage_prop_trt_adj_km,
    avg_te_trt_cens_adj_km,
    avg_se_trt_cens_adj_km,
    relative_bias_trt_cens_adj_km,
    coverage_prop_trt_cens_adj_km,
    avg_te_cens_adj_km,
    avg_se_cens_adj_km,
    relative_bias_cens_adj_km,
    coverage_prop_cens_adj_km,
    avg_te_cens_trt_adj_km,
    avg_se_cens_trt_adj_km,
    relative_bias_cens_trt_adj_km,
    coverage_prop_cens_trt_adj_km
  )

write_csv(file = "fin_asd_sum.csv", wider_asd_sum)
write_csv(file = "fin_res_sum.csv", wider_res_sum)
















