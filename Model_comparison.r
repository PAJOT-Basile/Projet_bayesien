
DIC1 <- M1$BUGSoutput$DIC
DIC2 <- M2$BUGSoutput$DIC
DIC3 <- M3$BUGSoutput$DIC
DIC4 <- M4$BUGSoutput$DIC
DIC5 <- M5$BUGSoutput$DIC
DIC6 <- M6$BUGSoutput$DIC
DIC7 <- M7$BUGSoutput$DIC
DIC8 <- M8$BUGSoutput$DIC
DIC9 <- M9$BUGSoutput$DIC

DICs <- cbind(DIC1, DIC2,
              DIC3, DIC4,
              DIC5, DIC6,
              DIC7, DIC8,
              DIC9) %>% 
  as.data.frame()


Best_model <- M6
Best_model

res <- Best_model$BUGSoutput$sims.matrix %>% 
  as.data.frame()
exp(mean(res$mu.0)) * mean(gopher$Area)

exp(mean(res$mu.0) + mean(res$b.prev)*mean(gopher$prev)) * mean(gopher$Area)

mean(res$sd.s)


res_dis <- M9$BUGSoutput$sims.matrix
