

# Flat vs reboot ----------------------------------------------------------

model <- model_RxODE <-  RxODE({
  d/dt(Central) <- -ke *  max(Central,0)
  Conc <- Central/Vd

  tumVol <- X1 + X2 + X3 + X4
  growth <- lambda0 *  max(X1,0)/((1 + (lambda0 * max(tumVol,0)/lambda1)^psi)^(1/psi))

  X1(0) <- w0

  d/dt(X1) <- growth - X1 * max(Conc,0) * k2
  d/dt(X2) <- X1 * max(Conc,0) * k2 - k1 * X2
  d/dt(X3) <- k1 * (X2 - X3)
  d/dt(X4) <- k1 * (X3 - X4)

})




proto <- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,0.01), amt = 0, evid = 0)) # change the step, (0.01,0.1,1), you will see !

id <- tibble(ke = 1, w0 = 50, k2 = 3.4, lambda0 = 0.859, Vd = 4, lambda1 = 63.4, k1 = 0.5, psi = 20)

id2 <-  tibble(ke = 1, w0 = 50, k2 = 2.93, lambda0 = 1, Vd = 4, lambda1 = 86.1, k1 = 0.5, psi = 20)


res <- model$solve(id2, proto, c(X1 = 50)) %>%   as.tibble() %>% mutate(supposed = "above") %>%
  bind_rows(model$solve(id, proto, c(X1 = 50))%>% as.tibble() %>%  mutate(supposed = "below"))


res %>%
  ggplot()+
  geom_line(aes(time, tumVol, col = supposed))+
  scale_y_log10()


# Can we avoid that? ------------------------------------------------------


model <- model_RxODE <-  RxODE({
  d/dt(Central) <- -ke *  max(Central,0)
  Conc <- Central/Vd

  tumVol <- X1 + X2 + X3 + X4
  growth <- lambda0 *  max(X1,0)/((1 + (lambda0 * max(tumVol,0)/lambda1)^psi)^(1/psi))

  X1(0) <- w0

  d/dt(X1) <- growth - X1 * max(Conc,0) * k2
  d/dt(X2) <- X1 * max(Conc,0) * k2 - k1 * X2
  d/dt(X3) <- k1 * (X2 - X3)
  d/dt(X4) <- k1 * (X3 - X4)

})


proto <- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,1), amt = 0, evid = 0)) # change the step, (0.01,0.1,1), you will see !


proto2 <- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,0.01), amt = 0, evid = 0)) # change the step, (0.01,0.1,1), you will see !
#
# id <- tibble(ke = 1, w0 = 50, k2 = 3.4, lambda0 = 0.859, Vd = 4, lambda1 = 63.4, k1 = 0.5, psi = 20)
#
# id2 <-  tibble(ke = 1, w0 = 50, k2 = 2.93, lambda0 = 1, Vd = 4, lambda1 = 86.1, k1 = 0.5, psi = 20)


res2 <- model$solve(id2, proto, c(X1 = 50)) %>%   as.tibble() %>% mutate(supposed = "step1") %>%
  bind_rows(model$solve(id2, proto2, c(X1 = 50))%>% as.tibble() %>%  mutate(supposed = "step0.01"))

res2 %>%
  filter(supposed == "step1") %>%


res2 %>%
  ggplot()+
  geom_line(aes(time, tumVol, col = supposed))+
  scale_y_log10()


res2 %>%
  select(time, tumVol, supposed) %>%
  spread(supposed, tumVol) %>%
  filter(!is.na(step1)) %>%
  filter(round(step0.01,1) != round(step1,1))

res2 %>%
  filter(time%in%(res2$time[res2$supposed == "step1" ])) %>%
  filter(time >25) %>%
  arrange(time)

res2 %>%
  filter(supposed == "1")

# Can we avoid that? (same modif model)------------------------------------------------------


model <- model_RxODE <-  RxODE({
  d/dt(Central) <- -ke *  max(Central,0)
  Conc <- Central/Vd

  tumVol <- X1 + X2 + X3 + X4
  growth <- lambda0 *  max(X1,0)/((1 + (lambda0 * max(tumVol,0)/lambda1)^psi)^(1/psi))

  X1(0) <- w0

  d/dt(X1) <- max(growth,0) - X1 * max(Conc,1E-10) * k2
  d/dt(X2) <- X1 * max(Conc,1E-10) * k2 - k1 * X2
  d/dt(X3) <- k1 * (X2 - X3)
  d/dt(X4) <- k1 * (X3 - X4)

})


proto <- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,1), amt = 0, evid = 0)) # change the step, (0.01,0.1,1), you will see !


proto2 <- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,0.01), amt = 0, evid = 0)) # change the step, (0.01,0.1,1), you will see !
#
# id <- tibble(ke = 1, w0 = 50, k2 = 3.4, lambda0 = 0.859, Vd = 4, lambda1 = 63.4, k1 = 0.5, psi = 20)
#
# id2 <-  tibble(ke = 1, w0 = 50, k2 = 2.93, lambda0 = 1, Vd = 4, lambda1 = 86.1, k1 = 0.5, psi = 20)


res2 <- model$solve(id2, proto, c(X1 = 50)) %>%   as.tibble() %>% mutate(supposed = "step1") %>%
  bind_rows(model$solve(id2, proto2, c(X1 = 50))%>% as.tibble() %>%  mutate(supposed = "step0.01"))



  res2 %>%
  ggplot()+
  geom_line(aes(time, tumVol, col = supposed))+
  scale_y_log10()


res2 %>%
  select(time, tumVol, supposed) %>%
  spread(supposed, tumVol) %>%
  filter(!is.na(step1)) %>%
  filter(round(step0.01,1) != round(step1,1))

res2 %>%
  filter(time%in%(res2$time[res2$supposed == "step1" ])) %>%
  filter(time >25) %>%
  arrange(time)

res2 %>%
  filter(supposed == "step1")
