library(RxODE)

# Note: Maybe this plot is not reproducible on every machine, as its purpose is to show a computational error
# that might depends on computing power (?)

# Flat vs reboot ----------------------------------------------------------

model <- model_RxODE <-  RxODE({
  d/dt(Central) <- -ke *  max(Central,0)
  Conc <- Central/Vd

  tumVol <- X1 + X2 + X3 + X4
  growth <- lambda0 *  max(X1,0)/((1 + (lambda0 * tumVol/lambda1)^psi)^(1/psi))

  X1(0) <- w0

  d/dt(X1) <- growth - X1 * Conc * k2
  d/dt(X2) <- X1 * Conc * k2 - k1 * X2
  d/dt(X3) <- k1 * (X2 - X3)
  d/dt(X4) <- k1 * (X3 - X4)

})




proto <- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,1), amt = 0, evid = 0)) # change the step, (0.01,0.1,1), you will see !


proto2 <- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,0.1), amt = 0, evid = 0)) # change the step, (0.01,0.1,1), you will see !

proto3<- tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
  bind_rows(crossing(cmt = "Central", time = seq(0,52,0.01), amt = 0, evid = 0))

below <- tibble(ke = 1, w0 = 50, k2 = 3.4, lambda0 = 0.859, Vd = 4, lambda1 = 63.4, k1 = 0.5, psi = 20) #below

above <-  tibble(ke = 1, w0 = 50, k2 = 2.93, lambda0 = 1, Vd = 4, lambda1 = 86.1, k1 = 0.5, psi = 20) #above


res <- model$solve(above, proto, c(X1 = 50)) %>%   as.tibble() %>% mutate( step = "1", pos = "above") %>%
  bind_rows(model$solve(below, proto, c(X1 = 50))%>% as.tibble() %>%  mutate( step = "1", pos = "below")) %>%
bind_rows(model$solve(below, proto2, c(X1 = 50))%>% as.tibble() %>%  mutate( step = "0.1", pos = "below")) %>%
bind_rows(model$solve(above, proto2, c(X1 = 50))%>% as.tibble() %>%  mutate( step = "0.1", pos = "above")) %>%
bind_rows(model$solve(below, proto3, c(X1 = 50))%>% as.tibble() %>%  mutate( step = "0.01", pos = "below")) %>%
  bind_rows(model$solve(above, proto3, c(X1 = 50))%>% as.tibble() %>%  mutate( step = "0.01", pos = "above")) %>%
  mutate(pos = if_else(pos == "above", "\nAbove\nk2 = 2.93\nlbd0 = 1\nlbd1 = 86.1", "\nBelow\nk2 = 3.4\nlbd0 = 0.85\nlbd1 = 63.4\n"))


res %>%
  filter(step != "0.1") %>%
   ggplot()+
  geom_segment(aes(x = 40, xend = 40, y = 10, yend = 1000), size = 2)+
  geom_line(aes(time, tumVol, col = pos), size = 2)+
  geom_line(data = res %>% filter(step == "0.01") %>% select(-step), aes(time, tumVol, col = pos),alpha = 1, lty = 3)+
  scale_y_log10()+
  facet_wrap(~paste0("step = ", step))+
  theme_bw()+
  labs(x = "Time (days)", y = "Tumor Volume (mm3)", col = "Theoretical\nPosition")

