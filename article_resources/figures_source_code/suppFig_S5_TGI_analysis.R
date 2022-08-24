library(PaSM)


model  <-  RxODE({
  d/dt(Central) <- -ke *  max(Central,0)
  Conc <- Central/Vd

  tumVol <- X1 + X2 + X3 + X4
  psi <- 20
  growth <- lambda0 *  max(X1,0)/((1 + (lambda0 * max(tumVol,0)/lambda1)^psi)^(1/psi))

  X1(0) <- w0

  d/dt(X1) <- max(growth,0) - X1 * max(Conc,0) * k2
  d/dt(X2) <- X1 * max(Conc,1E-10) * k2 - k1 * X2
  d/dt(X3) <- k1 * (X2 - X3)
  d/dt(X4) <- k1 * (X3 - X4)

})






# parameters impact -------------------------------------------------------


defaults_values <- tibble(ke = 1, Vd = 5, lambda0 = 0.05, lambda1 = 10, w0 = 50, k2 = 1, k1 = 0.5)

map(parse_exprs(names(defaults_values)), function(x){

  defaults_values %>%
    bind_rows(defaults_values %>% mutate(!!enexpr(x) := !!enexpr(x) * 5)) %>%
    bind_rows(defaults_values %>% mutate(!!enexpr(x) := !!enexpr(x) / 5)) %>%
    mutate(param = deparse(x), values = c("x1", "x5", "x0.2"))

}) %>%
  bind_rows() %>%
  rowid_to_column("id") -> idsforsimuls



protocol <- tibble(time = 0, cmt = "Central", amt = 10, evid = 1) %>%
  bind_rows(crossing(time = 0:120, cmt = "Central", amt = 0, evid = 0)) %>%
  crossing(id = unique(idsforsimuls$id)) %>%
  arrange(id, time)


simulations <- model$solve(idsforsimuls %>% select(-param, -values), protocol) %>%
  as_tibble()

plot1 <- simulations %>%
  left_join(idsforsimuls %>%  select(id, param, values), by = "id") %>%
  left_join(defaults_values %>% gather("param", "value")) %>%
  mutate(param = paste0(param, " (", value, ")")) %>%
  ggplot()+
  labs(x = "Time", y = "Tumor Volume", col = "Value\n(Ref)")+
  geom_line(aes(time, tumVol, col = values) , size = 2 )+
  facet_wrap(~param, scales = "free")+
  scale_y_log10()+
  theme_bw(); plot1




# Exemple W0 chage --------------------------------------------------------



idsw0 <- tibble(k1 = c(0, 0), k2 = c(7.3, 7.3), ke = c(0.3, 0.3), lambda0 = c(0.64, 0.64), lambda1 = c(87.6, 87.6), psi = c(20, 20), Vd = c(184, 184), w0 = c(320, 500), id = 1:2)

idsw0 <- tibble(k1 = 0, k2 = 7.3, ke = 1, lambda0 = 0.64, lambda1 = 87.6, Vd = 184, w0 = c(320,500))

protocolw0 <- tibble(time = 0, cmt = "Central", amt = 50, evid = 1) %>%
  bind_rows(crossing(time = 0:60, cmt = "Central", amt = 0, evid = 0)) %>%
  crossing(id =1:2) %>%
  arrange(id, time)

plot2 <- model$solve(idsw0, protocolw0 ) %>%
  as_tibble() %>%
   mutate(X1_on_tumVol = X1/tumVol) %>%
   left_join(idsw0) %>%
   gather("cmt", "value",tumVol,X1_on_tumVol) %>%
   ggplot()+
  geom_line(aes(time, value, col = factor(w0)) , size = 2 )+
  scale_y_log10()+
  labs(x = "Time", y = "Tumor Volume (mm3)", col = "w0")+
   # geom_hline(data = tibble(x = 78, cmt = "tumVol"), aes(yintercept = x), lty = 2)+
   # geom_text(data = tibble(x = 78, cmt = "tumVol"), aes(x = 30, y = 60, label = "Exponential Growth"), lty = 2)+
   # geom_text(data = tibble(x = 78, cmt = "tumVol"), aes(x = 30, y = 200, label = "Linear Growth"), lty = 2)+
   # geom_hline(data = tibble(x = 78, cmt = "tumVol"), aes(yintercept = x), lty = 2)+
   facet_wrap(~cmt, scales = "free", ncol = 1)+
  theme_bw();plot2


# Exemple k1 change

idsk1 <- tibble(k1 = c(0.1,5), k2 = 2, ke = 1, lambda0 = 4, lambda1 = 200, Vd = 5, w0 = 50, id = 1:2)


plot3 <- model$solve(idsk1, protocolw0 ) %>%

  as_tibble() %>%
   filter(time < 20) %>%
   left_join(idsk1) %>%
   mutate(X1_on_tumVol = X1/tumVol) %>%
   left_join(idsw0) %>%
   gather("cmt", "value",tumVol,X1_on_tumVol) %>%
   ggplot()+
   geom_line(aes(time, value, col = factor(k1)) , size = 2 )+
   scale_y_log10()+
  labs(x = "Time", y = "Tumor Volume (mm3)", col = "k1")+
     facet_wrap(~cmt, scales = "free", ncol = 1)+
   theme_bw();plot3



# Final plot --------------------------------------------------------------

cowplot::plot_grid(plot1, plot2, plot3, nrow =1, labels = LETTERS, rel_widths = c(2,1,1))

# Optimizing number of patients per salve



library(microbenchmark)

allids <- crossing(ke = 1, Vd = 5, lambda0 = seq(0.05,0.15,0.01), lambda1 = 10:15, w0 = 40:60, k2 = seq(0.6,1.4,0.1), k1 = 0.5) %>%
  rowid_to_column("id")


protocols <- tibble(time = 0, cmt = "Central", amt = 50, evid = 1) %>%
  bind_rows(crossing(time = 0:60, cmt = "Central", amt = 0, evid = 0)) %>%
  crossing(id = unique(allids$id)) %>%
  arrange(id, time)

alltimes <- tibble(nsimul = c(1,10,100,250,500,1000,2000,5000,10000)) %>%
  mutate(time_dbl = map_dbl(nsimul, function(x){

    idtemps <- allids %>% filter(id %in% 1:x)
    protocoltemp <- protocols  %>% filter(id %in% 1:x)

    test <-  microbenchmark(model$solve(idtemps, protocoltemp ), times = 5, unit = "ms")

    # time in ms
    median(test$time / 1E6)

  }))

alltimes %>%
  mutate(time_per_id = time_dbl / nsimul ) %>%
  ggplot()+
  geom_point(aes(nsimul, time_per_id))+
  geom_line(aes(nsimul, time_per_id))+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()


