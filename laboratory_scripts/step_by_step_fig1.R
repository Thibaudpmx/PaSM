
self <- VP_proj_creator$new()


self$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,45))


self$targets$max <- c(62,200)

param <- crossing(k1 = 0.5, w0 = 50, k2 = seq(2,6,2), ke = 1, lambda0 = seq(0.02,0.06,0.01), lambda1 = 12, Vd = 40, psi = 20) %>%
  rowid_to_column("id")

param %>% arrange(k2, desc(lambda0)) %>% slice(1) %>% pull(id) -> idlowest

events <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows( self$protocols$dose50 %>% mutate(evid = 0) %>% select(-time) %>% crossing(time = 0:50)) %>%
  arrange(time) %>%
  crossing(id = 1:nrow(param))

simul <- self$model$solve(param, events, c(X1 = 50) ) %>%
  as_tibble() %>%
  left_join(param) %>%
  mutate(l0 = if_else(lambda0 == 0.06 , " 0.06", "< 0.06"))


# same step by step
self$data %>%
  filter(cmt %in%  unique( self$targets$cmt) & protocol %in% unique( self$targets$protocol)) %>%
  ggplot()+
  geom_segment(data =  self$targets,
               aes(x = time, xend = time, y = min, yend = max), col ="black", size = 2)+
  scale_y_log10()+
  theme_bw()+

  labs(x = "Time (days)", y = "Tumor Volume (mm3)" , col = "k2", alpha = "lambda0", size = "lambda0") +
  #Patient inside targets
  geom_line(data = simul3, aes(time, tumVol), col = "darkgrey")+
  geom_segment(aes(x =  12, xend = 25, y = 40, yend = 120), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))+
  geom_segment(aes(x = 25, xend = 45, y = 120, yend = 100), arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"))+
  geom_label(aes(x = 25, y = 120, label = "Targets"))+
  # just the ref, maybe to remove after
  # geom_line(data=simul %>% filter(l0 == " 0.06" & k2 ==2), aes(time, tumVol, group = id), col = "red")+
  # # geom_line(data=simul %>% filter(l0 == " 0.06"), aes(time, tumVol, group = id, col = factor(k2)), size = 2)+
  geom_line(data=simul, aes(time, tumVol, group = id, size = l0, alpha = l0, col = factor(k2)))+
  scale_alpha_manual(values = c(1,0.3))+
  scale_size_manual(values = c(2,1))

