## ---------------------------
## Script name: Figure R article
##
## Purpose of script: production of figure 2
##c
## Author: Thibaud Derippe
##
## Date Created: 2022-04-04 (day of code repartition from R6object.R)
##
## Under GPL-3 License
## Email: thibaud.derippe@gmail.com
## GitHub: https://github.com/Thibaudpmx/QSPVP
## ---------------------------
##
## Notes: Stochasticity in the algorithm, I put a seed
## but results might vary on another machine
##
##
## ---------------------------


# Figure 2 ----------------------------------------------------------------

# Use the same target as in Figure 1



## One increase one decrease (k2 and lambda0)

self <- VP_proj_creator$new()


self$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,45))


param <- crossing(k1 = 0.5, k2 = seq(0,6,0.2), ke = 1, lambda0 = seq(0,0.7,0.02), lambda1 = 12, Vd = 40, psi = 20, w0 = 50)
param$k2 <- round(param$k2, 1)
# bind_rows(crossing(k1 = 0.5, k2 = c(0.8,0.9), ke = 1, lambda0 = c(0.05,0.03), lambda1 = 12, Vd = 40, psi = 20)) %>%

# param <- crossing(k1 = 0.5, k2 = seq(2,6,2), ke = 1, lambda0 = seq(0.02,0.06,0.01), lambda1 = 12, Vd = 40, psi = 20) %>%
#   bind_rows(crossing(k1 = 0.5, k2 = c(0.8,0.9), ke = 1, lambda0 = c(0.05,0.03), lambda1 = 12, Vd = 40, psi = 20)) %>%
#   rowid_to_column("id")

# param %>%
# ggplot()+
#   geom_point(aes(x = k2, y = lambda0))

self$add_VP(param, use_green_filter = T)

self$n_filter_reduc()
self$plot_2D(k2, lambda0, plotMain = T, plotoreturn = 3)
self$fill_simul()
self$plot_VP()

set.seed(1653)
param20 <-   param %>%
  sample_n(20) %>%
  filter(!(k2 == 0.2 & lambda0 == 0.46)) %>%
  bind_rows(

    param %>% filter(k2 == 3.6 & lambda0 == 0.44)

  )

plot1 <-  param20 %>%
  ggplot()+
  geom_point(aes(x = k2, y = lambda0, shape = "VP\nsampled\namong\ncohort"))+
  labs(shape = "")+
  theme_bw();plot1

ids <-   param20 %>%
  rowid_to_column("id")

events <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows(tibble(time = self$times, evid = 0, amt = 0, cmt = "Conc")) %>%
  crossing(id = 1:20) %>%
  arrange(id, time)

simulations <- self$model$solve(ids, events, c(X1 = 50)) %>% as_tibble()

analysis <- simulations %>%
  filter(time %in% self$targets$time) %>%
  mutate(cmt = "tumVol") %>%
  left_join( self$targets) %>%
  mutate(below = tumVol > min, above = tumVol < max) %>%
  mutate(type = case_when(below & above ~ "Accepted",
                          below & !above ~ "Above",
                          !below & above ~ "Below"))

analysis <- analysis %>% distinct(id, type) %>%
  group_by(id) %>%
  nest() %>%
  mutate(label = map_chr(data, function(x){

    if(nrow(x) == 2 & "Accepted" %in% x$type ) x <- x %>% filter(x != "Accepted")

    paste0(x$type, collapse = "&")
  } ))


plot2 <- simulations %>%
  left_join(analysis %>%  select(id, label)) %>%
  ggplot()+
  geom_line(aes(time, tumVol, group = id))+
  geom_segment(data = self$targets, aes(x = time, xend = time, y = min, yend = max), col = "red", size = 2)+
  scale_y_log10()+
  facet_wrap(~label) +
  labs(x = "Time (days)", y = "Tumor volume (mm3)")+
  theme_bw(); plot2


param20 <- param20 %>%
  rowid_to_column("id") %>%
  left_join(analysis %>%  select(id, label))

forsquare <- param20 %>%
  # filter(label) %>%
  mutate(k2min = if_else(label == "Above", 0, k2),
         k2max = if_else(label == "Above", k2, Inf),
         lambda0min = if_else(label == "Below", 0, lambda0),
         lambda0max = if_else(label == "Below", lambda0, Inf))

plot3 <-
  param20 %>%
  ggplot()+
  geom_point(aes(x = k2, y = lambda0, col = label))+
  geom_rect(data = forsquare %>% filter(label != "Accepted"),
            aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = label ),
            alpha = 0.2, col = "black")+
  theme_bw()+
  scale_color_manual(values = c("red", "darkgreen", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate"))+
  labs(col = "VPs", fill = "Filters"); plot3

filter_reduc(forsquare %>% filter(label == "Above"), filtre = self$make_filters()[1]) %>%
  bind_rows(

    filter_reduc(forsquare %>% filter(label == "Below"), filtre = self$make_filters()[2])

  ) ->filtersreduc

plot4 <-
  param20 %>%
  ggplot()+
  geom_point(aes(x = k2, y = lambda0, col = label))+
  geom_rect(data = filtersreduc %>% filter(label != "Accepted"),
            aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = label ),
            alpha = 0.2, col = "black")+
  theme_bw()+
  scale_color_manual(values = c("red", "darkgreen", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate"))+
  labs(col = "VPs", fill = "Filters"); plot4

filtersreduc %>%
  filter(k2max < Inf) %>%
  mutate(filter = paste0("k2 <=", k2, " & lambda0 >= ", lambda0 )) %>%
  pull(filter) -> pointsabove

pointsabove <- paste0("(", pointsabove, ")") %>% paste0(collapse = "|")

filtersreduc %>%
  filter(k2max == Inf) %>%
  mutate(filter = paste0("k2 >=", k2, " & lambda0 <= ", lambda0 )) %>%
  pull(filter)  -> pointsbelow

pointsbelow <- paste0("(", pointsbelow, ")") %>% paste0(collapse = "|")


plot5 <-
  param20 %>%
  ggplot()+
  geom_point(data = param %>% filter(!!parse_expr(pointsabove)), aes(k2, lambda0,  col = "Above"))+
  geom_point(data = param %>% filter(!!parse_expr(pointsbelow)), aes(k2, lambda0,  col = "Below"))+
  geom_point(aes(x = k2, y = lambda0, col = label))+
  geom_rect(data = filtersreduc %>% filter(label != "Accepted"),
            aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = label ),
            alpha = 0.2)+
  theme_bw()+
  scale_color_manual(values = c("red", "darkgreen", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate"))+
  labs(col = "VPs", fill = "Filters"); plot5




forsquaregreen <-

  param20 %>%
  filter(label == "Accepted") %>%
  mutate(label = "Above") %>%
  bind_rows(
    param20 %>%
      filter(label == "Accepted") %>%
      mutate(label = "Below")
  ) %>%
  # filter(label) %>%
  mutate(k2min = if_else(label == "Above", 0, k2),
         k2max = if_else(label == "Above", k2, Inf),
         lambda0min = if_else(label == "Below", 0, lambda0),
         lambda0max = if_else(label == "Below", lambda0, Inf)) %>%
  mutate(label2 = if_else(label == "Above", "Above\nLower Limit", "Below\nUpper Limit"))


forsquaregreen

plot6 <- plot5 +
  geom_rect(data = forsquaregreen %>% slice(-5, -8), aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Green"), alpha = 0.2)+
  geom_rect(data = forsquaregreen %>% slice(10), aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Green"), alpha = 0, col = "darkgreen",  lty= 3)+
  geom_rect(data = forsquaregreen %>% slice(3), aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Green"), alpha = 0, col = "darkgreen",  lty= 2)+
  # forsquaregreen %>% slice(-5, -8)+
  scale_fill_manual(values = c("red", "chocolate", "darkgreen"))+
  # geom_rect(aes(xmin = 3.6, xmax = 3.8, ymin = 0.4, ymax = 0.44, fill = "Green"); )+
  geom_segment(aes(x = 3.6, xend = 3.6, y= 0.44, yend = 0.4 ), col = "darkgreen")+
  geom_segment(aes(x = 3.8 , xend = 3.8 , y= 0.44, yend = 0.4 ), col = "darkgreen")+
  geom_segment(aes(x = 3.6, xend = 3.8 , y= 0.44, yend = 0.44 ), col = "darkgreen")+
  geom_segment(aes(x = 3.6, xend = 3.8, y= 0.4, yend = 0.4 ), col = "darkgreen")+
  geom_point(data = param %>% filter( k2>= 3.6 & k2<=3.8 & lambda0>= 0.4 &lambda0<=0.44), aes(k2, lambda0), col = "darkgreen", alpha = 0.2); plot6

param %>%
  left_join( param %>% filter(!!parse_expr(pointsabove)) %>% mutate(nop = T)) %>%
  filter(is.na(nop)) %>%
  select(-nop) %>%
  left_join( param %>% filter(!!parse_expr(pointsbelow)) %>% mutate(nop = T)) %>%
  filter(is.na(nop)) %>%
  select(-nop) %>%
  filter(!( k2>= 3.6 & k2<=3.8 & lambda0>= 0.4 &lambda0<=0.44)) %>%
  sample_n(20) -> news

plot7 <- param20 %>%
  ggplot()+
  # geom_point(aes(x = k2, y = lambda0, col = label))+
  geom_point(data = param %>% filter(!!parse_expr(pointsabove)), aes(k2, lambda0,  col = "Extrap. Above"))+
  geom_point(data = param %>% filter(!!parse_expr(pointsbelow)), aes(k2, lambda0,  col = "Extrap. Below"))+
  geom_point(data = param %>% filter( k2>= 3.6 & k2<=3.8 & lambda0>= 0.4 &lambda0<=0.44), aes(k2, lambda0, col = "Extrap. Accepted"))+
  geom_point(data = param20, aes(x = k2, y = lambda0, col = "From RxODE"))+
  scale_color_manual(values = c("red", "darkgreen", "chocolate", "black"))+
  labs(col = "   VPs")+
  theme_bw() ; plot7



plot8 <- param20 %>%
  ggplot()+
  geom_point(aes(x = k2, y = lambda0, col = label))+
  geom_point(data = param %>% filter(!!parse_expr(pointsabove)) , aes(k2, lambda0,  col = "Above"))+
  geom_point(data = param %>% filter(!!parse_expr(pointsbelow)), aes(k2, lambda0,  col = "Below"))+
  geom_point(data = param %>% filter( k2>= 3.6 & k2<=3.8 & lambda0>= 0.4 &lambda0<=0.44), aes(k2, lambda0, col = "Accepted"))+
  # geom_point(data = param20, aes(x = k2, y = lambda0))+
  scale_color_manual(values = c("red", "darkgreen", "chocolate"))+
  labs(col = "VPs")+
  theme_bw()+
  geom_point(data = news, aes(k2, lambda0))

# plot7 <-
#   plot5 +
#   geom_point(data = param %>% filter( k2>= 3.6 & k2<=3.8 & lambda0>= 0.4 &lambda0<=0.44), aes(k2, lambda0), col = "darkgreen")
#              aes(xmin = 3.6, xmax = 3.8, ymin = 0.4, ymax = 0.44, fill = "Green"), alpha = 0.05)+
#     scale_fill_manual(values = c("red", "chocolate", "darkgreen")))
#   geom_
# filtersreduc


self$n_filter_reduc()


plot9 <- ggplot()+
  geom_point(data = self$filters_neg_above, aes(k2, lambda0, col = "Above"))+
  geom_rect(data = self$filters_neg_above,aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above" ), alpha = 0.1)+
  geom_point(data = self$filters_neg_below, aes(k2, lambda0, col = "Below"))+
  geom_point(data = self$poolVP, aes(k2, lambda0, col = "Accepted"))+
  geom_rect(data = self$filters_neg_below,aes(xmin = k2 , xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below" ), alpha = 0.1)+

  theme_bw()+
  scale_color_manual(values = c("red", "darkgreen", "chocolate"))+
  scale_fill_manual(values = c("red",  "chocolate"))+
  labs(col = "VPs", fill = "Filters"); plot9


plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,plot7, plot8,plot9, nrow = 3, labels = LETTERS)


