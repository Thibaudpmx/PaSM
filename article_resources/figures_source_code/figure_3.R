## ---------------------------
## Script name: Figure R article
##
## Purpose of script: production of figure 3
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

library(ggforce)


# prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100.184, 431.005),
# max = c(100.185, 431.006))

prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100, 431),
                    max = c(100.5, 431.5))

self <- VP_proj_creator$new()

self$set_targets(manual = prototiny)
self$targets <- prototiny




k2_1 <- seq(0,3,0.4)
lambda0_1 <- seq(0,1.4,0.2)
points <- crossing(k2_1, lambda0_1)

plot1 <- ggplot()+
  geom_point(data = points, aes(k2_1, lambda0_1))+
  geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 2)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1), lty = 2)+
  theme_bw()+
  theme(line = element_blank())+
  scale_x_continuous(breaks = k2_1)+
  scale_y_continuous(breaks = lambda0_1)+
  labs(x = "K2", y = "lambda0"); plot1


VP_df <- crossing(k1 = c(0.5),
                  k2 = k2_1,
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =lambda0_1,
                  lambda1 = c(12),
                  w0 = 50,
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } ) %>% mutate(psi = 20)


ids <-   VP_df %>%
  rowid_to_column("id")

events <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows(tibble(time = self$times, evid = 0, amt = 0, cmt = "Conc")) %>%
  crossing(id = 1:nrow(VP_df)) %>%
  arrange(id, time)

simulations <- self$model$solve(ids, events, c(X1 = 50)) %>% as_tibble()

simulations %>%
  filter(time %in% prototiny$time) %>%
  left_join(prototiny) %>%
  mutate(test = case_when(tumVol > max  ~ "Above",
                          tumVol < max ~ "Below",
                          T ~ "yes")) %>%
  distinct(id, test) -> idoutput


self$add_VP(VP_df)
self$compute_zone_maybe()
maybe <- self$zone_maybe

self$plot_2D(lambda0, k2, plotMain = F)

plot2 <-
  mtcars %>%
  ggplot()+
  geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 2, alpha  = 0.3)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1), lty = 2, alpha  = 0.3)+
  # geom_point(data = ids %>% left_join(idoutput), aes(k2, lambda0))+

  geom_rect(data= self$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "To explore" ), alpha = 0.6,  col = "black")+
  scale_x_continuous(breaks = k2_1)+
  scale_y_continuous(breaks = lambda0_1)+
  geom_text(data = maybe %>% rowid_to_column("id"), aes(x = (k2min + k2max) /2, y = (lambda0min + lambda0max)/2, label = id) )+
  geom_segment(aes(x = 0, y = 0, xend = Inf, yend = 0, col = "Below"))+
  geom_point(data= self$filters_neg_below, aes(k2, lambda0, col = "Below"))+
  geom_point(data= self$filters_neg_above, aes(k2, lambda0, col = "Above"))+
  scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey"))+
  labs(col = "", fill = "", x = "K2", y = "lambda0")+
  theme_bw()+
  guides(col = F)+
  theme(line = element_blank()); plot2

plot_grid(plot1, plot2)

digit <- 3

maybe %>% slice(1) -> x


maybe %>%
  rowid_to_column("id") %>%
  group_split(id) %>%
  map(function(x){

    difx <-  x$k2max -   x$k2min
    dify <-  x$lambda0max -   x$lambda0min

    crossing(k2 = seq(x$k2min, x$k2max, difx/6) %>% round(3),
             lambda0 = seq(x$lambda0min, x$lambda0max, dify/6) %>% round(3)) %>%
      mutate(bloc = x$id)

  }) %>%
  bind_rows() -> newVPs

newVPs %>%
  left_join(

    newVPs %>%
      group_by(bloc) %>%
      summarise(k2min = min(k2), k2max = max(k2), lambda0min = min(lambda0), lambda0max = max(lambda0))

  ) -> forseg


plot3 <- ggplot()+
  geom_segment(data = forseg, aes(x = k2, xend = k2, y = lambda0, yend = lambda0max), lty = 3)+
  geom_segment(data = forseg, aes(x = k2, xend = k2max, y = lambda0, yend = lambda0), lty = 3)+
  # geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 2, alpha  = 0.3)+
  # geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1), lty = 2, alpha  = 0.3)+
  # geom_point(data = ids %>% left_join(idoutput), aes(k2, lambda0))+

  geom_rect(data= self$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "To explore"), alpha = 0.6,  col = "black")+
  # scale_x_continuous(breaks = k2_1)+
  # scale_y_continuous(breaks = lambda0_1)+
  # geom_segment(aes(x = 0, y = 0, xend = Inf, yend = 0, col = "Below"))+
  # geom_point(data= self$filters_neg_below, aes(k2, lambda0, col = "Below"))+
  # geom_point(data= self$filters_neg_above, aes(k2, lambda0, col = "Above"))+
  scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey"))+
  labs(col = "", fill = "", x = "K2", y = "lambda0")+
  theme_bw()+
  guides(col = F)+

  theme(line = element_blank())+
  geom_point(data = newVPs, aes(k2, lambda0)); plot3
# geom_label(data = maybe %>% rowid_to_column("id"), aes(x = (k2min + k2max) /2, y = (lambda0min + lambda0max)/2, label = id), fill = "grey" )



self2 <- self$clone(deep = T)

newVPs %>%
  mutate(forjoin = 1) %>%
  left_join(
    VP_df %>% select(-k2, -lambda0) %>% mutate(forjoin = 1) %>% distinct()
  ) -> VP_df2


self2$add_VP(VP_df2)
self2$n_filter_reduc()

# self2$add_VP()


ggplot()+
  # geom_segment(data = forseg, aes(x = k2, xend = k2, y = lambda0, yend = lambda0max), lty = 3)+
  # geom_segment(data = forseg, aes(x = k2, xend = k2max, y = lambda0, yend = lambda0), lty = 3)+
  # geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 2, alpha  = 0.3)+
  # geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1), lty = 2, alpha  = 0.3)+
  # geom_point(data = ids %>% left_join(idoutput), aes(k2, lambda0))+

  geom_rect(data= self2$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self2$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Explored"), alpha = 0.6,  col = "black")+
  # scale_x_continuous(breaks = k2_1)+
  # scale_y_continuous(breaks = lambda0_1)+
  # geom_segment(aes(x = 0, y = 0, xend = Inf, yend = 0, col = "Below"))+
  # geom_point(data= self$filters_neg_below, aes(k2, lambda0, col = "Below"))+
  # geom_point(data= self$filters_neg_above, aes(k2, lambda0, col = "Above"))+
  scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey"))+
  labs(col = "", fill = "", x = "K2", y = "lambda0")+
  theme_bw()+
  guides(col = F)+
  geom_point(data = newVPs, aes(k2, lambda0))+
  theme(line = element_blank())  -> plot4;plot4



# plot_grid(plot1, plot2, plot3, plot4)



self2$compute_zone_maybe()
maybe2 <- self2$zone_maybe
maybe2 %>%
  rowid_to_column("id") %>%
  group_split(id) %>%
  map(function(x){

    difx <-  x$k2max -   x$k2min
    dify <-  x$lambda0max -   x$lambda0min

    crossing(k2 = seq(x$k2min, x$k2max, difx/40) %>% round(5),
             lambda0 = seq(x$lambda0min, x$lambda0max, dify/40) %>% round(5)) %>%
      mutate(bloc = x$id)

  }) %>%
  bind_rows() -> newVPs2;newVPs2

self3 <- self2$clone(deep = T)

newVPs2 %>%
  mutate(forjoin = 1) %>%
  left_join(
    VP_df2 %>% select(-k2, -lambda0,-bloc ) %>% mutate(forjoin = 1) %>% distinct()
  ) -> VP_df3


self3$add_VP(VP_df3)

self3$n_filter_reduc()
self3$compute_zone_sure()



plot5 <- ggplot()+

  geom_rect(data= self3$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self3$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe2,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Explored"), alpha = 0.6,  col = "black")+
  # scale_x_continuous(breaks = k2_1)+
  # scale_y_continuous(breaks = lambda0_1)+
  # geom_segment(aes(x = 0, y = 0, xend = Inf, yend = 0, col = "Below"))+
  # geom_point(data= self$filters_neg_below, aes(k2, lambda0, col = "Below"))+
  # geom_point(data= self$filters_neg_above, aes(k2, lambda0, col = "Above"))+
  geom_rect(data = self3$zone_sure,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "VPs"), alpha = 0.6,  col = "black")+
  scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey", "darkgreen"))+
  labs(col = "", fill = "", x = "K2", y = "lambda0")+
  theme_bw()+
  guides(col = F)+
  theme(line = element_blank())+
  geom_point(data = self3$poolVP, aes(k2, lambda0), col = "darkgreen")+
  facet_zoom(xlim = c(1.15,1.26), ylim = c(0.19,0.225), split = F, zoom.size = 1);plot5

VTplot <- self3$plot_VP()+
  geom_point(data = self$targets %>% mutate(label = paste0("\n\nTime = ", time, "\n[", min, "-", max, "]\n\n")),
             aes(time, min,col = label), size = 4)+
  theme_bw()+
  facet_wrap(~"Accepted VPs")+
  labs(x = "Time (days)", y = "Tumor Volume (mm3)", col = "Targets"); VTplot


prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100, 431),
                    max = c(100.5, 431.5))


plothm <- ggplot()+
  geom_rect(data = self3$zone_maybe, aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max))+
  # geom_segment(data = self3$zone_maybe, aes(x = k2min, xend = k2max, y = 1, yend = 1))+
  # geom_rect(data = self3$zone_sure, aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max), col = "darkgreen")+
  theme_bw()+
  facet_wrap(~"Possible zone of VPs kept in memory")


plot_grid(plot_grid(plot1, plot2, plot3, plot4, labels = c("A", "B", "C", "D")), plot5, plot_grid(VTplot, ggplot()+geom_text(aes(x = 1, y=1, label = "?")), labels = c("F", "G")), ncol = 1, rel_heights = c(2,1,1), labels = c("", "E") )




# end so far of the figure ------------------------------------------------
# Beneath are various tries (to delete? )

#
# seq(0,3,0.0015) %>% length() *
#   seq(0,1.4,0.000820) %>% length()
#
# #
# # self4 <- self3$clone(deep = T)
# # self3$compute_zone_maybe()
# # maybe4 <- self3$zone_maybe
# #
# # maybe4 %>%
# #   rowid_to_column("id") %>%
# #   group_split(id) %>%
# #   map(function(x){
# #
# #     difx <-  x$k2max -   x$k2min
# #     dify <-  x$lambda0max -   x$lambda0min
# #
# #     crossing(k2 = seq(x$k2min, x$k2max, 0.0015) %>% round(5),
# #              lambda0 = seq(x$lambda0min, x$lambda0max, 0.000820) %>% round(5)) %>%
# #       mutate(bloc = x$id)
# #
# #   }) %>%
# #   bind_rows() -> final_VPs
# #
# #
# # self4$poolVP <- self4$poolVP %>% slice(0)
# #
# #
# #
# # final_VPs %>%
# #   mutate(forjoin = 1) %>%
# #   left_join(
# #     VP_df2 %>% select(-k2, -lambda0,-bloc ) %>% mutate(forjoin = 1) %>% distinct()
# #   ) -> VP_dffinal
# #
# # self4$add_VP(VP_dffinal)
# #
# # self4$compute_zone_sure()
# # self4$compute_zone_maybe()
# #
# #
# #
# # plot6 <- ggplot()+
# #
# #   geom_rect(data= self4$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
# #   geom_rect(data= self4$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
# #   geom_rect(data = maybe4,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Explored"), alpha = 0.6)+
# #   # scale_x_continuous(breaks = k2_1)+
# #   # scale_y_continuous(breaks = lambda0_1)+
# #   # geom_segment(aes(x = 0, y = 0, xend = Inf, yend = 0, col = "Below"))+
# #   # geom_point(data= self$filters_neg_below, aes(k2, lambda0, col = "Below"))+
# #   # geom_point(data= self$filters_neg_above, aes(k2, lambda0, col = "Above"))+
# #   geom_rect(data = self4$zone_sure,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "VPs"), alpha = 0.6)+
# #   scale_color_manual(values = c("red", "chocolate"))+
# #   scale_fill_manual(values = c("red", "chocolate", "grey", "darkgreen"))+
# #   labs(col = "", fill = "", x = "K2", y = "lambda0")+
# #   theme_bw()+
# #   guides(col = F)+
# #   theme(line = element_blank())+
# #   geom_point(data = self4$poolVP, aes(k2, lambda0), col = "darkgreen")+
# #   facet_zoom(xlim = c(1.15,1.26), ylim = c(0.15,0.27), split = F, zoom.size = 1)
# #
# #
# # plot_grid(plot_grid(plot1, plot2, plot3, plot4), plot5, plot6, ncol = 1, rel_heights = c(2,1,1) )
#
#   # prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100, 431),
#   #                     max = c(100.5, 431.5))
#
#   self2$set_targets(manual = prototiny)
#
# k2_2 <- seq(1.2,1.6,0.002)
# lambda0_2 <- seq(0.8,1,0.005)
#
# # k2_2 <- seq(1.2,1.6,0.00001)
# # lambda0_2 <-  1
# # k2_2 <- seq(0.8,1.2,0.04)
# # lambda0_2 <- seq(0.2,0.4,0.05)
#
# VP_df2 <- crossing(k1 = c(0.5),
#                    k2 = k2_2,
#                    ke = 1 ,#*  seq(0.6,1.4,0.2),
#                    lambda0 =lambda0_2,
#                    lambda1 = c(12),
#                    Vd =  40) %>% #c(0.8,1,1.2)) %>%
#   map_df(function(x){
#
#     if(is.character(x)) return(x)
#     round(x,3)
#
#   } ) %>% mutate(psi = 20); nrow(VP_df2)
#
# self2$add_VP(VP_df2)
# #
# # ids2 <-   VP_df2 %>%
# #   rowid_to_column("id")
# #
# # events2 <- self2$protocols$dose50 %>% mutate(evid = 1) %>%
# #   bind_rows(tibble(time = self$times, evid = 0, amt = 0, cmt = "Conc")) %>%
# #   crossing(id = 1:nrow(VP_df2)) %>%
# #   arrange(id, time)
# #
# # simulations2 <- self2$model$solve(ids2, events2, c(X1 = 50)) %>% as_tibble()
# #
# # simulations2 %>%
# #   filter(time %in% prototiny$time) %>%
# #   left_join(prototiny) %>%
# #   mutate(test = case_when(tumVol > max  ~ "Above",
# #                           tumVol < max ~ "Below",
# #                           T ~ "yes")) %>%
# #   distinct(id, test) -> idoutput2
#
#
# # self2$add_VP(VP_df2)
# self2$plot_2D(k2, lambda0)
# self2$compute_zone_maybe()
# self2$compute_zone_sure()
# self2$n_filter_reduc()
# maybe2 <- self2$zone_maybe
# self2$zone_sure
#
#
#
#
# plot3 <- plot2+
#   geom_segment(aes(x = 1.2, xend = 1.6, y = 0.8, yend = 0.8))+
#   geom_segment(aes(x = 1.2, xend = 1.6, y = 1, yend = 1))+
#   geom_segment(aes(x = 1.6, xend = 1.6, y = 0.8, yend = 1))+
#   geom_segment(aes(x = 1.2, xend = 1.2, y = 0.8, yend = 1))+
#   # geom_point(data = self2$poolVP, aes(k2, lambda0), col  = "darkgreen")+
#   geom_rect(data= self2$zone_sure, aes(xmin =k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "VPs"), alpha = 0.4)+
#   # geom_segment(data = tibble(k2_2), aes(x = k2_2, xend = k2_2, y = min(lambda0_2), yend = max(lambda0_2)), lty = 3)+
#   # geom_segment(data = tibble(lambda0_2), aes(x = min(k2_2), xend = max(k2_2), y = lambda0_2, yend = lambda0_2), lty = 3)+
#   geom_rect(data= self2$filters_neg_above, aes(xmin = min(k2_2), xmax = k2, ymin = lambda0, ymax = max(lambda0_2), fill = "Above"), alpha = 0.4)+
#   geom_rect(data= self2$filters_neg_below, aes(xmin = k2, xmax = max(k2_2), ymin = min(lambda0_2), ymax = lambda0, fill = "Below"), alpha = 0.4)+
#   facet_zoom(xlim = c(1.2,1.6), ylim = c(0.8,1), zoom.size = 1, show.area = T)+
#   scale_fill_manual(values = c("red", "chocolate", "grey", "darkgreen"))
#
# # plot3 <- maybe %>%
# #   rowid_to_column("bloc") %>%
# #   ggplot()+
# #   geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 2, alpha  = 0.3)+
# #   geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1), lty = 2, alpha  = 0.3)+
# #   geom_rect(aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = factor(bloc) ), alpha = 0.6, col = "black")+
# #   theme_bw()+
# #   # geom_segment(aes(x = k2min, xend = k2min, y =lambda0min, yend = lambda0max))+
# #   # geom_segment(aes(x = k2min, xend = k2min, y =lambda0min, yend = lambda0max))+
# #   scale_x_continuous(breaks = k2_1)+
# #   scale_y_continuous(breaks = lambda0_1)+
# #   labs(col = "", fill = "Blocs", x = "K2", y = "lambda0")+
# #   guides(fill= F)+
# #   theme(line = element_blank()); plot3
#
# plot_grid(plot_grid(plot1, plot2), plot3,ncol = 1 )
#
#
#
#
#
#
#
# self$plot_2D()
# # crossing( k2 = seq(0,3,0.4),  lambda0 =seq(0,1.4,0.2)) %>%
# #   ggplot()+
# #   geom_vline(aes(xintercept = k2))+
# #   geom_hline(aes(yintercept = lambda0))+
# #   geom_point(aes(k2, lambda0))+
# #   geom_point(data = crossing( k2 = seq(2,2.4,0.05),  lambda0 =seq(0.6,0.8,0.01)),aes(k2, lambda0) )+
# #   geom_segment(data = tibble(k2 =  seq(2,2.4,0.05)), aes(x = k2, xend = k2, y = 0.6, yend = 0.8))+
# #   geom_segment(data = tibble(lambda0 =seq(0.6,0.8,0.01)), aes(x = 2, xend = 2.4, y = lambda0, yend = lambda0))+
# #   theme_bw()+
# #   facet_zoom(xlim = c(2,2.4), ylim = c(0.6,0.8), zoom.size = 1, split = F,show.area =  )
#
#
#
#
# # VP_df <- crossing(k1 = c(0.5),
# #                           k2 = seq(0,3,0.005),
# #                           ke = 1 ,#*  seq(0.6,1.4,0.2),
# #                           lambda0 =seq(0,1,0.005),
# #                           lambda1 = c(12),
# #                           Vd =  40) %>% #c(0.8,1,1.2)) %>%
# #   map_df(function(x){
# #
# #     if(is.character(x)) return(x)
# #     round(x,3)
# #
# #   } )

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,3,0.005),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,1,0.005),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

self$poolVP
self$plot_2D(k2, lambda0, plotoreturn = 1)


# analyse above
lamda0unique<-  c(self$filters_neg_above$lambda0, self$filters_neg_below$lambda0 )  %>% unique

lambda0seg <-  tibble(lambda0max = lamda0unique) %>%
  add_row(lambda0max = 0) %>%
  add_row(lambda0max = 1) %>%
  arrange(lambda0max) %>%
  mutate(lambda0min = lag(lambda0max)) %>%
  slice(-1)

# lambda0seg %>%
#   rowid_to_column() %>%
#   ggplot()+
#   geom_rect(aes(xmin = 0, xmax = 1.3, ymin = lambda0min, ymax =lambda0max, fill = factor(rowid)), col = "black")+
#   guides(fill = F)

k2unique<-  c(self$filters_neg_above$k2, self$filters_neg_below$k2 )  %>% unique

k2seg <- tibble(k2max = k2unique) %>%
  add_row(k2max = 0) %>%
  add_row(k2max = 1.2) %>%
  arrange(k2max) %>%
  mutate(k2min = lag(k2max)) %>%
  slice(-1)

# k2seg %>%
#  rowid_to_column() %>%
#  ggplot()+
#  geom_rect(aes(xmin = k2min, xmax = k2max, ymin = 0, ymax =1, fill = factor(rowid)), col = "black")+
#  guides(fill = F)

allsquares <- crossing(lambda0seg, k2seg) %>%
  filter(k2min != k2max) %>%
  filter(lambda0min != lambda0max)


allsquares0 <- allsquares

# allsquares0 %>%
#  rowid_to_column() %>%
#  ggplot()+
#  geom_rect(aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax =lambda0max, fill = factor(rowid)),col = "black")+
#  guides(fill = F)

filters <- self$make_filters() %>%
  map_chr(~ gsub("line\\$", "", .x))

filter_above <- "!(lambda0min >= ref$lambda0 &  k2max <= ref$k2)"
for(a in 1:nrow(self$filters_neg_above)){

  ref <- self$filters_neg_above %>% slice(a)

  allsquares <- allsquares %>%
    filter(!!parse_expr(filter_above))

}


filter_below<- "!(lambda0max <= ref$lambda0 &  k2min >= ref$k2)"
for(a in 1:nrow(self$filters_neg_below)){

  ref <- self$filters_neg_below %>% slice(a)

  allsquares <- allsquares %>%
    filter(!!parse_expr(filter_below))

}


allsquares %>%
  rowid_to_column() %>%
  ggplot()+
  geom_rect(aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax =lambda0max, fill = factor(rowid)),col = "black")+
  labs(fill = "bloc")

# Now lets use those square

VP_df2 <- c %>%
  mutate(data = pmap(list(lambda0max, lambda0min, k2max, k2min), function(lambda0max, lambda0min, k2max, k2min){

    crossing(k2 = seq(k2min, k2max, 0.0002), lambda0 = seq(lambda0min, lambda0max, 0.0002) )

  })) %>%
  unnest()%>%
  select(lambda0, k2) %>%
  mutate(k1 = 0.5, ke = 1, lambda1 = 12, Vd = 40)




self$add_VP(VP_df2, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

self$plot_VP(nmax = 20)+
  geom_point(aes(x = 12, y = 100.1846), col = "red")+
  geom_point(aes(x = 40, y = 431.0057), col = "red")

self$plot_VP(nmax = 10)
self$n_filter_reduc()

length(seq(0,3,0.0002)) * length(seq(seq(0,1,0.0002))) * 0.6 / 2000/3600 # temps en heure

self$plot_2D(k2, lambda0, plotoreturn = 1)


# analyse above
lamda0unique<-  c(self$filters_neg_above$lambda0, self$filters_neg_below$lambda0 )  %>% unique

lambda0seg <-  tibble(lambda0max = lamda0unique) %>%
  add_row(lambda0max = 0) %>%
  add_row(lambda0max = 1) %>%
  arrange(lambda0max) %>%
  mutate(lambda0min = lag(lambda0max)) %>%
  slice(-1)

# lambda0seg %>%
#   rowid_to_column() %>%
#   ggplot()+
#   geom_rect(aes(xmin = 0, xmax = 1.3, ymin = lambda0min, ymax =lambda0max, fill = factor(rowid)), col = "black")+
#   guides(fill = F)

k2unique<-  c(self$filters_neg_above$k2, self$filters_neg_below$k2 )  %>% unique

k2seg <- tibble(k2max = k2unique) %>%
  add_row(k2max = 0) %>%
  add_row(k2max = 1.2) %>%
  arrange(k2max) %>%
  mutate(k2min = lag(k2max)) %>%
  slice(-1)

# k2seg %>%
#  rowid_to_column() %>%
#  ggplot()+
#  geom_rect(aes(xmin = k2min, xmax = k2max, ymin = 0, ymax =1, fill = factor(rowid)), col = "black")+
#  guides(fill = F)

allsquares <- crossing(lambda0seg, k2seg) %>%
  filter(k2min != k2max) %>%
  filter(lambda0min != lambda0max)


allsquares0 <- allsquares

# allsquares0 %>%
#  rowid_to_column() %>%
#  ggplot()+
#  geom_rect(aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax =lambda0max, fill = factor(rowid)),col = "black")+
#  guides(fill = F)

filters <- self$make_filters() %>%
  map_chr(~ gsub("line\\$", "", .x))

filter_above <- "!(lambda0min >= ref$lambda0 &  k2max <= ref$k2)"
for(a in 1:nrow(self$filters_neg_above)){

  ref <- self$filters_neg_above %>% slice(a)

  allsquares <- allsquares %>%
    filter(!!parse_expr(filter_above))

}


filter_below<- "!(lambda0max <= ref$lambda0 &  k2min >= ref$k2)"
for(a in 1:nrow(self$filters_neg_below)){

  ref <- self$filters_neg_below %>% slice(a)

  allsquares <- allsquares %>%
    filter(!!parse_expr(filter_below))

}


allsquares %>%
  rowid_to_column() %>%
  ggplot()+
  geom_rect(aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax =lambda0max, fill = factor(rowid)),col = "black")+
  labs(fill = "bloc")+ guides(fill = F)+
  geom_point(data = self$poolVP, aes(k2, lambda0), col="red")

length(seq(0,3,0.0002)) * length(seq(seq(0,1,0.0002))) * 0.35 / 2000/3600 # temps en heure


