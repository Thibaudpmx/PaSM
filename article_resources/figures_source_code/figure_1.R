## ---------------------------
## Script name: Figure 1 article
##
## Purpose of script: production of figure 1
##c
## Author: Thibaud Derippe
##
## Date Created: 2022-04-04 (day of code repartition from R6object.R)
##
## Under GPL-3 License
## Email: thibaud.derippe@gmail.com
## GitHub: https://github.com/Thibaudpmx/PaSM
## ---------------------------
##
## Notes: No intermediate results; script can be launched to
##  directly produce the figure with no long computational time
##
##
## ---------------------------


library(PaSM)
library(cowplot)
# Demo rejection ----------------------------------------------------------


# Create a project
self <- VP_proj_creator$new() # note: by default, the path of the source file read is
                              # file.path( find.package("PaSM"), "config_Simeoni.r")

# Determine some targets
self$set_targets(manual = tibble(protocol = "dose50", cmt = "tumVol", time = c(12,45), min = c(21,50), max = c(62,200)))

# Set the values of the potential VP of interest, with just different K2 and lambda 0
param <- crossing(k1 = 0.5, w0 = 50, k2 = seq(2,6,2), ke = 1, lambda0 = seq(0.02,0.06,0.01), lambda1 = 12, Vd = 40, psi = 20) %>%
  rowid_to_column("id")

# Take the VP with the loxest values (high)
# param %>% arrange(k2, desc(lambda0)) %>% slice(1) %>% pull(id) -> idlowest

# Create the events table
events <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows( self$protocols$dose50 %>% mutate(evid = 0) %>% select(-time) %>% crossing(time = 0:50)) %>%
  arrange(time) %>%
  crossing(id = 1:nrow(param))

# And perform the simulations
simul <- self$model$solve(param, events, c(X1 = 50) ) %>%
  as_tibble() %>%
  left_join(param) %>%
  mutate(l0 = if_else(lambda0 == 0.06 , " 0.06", "< 0.06"))



# Simul 3 correspond to a patient inside the target
param3 <- crossing(k1 = 0.5,w0 = 50, k2 = c(0.8), ke = 1, lambda0 = c(0.05), lambda1 = 12, Vd = 40, psi = 20) %>%
  rowid_to_column("id")

events3 <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows( self$protocols$dose50 %>% mutate(evid = 0) %>% select(-time) %>% crossing(time = 0:50)) %>%
  arrange(time) %>%
  crossing(id = 1:nrow(param3))

simul3 <- self$model$solve(param3, events3, c(X1 = 50) ) %>%
  as_tibble()

# Perform the main plot 1
plot1 <-
  ggplot()+
  # Targets as black segments
  geom_segment(data =  self$targets,
               aes(x = time, xend = time, y = min, yend = max), col ="black", size = 2)+
  geom_segment(aes(x =  12, xend = 25, y = 40, yend = 120), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))+
  geom_segment(aes(x = 30, xend = 45, y = 120, yend = 100), arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"))+
  geom_label(aes(x = 25, y = 120, label = "Targets"))+
  # Accepted profil (simul3)
  geom_line(data = simul3, aes(time, tumVol), col = "darkgrey")+
  geom_segment(aes(x =  6, xend = 8, y =  45 , yend = 70), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), col = "darkgrey", size = 0.5)+
  geom_text(aes(7.5, 90, label = "(VP accepted)"), col = "darkgrey", size = 3.5)+
  # Simulations of all rejected ID (simul)
  geom_line(data=simul, aes(time, tumVol, group = id, size = l0, alpha = l0, col = factor(k2)))+
  scale_size_manual(values = c(2,1))+ #  with thicker profiles for lambda0 == 0.06
  scale_alpha_manual(values = c(1,0.3))+ # in addition to an higher alpha
  labs(col = "k2", alpha = "lambda0", size= "lambda0", x = "Time (days)", y = "Tumor Volume (mm3)")+
  # Y log scale
  scale_y_log10()+
  # Add the lambda0 values as texxt
  geom_text(data = simul%>% filter(time == 50), aes(x = 53, y = tumVol, label = lambda0, col = factor(k2)))+
  # Bracket  lower lambda0
  geom_segment(aes(x = 56, xend = 56, y = 9, yend = 60))+
  geom_segment(aes(x = 51, xend = 56, y = 9, yend = 9))+
  geom_segment(aes(x = 51, xend = 56, y = 60, yend = 60))+
  geom_segment(aes(x = 56, xend = 58, y = 20, yend = 20))+
  geom_text(aes(x = 70, y = 20, label = "lower lambda0"))+
  # Text higher k2
  geom_text(aes(x = 70, y = 7.3, label = "higher k2", col = "4"))+
  geom_text(aes(x = 70, y = 0.6, label = "higher k2", col = "6"))+
  #  Bracket  lower lambda0 and higher k2 (1/2)
  geom_segment(aes(x = 56, xend = 56, y = 0.7, yend = 5))+
  geom_segment(aes(x = 51, xend = 56, y = 5, yend = 5))+
  geom_segment(aes(x = 51, xend = 56, y = 0.7, yend = 0.7))+
  geom_segment(aes(x = 56, xend = 58, y = 2, yend = 2)) +
  geom_text(aes(x = 70, y = c(2), label = "lower lambda0\nand higher k2"))+
  #  Bracket  lower lambda0 and higher k2 (2/2)
  geom_segment(aes(x = 56, xend = 56, y = 0.04, yend = 0.4))+
  geom_segment(aes(x = 51, xend = 56, y = 0.04, yend = 0.04))+
  geom_segment(aes(x = 51, xend = 56, y = 0.4, yend = 0.4))+
  geom_segment(aes(x = 56, xend = 58, y = 0.15, yend = 0.15))+
  geom_text(aes(x = 70, y = c(0.15), label = "lower lambda0\nand higher k2"))+
  # Labels rejections and other
  coord_cartesian(xlim = c(0,80)) +
  geom_label(aes(x = 70, y = 80, label = "Step II: Rejection by\nextrapolation"))+
  geom_segment(aes(x = 70, xend = 70, y = 50, yend = 30), arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"))+
  geom_label(aes(x = 67, y = 300, label = "Step I: Rejection after RxODE"), col = "red")+
  geom_segment(aes(x = 70, xend = 50, y = 200, yend = 100), arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"), col = "red")+

  geom_text(aes(x = 70, y = c(0.07), label = "(all other parameters strictly equal)"), size = 3)+
  theme_bw();plot1




# Demo extrapolation acceptation ------------------------------------------



# Creation of the cohort (extreme patient)
paramaccext <- crossing(k1 = 0.5, k2 = c(0.8,0.9), ke = 1, lambda0 = c(0.05,0.03), lambda1 = 12, Vd = 40, psi = 20, w0 = 50) %>%
  rowid_to_column("id")


eventsaccext <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows( self$protocols$dose50 %>% mutate(evid = 0) %>% select(-time) %>% crossing(time = 0:50)) %>%
  arrange(time) %>%
  crossing(id = 1:nrow(paramaccext))

simulaccext <- self$model$solve(paramaccext, eventsaccext, c(X1 = 50) ) %>%
  as_tibble()%>%
  left_join(paramaccext) %>%
  mutate(name = paste0("\nlambda0 = ", lambda0,"\nk2 = ", k2))


# Simul patients extrapolated as accepted
paramaSure<- crossing(k1 = 0.5, k2 = seq(0.8,0.9,0.02), ke = 1, lambda0 = seq(0.03,0.05,0.005), lambda1 = 12, Vd = 40, psi = 20, w0 = 50) %>%
  rowid_to_column("id")

eventsSure<- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows( self$protocols$dose50 %>% mutate(evid = 0) %>% select(-time) %>% crossing(time = 0:50)) %>%
  arrange(time) %>%
  crossing(id = 1:nrow(paramaSure))

simulSure<- self$model$solve(paramaSure, eventsSure, c(X1 = 50)) %>%
  as_tibble() %>%
  left_join(paramaSure) %>%
  mutate(name = paste0("\nlambda0 = ", lambda0,"\nk2 = ", k2))

# plot generation
plot2 <-
  ggplot()+
  # targets, y scale and black and white theme
  geom_segment(data =  self$targets,
               aes(x = time, xend = time, y = min, yend = max), col ="black", size = 2)+
  scale_y_log10()+
  theme_bw()+
  # Create the ribbon of the accepted zone
  geom_ribbon(data = simulaccext %>% select(time, id, tumVol) %>% spread(id, tumVol), aes(x = time, ymin = `2`, ymax = `3`),
              alpha = 0.2, fill = "darkgreen")+
  # Add the two extreme patients
  geom_line(data = simulaccext %>% filter(id %in% 2:3),aes(time, tumVol, group = id, col = name), size = 2 )+
  # And the extrapolated ones
  geom_line(data = simulSure ,aes(time, tumVol, group = id), col = "darkgreen" )+
  # arrow below upper limits
  geom_segment(aes(x = 45, xend = 45, y = 120, yend  =  175), col = "#00BFC4", size = 2, arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))+
  geom_segment(aes(x = 45, xend = 45, y = 62.6, yend  =  90), size = 2, col = "#F8766D", arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"))+
  geom_text(aes(x = 67, y = 200, label ="VPs with lambda0 < 0.05 and\n k2 > 0.8 below upper limits" ), col = "#00BFC4")+
  # arrow above lower limits
  geom_segment(aes(x = 12, xend = 12, y = 31, yend  =   36.4), size = 2, col = "#00BFC4", arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))+
  geom_segment(aes(x = 12, xend = 12, y =  26.2 , yend  =  30), size = 2, col = "#F8766D", arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"))+
  geom_text(aes(x = 67, y = 75, label ="VPs with lambda0 > 0.03 and\n k2 < 0.9 above lower limits" ), col = "#F8766D")+
  # Description zone
  geom_label(aes(65, y = 120, label = "All VPs with\n0.03 < lambda0 < 0.05\n and 0.8 < k2 < 0.9\nwill be accepted\nin this green zone"))+
  # various
  coord_cartesian(xlim = c(0,80))+
  geom_text(aes(x = 65, y = 50, label = "(all other parameters strictly equal)"), size = 3.5)+
  labs(col = "VP parameters", x = "Time (days)", y = "Tumor Volume (mm3)");plot2



# Final grid plot ---------------------------------------------------------

plot_grid(plot1, plot2, labels = c("a", "b"))

# Save 300 dip for article

tiff(width = 4500, height = 1500,filename = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig1.tiff", res = 300)
plot_grid(plot1, plot2, labels = c("a", "b"))
dev.off()
shell.exec( "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig1.tiff")
