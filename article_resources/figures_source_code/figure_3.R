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
## GitHub: https://github.com/Thibaudpmx/PaSM
## ---------------------------
##
## Description of the second algorithm
##
##
##
## ---------------------------
library(PaSM)
library(cowplot)
library(ggforce)


# define the tiny targets for this analysis

prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100, 431),
                    max = c(100.5, 431.5))

# create a new project (still using the default Simeoni file)
self <- VP_proj_creator$new()
# Assignet the targets
self$set_targets(manual = prototiny)


# plotA -------------------------------------------------------------------

# Define the k2 and lambda0 values for the first iteration (step 0.4 and 0.2)

k2_1 <- seq(0,3,0.4)
lambda0_1 <- seq(0,1.4,0.2)

k2_1edgeless <- k2_1[2:(length(k2_1)-1)] # remove the two boundaries values for each
lambda0_1edgeless <-  lambda0_1[2:(length(lambda0_1)-1)]

points <- crossing(k2_1 = k2_1edgeless, lambda0_1 = lambda0_1edgeless) # cross the edgeless param values

# Perform the plot
plot1 <- ggplot()+
  geom_point(data = points, aes(alpha = "Equidistant\nfirst VPs",k2_1, lambda0_1))+
  geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 3)+
  geom_segment(data = tibble(lambda0_1), aes(lty = "Blocs of VPs\n", x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1))+
  geom_segment(data = tibble(lambda0_1), aes(lty = " Whole\nParameter\nspace\n", x = min(k2_1), xend = max(k2_1), y = min(lambda0_1), yend = min(lambda0_1)), size = 1.2)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = max(lambda0_1), yend = max(lambda0_1)), lty = 1, size = 1.2)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = min(k2_1), y = min(lambda0_1), yend = max(lambda0_1)), lty = 1, size = 1.2)+
  geom_segment(data = tibble(lambda0_1), aes(x = max(k2_1), xend = max(k2_1), y = min(lambda0_1), yend = max(lambda0_1)), lty = 1, size = 1.2)+
  theme_bw()+
  scale_linetype_manual(values = c(1,3))+
  theme(line = element_blank())+
  scale_x_continuous(breaks = k2_1)+
  scale_y_continuous(breaks = lambda0_1)+
  scale_alpha_manual(values = 1)+
  labs(x = "K2", y = "lambda0", lty = "", alpha = ""); plot1



# Plot B ------------------------------------------------------------------


# Let's create the cohort of VPs with other param having a single values  (arbitrarly picked)

VP_df <- crossing(k1 = c(0.5),
                  k2 = k2_1edgeless,
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =lambda0_1edgeless,
                  lambda1 = c(12),
                  w0 = 50,
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } ) %>% mutate(psi = 20)



# And use algorithm 1 to have the results
self$add_VP(VP_df)
# Then compute the maybe zone
self$compute_zone_maybe()
maybe <- self$zone_maybe

# And replace the inf by the max values of lambda0 and k2
maybe$lambda0max[maybe$lambda0max == Inf] <- max(lambda0_1)
maybe$k2max[maybe$k2max == Inf] <- max(k2_1)

# Before finally perform plot number 2

plot2 <-
  ggplot()+
    theme_bw()+

  geom_rect(data= self$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "To explore" ), alpha = 0.6,  col = "black", size = 1.2)+
  scale_x_continuous(breaks = k2_1)+
  scale_y_continuous(breaks = lambda0_1)+

  geom_text(data = maybe %>% rowid_to_column("id"), aes(x = (k2min + k2max) /2, y = (lambda0min + lambda0max)/2, label = id) )+
  geom_segment(aes(x = 0, y = 0, xend = Inf, yend = 0, col = "Below"))+
   scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey"))+
  geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 3)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1), lty = 3)+
  geom_point(data= self$filters_neg_below, aes(k2, lambda0, fill = "Below"), shape = 21, size = 2)+
  geom_point(data= self$filters_neg_above, aes(k2, lambda0, fill = "Above"), shape = 21, size = 2)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = min(lambda0_1), yend = min(lambda0_1)), lty = 1)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = max(k2_1), y = max(lambda0_1), yend = max(lambda0_1)), lty = 1)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_1), xend = min(k2_1), y = min(lambda0_1), yend = max(lambda0_1)), lty = 1)+
  geom_segment(data = tibble(lambda0_1), aes(x = max(k2_1), xend = max(k2_1), y = min(lambda0_1), yend = max(lambda0_1)), lty = 1)+
  labs(col = "", fill = "", x = "K2", y = "lambda0")+
  theme_bw()+
  guides(col = F)+
  theme(line = element_blank()); plot2

# plot_grid(plot1, plot2)



# Plot C ------------------------------------------------------------------

# Compute manually the extreme VP in each blocs
# Good way to see how the algorithm to reduce the blocs works

testabove <- self$clone(deep = T) # copy above
testbelow <- self$clone(deep = T) # copy below

testabove$targets$max <- Inf  # target from min to Inf
testbelow$targets$min <-- Inf # target from -Inf to max

maybeabove <- maybe
maybebelow <- maybe

# Copy paste slighlty adjust of the real algorithm
# goal is to select only the lowest and highest patient of each blocs
for(a in self$param){

  if(a %in% self$param_increase[[1]]){

    names(maybeabove)[names(maybeabove) == paste0(a, "max")] <- a
    names(maybebelow)[names(maybebelow) == paste0(a, "min")] <- a
  }else if(a %in% self$param_reduce[[1]] ){

    names(maybeabove)[names(maybeabove) == paste0(a, "min")] <- a
    names(maybebelow)[names(maybebelow) == paste0(a, "max")] <- a
  }else{

    #to handle the rest...
    names(maybeabove)[names(maybeabove) == paste0(a, "min")] <- a
    names(maybebelow)[names(maybebelow) == paste0(a, "max")] <- a
  }
}

# gather all ids
ids <-   bind_rows(maybeabove %>% rowid_to_column("bloc"),
                   maybebelow %>% rowid_to_column("bloc")) %>%
  rowid_to_column("id") %>%
  select(!!!parse_exprs(self$param), bloc ,id) %>%
  mutate(psi = 20)

events <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows(tibble(time = self$times, evid = 0, amt = 0, cmt = "Conc")) %>%
  crossing(id = 1:nrow(ids)) %>%
  arrange(id, time)

# and perform their simulations
simulations_blocs <- self$model$solve(ids, events, c(X1 = 50)) %>% as_tibble()

# Compute the distance (not really used in practice...)
simulations_blocs %>%
  left_join(self$targets) %>%
  filter(!is.na(min)) %>%
  mutate(distance = abs(tumVol - min)) %>%
  left_join(ids %>% select(id, bloc)) %>%
  group_by(bloc) %>%
  summarise(distance = sum(distance)) %>%
filter(bloc %in% c(1,3,5))-> distance

# and perform the plots!
plot3 <- simulations_blocs %>%
  left_join(ids %>% select(id, bloc)) %>%
  filter(bloc %in% c(1,3,5,7)) %>%
  ggplot()+

  geom_line(aes(time, tumVol, group = id))+
  geom_point(data = self$targets, aes(time, min, col = "targets"))+
  facet_wrap(~paste0("Bloc ", bloc))+
  scale_color_manual(values = "red")+
  geom_segment(data = tibble(bloc = 7, xmin = c(30,30 ), xmax = c(45,45), ymin = c(50,100), ymax = c(100,50)),
               aes(x = xmin, xend = xmax, y = ymin, yend=ymax), col = "red", size = 2)+
  geom_segment(data = tibble(bloc =c(1,3,5), xmin = c(30), xmax = c(35), ymin = c(70), ymax = c(50)),
               aes(x = xmin, xend = xmax, y = ymin, yend=ymax), col = "darkgreen", size = 2)+
  geom_segment(data = tibble(bloc = c(1,3,5), xmin = c(35), xmax = c(45), ymin = c(50), ymax = c(100)),
               aes(x = xmin, xend = xmax, y = ymin, yend=ymax), col = "darkgreen", size = 2)+
  geom_text(data = distance, aes(x = 30, y = 30, label = paste0("distance: ", round(distance))))+
  theme_bw()+
  labs(x = "Time (days)", y = "Tumor Volume (mm3)", col = "")+
  scale_y_log10(); plot3

# plot_grid(plot1, plot2, plot3)






# Plot D ------------------------------------------------------------------

# Perform the real maybe zone reduction
maybe <- reduce_maybe2(maybe) # the seventh zones as dissapeared !

# maybe %>% slice(1) -> x

# In each zone, compute the new equidistant VPs, removing the borders
maybe %>%
  rowid_to_column("id") %>%
  group_split(id) %>%
  map(function(x){

    difx <-  x$k2max -   x$k2min
    dify <-  x$lambda0max -   x$lambda0min

    crossing(k2 = (seq(x$k2min, x$k2max, difx/7) %>% round(3))[2:7],
             lambda0 = (seq(x$lambda0min, x$lambda0max, dify/7)%>% round(3))[2:7] ) %>%
      mutate(bloc = x$id)

  }) %>%
  bind_rows() -> newVPs

# same but with bords for plot
maybe %>%
  rowid_to_column("id") %>%
  group_split(id) %>%
  map(function(x){

    difx <-  x$k2max -   x$k2min
    dify <-  x$lambda0max -   x$lambda0min

    crossing(k2 = (seq(x$k2min, x$k2max, difx/7) %>% round(3)),
             lambda0 = (seq(x$lambda0min, x$lambda0max, dify/7)%>% round(3)) ) %>%
      mutate(bloc = x$id)

  }) %>%
  bind_rows() -> newVPswithbord

# create the data for producing the segments (the borders)
newVPswithbord %>%
  left_join(

    newVPs %>%
      group_by(bloc) %>%
      summarise(k2min = min(k2), k2max = max(k2), lambda0min = min(lambda0), lambda0max = max(lambda0))

  ) -> forseg

# reduction of the domain
k2_2 <- k2_1[k2_1 <=1.6]


# plot production
plot4 <- ggplot()+
  geom_segment(data = forseg, aes(x = k2, xend = k2, y = lambda0, yend = lambda0max), lty = 1)+
  geom_segment(data = forseg, aes(x = k2, xend = k2max, y = lambda0, yend = lambda0), lty = 1)+
  geom_rect(data= self$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "To explore"), alpha = 0.6,  col = "black", size = 1.2)+
  scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey"))+
  labs(col = "", fill = "", x = "K2", y = "lambda0")+
  geom_segment(data = tibble(k2_2), aes(x = k2_2, xend = k2_2, y = min(lambda0_1), yend = max(lambda0_1)), lty = 3)+
  geom_segment(data = tibble(lambda0_1), aes( x = min(k2_2), xend = max(k2_2), y = lambda0_1, yend = lambda0_1), lty = 3)+
  theme_bw()+
   geom_segment(data = tibble(lambda0_1), aes(x = min(k2_2), xend = max(k2_2), y = min(lambda0_1), yend = min(lambda0_1)), lty = 1)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_2), xend = max(k2_2), y = max(lambda0_1), yend = max(lambda0_1)), lty = 1)+
  geom_segment(data = tibble(lambda0_1), aes(x = min(k2_2), xend = min(k2_2), y = min(lambda0_1), yend = max(lambda0_1)), lty = 1)+
  geom_segment(data = tibble(lambda0_1), aes(x = max(k2_2), xend = max(k2_2), y = min(lambda0_1), yend = max(lambda0_1)), lty = 1)+
  guides(col = F)+
  geom_point(data= self$filters_neg_below, aes(k2, lambda0, fill = "Below"), shape = 21, size = 2)+
  geom_point(data= self$filters_neg_above, aes(k2, lambda0, fill = "Above"), shape = 21, size = 2)+
  theme(line = element_blank())+
  geom_point(data = newVPs, aes(k2, lambda0)); plot4

# plot_grid(plot1, plot2, plot3, plot4, labels = LETTERS)



# Plot E ------------------------------------------------------------------

# ok let's analyse all the same blocs (here all pooled, in the real algorihtm within a loop)

self2 <- self$clone(deep = T)

newVPs %>%
  mutate(forjoin = 1) %>%
  left_join(
    VP_df %>% select(-k2, -lambda0) %>% mutate(forjoin = 1) %>% distinct()
  ) -> VP_df2


self2$add_VP(VP_df2)
self2$n_filter_reduc()


# Compute the new maybe zone,  same steps as before
self2$compute_zone_maybe()

maybe2 <- self2$zone_maybe
maybe2$lambda0max[maybe2$lambda0max == Inf] <- max(lambda0_1)
maybe2$k2max[maybe2$k2max == Inf] <- max(k2_2)

maybe2 <-  reduce_maybe2(maybe2, obj = self2)

# and perform the plot
ggplot()+
  geom_rect(data= self2$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self2$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Explored"), alpha = 0.6,  col = "black")+
  geom_rect(data = maybe2,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Explored"), alpha = 0.6,  col = "black", size = 1.2)+
scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey"))+
  labs(col = "", fill = "", x = "K2", y = "lambda0")+
  theme_bw()+
  guides(col = F)+
  geom_point(data= self2$filters_neg_below, aes(k2, lambda0, fill = "Below"), shape = 21, size = 2)+
  geom_point(data= self2$filters_neg_above, aes(k2, lambda0, fill = "Above"), shape = 21, size = 2)+
  geom_segment(data = tibble(k2_2), aes(x = k2_2, xend = k2_2, y = min(lambda0_1), yend = max(lambda0_1)), lty = 3)+
  geom_segment(data = tibble(lambda0_1), aes( x = min(k2_2), xend = max(k2_2), y = lambda0_1, yend = lambda0_1), lty = 3)+
  theme(line = element_blank())  -> plot5;plot5




# plot F ------------------------------------------------------------------

# Compute the new VPs, basically same steps as before
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

# Of not, not every possible VPS have been computed here but for grpahical explenation it is largely enough.

plot6 <- ggplot()+

  geom_rect(data= self3$filters_neg_above, aes(xmin = 0, xmax = k2, ymin = lambda0, ymax = Inf, fill = "Above", col = "Above"), alpha = 1)+
  geom_rect(data= self3$filters_neg_below, aes(xmin = k2, xmax = Inf, ymin = 0, ymax = lambda0, fill = "Below", col = "Below"), alpha = 1)+
  geom_rect(data = maybe2,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = "Explored"), alpha = 0.6,  col = "black")+
  scale_color_manual(values = c("red", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate", "grey", "darkgreen"))+
  labs(col = "", fill = "", x = "K2", y = "lambda0", alpha = "")+
  theme_bw()+
  guides(col = F)+
  theme(line = element_blank())+
  geom_point(data = self3$poolVP, aes(k2, lambda0, alpha= "Accepted\nVPs"), col = "darkgreen")+
  scale_alpha_manual(values = 1)+
  geom_segment(data = tibble(k2_1), aes(x = k2_1, xend = k2_1, y = min(lambda0_1), yend = max(lambda0_1)), lty = 2)+
  geom_segment(data = tibble(lambda0_1), aes( x = min(k2_1), xend = max(k2_1), y = lambda0_1, yend = lambda0_1), lty = 2)+
  geom_rect(data = maybe,aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max), alpha = 0.0,lty = 1,size = 1.2,  col = "black")+
  facet_zoom(xlim = c(1.15,1.26), ylim = c(0.19,0.225), split = F, zoom.size = 1);plot6



# plot G --------------------------------------------------------------


VTplot <- self3$plot_VP(nmax = Inf)+
  geom_point(data = self$targets %>% mutate(label = paste0("\n\nTime = ", time, "\n[", min, "-", max, "]\n\n")),
             aes(time, min,col = label), size = 4)+
  theme_bw()+
  facet_wrap(~"Accepted VPs")+
  labs(x = "Time (days)", y = "Tumor Volume (mm3)", col = "Targets"); VTplot


# Final Grid
plot_grid(plot_grid(plot1, plot2, plot3,plot4, nrow = 1, labels = c("a", "b", "c", "d")),
          plot_grid(plot5, plot6, VTplot, nrow = 1, rel_widths = c(1,2,1),labels =  c("e", "f", "g")), ncol = 1)



# 300 dpi image -----------------------------------------------------------


them <- theme(plot.title = element_text(hjust = 0.5))
tiff(width = 4700, height = 2000,filename = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig3.tiff", res = 300)


plot_grid(plot_grid(plot1+ ggtitle("Parameter space division by\ngenerating equidistant VPs")+them,
                    plot2+ ggtitle("First algorithm results and \nplausibility zones computation")+them,
                    plot3+ ggtitle("Verification of plausibility zones")+them,
                    plot4+ ggtitle("New parameter space division\nfor each plausibility zone")+them, nrow = 1, labels = c("a", "b", "c", "d")),
          plot_grid(plot5+ ggtitle("Zoom in on each plausibility zone")+them,
                    plot6+ ggtitle("Final parameter space mapping")+them,
                    VTplot+ ggtitle("Simulations of accepted VPs")+them,
                    nrow = 1, rel_widths = c(1,2,1),labels =  c("e", "f", "g")), ncol = 1)


dev.off()
shell.exec( "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig3.tiff")

