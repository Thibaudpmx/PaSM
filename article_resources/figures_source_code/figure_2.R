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
## GitHub: https://github.com/Thibaudpmx/PaSM
## ---------------------------
##
## Notes: 1) Stochasticity in the algorithm, I put a seed for full
## reproductibility in my machine, but results might vary on another one
## (initial sampling of VP to perform RxODE)
##
## 2) Figure are linked one to another (progress throughout the algorithm)
## such as it is necessar to launch the full script (no direct acces to plot H for instance)
## ---------------------------


library(PaSM)
library(cowplot)

# First, lets create an object for the analysis
# Creation of the project
self <- VP_proj_creator$new()

# Taret determine with data
self$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,45))

# Plot A ------------------------------------------------------------------

# Productoin of the cohort of VP to analyse
param <- crossing(k1 = 0.5, k2 = seq(0,6,0.2), ke = 1, lambda0 = seq(0,0.7,0.02), lambda1 = 12, Vd = 40, psi = 20, w0 = 50)

# param %>% filter(k2 == 2.4) # Strange behavior, no line found...
param$k2 <- round(param$k2, 1)
# param %>% filter(k2 == 2.4) # solved by the above line. (impact plot E)

set.seed(1653) ### ATTENTION, HERE IS THE SEED, SO RESULT CAN BE DIFFERENT ON YOUR COMPUTER


# Sample 20 patients among param
param20 <-   param %>%
  sample_n(20) %>%
  # Ok so this initial sampling did not create a full accepted extrapolated zone
  # Therefore, I forced the creation of one by manually
  # removing a useless profile ##"
  filter(!(k2 == 0.2 & lambda0 == 0.46)) %>%
  bind_rows(  param %>% filter(k2 == 3.6 & lambda0 == 0.44) ) # and adding a decisive one for having a green zone

# Note: this manual modification is dependent on my sampling / seed, therefore you can remove the two last line
# as your seed will most proably not produce the same sampling anyway

# Creation of the plot A
plot1 <-  param20 %>%
  ggplot()+
  geom_point(aes(x = k2, y = lambda0, shape = "VP\nsampled\namong\ncohort"))+
  labs(shape = "")+
  theme_bw();plot1


# Plot B ------------------------------------------------------------------

# Manually compute all ids

ids <-   param20 %>%
  rowid_to_column("id")

events <- self$protocols$dose50 %>% mutate(evid = 1) %>%
  bind_rows(tibble(time = self$times, evid = 0, amt = 0, cmt = "Conc")) %>%
  crossing(id = 1:20) %>%
  arrange(id, time)

simulations <- self$model$solve(ids, events, c(X1 = 50)) %>% as_tibble()

# For each target, set values to  "accepted", "above" or "below".

analysis <- simulations %>%
  filter(time %in% self$targets$time) %>%
  mutate(cmt = "tumVol") %>%
  left_join( self$targets) %>%
  mutate(below = tumVol > min, above = tumVol < max) %>%
  mutate(type = case_when(below & above ~ "Accepted",
                          below & !above ~ "Above",
                          !below & above ~ "Below"))

# Label patient as accepted if their two rows (targets) are accepted only
# Of note, no patient were above and below, otherwise an additional
# step would have been need to allocate them
analysis <- analysis %>% distinct(id, type) %>%
  group_by(id) %>%
  nest() %>%
  mutate(label = map_chr(data, function(x){

    if(nrow(x) == 2 & "Accepted" %in% x$type ) x <- x %>% filter(x != "Accepted")

    paste0(x$type, collapse = "&")
  } ))



# Produce the plot
plot2 <- simulations %>%
  left_join(analysis %>%  select(id, label)) %>%
  ggplot()+
  geom_line(aes(time, tumVol, group = id))+
  geom_segment(data = self$targets, aes(x = time, xend = time, y = min, yend = max), col = "red", size = 2)+
  scale_y_log10()+
  facet_wrap(~label) +
  labs(x = "Time (days)", y = "Tumor volume (mm3)")+
  theme_bw(); plot2



# Plot C ------------------------------------------------------------------

# Lets add the label "below", "above" etc to the param values

param20 <- param20 %>%
  rowid_to_column("id") %>%
  left_join(analysis %>%  select(id, label))


# Compute the square with 0 or Inf limit according the parameter
forsquare <- param20 %>%
  # filter(label) %>%
  mutate(k2min = if_else(label == "Above", 0, k2),
         k2max = if_else(label == "Above", k2, Inf),
         lambda0min = if_else(label == "Below", 0, lambda0),
         lambda0max = if_else(label == "Below", lambda0, Inf))

# And compute the plot
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



# plot D ------------------------------------------------------------------

# Let's directly use the package fonction to reduce the filters before making the same plot

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




# Plot E ------------------------------------------------------------------

# Lets compute the filters (right place to see how the main algorithm works !)

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
  # Add point above thanks to the filter created above
  geom_point(data = param %>% filter(!!parse_expr(pointsabove)), aes(k2, lambda0,  col = "Above"))+
  # same for point below
  geom_point(data = param %>% filter(!!parse_expr(pointsbelow)), aes(k2, lambda0,  col = "Below"))+
  geom_point(aes(x = k2, y = lambda0, col = label))+

  geom_rect(data = filtersreduc %>% filter(label != "Accepted"),
            aes(xmin = k2min, xmax = k2max, ymin = lambda0min, ymax = lambda0max, fill = label ),
            alpha = 0.2)+
  theme_bw()+
  scale_color_manual(values = c("red", "darkgreen", "chocolate"))+
  scale_fill_manual(values = c("red", "chocolate"))+
  labs(col = "VPs", fill = "Filters"); plot5


# Plot F ------------------------------------------------------------------


# Compute the square

forsquaregreen <-
  # Consider accepted patient as above
  param20 %>%
  filter(label == "Accepted") %>%
  mutate(label = "Above") %>%
  bind_rows(
    # Consider accepted patient as below
    param20 %>%
      filter(label == "Accepted") %>%
      mutate(label = "Below")
  ) %>%
  # create the square as previously
  mutate(k2min = if_else(label == "Above", 0, k2),
         k2max = if_else(label == "Above", k2, Inf),
         lambda0min = if_else(label == "Below", 0, lambda0),
         lambda0max = if_else(label == "Below", lambda0, Inf)) %>%
  mutate(label2 = if_else(label == "Above", "Above\nLower Limit", "Below\nUpper Limit"))


# perform te plot base of plot E
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


# Plot G ------------------------------------------------------------------

# nothing  need to be computed before, just a summary of the end of first algorithm loop

plot7 <- param20 %>%
  ggplot()+
  # geom_point(aes(x = k2, y = lambda0, col = label))+
  geom_point(data = param %>% filter(!!parse_expr(pointsabove)), aes(k2, lambda0,  col = "Extrap. Above"))+
  geom_point(data = param %>% filter(!!parse_expr(pointsbelow)), aes(k2, lambda0,  col = "Extrap. Below"))+
  geom_point(data = param %>% filter( k2>= 3.6 & k2<=3.8 & lambda0>= 0.4 &lambda0<=0.44), aes(k2, lambda0, col = "Interp. Accepted"))+
  geom_point(data = param20, aes(x = k2, y = lambda0, col = "Simulated\n(original sampling)"))+
  scale_color_manual(values = c("red", "chocolate",   "darkgreen", "black"))+
  labs(col = "   VPs")+
  theme_bw() ; plot7



# plot H ------------------------------------------------------------------

# Compute the new patient to sample
param %>%
  # remove patients already rejected
  left_join( param %>% filter(!!parse_expr(pointsabove)) %>% mutate(nop = T)) %>%
  filter(is.na(nop)) %>%
  select(-nop) %>%
  left_join( param %>% filter(!!parse_expr(pointsbelow)) %>% mutate(nop = T)) %>%
  filter(is.na(nop)) %>%
  select(-nop) %>%
  # And already accepted
  filter(!( k2>= 3.6 & k2<=3.8 & lambda0>= 0.4 &lambda0<=0.44)) %>%
  # And then sampled them (seed is probably still one, but not important here anyway)
  sample_n(20) -> news


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



# Plot I ------------------------------------------------------------------


# Perform the full analysis with our algorithm
self$add_VP(param, use_green_filter = T)

# Reduce the filter
self$n_filter_reduc()


# and perform the final plot
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



# Final grid --------------------------------------------------------------


plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,plot7, plot8,plot9, nrow = 3, labels = letters)



# Save 300 dip for article
them <- theme(plot.title = element_text(hjust = 0.5))
tiff(width = 3500, height = 2300,filename = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig2.tiff", res = 300)
plot_grid(plot1 + ggtitle("Parameter sampling")+them,
          plot2 + ggtitle("Simulations")+them,
          plot3 + ggtitle("Computation of rejection zones")+them,
          plot4 + ggtitle("Reduction of rejection zones")+them,
          plot5 + ggtitle("Rejection extrapolations")+them,
          plot6 + ggtitle("Computation and reduction\n of acceptance zones")+them,
          plot7 + ggtitle("Acceptance interpolations\nand first iteration summary")+them,
          plot8 + ggtitle("New parameter sampling")+them,
          plot9 + ggtitle("Final parameter space mapping"), labels = letters[1:9])
dev.off()
shell.exec( "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig2.tiff")

