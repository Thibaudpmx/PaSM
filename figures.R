library(peccary)
library(QSPVP)
library(RxODE)
library(progress)
library(R6)
library(crayon)
# Step 1 create or load project ---------------------------------------------------
create_VT_Project("D:/these/Second_project/QSP/modeling_work/VT_simeoni")
# Two compartments --------------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


## One increase one decrease (k2 and lambda0)

in_dec <- VP_proj_creator$new()


in_dec$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,19, 30,45))

VP_df_in_dec <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.2),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.16,0.025),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )



in_dec$add_VP(VP_df_in_dec, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)


## Two increases

in_in <- VP_proj_creator$new()

in_in$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,19, 30,45))

VP_df_in_in <- crossing(k1 = c(0.5),
                         k2 = 2,
                         ke = 1 ,#*  seq(0.6,1.4,0.2),
                         lambda0 = seq(0.05,1,0.05),
                         lambda1 = c(12),
                         Vd =  seq(0,80,1)) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )



in_in$add_VP(VP_df_in_in, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
in_in$n_filter_reduc()




plot_grid(nrow = 1,
# data_segment_plot(data = in_dec$data, targets = in_dec$targets) + theme_bw(),
in_dec$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = F) + ggtitle("", "one increaser - one decreaser"),
# in_dec$plot_2D(x = lambda0, y = k2, plotoreturn = 1,add_point = F) + ggtitle("", "one increaser - one decreaser"),
in_in$plot_2D(x = Vd, y = lambda0, plotoreturn = 1,add_point = F) + ggtitle("", "two increasers")
)



# Several time points --------------------------------------------------------

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
twozone <- VP_proj_creator$new()
above_or_below <- VP_proj_creator$new()
above_or_below2 <- VP_proj_creator$new()

twozone$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,45))
above_or_below$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(45))
above_or_below2$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12))
twozone$targets <- twozone$targets  %>% mutate(max = c(93,300))
above_or_below$targets <- twozone$targets  %>% slice(2)
above_or_below2$targets <- twozone$targets  %>% slice(1)

VP_df_twozone <- crossing(k1 = c(0.5),
                         k2 = seq(0,3,0.05),
                         ke = 1 ,#*  seq(0.6,1.4,0.2),
                         lambda0 =seq(0,1,0.02),
                         lambda1 = c(12),
                         Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

twozone$add_VP(VP_df_twozone, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
above_or_below$add_VP(VP_df_twozone, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
above_or_below2$add_VP(VP_df_twozone, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)

twozone$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T)
above_or_below$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T)

above_or_below2$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T) +
  scale_y_continuous(limits = c(0,0.2))


# Computa patient



line <- VP_df_twozone %>%
  rowid_to_column("id") %>%
  mutate(protocol = "dose50")


protocol <-  line %>%
  mutate(protocol2 = map(protocol, ~ self$protocols[[.x]])) %>%
  select(id, protocol2) %>%
  unnest(protocol2)

res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                    icmt = twozone$initial_cmt_values, time_vec =twozone$times,
                    pardf = twozone$parameters_default_values, model = twozone$model) %>%
  as_tibble()#;res

# patient both above and below

res %>%
  filter(time == 45 & tumVol > 300 ) %>% pull(id) -> idtem

res %>%
  filter(id %in% idtem) %>%
  filter(time == 12 & tumVol < 21.3) %>% pull(id) %>% unique() -> idbothbelowabove

res %>%
  filter(id %in% idbothbelowabove) %>%
  ggplot()+
  geom_line(aes(time, tumVol, group = factor(id)))+
  scale_y_log10()+
  geom_segment(data = twozone$targets, aes(x = time, xend = time, y = min, yend = max), col ="red", size = 2)+
  theme_bw()+
  facet_wrap(~"Above AND Below")

# patient above or below with only time = 42

res %>%
  filter(time == 45 & (tumVol > 300 | tumVol <  50.3)) %>%
  mutate(forwrap = if_else(tumVol > 300, "Strictly above", "Strictly below")) %>%
  select(id, forwrap) -> idtemp

res %>%
  # filter(id %in% idtemp$id) %>%
  left_join(idtemp) %>%
  mutate(forwrap = if_else(is.na(forwrap), "In target", forwrap)) %>%
  ggplot()+
  geom_line(aes(time, tumVol, group = factor(id)))+
  scale_y_log10()+
  geom_segment(data = twozone$targets %>% slice(2), aes(x = time, xend = time, y = min, yend = max), col ="red", size = 2)+
  theme_bw()+
  # scale_color_manual(values = c("red", "chocolate"))+
  facet_wrap(~forwrap)

# patient above or below with only time = 12

res %>%
  filter(time == 12 & (tumVol > 92.98855 | tumVol <   21.34172)) %>%
  mutate(forwrap = if_else(tumVol > 300, "Strictly above", "Strictly below")) %>%
  select(id, forwrap) -> idtemp

res %>%
  # filter(id %in% idtemp$id) %>%
  left_join(idtemp) %>%
  mutate(forwrap = if_else(is.na(forwrap), "In target", "Rejected")) %>%
  ggplot()+
  geom_line(aes(time, tumVol, group = factor(id)))+
  scale_y_log10()+
  geom_segment(data = twozone$targets %>% slice(1), aes(x = time, xend = time, y = min, yend = max), col ="red", size = 2)+
  theme_bw()+
  # scale_color_manual(values = c("red", "chocolate"))+
  facet_wrap(~forwrap)





# Several YTYPE --------------------------------------------------------

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
pk_only <- VP_proj_creator$new()
pd_only <- VP_proj_creator$new()
pk_pd <- VP_proj_creator$new()

pk_pd$set_targets(filter = Dose==50 , timeforce = c(2, 45))
pk_pd$targets <- pk_pd$targets %>% ungroup %>%  slice(2,3 )
pk_pd$targets$min[[2]] <- 0.02
pk_pd$targets$max[[1]] <- 100

pk_only$set_targets(filter = Dose==50 & cmt != "tumVol"  , timeforce = c(2))
pk_only$targets$min <- 0.02

pd_only$set_targets(filter = Dose==50  & cmt == "tumVol" , timeforce = c(45))
pd_only$targets <- pk_pd$targets %>% slice(1)

VP_df_pkpd <- crossing(k1 = c(0.5),
                          k2 = seq(0,3,0.2),
                          ke = seq(0,4,0.2) ,#*  seq(0.6,1.4,0.2),
                          lambda0 = 0.05,
                          lambda1 = c(12),
                          Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )



pk_only$add_VP(VP_df_pkpd, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
pd_only$add_VP(VP_df_pkpd, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
pk_pd$add_VP(VP_df_pkpd, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)

pk_only$plot_VP()+ scale_x_continuous(limits = c(0,5)) + scale_y_continuous() + theme_bw()
pk_only$plot_2D(x = k2, y = ke, plotoreturn = 1,add_point = T)

pd_only$plot_VP()+ theme_bw()
pd_only$plot_2D(x = k2, y = ke, plotoreturn = 1,add_point = T)


# Rare value --------------------------------------------------------------

# demo <- tibble( )
#
# ggplot()+
#   geom_line

# determine a tiny tiny target

line <- tibble(k1 = 0.5,
                 k2 = 1.2344,
                 ke = 1 ,#*  seq(0.6,1.4,0.2),
                 lambda0 = 0.4644,
                 lambda1 = c(12),
                 Vd =  40) %>%
  rowid_to_column("id") %>%
  mutate(protocol = "dose50")


protocol <-  line %>%
  mutate(protocol2 = map(protocol, ~ self$protocols[[.x]])) %>%
  select(id, protocol2) %>%
  unnest(protocol2)

res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                    icmt = self$initial_cmt_values, time_vec =self$times,
                    pardf = self$parameters_default_values, model = self$model) %>%
  as_tibble()#;res

res %>%
  ggplot()+
  geom_line(aes(time, tumVol))

res %>%
  filter(time == 12 | time == 40) %>%
  pull(tumVol)
## end determine the value
source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100.1845, 431.0056),
                    max = c(100.1847, 431.0058))
self <- VP_proj_creator$new()

self$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,45))
self$targets <- prototiny


# VP_df <- crossing(k1 = c(0.5),
#                           k2 = seq(0,3,0.005),
#                           ke = 1 ,#*  seq(0.6,1.4,0.2),
#                           lambda0 =seq(0,1,0.005),
#                           lambda1 = c(12),
#                           Vd =  40) %>% #c(0.8,1,1.2)) %>%
#   map_df(function(x){
#
#     if(is.character(x)) return(x)
#     round(x,3)
#
#   } )

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,3,0.05),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,1,0.05),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)


self$filters_neg_below %>%
  filter(iddummy == 107)
#
# self$plot_VP()
# self$n_filter_reduc()
# self$poolVP
# self$plot_2D(k2, lambda0, plotoreturn = 1)+
#   geom_point(aes(x = 1.2344, y = 0.4644), col = "darkgreen") +
#   geom_point(aes(x = 0.15, y = 0.05)) +
#   guides(shape = F)+
#   # geom_point(data = VP_df, aes(k2, lambda0))
#   geom_point(data = self$poolVP, aes(x = k2, y = lambda0), col = "darkgreen")

self$plot_2D(k2, lambda0, plotoreturn = 1)+
  geom_point(aes(x = 1.2344, y = 0.4644), col = "darkgreen") +
  # geom_point(aes(x = 0.15, y = 0.05)) +
  guides(shape = F)+
  # geom_point(data = VP_df, aes(k2, lambda0))
  geom_point(data = self$poolVP, aes(x = k2, y = lambda0), col = "darkgreen")+
  geom_hline(yintercept = unique(VP_df$lambda0))+
  geom_vline(xintercept = unique(VP_df$k2)) -> plottemp; plottemp


plottemp +
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))
# self$filters_neg_below %>%
#   filter(k2 == 0)

# round 2: capture the white space

demo <- self$filters_neg_above %>%
  slice(1)

# how to be sure there is no place?
VP_df2 <- crossing(k1 = c(0.5),
                  k2 = seq(0.4-0.05,0.4+0.05,0.0002),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0.1-0.05,0.1 + 0.05,0.0002),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } );VP_df2





self2 <- self
self2$filters_neg_above <- tibble()
self2$filters_neg_below <- tibble()

self$add_VP(VP_df2, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

self2$plot_VP()
self$plot_2D(k2, lambda0, plotoreturn = 1)

self2$plot_2D(k2, lambda0, plotoreturn = 1)+
  geom_hline(yintercept = unique(VP_df2$lambda0))+
  geom_vline(xintercept = unique(VP_df2$k2))


  scale_x_continuous(limits = c(1.15,1.15+0.05))+
  scale_y_continuous(limits = c(0.2-0.05,0.2))

  VP_df3<- crossing(k1 = c(0.5),
                    k2 = seq(0 ,0 +0.05,0.0001),
                    ke = 1 ,#*  seq(0.6,1.4,0.2),
                    lambda0 =seq(0,0.05 + 0.05,0.0001),
                    lambda1 = c(12),
                    Vd =  40) %>% #c(0.8,1,1.2)) %>%
    map_df(function(x){

      if(is.character(x)) return(x)
      round(x,3)

    } );VP_df3

self$add_VP(VP_df3, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

self$plot_2D(k2, lambda0, plotoreturn = 1)+
geom_hline(yintercept = unique(VP_df$lambda0))+
  geom_vline(xintercept = unique(VP_df$k2))


# Rare value loop ---------------------------------------------------------
source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100.1845, 431.0056),
                    max = c(100.1847, 431.0058))
self <- VP_proj_creator$new()

self$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,45))
self$targets <- prototiny

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,3,0.05),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,1,0.05),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)


self$plot_2D(k2, lambda0, plotoreturn = 1)+
  geom_point(aes(x = 1.2344, y = 0.4644), col = "darkgreen") +
  # geom_point(aes(x = 0.15, y = 0.05)) +
  guides(shape = F)+
  # geom_point(data = VP_df, aes(k2, lambda0))
  geom_point(data = self$poolVP, aes(x = k2, y = lambda0), col = "darkgreen")+
  geom_hline(yintercept = unique(VP_df$lambda0))+
  geom_vline(xintercept = unique(VP_df$k2)) -> plottemp; plottemp

demo <- c(1:3)


# round 2: capture the white space
self$filters_neg_above <-  self$filters_neg_above %>%
  mutate(level = 1)

self$filters_neg_below<-  self$filters_neg_below %>%
  mutate(level = 1)



filter1 <- self$filters_neg_above
filter2 <- self$filters_neg_below
for(a in 1:nrow(filter1) ){

  line <- filter1 %>%
    slice(a)

  VP_df_temp <- crossing(k1 = c(0.5),
                     k2 = seq(min(0,line$k2-0.05),line$k2+0.05,0.001),
                     ke = 1 ,#*  seq(0.6,1.4,0.2),
                     lambda0 =seq(min(0,line$lambda0-0.05),line$lambda0 + 0.05,0.001),
                     lambda1 = c(12),
                     Vd =  40) %>% #c(0.8,1,1.2)) %>%
    map_df(function(x){

      if(is.character(x)) return(x)
      round(x,3)

    } );VP_df_temp

  self$add_VP(VP_df_temp, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

}


self$plot_2D(k2, lambda0, plotoreturn = 1)
self$n_filter_reduc()

for(a in 1:nrow(filter2) ){

  line <- filter2 %>%
    slice(a)

  VP_df_temp <- crossing(k1 = c(0.5),
                         k2 = seq(min(0,line$k2-0.05),line$k2+0.05,0.001),
                         ke = 1 ,#*  seq(0.6,1.4,0.2),
                         lambda0 =seq(min(0,line$lambda0-0.05),line$lambda0 + 0.05,0.001),
                         lambda1 = c(12),
                         Vd =  40) %>% #c(0.8,1,1.2)) %>%
    map_df(function(x){

      if(is.character(x)) return(x)
      round(x,3)

    } );VP_df_temp

  self$add_VP(VP_df_temp, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

}

self$n_filter_reduc()

self$poolVP

self$plot_2D(k2, lambda0, plotoreturn = 1)


self$plot_2D(k2, lambda0, plotoreturn = 1)+
  coord_cartesian(xlim =  c(0.4,0.4001), ylim=c(0.08,0.1) )



self$plot_2D(k2, lambda0, plotoreturn = 1)+
  coord_cartesian(xlim =  c(1.233,1.236), ylim=c(0.464,0.465) )

self2 <- self
self2$filters_neg_above <- tibble()
self2$filters_neg_below <- tibble()

self$add_VP(VP_df2, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

self2$plot_VP()
self$plot_2D(k2, lambda0, plotoreturn = 1)

self2$plot_2D(k2, lambda0, plotoreturn = 1)+
  geom_hline(yintercept = unique(VP_df2$lambda0))+
  geom_vline(xintercept = unique(VP_df2$k2))


scale_x_continuous(limits = c(1.15,1.15+0.05))+
  scale_y_continuous(limits = c(0.2-0.05,0.2))

VP_df3<- crossing(k1 = c(0.5),
                  k2 = seq(0 ,0 +0.05,0.0001),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.05 + 0.05,0.0001),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } );VP_df3

self$add_VP(VP_df3, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

self$plot_2D(k2, lambda0, plotoreturn = 1)+
  geom_hline(yintercept = unique(VP_df$lambda0))+
  geom_vline(xintercept = unique(VP_df$k2))


# various -----------------------------------------------------------------

a <- seq(0,1000,142.8571)

crossing(d = a, b = a) %>%
ggplot()+
  geom_point(aes(d, b))+
  geom_hline(yintercept = a)+
  geom_vline(xintercept = a)+
  theme_bw()+
  labs(x = "Bcl2", x= "Mcl1")

allpoints <- crossing(Mcl1 = a, Bcl2 = a, BIM = a)

plot_ly()%>%
  add_markers(type = "scatter3d",
              mode = "markers",
              data = allpoints,
              x = ~Mcl1,
              y = ~Bcl2,
              z = ~BIM,
              opacity = 1)



plot_grid(
ggplot()+
  geom_rect(aes(xmin = 2, xmax = 5, ymin = 0, ymax = 4, fill = "below"), alpha = 0.3)+
  geom_rect(aes(xmin = 0, xmax = 3, ymin = 3, ymax = 5, fill = "above"), alpha = 0.3)+
  geom_point(aes(x = 3,3), col = "red")+
  geom_point(aes(x = 2,4), col = "chocolate")+
  scale_fill_manual(values = c("red", "chocolate"))+
  theme_bw()+
  labs(x = "k2", y = "lambda0")+
  ggtitle("Initial"),

ggplot()+
  geom_rect(aes(xmin = 2, xmax = 5, ymin = 0, ymax = 4, fill = "below"), alpha = 0.3)+
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 3, ymax = 5, fill = "above"), alpha = 0.3)+
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 3, fill = "desired"), alpha = 0.3)+
  geom_point(aes(x = 2,3), col = "red")+
  geom_point(aes(x = 2,4), col = "chocolate")+
  scale_fill_manual(values = c("red", "chocolate", "darkgreen"))+
  theme_bw()+
  labs(x = "k2", y = "lambda0")+
  ggtitle("k2 modif"),

ggplot()+
  geom_rect(aes(xmin = 2, xmax = 5, ymin = 0, ymax = 4, fill = "below"), alpha = 0.3)+
  geom_rect(aes(xmin = 0, xmax = 3, ymin = 4, ymax = 5, fill = "above"), alpha = 0.3)+
  geom_rect(aes(xmin = 3, xmax = 5, ymin = 4, ymax = 5, fill = "desired"), alpha = 0.3)+
  geom_point(aes(x = 3,4), col = "red")+
  geom_point(aes(x = 2,4), col = "chocolate")+
  scale_fill_manual(values = c("red", "chocolate", "darkgreen"))+
  theme_bw()+
  labs(x = "k2", y = "lambda0")+
  ggtitle("lambda modif"), nrow = 1
)


plot_grid(
  ggplot()+
    geom_rect(aes(xmin = 2, xmax = 5, ymin = 0, ymax = 4, fill = "below"), alpha = 0.3)+
    geom_rect(aes(xmin = 0, xmax = 3, ymin = 4.5, ymax = 5, fill = "above"), alpha = 0.3)+
    geom_point(aes(x = 3,4.5), col = "red")+
    geom_point(aes(x = 2,4), col = "chocolate")+
    scale_fill_manual(values = c("red", "chocolate"))+
    theme_bw()+
    labs(x = "k2", y = "lambda0"),

  ggplot()+
    geom_rect(aes(xmin = 2, xmax = 5, ymin = 0, ymax = 4, fill = "below"), alpha = 0.3)+
    geom_rect(aes(xmin = 0, xmax = 2, ymin = 4.5, ymax = 5, fill = "above"), alpha = 0.3)+
    geom_point(aes(x = 2,4.5), col = "red")+
    geom_point(aes(x = 2,4), col = "chocolate")+
    scale_fill_manual(values = c("red", "chocolate"))+
    theme_bw()+
    labs(x = "k2", y = "lambda0"),

  ggplot()+
    geom_rect(aes(xmin = 2, xmax = 5, ymin = 0, ymax = 4, fill = "below"), alpha = 0.3)+
    geom_rect(aes(xmin = 0, xmax = 3, ymin = 4.5, ymax = 5, fill = "above"), alpha = 0.3)+
    geom_point(aes(x = 3,4.5), col = "red")+
    geom_point(aes(x = 2,4), col = "chocolate")+
    scale_fill_manual(values = c("red", "chocolate"))+
    theme_bw()+
    labs(x = "k2", y = "lambda0"), nrow = 1
)


self$filters_neg_above %>%
  ggplot()+
  geom_line(aes(k2, lambda0,col = "above"))+
  geom_line(data = self$filters_neg_below, aes(k2, lambda0, col = "below"))+
  scale_y_log10()


# Capture white square ----------------------------------------------------

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")




prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100.1845, 431.0056),
                    max = c(100.1847, 431.0058))
self <- VP_proj_creator$new()

self$set_targets(manual = prototiny)
self$targets <- prototiny


# VP_df <- crossing(k1 = c(0.5),
#                           k2 = seq(0,3,0.005),
#                           ke = 1 ,#*  seq(0.6,1.4,0.2),
#                           lambda0 =seq(0,1,0.005),
#                           lambda1 = c(12),
#                           Vd =  40) %>% #c(0.8,1,1.2)) %>%
#   map_df(function(x){
#
#     if(is.character(x)) return(x)
#     round(x,3)
#
#   } )

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

VP_df2 <- allsquares %>%
  mutate(data = pmap(list(lambda0max, lambda0min, k2max, k2min), function(lambda0max, lambda0min, k2max, k2min){

      crossing(k2 = seq(k2min, k2max, 0.0002), lambda0 = seq(lambda0min, lambda0max, 0.0002) )

  })) %>%
  unnest()%>%
  select(lambda0, k2) %>%
  mutate(k1 = 0.5, ke = 1, lambda1 = 12, Vd = 40)




self$add_VP(VP_df2, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F,keepRefFiltaftDis = T)

self$poolVP %>%
  slice(1:10) %>%
  unnest() %>%
  ggplot()+
  geom_line(aes(time, tumVol))+
  geom_point(aes(x = 12, y = 100.1846), col = "red")+
  geom_point(aes(x = 40, y = 431.0057), col = "red")

self$n_filter_reduc()


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



# Time vs  pct vert -------------------------------------------------------

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.0025),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.16,0.0025),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")







protoinf <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 0, max = Inf)
infi <- VP_proj_creator$new()
infi$set_targets(manual = protoinf)
infi$targets
self <- infi
tinf <- Sys.time()
infi$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
tinf <- difftime(Sys.time(), tinf, units = "s")

protobase00 <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 0.1, max = 1500)
base00 <- VP_proj_creator$new()
base00$set_targets(manual = protobase00)
base00$targets


tbase00 <- Sys.time()
base00$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
tbase00 <- difftime(Sys.time(), tbase00, units = "s")


protobase0 <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 25, max = 700)
base0 <- VP_proj_creator$new()
base0$set_targets(manual = protobase0)
base0$targets


tbase0 <- Sys.time()
base0$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
tbase0 <- difftime(Sys.time(), tbase0, units = "s")



protobase <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 50, max = 400)
base <- VP_proj_creator$new()
base$set_targets(manual = protobase)
base$targets


tbase <- Sys.time()
base$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
tbase <- difftime(Sys.time(), tbase, units = "s")



protobase2 <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 100, max = 200)
base2 <- VP_proj_creator$new()
base2$set_targets(manual = protobase2)
base2$targets


tbase2 <- Sys.time()
base2$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
tbase2 <- difftime(Sys.time(), tbase2, units = "s")


protobase3 <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 150, max = 175)
base3 <- VP_proj_creator$new()
base3$set_targets(manual = protobase3)
base3$targets


tbase3 <- Sys.time()
base3$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
tbase3 <- difftime(Sys.time(), tbase3, units = "s")


protonone <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 150, max = 150.000000001)
none <- VP_proj_creator$new()
none$set_targets(manual = protonone)
none$targets


tnone <- Sys.time()
none$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)
tnone <- difftime(Sys.time(), tnone, units = "s")

library(ggpubr)


restime <- tibble(time = c(tinf,tbase0,tbase00, tbase, tbase2, tbase3, tnone)) %>%
  mutate(time = as.double(time),
         nVP =  list(infi, base0, base00, base, base2, base3, none) %>%
           map_dbl(~ nrow(.x$poolVP))) %>%
  mutate(pctVP = nVP  / nrow(VP_df)) %>%
  mutate(x = 1-pctVP)

restime <-  tibble(x = c(0,0.571,0.101,0.733,0.925,0.982,1), time = c(83.6,35.2,77.9,24.6,12.8,6.33,5.86))
restime %>%
  ggscatter( x = "x", y = "time",
            add = "reg.line",                                 # Add regression line
            conf.int = F,                                  # Add confidence interval
            add.params = list(color = "chocolate",
                              fill = "chocolate")
  )+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 30)  +
  labs(x = "Fraction VP rejection", y = "Time to compute all VPs (second)")+ # Add correlation coefficient
  geom_hline(aes(yintercept = 70, col = "Reference time brut way (64-70 sec)"), lty = 2)+
  geom_hline(aes(yintercept = 64, col = "Reference time brut way (64-70 sec)"), lty = 2)+
  labs(col = "")

# tibble(time = c(tinf,tbase0,tbase00, tbase, tbase2, tbase3, tnone)) %>%
#   mutate(time = as.double(time),
#          nVP =  list(infi, base0, base00, base, base2, base3, none) %>%
#            map_dbl(~ nrow(.x$poolVP))) %>%
#   mutate(pctVP = nVP  / nrow(VP_df)) %>%
#   mutate(x = 1-pctVP) %>%
#   ggplot()+
#   geom_point(aes(1-pctVP, time ))+
#   # geom_point(aes(1-pctVP, time ))+
#   theme_bw() +
#   labs(x = "Percentage VP rejection", y = "Time to compute all VPs (second)")

tibble(time = c(tinf,tbase0,tbase00, tbase, tbase2, tbase3, tnone)) %>%
  mutate(time = as.double(time),
         nVP =  list(infi, base0, base00, base, base2, base3, none) %>%
           map_dbl(~ nrow(.x$poolVP))) %>%
  mutate(pctVP = nVP  / nrow(VP_df)) %>%
  mutate(pct_Rejec = 1-pctVP) %>%
  select(pct_Rejec, time)




# do it manually ----------------------------------------------------------

data_VT <- read.table("D:/Peccary_Annexe/Exemple_demo/DATA/Simeoni.txt", header = T, sep = ";", na.strings = ".") %>%
  mutate(protocol = paste0("dose", Dose), cmt = if_else(YTYPE == 2, "tumVol", "Conc")) %>%
  as_tibble %>%
  filter(!is.na(cmt))

model_RxODE <-  RxODE({
  d/dt(Central) <- -ke * Central
  Conc <- Central/Vd

  tumVol <- X1 + X2 + X3 + X4
  growth <- lambda0 * X1/((1 + (lambda0 * tumVol/lambda1)^psi)^(1/psi))
  d/dt(X1) <- growth - X1 * Conc * k2
  d/dt(X2) <- X1 * Conc * k2 - k1 * X2
  d/dt(X3) <- k1 * (X2 - X3)
  d/dt(X4) <- k1 * (X3 - X4)
})


targetsformanual <-  tribble(~protocols, ~time, ~cmt, ~ min, ~max,
                    "dose50",45 ,"tumVol", 150,175)
initial_cmt_values <- c(X1 = 50) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,52, 1) # times you want to see your observations

t00 <- Sys.time()
demo <- VP_df %>%
  rowid_to_column("id") %>%
  mutate(group = ceiling(id/2000)) %>%
  # filter(group == 2) -> x
  group_split(group) %>%
  map(function(x){


    events <-
      x %>%
      select(id) %>%
      mutate(cmt = "Central", amt = 50, evid = 1, time = 0) %>%
      bind_rows(

        x %>%
          select(id) %>%
          crossing(time = times) %>%
          mutate(amt = 0, evid = 0,cmt = "Central")
          # mutate(cmt = "Central", amt = 50, evid = 1)

      ) %>%
      arrange(id)



res <- model_RxODE$solve(x %>% mutate(psi = 20), events, initial_cmt_values) %>%
  as_tibble

res %>%
  filter(time == 40) %>%
  filter(tumVol > 150 | tumVol < 175) %>% pull(id) -> idtokeep

x %>%
  rowid_to_column() %>%
  filter(rowid %in% idtokeep) %>%
  left_join(res %>%   filter(id %in% idtokeep) %>%  rename(rowid = id) %>% group_by(rowid) %>% nest())
  }) %>%
  invoke(.fn = bind_rows)
difftime(Sys.time(),t00, units = "s")

# testf <- function(){
#
#
#   res %>%
#     filter(time == 40) %>%
#     filter(tumVol < 150 | tumVol > 175) %>%
#     pull(id) -> idtemp
#
#   VP_df2 %>%
#     left_join(line %>% distinct(id, rowid)) %>%
#     filter(!id %in% idtemp )
#
# }
# library(microbenchmark)
# microbenchmark(testf(), unit = "s")

# Influence number of param -------------------------------------------------------

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.7),
                  ke = seq(0,8,0.7) ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.16,0.014),
                  lambda1 = seq(9,14.5,0.5),
                  Vd =  31:42) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } );VP_df

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

protohere <- tibble(protocol = "dose50", cmt = "tumVol", time = c(45), min = 150, max = 175)


# 0 param
noimpact0 <- VP_proj_creator$new()
noimpact0$set_targets(manual = protohere)
tnoimpact0 <- Sys.time()
noimpact0$add_VP(npersalve = 2000, VP_df, reducefilteratend = F)
tnoimpact0 <- difftime(Sys.time(), tnoimpact0, units = "s")
rm(noimpact0)

# 1 param

noimpact1 <- VP_proj_creator$new()
noimpact1$set_targets(manual = protohere)
noimpact1$param_increase$tumVol <- noimpact1$param_increase$tumVol[-1]
tnoimpact1 <- Sys.time()
noimpact1$add_VP(npersalve = 2000, VP_df, reducefilteratend = F)
tnoimpact1 <- difftime(Sys.time(), tnoimpact1, units = "s")
rm(noimpact1)

# 2 param

noimpact2 <- VP_proj_creator$new()
noimpact2$set_targets(manual = protohere)
noimpact2$param_increase$tumVol <- noimpact2$param_increase$tumVol[-c(1,2)]
tnoimpact2 <- Sys.time()
noimpact2$add_VP(npersalve = 2000,VP_df, reducefilteratend = F)
tnoimpact2 <- difftime(Sys.time(), tnoimpact2, units = "s")
rm(noimpact2)

# 3 param

noimpact3 <- VP_proj_creator$new()
noimpact3$set_targets(manual = protohere)
noimpact3$param_increase$tumVol <- noimpact3$param_increase$tumVol[-c(1,2,3)]
tnoimpact3 <- Sys.time()
noimpact3$add_VP(npersalve = 2000,VP_df, reducefilteratend = F)
tnoimpact3 <- difftime(Sys.time(), tnoimpact3, units = "s")
rm(noimpact3)

# 4 param

noimpact4 <- VP_proj_creator$new()
noimpact4$set_targets(manual = protohere)
noimpact4$param_increase$tumVol <- noimpact4$param_increase$tumVol[-c(1:4)]
tnoimpact4 <- Sys.time()
noimpact4$add_VP(npersalve = 2000,VP_df, reducefilteratend = F)
tnoimpact4 <- difftime(Sys.time(), tnoimpact4, units = "s")
rm(noimpact4)


# 5 param

noimpact5 <- VP_proj_creator$new()
noimpact5$set_targets(manual = protohere)
noimpact5$param_increase$tumVol <- noimpact5$param_increase$tumVol[-c(1:4)]
noimpact5$param_reduce$tumVol <- character()
tnoimpact5 <- Sys.time()
noimpact5$add_VP(npersalve = 2000,VP_df, reducefilteratend = F)
tnoimpact5 <- difftime(Sys.time(), tnoimpact5, units = "s")
rm(noimpact5)


tibble(time = c(tnoimpact0, tnoimpact1, tnoimpact2, tnoimpact3, tnoimpact4, tnoimpact5), paramnoinflu = c(1:6))
tibble(time = c(22,30,65.8,97.85,110,109, 21, 26, 66, 101,120,122), paramnoinflu = c(0:5, 0:5),
       nsim = c(rep("1000", 6), rep("2000",6))) %>%
  mutate(paramnoinflu = 5 - paramnoinflu) %>%
  mutate(y = paste0(paramnoinflu, "/5")) %>%
  ggplot()+
  geom_hline(aes(yintercept = 66, lty = "Ref\n66\nsec"))+
  scale_linetype_manual(values = 2)+
  geom_point(aes(y, time, col = nsim)) +
  geom_line(aes(paramnoinflu + 1, time, col = nsim)) +
  theme_bw()+
  labs(col= "nsim", x = "Number of parameter involved area",y = "Time to perform all VPs (in seconds) ", lty = "")
