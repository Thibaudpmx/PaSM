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

