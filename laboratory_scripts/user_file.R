library(QSPVP)
library(profvis)
library(microbenchmark)
# Step 1 create or load project ---------------------------------------------------
create_VT_Project("D:/these/Second_project/QSP/modeling_work/VT_simeoni")

# Step 2: fulfil the file

shell.exec("D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config.r")
shell.exec("D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/2_verif.r")

# test
ind_param = tibble(k1 = 3, k2 = 0.4, ke = 27.5/52.8, lambda0 = 0.06, lambda1 = 6, Vd = 52.8)


celltheque


dose = 100
add_events <- tibble(cmt = c("Central")) %>%  # Administration sampling, with "concX" being replaced by concentration
  mutate(time = 0, amt = dose)

#
res <- simulations(ind_param = ind_param, add_events = add_events,returnSim = T)
simulations(ind_param = ind_param, add_events = add_events,returnSim = T) %>%
  filter(time >23)
simulations(ind_param = ind_param, add_events = add_events,returnSim = F)
#
# res %>%
#   gather("key", "value", Conc, tumVol) %>%
#   ggplot()+
#   geom_line(aes(time, value))+
#   facet_wrap(~key, scales = "free")+
#   scale_y_log10()
read.table("D:/Peccary_Annexe/Exemple_demo/Simeoni/closeIV/IndividualParameters/estimatedIndividualParameters.txt", header = T, sep = ",") %>%
  as_tibble %>%
  filter(id == 218)

# Step3: create celltheque
 #%>%
  # slice(1:(3*6000))# what you want to add in your celltheque as individuals


# targets <-  tribble(~protocols, ~time, ~cmt, ~ min, ~max,
#                     "dose0",10,"tumVol", 50, 120,
#                     "dose0",25,"tumVol", 70,220,
#                     "dose0",0,"Conc",0,1,
#                     "dose50",10,"tumVol", 40,125,
#                     "dose50",25,"tumVol", 30,70,
#                     "dose50",10,"Conc", 1e-6,1,
#                     "dose100",25,"tumVol", 3,167)


# Two compart ---------------------------------------------------------------


source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()

self$param_increase
self$make_filters()
# self$set_targets(filter = cmt == "tumVol", ntime = 6)
self$set_targets(filter = Dose ==50  & cmt == "tumVol",timeforce = c(12,19, 30,45))

VP_df <- crossing(k1 = c(0.5),
         k2 = seq(0,8,0.2),
         ke = 1 ,#*  seq(0.6,1.4,0.2),
         lambda0 =seq(0,0.16,0.025),
         lambda1 = c(12),
         Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = T)



self$plot_VP()

self$plot_2D(x = k2, y = lambda0,plotoreturn = 1, add_point = T)

self$plot_2D(x = k2, y = lambda0,toaddneg = VP_df, plotMain = T)
self$n_filter_reduc()
# self$fill_simul()
self$filters_neg_above
self$filters_neg_below
self$filters_pos_below
self$filters_pos_above


# 2cm big -----------------------------------------------------------------

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter = cmt == "tumVol"  ,timeforce = c(12,19, 30,45))

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.0025),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.16,0.0025),
                  lambda1 = c(12),
                  w0 = 50,
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

# self$add_VP(VP_df, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F, methodFilter = 1)

self$add_VP(VP_df, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F, methodFilter = 2, use_green_filter = F)


self$n_filter_reduc()
self$plot_VP(nmax = 2000)

self$targets

self$poolVP %>%
  unnest(simul) %>%
  filter(time == 12 & tumVol < 21) %>%
  pull(cellid) -> idtarget

self$plot_2D(x = k2, y = lambda0, plotMain = F, plotoreturn = 1, add_point = T)

self$n_filter_reduc()

ti <- Sys.time()
self$fill_simul()

difftime( Sys.time(), ti, units = "sec")
self$filters_neg_above
self$filters_neg_below
self$filters_pos_above
self$filters_pos_below



# three compart ---------------------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter = Dose == 50 & cmt == "tumVol",timeforce = c(12,19, 30,45))



VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.5),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  w0 = 50,
                  lambda0 =seq(0.04,0.16,0.025),
                  lambda1 = c(12),
                  Vd =  c(20:40)) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

# VP_df <- crossing(k1 = c(0.5),
#                   k2 = seq(0,8,0.1),
#                   ke = 1 ,#*  seq(0.6,1.4,0.2),
#                   lambda0 =seq(0,0.16,0.01),
#                   lambda1 = c(12),
#                   Vd =  c(0:40)) %>% #c(0.8,1,1.2)) %>%
#   map_df(function(x){
#
#     if(is.character(x)) return(x)
#     round(x,3)
#
#   } )

self$add_VP(VP_df, fillatend =F, reducefilteratend = F)

self$plot_VP()
self$n_filter_reduc()
# self$plot_2D(k2, lambda0, plotMain = T)
self$plot_3D(k2, lambda0, Vd)

self$filters_neg_above
self$filters_neg_below
self$filters_pos_above
self$filters_pos_below


# 4 dimensions and more

# Infinite dimension ------------------------------------------------------


source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()
#
# saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/longrun.RDS")

# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter = Dose == 50 & cmt == "tumVol",timeforce = c(12,19, 30,45))

VP_df <- crossing(k1 = c(0.5),
                  w0 = 50,
                  k2 = seq(0,8,0.1),
                  ke = seq(0.6,1.4,0.4),
                  lambda0 =seq(0,0.16,0.03),
                  lambda1 = c(10,12,14,33),
                  Vd =  c(0:40)) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

t0 <- Sys.time()
self$add_VP(VP_df, fillatend = F, reducefilteratend = F,use_green_filter = T, npersalve = 2000, time_compteur = T, pctActivGreen = 0.75)
Sys.time() - t0


self$timeTrack

self$plot_VP(nmax = 2000)
self$n_filter_reduc()

self$compute_zone_maybe()
self$compute_zone_sure()

self$zone_maybe
self$zone_sure


# verif
self <- VP_proj_creator$new()
#
# saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/longrun.RDS")

# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter = Dose == 50 & cmt == "tumVol",timeforce = c(12,19, 30,45))

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0.1,0.2,0.05),
                  ke = seq(0.6,1,0.1),
                  lambda0 =seq(0.03,0.06,0.01),
                  lambda1 = seq(10,12,0.5),
                  Vd =  seq(13,14,0.1)) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = F,use_green_filter = F, npersalve = 2000, time_compteur = F, pctActivGreen = 0.75)

self$plot_VP()


# w0 wrong ----------------------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(manual = tibble(protocol = "dose50", cmt = "tumVol", time = 40, min = 80, max = 85))

VP_df <- crossing(k1 = c(0),
                  k2 = 0:8,
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 = 50,
                  lambda1 = c(20,100,1),
                  w0 = seq(20,400,1),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )



# self$add_VP(VP_df, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F, methodFilter = 1)

self$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F, methodFilter = 2, use_green_filter = F)



# PK AND PD ---------------------------------------------------------------





source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



# self$set_targets(filter = cmt == "tumVol", ntime = 6)
self$set_targets(filter = Dose ==50 ,timeforce = c(3,19, 30,45))

tar <- self$targets ; tar$min[[5]] <- 0.005; tar$max[[2]] <-100; tar$max[[3]] <- 200; tar$max[[5]] <- 300
self$targets  <- tar

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.1),
                  ke = seq(0,3,0.2) ,#*  seq(0.6,1.4,0.2),
                  lambda0 = c(0.025),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend =F, use_green_filter = F)

self$plot_VP()
self$poolVP


# Lindner !  --------------------------------------------------------------


# Defining targets


source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner.r")

self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",740,"Pore", 20, 23
))


VP_df <- crossing(Bcl20 = seq(100,1000,100),
                  Bclxl0 = seq(100,1000,100),
                  Mcl10 = 50 ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(100,1000,100),
                  PUMA0 =seq(100,1000,100),
                  NOXA0 =  seq(100,1000,100),
                  BAXc0 = 1000,
                  BAK0 = 1000) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

lindner$add_VP(VP_df, fillatend = F, reducefilteratend = F)



self$plot_VP()



# Massive screening -------------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner.r")

self$targets <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                        "unique",740,"Pore", 20, 20.0001
)


domain <- tribble(~param, ~from, ~to, ~step,
                  "Bcl20", 0, 1000, 1 ,
                  "Bclxl0", 0, 1000, 1,
                  "Mcl10", 0,1000,1,
                  "BIM0", 0, 1000, 1 ,
                  "PUMA0", 0, 1000, 1,
                  "NOXA0", 0,1000,1,
                  "BAXc0", 1000,1000,0,
                  "BAK0" , 1000,1000,0

)

demo <- self$big_screening(domain)


demo2 <- self$big_screening(domain = demo[[1]])

demo3 <- self$big_screening(domain = demo2[[1]])

demo3 %>%
  unnest() %>%
  filter(time == 740) %>%
  pull(Pore)


# Test IF -----------------------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
self <- VP_proj_creator$new(sourcefile = "file:///D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_deux_elim.r")

self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "dose50",20,"Conc", 0.1, 23
))


self$param_reduce

self$make_filters(cmt = "Conc")

VP_df <- crossing(  k21  = seq(0,1,0.1),
                    k12  = seq(0,1,0.1),
                    ke = seq(0,1,0.1),
                    ke2 = seq(0,1,0.1),
                    V1 = seq(1,10,1)
) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = F)



self$plot_VP(nmax = 1000)



# Non usable param Handle -------------------------------------------------

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
self <- VP_proj_creator$new()
#
# saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/longrun.RDS")

# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter =  cmt == "tumVol",timeforce = c(12,19, 30,45))

self$targets$min <- map2_dbl(self$targets$min, self$targets$max, ~ mean(c(.x, .y))) %>% round

VP_df <- crossing(k2 = seq(0,8,0.1),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.16,0.01),
                  lambda1 = c(10,12,14,33),
                  Vd =  c(0:40)) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

fix_df <- tibble(k1 = seq(0,2,0.5))

t0 <- Sys.time()
self$add_VP(VP_df,fix_df =  fix_df, fillatend = F, reducefilteratend = F,use_green_filter = F, npersalve = 2000, time_compteur = F, pctActivGreen = 0.75)
Sys.time() - t0

self$plot_VP(nmax = 2000)



# Test Lindner normal -----------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",90,"TimeAbove", 1E-1, Inf
))

below <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")


below$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",90,"TimeAbove", -Inf, 1E-1
))

VP_df <- crossing(Bcl20 = seq(100,1000,200),
                  Bclxl0 = seq(100,1000,200),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(100,1000,200),
                  PUMA0 =seq(100,1000,200),
                  NOXA0 =  seq(100,1000,200),
                  BAXc0 = 1000,
                  BAK0 = 0) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )


# VP_df <- crossing(Bcl20 = 500,
#                   Bclxl0 =500,
#                   Mcl10 = 50 ,#*  seq(0.6,1.4,0.2),
#                   BIM0 = 500,
#                   PUMA0 = 500,
#                   NOXA0 =  seq(100,1000,100),
#                   BAXc0 = 1000,
#                   BAK0 = 1000) %>% #c(0.8,1,1.2)) %>%
#   map_df(function(x){
#
#     if(is.character(x)) return(x)
#     round(x,3)
#
#   } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore")
saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/death.RDS")

self$poolVP

below$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore")


self$p

self$poolVP %>%
  filter(TimeAbove_BU !=0) %>%
  slice(1:500) %>%
    unnest() %>%
  ggplot()+
  geom_line(aes(time, Pore, group =  factor(id)))+
  geom_hline(yintercept = 10, col = "red")

# Test Lindner normal vs 1 drug -----------------------------------------------------



  source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 50))


self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$protocols <- protocols



self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",300,"TimeAbove", -Inf, 1E-1,
                                  "apog",300,"TimeAbove", 1E-1, Inf
))


VP_df <- crossing(Bcl20 = seq(100,1000,200),
                  Bclxl0 = seq(100,1000,200),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(100,1000,200),
                  PUMA0 =seq(100,1000,200),
                  NOXA0 =  seq(100,1000,200),
                  BAXc0 = 1000,
                  BAK0 = 1000) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )


self$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore")

self$poolVP %>%
  arrange(id) %>%
  unnest(simul) %>%
  ggplot()+
  geom_line(aes(time, Pore, group = id))+
  facet_wrap(~protocol)+
  geom_hline(yintercept = 10, col = "red")+
  coord_cartesian(xlim = c(62,72))


self$poolVP %>%
  unnest(simul)

# Test Lindner normal vs 1 drug -----------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 50),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 25))


self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$protocols <- protocols



self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",300,"TimeAbove", -Inf, 1E-1,
                                  "apog",300,"TimeAbove", 1E-1, Inf,
                                  "ABT737",300,"TimeAbove", -Inf, 1E-1
))


VP_df <- crossing(Bcl20 = seq(100,1000,200),
                  Bclxl0 = seq(100,1000,200),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(100,1000,200),
                  PUMA0 =seq(100,1000,200),
                  NOXA0 =  seq(100,1000,200),
                  BAXc0 = 1000,
                  BAK0 = 1000) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )


self$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore")

self$poolVP %>%
  arrange(id) %>%
  unnest(simul) %>%
  mutate(protocol = if_else(protocol == "unique", "No Drug" , protocol)) %>%
  ggplot()+
  geom_line(aes(time, Pore, group = id))+
  facet_wrap(~protocol)+
  geom_hline(yintercept = 10, col = "red")+
  coord_cartesian(xlim = c(62,72))


# Test Lindner normal vs 1 drug -----------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 50),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 25))


self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$protocols <- protocols



self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",300,"TimeAbove", -Inf, 1E-1,
                                  "apog",300,"TimeAbove", -Inf, 1E-1,
                                  "ABT737",300,"TimeAbove", 1E-1, Inf
))


VP_df <- crossing(Bcl20 = seq(100,1000,200),
                  Bclxl0 = seq(100,1000,200),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(100,1000,200),
                  PUMA0 =seq(100,1000,200),
                  NOXA0 =  seq(100,1000,200),
                  BAXc0 = 1000,
                  BAK0 = 1000) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )


self$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore")

self$compute_zone_maybe()

self$poolVP %>%
  arrange(id) %>%
  unnest(simul) %>%
  ggplot()+
  geom_line(aes(time, Pore, group = id))+
  facet_wrap(~protocol)+
  geom_hline(yintercept = 10, col = "red")+
  coord_cartesian(xlim = c(62,72))


