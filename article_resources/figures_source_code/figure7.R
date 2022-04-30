## ---------------------------
## Script name: Figure R article
##
## Purpose of script: production of figure 7
##c
## Author: Thibaud Derippe
##
## Date Created: 2022-04-18 (trial codes previously in "eval_performance" file)
##
## Under GPL-3 License
## Email: thibaud.derippe@gmail.com
## GitHub: https://github.com/Thibaudpmx/QSPVP
## ---------------------------
##
##
##
##
##
## ---------------------------

library(QSPVP)

cohort_creator_Lindner <- function(nmodif){

  base <- crossing( Bcl20 = 500,
                    Bclxl0 = 200,
                    Mcl10 = 50,
                    BIM0 =  200,
                    PUMA0 = 200,
                    NOXA0 = 200 #*  seq(0.6,1.4,0.2),
  ) %>% #c(0.8,1,1.2)) %>%
    map_df(function(x){

      if(is.character(x)) return(x)
      round(x,3)

    } )


  nperparam <- ceiling(200000^(1/nmodif))

  list <-   map(1:nmodif, function(x){

    min <- base[[x]]/10
    max <- base[[x]]*10
    step <- (max - min)/(nperparam-1)

    expr(seq( !!min,  !!max, !!step  ))
  }
  )

  names(list) <- names(base[1:nmodif])

  if(nmodif < length(base)){

    output <- crossing(base[-c(1:nmodif)], !!!list)

  }else{


    output <- crossing(!!!list)
  }


  output %>%
    slice(1:200000) %>%
    mutate(BAXc0 = 1000, BAK0 = 1000)


}

self <-   VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",90,"Pore", 10, 20

))

self$times <- c(0, 50:90)
cohort <-  cohort_creator_Lindner(6)

self$add_VP(VP_df = cohort, use_green_filter = T, npersalve = 2000, time_compteur = T, keep = "Pore", timeSave = 50:80)




plotLindner <- function(obj, npatient = 2000, title = ""){


  obj$poolVP %>%
  arrange(id) %>%
    mutate(test = map_lgl(simul, ~class(.x)[[1]] == "tbl_df")) %>%
    filter(test) %>%
    mutate(protocol = case_when(protocol == "ABT737" ~ "ABT-737",
                                protocol == "apog" ~ "ApoG2",
                                T~ protocol))-> temp


  temp <- temp %>% left_join(temp %>%
                               group_by(id) %>%
                               tally )

  idtosample <- temp$id[temp$n == max(temp$n)] %>% unique


  temp %>%
    filter(id%in% sample(idtosample, min(npatient, length(idtosample))) ) %>%
    # slice(1:(min(npatient * obj$targets$protocol %>% unique %>% length(), nrow(temp)))) %>%
    unnest(simul) %>%
    mutate(protocol = if_else(protocol == "unique", "No Drug" , protocol)) %>%
    ggplot()+
    geom_line(aes(time, Pore, group = id))+
    facet_wrap(~protocol)+
    geom_hline(yintercept = 10, col = "red", size = 1.5)+
    coord_cartesian(xlim = c(50,100))+
    theme_bw()+
    labs(x = "Time (hours)", y = "Percentage Pore", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}
# Test Lindner normal ABT737-----------------------------------------------------

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Lindner_article")



protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 50),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 50))


protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 40),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 4))



self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$protocols <- protocols


self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",100,"TimeAbove", -Inf, 1E-1,
                                  "apog",100,"TimeAbove", -Inf, 1E-1,
                                  "ABT737",100,"TimeAbove", 1E-1, Inf
))



VP_df <- crossing(Bcl20 = seq(0,900,100),
                  Bclxl0 = seq(0,900,100),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(0,900,100),
                  PUMA0 =seq(0,900,100),
                  NOXA0 =  seq(0,900,200),
                  BAXc0 = 500,
                  BAK0 = 500) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )


self$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore")

self$n_filter_reduc()

self$poolVP %>%
  arrange(id) %>%
  unnest(simul) %>%
  mutate(protocol = if_else(protocol == "unique", "No Drug" , protocol)) %>%
  ggplot()+
  geom_line(aes(time, Pore, group = id))+
  facet_wrap(~protocol)+
  geom_hline(yintercept = 10, col = "red")+
  coord_cartesian(xlim = c(62,72))

saveRDS(self,"Apo50_ABT25_APOGW.RDS")


# Test Lindner normal APOG winner -----------------------------------------------------



protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 40),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 4))


self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$protocols <- protocols



self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",100,"TimeAbove", -Inf, 1E-1,
                                  "apog",100,"TimeAbove", 1E-1, Inf,
                                  "ABT737",100,"TimeAbove", -Inf, 1E-1
))


VP_df <- crossing(Bcl20 = seq(0,900,100),
                  Bclxl0 = seq(0,900,100),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(0,900,100),
                  PUMA0 =seq(0,900,100),
                  NOXA0 =  seq(0,900,200),
                  BAXc0 = 500,
                  BAK0 = 500) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )


self$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore", time_compteur = T)
selfAPOG <- self

selfAPOG$poolVP %>%
  arrange(id) %>%
  unnest(simul) %>%
  mutate(protocol = if_else(protocol == "unique", "No Drug" , protocol)) %>%
  ggplot()+
  geom_line(aes(time, Pore, group = id))+
  facet_wrap(~protocol)+
  geom_hline(yintercept = 10, col = "red")+
  coord_cartesian(xlim = c(62,72), ylim = c(7.5,12))

self$timeTrack$poolVP_compteur[grepl("green", tolower(names(self$timeTrack$poolVP_compteur)))]

timetable(self)



saveRDS(selfAPOG, "80_2_APOGW.RDS")



# Full analysis -----------------------------------------------------------



VP_df <- crossing(Bcl20 = seq(0,900,100),
                  Bclxl0 = seq(0,900,100),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(0,900,100),
                  PUMA0 =seq(0,900,100),
                  NOXA0 =  seq(0,900,200),
                  BAXc0 = 500,
                  BAK0 = 500) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 80),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 40))


# Step1: Die anyways ------------------------------------------------------


# Step 1: remove those who die anyway






ApoptoWODrug <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

ApoptoWODrug$protocols <- protocols



ApoptoWODrug$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",100,"TimeAbove",  1E-1, Inf
))



ApoptoWODrug$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore", time_compteur = T, keepRedFiltaftDis = T)

saveRDS(ApoptoWODrug, "ApoptoWODrug")



VP_df %>%
  left_join(ApoptoWODrug$poolVP %>% select(-simul)) %>%
  filter(is.na(id)) %>%
  select(1:8) -> VP_df2



# Step2: Dead Both drugs ------------------------------------------------------
protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 2000),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 2000))


BothDead <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

BothDead$protocols <- protocols


BothDead$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  # "unique",100,"TimeAbove", -Inf, 1E-1,
                                  "apog",100,"TimeAbove", 1E-1, Inf,
                                  "ABT737",100,"TimeAbove",  1E-1, Inf
))

BothDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore", time_compteur = T, keepRedFiltaftDis = T)

plotLindner(obj = BothDead)

saveRDS(BothDead, "BothDead")

PrevRej <-

  BothDead$VP_rejected %>%
  filter(!is.na(id)) %>%
    rename(TimeAbove = value) %>%
    select(time, protocol, TimeAbove, !!!BothDead$param)

PrevAcc <- BothDead$poolVP %>%
  unnest() %>%
  filter(time == 100) %>%
  select(time, protocol, TimeAbove, !!!BothDead$param)



Previous <- PrevAcc %>%
  bind_rows(
    PrevRej
  ) %>%
  rowid_to_column("id")



BothDead$timeTrack$poolVP_compteur$nsimul %>% sum() # why not equal? Oh, because the one accepted for a proto but not the other...



# Step3: Never good... ------------------------------------------------------

protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 2000),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 2000),
                   synergy = tibble(cmt = c("ABT737","ApoG2"), time = 50, amt = 2000))



NeverDead <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

NeverDead$protocols <- protocols


NeverDead$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                      # "unique",100,"TimeAbove", -Inf, 1E-1,
                                      "apog",100,"TimeAbove", -Inf, 1E-1,
                                      "ABT737",100,"TimeAbove", -Inf, 1E-1
))

#
# NeverDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
#                 time_compteur = T, keepRedFiltaftDis = T, PreviousResults = Previous)


NeverDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
                 time_compteur = T, keepRedFiltaftDis = T)




saveRDS(NeverDead, "NeverDead")

plotLindner(obj = NeverDead)



PrevRej <-

  NeverDead$VP_rejected %>%
  filter(!is.na(id)) %>%
  rename(TimeAbove = value) %>%
  select(time, protocol, TimeAbove, !!!NeverDead$param)

PrevAcc <- NeverDead$poolVP %>%
  unnest() %>%
  filter(time == 100) %>%
  select(time, protocol, TimeAbove, !!!NeverDead$param)


Previous <-  bind_rows(Previous, PrevAcc, PrevRej) %>%
    distinct()


# Step 4 Only ApoG ------------------------------------------------------------------



apoG <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

apoG$protocols <- protocols


apoG$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                       # "unique",100,"TimeAbove", -Inf, 1E-1,
                                       "apog",100,"TimeAbove",  1E-1, Inf,
                                       "ABT737",100,"TimeAbove", -Inf, 1E-1
))




apoG$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
            time_compteur = T, keepRedFiltaftDis = T)


apoG$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
                 time_compteur = T, keepRedFiltaftDis = T, PreviousResults = Previous)


apoG

saveRDS(apoG, "apoG")
# Step 5 Only ABT737 ------------------------------------------------------------------



ABT737 <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

ABT737$protocols <- protocols


ABT737$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  # "unique",100,"TimeAbove", -Inf, 1E-1,
                                  "apog",100,"TimeAbove",  -Inf, 1E-1,
                                  "ABT737",100,"TimeAbove",  1E-1, Inf
))


ABT737$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
              time_compteur = T, keepRedFiltaftDis = T)

ABT737$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
            time_compteur = T, keepRedFiltaftDis = T, PreviousResults = Previous)


plotLindner(obj = ABT737)
saveRDS(ABT737, "ABT737")
# Sum to verify



# Step 6 Synergi ------------------------------------------------------------------



synergy <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

synergy$protocols <- protocols


synergy$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  # "unique",100,"TimeAbove", -Inf, 1E-1,
                                  "apog",100,"TimeAbove",  -Inf, 1E-1,
                                  "ABT737",100,"TimeAbove", -Inf, 1E-1,
                                  "synergy",100,"TimeAbove", 1E-1, Inf
))




synergy$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
            time_compteur = T, keepRedFiltaftDis = T)



saveRDS(synergy, "synergy")

# Step 6 Synergi from start------------------------------------------------------------------



synergy <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

synergy$protocols <- protocols


synergy$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                     "unique",100,"TimeAbove", -Inf, 1E-1,
                                     "apog",100,"TimeAbove",  -Inf, 1E-1,
                                     "ABT737",100,"TimeAbove", -Inf, 1E-1,
                                     "synergy",100,"TimeAbove", 1E-1, Inf
))




synergy$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
               time_compteur = T, keepRedFiltaftDis = T)


plotLindner(synergy)
saveRDS(synergy, "synergy")



# Final -------------------------------------------------------------------



NApoptose <- ApoptoWODrug$poolVP %>% nrow
NBothDead <- BothDead$poolVP %>% nrow / 2
NNeverDead <- NeverDead$poolVP %>% nrow /2
NapoG <- apoG$poolVP %>% nrow /2
NABT737 <- ABT737$poolVP %>% nrow /2

NApoptose + NBothDead + NNeverDead + NapoG + NABT737


nsim <- 100

set.seed(23)
cowplot::plot_grid(

  plotLindner(ApoptoWODrug, npatient = nsim, title = paste0("Spontaneous Apoptosis (", nrow(ApoptoWODrug$poolVP)," VPs)")),
  plotLindner(NeverDead, npatient = nsim, title = paste0("Both ineffective drugs (", nrow(NeverDead$poolVP)/2," VPs)")),
  plotLindner(BothDead, npatient = nsim, title = paste0("Both effective drugs  (", nrow(BothDead$poolVP)/2," VPs)")),
  plotLindner(apoG, npatient = nsim, title = paste0("Effective ApoG2 only (", nrow(apoG$poolVP)/2," VPs)")) ,
  plotLindner(ABT737, npatient = nsim, title = paste0("Effective ABT-373 only (", nrow(ABT737$poolVP)/2," VPs)")),
  plotLindner(synergy, npatient = nsim, title = paste0("Effective synergy only (", nrow(synergy$poolVP)/3," VPs)")),
  labels = LETTERS

)

# Total time

(ApoptoWODrug$timeTrack$tTOTAL +
  BothDead$timeTrack$tTOTAL +
  ABT737$timeTrack$tTOTAL +
  NeverDead$timeTrack$tTOTAL +
  apoG$timeTrack$tTOTAL) / 60



