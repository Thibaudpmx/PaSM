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



plotLindner <- function(obj, npatient = 2000, title = ""){


  obj$poolVP %>%
  arrange(id) %>%
    mutate(test = map_lgl(simul, ~class(.x)[[1]] == "tbl_df")) %>%
    filter(test) %>%
    mutate(protocol = case_when(protocol == "ABT737" ~ "ABT-737",
                                protocol == "apog" ~ "ApoG2",
                                protocol == "synergy" ~ "Combination",
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
# original configuration-----------------------------------------------------

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Lindner_article")



protocols <- list( unique = tibble(cmt = "Bcl2", time = 0, amt = 0),
                   apog = tibble(cmt = "ApoG2", time = 50, amt = 2000),
                   ABT737 = tibble(cmt = "ABT737", time = 50, amt = 2000),
                   synergy = tibble(cmt = c("ABT737","ApoG2"), time = 50, amt = 2000))

VP_df <- crossing(Bcl20 = seq(0,900,100),
                  Bclxl0 = seq(0,900,100),
                  Mcl10 = c(0,50,100,150) ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(0,900,100),
                  PUMA0 =seq(0,900,100),
                  NOXA0 =  seq(0,900,200),
                  BAXc0 = 500,
                  BAK0 = 500)

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


BothDead <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

BothDead$protocols <- protocols

targetsBothDead <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
        "unique",100,"TimeAbove", -Inf, 1E-1,
        "apog",100,"TimeAbove", 1E-1, Inf,
        "ABT737",100,"TimeAbove",  1E-1, Inf
)

BothDeadFull <- BothDead$clone(deep = T)

BothDead$set_targets(manual = targetsBothDead %>% slice(-1))
BothDeadFull$set_targets(manual = targetsBothDead)

BothDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore", time_compteur = T, keepRedFiltaftDis = T)
if(nrow(BothDead$poolVP) >1) saveRDS(BothDead, "BothDead")
# BothDead <- readRDS("BothDead")

BothDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore", time_compteur = T, keepRedFiltaftDis = T)
if(nrow(BothDead$poolVP) >1) saveRDS(BothDead, "BothDead")


BothDeadFull$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore", time_compteur = T, keepRedFiltaftDis = T)
if(nrow(BothDeadFull$poolVP) >1) saveRDS(BothDeadFull, "BothDeadFull")

plotLindner(obj = BothDead)



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

saveRDS(Previous, "previousBothdead")

BothDead$timeTrack$poolVP_compteur$nsimul %>% sum() # why not equal? Oh, because the one accepted for a proto but not the other...


# Step 3 Only ApoG ------------------------------------------------------------------



apoG <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

apoG$protocols <- protocols
apoGFull <- apoG$clone()

protoApoG <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
        "unique",100,"TimeAbove", -Inf, 1E-1,
        "apog",100,"TimeAbove",  1E-1, Inf,
        "ABT737",100,"TimeAbove", -Inf, 1E-1
)

apoG$set_targets(manual = protoApoG %>% slice(-1))
apoGNoPrev <- apoG$clone(deep = T)
apoGFull$set_targets(manual = protoApoG%>% slice(-1))



apoG$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
            time_compteur = T, keepRedFiltaftDis = T, PreviousResults = Previous)



apoGNoPrev$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
            time_compteur = T, keepRedFiltaftDis = T)


if(nrow(apoG$poolVP) > 1 ) saveRDS(apoG, "apoG")

if(nrow(apoGNoPrev$poolVP) > 1 ) saveRDS(apoGNoPrev, "apoGNoPrev")



apoGFull$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
            time_compteur = T, keepRedFiltaftDis = T)

if(nrow(apoGFull$poolVP) > 1 ) saveRDS(apoGFull, "apoGFull")

apoGFull$timeTrack$tTOTAL/60



PrevRej <-

  apoG$VP_rejected %>%
  filter(!is.na(id)) %>%
  rename(TimeAbove = value) %>%
  select(time, protocol, TimeAbove, !!!NeverDead$param)

PrevAcc <- apoG$poolVP %>%
  unnest() %>%
  filter(time == 100) %>%
  select(time, protocol, TimeAbove, !!!NeverDead$param)


Previous <-  bind_rows(readRDS( "previousBothdead"), PrevAcc, PrevRej) %>%
  distinct()

saveRDS(Previous, "previousPostApoG")


# Step 4 Only ABT737 ------------------------------------------------------------------



ABT737 <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

ABT737$protocols <- protocols

ABT737Full <- ABT737$clone(deep = T)
ABT737targets <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                         "unique",100,"TimeAbove", -Inf, 1E-1,
                         "apog",100,"TimeAbove",  -Inf, 1E-1,
                         "ABT737",100,"TimeAbove",  1E-1, Inf
)

ABT737$set_targets(manual = ABT737targets %>% slice(-1) )
ABT737NoPrev <- ABT737$clone(deep = T)
ABT737Full$set_targets(manual = ABT737targets %>% slice(-1))



ABT737$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
              time_compteur = T, keepRedFiltaftDis = T, PreviousResults = Previous)


if(nrow(ABT737$poolVP) >1) saveRDS(ABT737, "ABT737")

ABT737NoPrev$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
              time_compteur = T, keepRedFiltaftDis = T)

if(nrow(ABT737NoPrev$poolVP) >1) saveRDS(ABT737NoPrev, "ABT737NoPrev")


ABT737Full$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
              time_compteur = T, keepRedFiltaftDis = T)

if(nrow(ABT737Full$poolVP) >1) saveRDS(ABT737Full, "ABT737Full")
ABT737Full$timeTrack$tTOTAL / 60

plotLindner(obj = ABT737)
# saveRDS(ABT737, "ABT737")
# Sum to verify

# Cmopute last df ---------------------------------------------------------

BothDead <- readRDS("BothDead")
OnlyApoG <- readRDS( "apoG")
onlyABT737 <- readRDS("ABT737")


ResistantsMono <- VP_df2 %>%
  left_join(BothDead$poolVP) %>%
  filter(is.na(protocol)) %>%
  select(1:8) %>%
  left_join(OnlyApoG$poolVP) %>%
  filter(is.na(protocol)) %>%
  select(1:8) %>%
  left_join(onlyABT737$poolVP) %>%
  filter(is.na(protocol)) %>%
  select(1:8)

saveRDS(ResistantsMono, "ResistantsMono")


# Step 6 Synergy ------------------------------------------------------------------



synergy <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

synergy$protocols <- protocols
synergyFull <- synergy$clone(deep = T)

protocolSynergy <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                           "unique",100,"TimeAbove", -Inf, 1E-1,
                           "apog",100,"TimeAbove",  -Inf, 1E-1,
                           "ABT737",100,"TimeAbove", -Inf, 1E-1,
                           "synergy",100,"TimeAbove", 1E-1, Inf
)

synergy$set_targets(manual = protocolSynergy %>% slice(-1))
synergyFull$set_targets(manual = protocolSynergy %>% slice(-1))



synergy$add_VP(ResistantsMono, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
               time_compteur = T, keepRedFiltaftDis = T)


if(nrow(synergy$poolVP) >1) saveRDS(synergy, "synergy")


synergyFull$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
               time_compteur = T, keepRedFiltaftDis = T)

if(nrow(synergyFull$poolVP) >1) saveRDS(synergyFull, "synergy")

synergy$timeTrack$tTOTAL / 60
synergyFull$timeTrack$tTOTAL / 60

# Step7: Single treatment therapy failure ------------------------------------------------------



NeverDead <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

NeverDead$protocols <- protocols
NeverDeadFull <- NeverDead$clone(deep = T)

ProtNeverDead <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                     "unique",100,"TimeAbove", -Inf, 1E-1,
                     "apog",100,"TimeAbove", -Inf, 1E-1,
                     "ABT737",100,"TimeAbove", -Inf, 1E-1
)

NeverDead$set_targets(manual = ProtNeverDead %>% slice(-1) )
NeverDeadFull$set_targets(manual = ProtNeverDead %>% slice(-1) )
#
NeverDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
                time_compteur = T, keepRedFiltaftDis = T, PreviousResults = Previous)


# NeverDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
#                  time_compteur = T, keepRedFiltaftDis = T)


NeverDeadFull$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
                 time_compteur = T, keepRedFiltaftDis = T)


if(nrow(NeverDead$poolVP) >1) saveRDS(NeverDead, "NeverDead")
if(nrow(NeverDeadFull$poolVP) >1) saveRDS(NeverDead, "NeverDeadFull")
NeverDeadFull$timeTrack$tTOTAL/60
plotLindner(obj = NeverDead)



# Step7: Combination treatment failure  ------------------------------------------------------



CombinationResist <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

CombinationResist$protocols <- protocols

protocolCombinationResist <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                           "unique",100,"TimeAbove", -Inf, 1E-1,
                           "apog",100,"TimeAbove",  -Inf, 1E-1,
                           "ABT737",100,"TimeAbove", -Inf, 1E-1,
                           "synergy",100,"TimeAbove",  -Inf, 1E-1
)

CombinationResist$set_targets(manual =protocolCombinationResist %>% slice(4))
#
CombinationResist$add_VP(VP_df, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
                 time_compteur = T, keepRedFiltaftDis = T)

saveRDS(CombinationResist, "CombinationResist")
# NeverDead$add_VP(VP_df2, fillatend = F, reducefilteratend = F, use_green_filter = T, keep = "Pore",
#                  time_compteur = T, keepRedFiltaftDis = T)

# Final Plot supplement-------------------------------------------------------------------
ApoptoWODrug <- readRDS("ApoptoWODrug")
ABT737 <- readRDS("ABT737")
apoG <- readRDS("apoG")
BothDead <- readRDS("BothDead")
Synergy <- readRDS("Synergy")
CombinationResist <- readRDS("CombinationResist")



nsim <- 100

set.seed(23)

# plotLindner(NeverDead, npatient = nsim, title = paste0("Both ineffective drugs (", nrow(NeverDead$poolVP)/2," VPs)")),

cowplot::plot_grid(

  plotLindner(ApoptoWODrug, npatient = nsim, title = paste0("Spontaneous Apoptosis (", length(unique(ApoptoWODrug$poolVP$id))," VPs)")),
   plotLindner(BothDead, npatient = nsim, title = paste0("Both effective drugs  (", length(unique(BothDead$poolVP$id))," VPs)")),
  plotLindner(apoG, npatient = nsim, title = paste0("Effective ApoG2 only (",length(unique(apoG$poolVP$id))  ," VPs)")) ,
  plotLindner(ABT737, npatient = nsim, title = paste0("Effective ABT-373 only (", length(unique(ABT737$poolVP$id))," VPs)")),
  plotLindner(Synergy, npatient = nsim, title = paste0("Effective combination only (", length(unique(Synergy$poolVP$id)) ," VPs)")),
  plotLindner(CombinationResist, npatient = nsim, title = paste0("Combination resistance (", length(unique(CombinationResist$poolVP$id)) ," VPs)")),
  labels = LETTERS

)



