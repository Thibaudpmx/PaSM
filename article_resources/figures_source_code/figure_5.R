## ---------------------------
## Script name: Figure R article
##
## Purpose of script: production of figure 6
##c
## Author: Thibaud Derippe
##
## Date Created: 2022-04-18 (trial codes previously in "eval_performance" file)
##
## Under GPL-3 License
## Email: thibaud.derippe@gmail.com
## GitHub: https://github.com/Thibaudpmx/PaSM
## ---------------------------
##
##
##
##
##
## ---------------------------

library(PaSM)



# Several time points --------------------------------------------------------

# source("D:/these/Second_project/QSP/PaSM/R/R6object.R")
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
                          w0 = 50,
                          Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

twozone$add_VP(VP_df_twozone, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
above_or_below$add_VP(VP_df_twozone, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
above_or_below2$add_VP(VP_df_twozone, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)

twozone$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T)


subplotA <- above_or_below2$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T)+ ggtitle("Time 12 only"); subplotA
subplotB <- above_or_below$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T) + ggtitle("Time 45 only"); subplotB
subplotC <-  twozone$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T) + ggtitle("Time 12 and 45"); subplotC
  scale_y_continuous(limits = c(0,0.2));



# Computa patient



line <- VP_df_twozone %>%
  rowid_to_column("id") %>%
  mutate(protocol = "dose50")


protocol <-  line %>%
  mutate(protocol2 = map(protocol, ~ above_or_below2$protocols[[.x]])) %>%
  select(id, protocol2) %>%
  unnest(protocol2)

res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                    icmt = twozone$initial_cmt_values, time_vec =twozone$times,
                    pardf = twozone$parameters_default_values, model = twozone$model) %>%
  as_tibble()#;res


# Id above and below

res %>%
  filter(time == 45 & tumVol > 300 ) %>% pull(id) -> idtem

res %>%
  filter(id %in% idtem) %>%
  filter(time == 12 & tumVol < 21.3) %>% pull(id) %>% unique() -> idbothbelowabove



res <- res %>%
  mutate(abovebelow = if_else(id == idbothbelowabove, T, F))


res %>%
  filter(time == 45 & (tumVol > 300 | tumVol <  50.3)) %>%
  mutate(forwrap = if_else(tumVol > 300, "Strictly above", "Strictly below")) %>%
  select(id, forwrap) -> idtemp45

res %>%
  filter(time == 12 & (tumVol > 92.98855 | tumVol <   21.34172)) %>%
  mutate(forwrap = if_else(tumVol > 92.98855, "Strictly above", "Strictly below")) %>%
  distinct(id, forwrap) -> idtemp12



res2 <- res %>%
  left_join(idtemp45) %>%
  rename(t45 = forwrap) %>%
  mutate(t45 = if_else(is.na(t45), "Accepted", t45)) %>%
  left_join(idtemp12) %>%
  rename(t12 = forwrap) %>%
  mutate(t12 = if_else(is.na(t12), "Accepted", t12))


res2 %>%
  mutate(finalwrap = case_when( t45 == "Strictly above" & t12 == "Strictly above" ~ "Both Above",
                                (t45 == "Strictly above" & t12 == "Strictly below")|(t12 == "Strictly above" & t45 == "Strictly below") ~ "Above and Below",
                                t45 == "Strictly below" & t12 == "Strictly below" ~ "Both Below",
                                t45 == "Accepted" & t12 == "Accepted" ~ "Accepted both",
                                t45 == "Accepted" & t12 != "Accepted" ~ "Accepted T45 only",
                                T ~"Accepted T12 only"
                                )) %>%
  # mutate(finalwrap = if_else(finalwrap %in% c("Both Below", "Both Above" ), "Both Above or Below", finalwrap)) %>%
  # filter(finalwrap == "Above and Below") %>% filter(time %in% c(12,45)) %>% filter(tumVol > 100 &  time == 12)
  ggplot()+
  # geom_point(aes(tume, tumVol, group = id))+
  geom_line(aes(time, tumVol, group = id))+
  scale_y_log10()+
  facet_wrap(~finalwrap)+
  geom_segment(data = above_or_below$targets, aes(x = time, xend =  time, y = min, yend = max), col = "red", size = 2)+
  geom_segment(data = above_or_below2$targets, aes(x = time, xend =  time, y = min, yend = max), col = "red", size = 2)+
  theme_bw() -> subplotD;subplotD

plotSeveralTimepoint <- plot_grid(subplotA, subplotB,subplotD,  subplotC, nrow = 1, labels = letters); plotSeveralTimepoint





# Several OoI -------------------------------------------------------------

# Several YTYPE --------------------------------------------------------

# source("D:/these/Second_project/QSP/PaSM/R/R6object.R")
pk_only <- VP_proj_creator$new()
pd_only <- VP_proj_creator$new()
pk_pd <- VP_proj_creator$new()

pk_pd$set_targets(filter = Dose==50 , timeforce = c(2, 20))
pk_pd$targets <- pk_pd$targets %>% ungroup %>%  slice(2,3 )
pk_pd$targets$min[[2]] <- 0.02
pk_pd$targets$max[[1]] <- 80

pk_only$targets <- pk_pd$targets %>% slice(2)

pd_only$targets <- pk_pd$targets %>% slice(1)

VP_df_pkpd <- crossing(k1 = c(0.5),
                       k2 = seq(0,3,0.2),
                       ke = seq(0,4,0.2) ,#*  seq(0.6,1.4,0.2),
                       lambda0 = 0.05,
                       lambda1 = c(12),
                       w0 = 50,
                       Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )



pk_only$add_VP(VP_df_pkpd, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
pd_only$add_VP(VP_df_pkpd, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)
pk_pd$add_VP(VP_df_pkpd, fillatend = F, reducefilteratend = T,  npersalve = 2000,  time_compteur = F)



line <- VP_df_pkpd %>%
  rowid_to_column("id") %>%
  mutate(protocol = "dose50")


protocol <-  line %>%
  mutate(protocol2 = map(protocol, ~ pk_only$protocols[[.x]])) %>%
  select(id, protocol2) %>%
  unnest(protocol2)

res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                    icmt = twozone$initial_cmt_values, time_vec =seq(0,25,0.1),
                    pardf = twozone$parameters_default_values, model = twozone$model) %>%
  as_tibble()#;res

res %>%
  filter(time == 20)
  ggplot()+
  geom_line(aes(time, tumVol, group = id))

res %>%
  filter(time ==  pk_only$targets$time ) %>%
  mutate(PK = case_when(Conc < pk_only$targets$min ~ "Below",
                      Conc > pk_only$targets$max ~ "Above",
                      T ~ "Accepted")) %>%
  distinct(id, PK)-> PKID


res %>%
  filter(time == pd_only$targets$time ) %>%
  mutate(PD = case_when(tumVol < pd_only$targets$min ~ "Below",
                        tumVol > pd_only$targets$max ~ "Above",
                        T ~ "Accepted")) %>%
  distinct(id, PD)-> PDID



res %>%
  left_join(PKID) %>%
  left_join(PDID) %>%
  mutate(PKPD = case_when(PK == "Accepted" & PD == "Accepted" ~ "Accepted",
                          PK != "Accepted" & PD == "Accepted" ~ "PK rejection",
                          PK == "Accepted" & PD != "Accepted" ~ "PD rejection",
                          T ~ "both rejection")) %>%
  mutate(straightforward = if_else(PK == "Accepted" & PD == "Accepted", "Accepted", "Rejected")) %>%
  gather("cmt", "value", Conc, tumVol) %>%
  filter(!(time > 5 & cmt == "Conc")) %>%
  ggplot()+
  geom_line(aes(time, value, group = id, col = straightforward), alpha = 0.3)+
  geom_segment(data = pd_only$targets, aes(x = 0, xend = 0, y = 0, yend= 0, col = "Accepted"))+ # just to have legend with alpha1
  geom_segment(data = pd_only$targets, aes(x = time, xend = time, y = min, yend= max), col = "red", size = 2)+
  geom_segment(data = pk_only$targets, aes(x = time, xend = time, y = min, yend= max), col = "red", size = 2)+

  facet_wrap(~cmt, scales = "free", ncol = 1)+
  scale_y_log10()+
  theme_bw() +
  labs(x = "Time", y = "Output of Interest", col = NULL) +
  theme(legend.position = c(0.2,0.2))-> pdplot;pdplot

plotSeveralYTYPE <- plot_grid(


  pk_only$plot_2D(x = k2, y = ke, plotoreturn = 1,add_point = T)+ ggtitle("PK only"),
  pd_only$plot_2D(x = k2, y = ke, plotoreturn = 1,add_point = T)+ ggtitle("PD only"),
  pdplot,
  pk_pd$plot_2D(x = k2, y = ke, plotoreturn = 1,add_point = T)+ ggtitle("PK-PD"), nrow = 1,labels = letters[-(1:4)]


)




# Save 300 dip for article
tiff(width = 4500, height = 2300,filename = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig5.tiff", res = 300)
plot_grid(plotSeveralTimepoint , plotSeveralYTYPE,  ncol = 1)
dev.off()
shell.exec( "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/fig5.tiff")


