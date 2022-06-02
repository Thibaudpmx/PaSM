## ---------------------------
## Script name: Figure R article
##
## Purpose of script: production of figure 3
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

# Data generation  VP to seek---------------------------------------------------------

# Please provide the folder where you want the generated data to be stored
root <- "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data"

new_wd <- file.path(root, "algo2")
if(!file.exists(new_wd))  dir.create(new_wd)


# Step 1 create the hide and seek VP

self <- VP_proj_creator$new()


VP <- tibble(k2 = 1.63, lambda0 = 0.85, ke = 0.7, Vd = 37, lambda1 = 72, w0 = 50, k1 = 0.5, id = 1:3, psi = 20)

events <- tibble(cmt = "Central", time = 0, amt = c(0,50,100), evid = 1, id = 1:3) %>%
          bind_rows(tibble(cmt = "tumVol", time = 0:52 , amt = 0, evid =  0) %>% crossing(id = 1:3)) %>%
  arrange(id, time)



## just to verify he is still there
# newVPs %>%
#   filter(k2 == 1.63, lambda0 == 0.85, ke == 0.7, Vd == 37, lambda1 == 72, w0 == 50, k1 == 0.5)
#
# maybe %>% rowid_to_column() %>%
#   filter(k2min<= 1.63, k2max >= 1.63, kemin <= 0.7 & kemax >= 0.7, lambda0min <= 0.85, lambda0max >= 0.85,
#          lambda1min <=72 & lambda1max >= 72, Vdmin <= 37, Vdmax >=37)

####

simul <- self$model$solve(VP, events) %>%
  as_tibble


# Plot 0 ------------------------------------------------------------------



plot0 <- simul %>%
  # filter(time %in% c(8,30)) %>%
  # select(id, time, tumVol) %>%
  mutate(protocol = case_when(id == 1 ~ 0,
                              id == 2 ~ 50,
                              id == 3 ~ 100)) %>%
                              {temppoint <<- .} %>%
  filter(time <=40) %>%
  ggplot()+
  geom_line(aes(time, tumVol, col = factor(protocol))) +
  scale_y_log10()+
  theme_bw()+
  labs( x = "Time (days)", y = "Tumor Volume (mm3)", col = "Dose", shape = "")+
  geom_vline(xintercept = c(8,30), lty = 2)+
  geom_point(data = temppoint %>% filter(time %in% c(8,30)), aes( time, tumVol, shape = "Target")); plot0

simul %>%
  filter(time %in% c(8,30)) %>%
  select(id, time, tumVol) %>%
  mutate(protocol = case_when(id == 1 ~ "dose0",
                              id == 2 ~ "dose50",
                              id == 3 ~ "dose100")) %>%
  select(-id) %>%
  mutate(cmt = "tumVol", min = ceiling(tumVol) - 1, max = floor(tumVol) + 1 ) %>%
  select(protocol, cmt, time, min, max) ->  prototiny


simul %>%
  filter(time %in% c(0, 1, 2)) %>%
  select(id, time, Conc) %>%
  filter(id != 1) %>%
  mutate(protocol = case_when(id == 1 ~ "dose0",
                              id == 2 ~ "dose50",
                              id == 3 ~ "dose100")) %>%
  select(-id) %>%
  mutate(cmt = "Conc", min = Conc *0.95, max = Conc *1.05 ) %>%
  # slice(-1) %>%
  select(protocol, cmt, time, min, max) ->  prototinyPK

  prototinypkpd <- bind_rows( prototiny, prototinyPK)


# Analyze PD only--------------------------------------------


self <- VP_proj_creator$new()
self$set_targets(manual = prototiny)


# selfmtd2 <- VP_proj_creator$new()

# selfmtd2$set_targets(manual = prototiny)


npersalve = 2E5
npersalveFinal = 1E6


domain <- tribble(~param, ~from, ~to, ~by,
                  "k2", 0, 3, 0.01 ,
                  "lambda0", 0, 1.4, 0.01,
                  "ke", 0, 2,0.1,
                  "Vd", 1,40,1,
                  "lambda1", 0,240,1
)
fix <-c(k1 = 0.5, w0 = 50)


filePD <- file.path(root, "algo2", "pdonly.RDS")

self$algo2(domain = domain, fix = fix,npersalve = 2E5, npersalveFinal = 1E6,save_every = 5, file = filePD,method = 1 )

# Note: useless to finish the algorithm, just the first step count
filePDmtd2 <- file.path(root, "algo2", "pdonly_mtd2.RDS")
self$algo2(domain = domain, fix = fix,npersalve = 2E5, npersalveFinal = 1E6,save_every = 5, file = filePDmtd2,method = 2 )

# Analyze PKPD ------------------------------------------------------------


self <- VP_proj_creator$new()
self$set_targets(manual = prototinypkpd)


selfwbdr <- VP_proj_creator$new()
selfwbdr$set_targets(manual = prototinypkpd)

domain <- tribble(~param, ~from, ~to, ~by,
                  "k2", 0, 6, 0.01 ,
                  "lambda0", 0, 3.4, 0.01,
                  "ke", 0, 5,0.1,
                  "Vd", 1,40,0.1,
                  "lambda1", 0,240,0.1
)

# ndomain(domain)




filePKPD <- file.path(root, "algo2", "pkpd.RDS")
self$algo2(domain = domain, fix = fix,npersalve =  npersalve,includeBorderZoom = F, npersalveFinal = npersalveFinal, method = 1, file = filePKPD, save_every = 2)


filePKPDwbdr <- file.path(root, "algo2", "pkpdwbdr.RDS")
self$algo2(domain = domain, fix = fix,npersalve =  npersalve,includeBorderZoom = T, npersalveFinal = npersalveFinal, method = 1, file = filePKPDwbdr, save_every = 2)


# Plot A -------------------------------------------------------------

  PDonly <-  readRDS(filePD)
  PDonlyMtd2 <-  readRDS(filePDmtd2)

  # Just verify it provides the same results
  PDonly$algo2list$tree %>% slice(1)
  PDonlyMtd2$algo2list$tree %>% slice(1)

  # and let's analyse
 PDtree <-  PDonly$algo2list$tree %>% # used for workflow
  mutate(round = map_dbl(Name,~ str_split(.x, "_")[[1]] %>% length))

timetrackFirstAna <- PDonly$algo2list$first_timeTrak



 # Labels
 labels <-  c("Algo1", "Filter\nReduc", "Comput\nzone maybe", "Blocs\nreduction", "other")
 stepname <- factor(x = labels, levels =
                      labels)
 # Time computations
 algo1t <- timetrackFirstAna$round1$FirstAlgo1 %>% as.double()
 filterreduc <- timetrackFirstAna$round1$FilterReduc %>% as.double()
 zonemaybet <- timetrackFirstAna$round1$ZoneMaybe %>% as.double()
 blocreducT <-  timetrackFirstAna$tMaybeReduce%>% as.double()
 diffT <- as.double(PDtree$time[[1]]) - algo1t  -   filterreduc - zonemaybet -  blocreducT

 timeMt2reducemaybe <- PDonlyMtd2$algo2list$first_timeTrak$tMaybeReduce %>% as.double()
 diffTMT2  <- as.double(PDonlyMtd2$algo2list$tree$time[[1]] - timeMt2reducemaybe)

 template <- tibble(step =stepname, time = c(algo1t,filterreduc,zonemaybet,blocreducT, diffT) ,method = "Full") %>%
   bind_rows(

     tibble(step = stepname, time = c(0,0,0,timeMt2reducemaybe, diffTMT2),method = "alt")

   )



 plotA <- template %>%
   ggplot()+
   geom_bar(aes(x = step,y= time, fill = method), stat="identity",col = "black", position = "dodge")+
   geom_text(aes(x = step,y= time +5, label = paste0(round(time,0), "s"), group =  method),position = position_dodge(width = 0.9))+
   labs(y = "Time(sec)")+
   theme_bw(); plotA


 # Plot B -------------------------------------------------------------

 tree <- PDonly$algo2list$tree
 plotB <- tree   %>%
   slice(-1) %>%
   mutate(time = as.double(time)) %>%
   ggplot()+
   geom_histogram(aes(time, fill = factor(after/3)),col = "black")+
   scale_x_log10()+
   theme_bw()+
   labs(x = "Time (sec)", fill = "VP\nfound"); plotB







# Analysis pkpd set -------------------------------------------------------


 PKPD <-  readRDS(filePKPD)




PKPDtree <-  PKPD$algo2list$tree %>% # used for workflow
    mutate(round = map_dbl(Name,~ str_split(.x, "_")[[1]] %>% length))

PKPDtree$time[PKPDtree$round == 2][-1] %>% sum() / 60
PKPDtree$time[PKPDtree$round == 2][-1] %>% min
PKPDtree$time[PKPDtree$round == 2][-1] %>% max

PKPD$poolVP$id %>% unique # why is there so many VPS. (oh I should have relaunch it.....)

PKPDtree %>%
  filter(round == 2) %>%
    # slice(-1) %>%
    filter(after > 0) %>%
    # nrow() # to know the number of blocs
    pull(after ) %>% sum

PKPDtree %>%
  filter(round == 3) %>%
  summarise(sum(time), min(time), max(time))
  pull(time ) %>% sum

PKPDtree$time %>% sum





self$algo2list$tree %>%
  summarise(time = sum(time))

2609.409 / 60

self$algo2list$first %>%
  group_by(blocsPool) %>%
  summarise(sum = sum(temp3)) %>%
  pull(sum) %>% median




plot0 <- simul %>%
  # filter(time %in% c(8,30)) %>%
  # select(id, time, tumVol) %>%
  mutate(protocol = case_when(id == 1 ~ 0,
                              id == 2 ~ 50,
                              id == 3 ~ 100)) %>%

  filter(time <=40) %>%

  gather("OoI", "value", tumVol, Conc) %>%
  filter(!(time >5 &  OoI == "Conc")) %>%
  {temppoint <<- .} %>%
  ggplot()+
  geom_line(aes(time, value, col = factor(protocol))) +
  scale_y_log10()+
  facet_wrap(~OoI, scales = "free", ncol =1)+
  theme_bw()+
  labs( x = "Time (days)", y = "Tumor Volume (mm3)", col = "Dose", shape = "")+
  geom_vline(data = tibble(OoI = "tumVol", x = c(8,30)), aes( xintercept = x), lty = 2)+
  geom_vline(data = tibble(OoI = "Conc", x = c(0,1,2)), aes( xintercept = x), lty = 2)+
  geom_point(data = temppoint %>% filter((time %in% c(8,30) & OoI == "tumVol") | (time %in% c(0,1,2) & OoI != "tumVol" & id != 1) ), aes( time, value, shape = "Target")); plot0


cowplot::plot_grid(plot0,"3",  plotA, plotB, labels = LETTERS)

cowplot::plot_grid(plot0,"3",  plotA, plotB,  "3","4", labels = c("A", "B", "C", "D", "", "E"))
# Backup ------------------------------------------------------------------

# file <- "D:/these/Second_project/QSP/modeling_work/VT_simeoni/fig5_data2.RDS"
#
# self$algo2list <- test$algo2list
# self$poolVP <- test$poolVP
# self$targets <- test$targets
#
# self <- test
# test <- readRDS(file)
# test$algo2list
# # Just to confirm after first analyse there is still a bloc containing our VPs
# maybeFinal %>%
#   filter(k2min <= VP$k2[[1]], k2max >= VP$k2[[1]],
#          kemin <= VP$ke[[1]], kemax >= VP$ke[[1]],
#          Vdmin <= VP$Vd[[1]], Vdmax >= VP$Vd[[1]],
#          lambda0min <= VP$lambda0[[1]], lambda0max >= VP$lambda0[[1]],
#          lambda1min <= VP$lambda1[[1]], lambda1max >= VP$lambda1[[1]]
#   )
#
#
# # It woooooooorks
#
# map(1:7, function(x){
#
#   self$algo2list[[paste0("first_",x)]]$poolVP
# }) %>%
#   bind_rows() -> allVPs
#
# allVPs %>%
#   arrange(id)
# unnest() %>%
#   ggplot()+
#   geom_line(aes(time, tumVol, col = factor(id)))+
#   facet_wrap(~protocol)+
#   scale_y_log10()
#
#
# self$algo2list$tree %>%
#   pull(time) %>%
  # sum() / 60
#


