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





# Step 1 create the hide and seek VP

self <- VP_proj_creator$new()


VP <- tibble(k2 = 1.63, lambda0 = 0.85, ke = 0.6, Vd = 36, lambda1 = 72, w0 = 50, k1 = 0.5, id = 1:3, psi = 20)

events <- tibble(cmt = "Central", time = 0, amt = c(0,50,100), evid = 1, id = 1:3) %>%
          bind_rows(tibble(cmt = "tumVol", time = 0:52 , amt = 0, evid =  0) %>% crossing(id = 1:3)) %>%
  arrange(id, time)


simul <- self$model$solve(VP, events) %>%
  as_tibble

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


# data generation lets find it --------------------------------------------


self <- VP_proj_creator$new()

self$set_targets(manual = prototiny)

npersalve = 2E5
npersalveFinal = 1E6


domain <- tribble(~param, ~from, ~to, ~digits,
                  "k2", 0, 3, 2 ,
                  "lambda0", 0, 1.4, 2,
                  "ke", 0, 2,1,
                  "Vd", 0,40,0,
                  "lambda1", 0,240,0
)


domain <- tribble(~param, ~from, ~to, ~by,
                  "k2", 0, 3, 0.01 ,
                  "lambda0", 0, 1.4, 0.01,
                  "ke", 0, 2,0.1,
                  "Vd", 0,40,1,
                  "lambda1", 0,240,1
)

fix <-c(k1 = 0.5, w0 = 50)
self$algo2(domain = domain, fix = fix,npersalve =  npersalve, npersalveFinal = npersalveFinal)

fix <-c(k1 = 0.5, w0 = 50)

file <- "D:/these/Second_project/QSP/modeling_work/VT_simeoni/fig5_data.RDS"


# Just to confirm after first analyse there is still a bloc containing our VPs
maybeFinal %>%
  filter(k2min <= VP$k2[[1]], k2max >= VP$k2[[1]],
         kemin <= VP$ke[[1]], kemax >= VP$ke[[1]],
         Vdmin <= VP$Vd[[1]], Vdmax >= VP$Vd[[1]],
         lambda0min <= VP$lambda0[[1]], lambda0max >= VP$lambda0[[1]],
         lambda1min <= VP$lambda1[[1]], lambda1max >= VP$lambda1[[1]]
  )


# It woooooooorks

map(1:7, function(x){

  self$algo2list[[paste0("first_",x)]]$poolVP
}) %>%
  bind_rows() -> allVPs

allVPs %>%
  arrange(id)
  unnest() %>%
  ggplot()+
  geom_line(aes(time, tumVol, col = factor(id)))+
  facet_wrap(~protocol)+
  scale_y_log10()


  self$algo2list$tree %>%
    pull(time) %>%
    sum() / 60



# Plot making -------------------------------------------------------------

  self <-  readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/fig5_data.RDS")


tree <- self$algo2list$tree

# Just a template, to replace by tracking time inside the function
template <- tibble(step = c("1_algo1", "2_filter", "3_maybe", "4_reduce", "5_other"), time = c(21.7,15.19,42.6,75.5, 7) ,method = "Full") %>%
  bind_rows(

    tibble(step = c("1_algo1", "2_filter", "3_maybe", "4_reduce", "5_other"), time = c(0,0,0,168, 7),method = "alt")


  )



plotA <- template %>%
  ggplot()+
  geom_bar(aes(x = step,y= time, fill = method), stat="identity",col = "black")+
  labs(y = "Time(sec)")+
  theme_bw(); plotA

template %>%
  ggplot()+
  geom_bar(aes(x = method,y= time, fill = step), stat="identity",col = "black")+
  labs(y = "Time(sec)")+
  theme_bw()

template %>%
  ggplot()+
  geom_bar(aes(x = fct_reorder(step, time, .desc = T),y= time), stat="identity")+
  theme_bw()

plotB <- tree   %>%
  slice(-1) %>%
  mutate(time = as.double(time)) %>%
  ggplot()+
  geom_histogram(aes(time, fill = factor(after/3)),col = "black")+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Time (sec)", fill = "VP\nfound"); plotB




domain <- tribble(~param, ~from, ~to, ~digits,
                  "k2", 0, 3, 2 ,
                  "lambda0", 0, 1.4, 2,
                  "ke", 0, 2,1,
                  "Vd", 0,40,0,
                  "lambda1", 0,240,0
)

fix <-c(k1 = 0.5, w0 = 50)


self$algo2list$tree %>%
  summarise(time = sum(time))

2609.409 / 60

self$algo2list$first %>%
  group_by(blocsPool) %>%
  summarise(sum = sum(temp3)) %>%
  pull(sum) %>% median

cowplot::plot_grid(plot0,"3",  plotA, plotB, labels = LETTERS)


