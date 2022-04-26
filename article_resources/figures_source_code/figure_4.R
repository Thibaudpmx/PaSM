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
## For plot A: first a script to generate the data, stored into various RDS files
## Then a second script to read all those RDS files and produce the plot.
##
##
## ---------------------------

library(QSPVP)



#  plot A - data generation: Simeoni analysis -----------------------------------------------





# Goal: vary percentage of acceptance/rejection and number of varying parameter
# Method: generate 200.000 VPs, compute them manually, then use the quantiles function
# to adjust the target, determining the pct of acceptance needed


# First, put where you want to save all the data (to heavy to put in GitHub, you will have to regenerate them yourself ! )
setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Simeoni")

# We first need a function to automatically create the 200.000 VPs according the number of varying parameters, so here it is


cohort_creator <- function(nmodif){

  base <- crossing( k2 = 0.5,
                    lambda0 = 0.1,
                    Vd =  40,
                    lambda1 = 12,
                    ke = 1 ,
                    w0 = 50
  )


  nperparam <- ceiling(200000^(1/nmodif)) # count how many parameter values per paramer, superior round


  list <-   map(1:nmodif, function(x){ # for each parameter, create the sequence function

    min <- base[[x]]/5
    max <- base[[x]]*5
    step <- (max - min)/(nperparam-1)

    expr(seq( !!min,  !!max, !!step  ))
  }
  )

  names(list) <- names(base[1:nmodif])

  if(nmodif < length(base)){ # create the data.frame of the cohort

    output <- crossing(base[-c(1:nmodif)], !!!list)

  }else{


    output <- crossing(!!!list)
  }

  set.seed(255661) # let's add a seed to have a deterministic behavior

  output %>%
    sample_n(200000) %>% # and randomly sample 200.000 VPs (due to the ceiling function, nrow(output) > 200K)
    mutate(k1 = 0.5, psi = 20) #%>%
    # map(~length(unique(.x)))


}


# Then we perform the main big loops.
# For each number of varying parameter (1 to 5)
#     + Simulate all VPs by having unlimited target and green filter disabled, save the result in a RDS file
#     + For all percentage of acceptance/rejection wanted
#             + Take all reference simulations and compute the needed targets to reach this acceptance rejetion
#             + Use our algorithm 5 times and, for each time, save the result in a RDS file


pcttargets <- c(1,0.875,0.99, 0.75,0.5,0.25,0.125,0.01,0)

a <- 5; b <- 1

for(a in 1:5){ # For each number of varying parameter to try

  cat(paste0("a = ", a))

  cohort <- cohort_creator(nmodif = a) # create the cohort of 200.000 VPs

  self <- VP_proj_creator$new() # create a new VP project
  target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min = -1E99, max = 1E99)
  self$set_targets(manual =target) # Note: can't use (-) Inf

  namepct <- paste0("Ref_",a, ".RDS") # Name of all simulations file


  if(!file.exists(namepct)){ # If it does not exist yet, compute it

    self2 <- self$clone(deep = T) # Create a clone
    self2$add_VP(cohort, use_green_filter = F, npersalve = 2000) # Compute, not using the green filter
    self2$poolVP %>%
      unnest(simul) %>%
      filter(time %in% target$time) %>% # take only the tumVol at the time of interest
      saveRDS(file = namepct )

    rm(self2)
  }

  alltumVol <- readRDS(namepct)  %>% pull(tumVol) # Read the Ref file and take all 200.000 tumVol values


  for(b in pcttargets){ # for every percentage

    for(d in 1:5){ # because each analyse repeated five time

    newnames <- paste0(a,"_", b, "_", d, ".RDS") # compute the name of the file (nparamvarying_pcttarget_iteration)

    if(!file.exists(newnames)){ # If the file does not exist yet
      # determine new targets
      quant <- ( 1-b)/2 # determine the probability for quantile function
      prototemp <- target # and apply the new min and max
      prototemp$min <- quantile(alltumVol, probs = quant)
      prototemp$max <- quantile(alltumVol, probs = 1 - quant)


      self2 <- self$clone(deep = T) # create a copy of self
      self2$set_targets(manual = prototemp) # use the new target


      # and do the computation ! With time_compteur activated for doing further analyses
      self2$add_VP(time_compteur = T, cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0, saveVPRej = F)


      saveRDS(self2, file = newnames)
    }

    }

  }# end for pcttargets

  if(a == 5){

  for(b in pcttargets){ # for every percentage

    for(d in 1:5){ # because each analyse repeated five time

      newnames <- paste0(a,"_", b, "_", d, "_alt.RDS") # compute the name of the file (nparamvarying_pcttarget_iteration)

      if(!file.exists(newnames)){ # If the file does not exist yet
        # determine new targets◘
        prototemp <- target # and apply the new min and max
        prototemp$min <- 0
        prototemp$max <- quantile(alltumVol, probs = b)

        if(b == 0)  prototemp$max <-  prototemp$max *0.8

                  self2 <- self$clone(deep = T) # create a copy of self
        self2$set_targets(manual = prototemp) # use the new target


        # and do the computation ! With time_compteur activated for doing further analyses
        self2$add_VP(time_compteur = T, cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0, saveVPRej = F)


        saveRDS(self2, file = newnames)
      }

    }

  }# end for pcttargets
  } # end for a == 6

}



# Now analyses the results


timetable <- function(self){

  tt <- self$timeTrack

  total <-tt$tTOTAL


  loop <- tt$poolVP_compteur %>%
    mutate( Treduc_filter_neg_both  =  Treduc_filter_neg_below +  Treduc_filter_neg_above) %>%
    mutate( RedFilter = TapplyRedFilterBoth + Treduc_filter_neg_both, RxODE = timemodel) %>%
    mutate(patientRemovedperSec = as.double(NremovedRedFilter) / as.double((RedFilter + timemodel )),
           refpatientRemovedperSec =  as.double(nsimul) / as.double((timemodel )))


  if("Treduc_filter_green_below_tumVol" %in% names( tt$poolVP_compteur)){

    loop <- loop %>%
      mutate( Treduc_filter_pos_both  =  Treduc_filter_green_below_tumVol +  Treduc_filter_green_above_tumVol) %>%
      mutate( TapplyGreenFilter  =  TapplyGreenFilterBelow_tumVol +  TapplyGreenFilterAbove_tumVol,
              GreenFilter =  Treduc_filter_pos_both + TapplyGreenFilter)



  }else{

    loop <- loop %>%
      mutate( Treduc_filter_pos_both  =  0) %>%
      mutate( TapplyGreenFilter  =  0,
              GreenFilter =  0)

  }

  # loop %>%
  #   select(timesampleline, protocoldf, timemodel, Treduc_filter_neg_both ,TapplyRedFilterBoth, TremoveNegDirect,
  #          Treduc_filter_pos_both, TapplyGreenFilter, timesimlandjoin,res2) %>%
  #   # loop %>%
  #   # select(timesampleline, protocoldf, timemodel, res2, filter_neg_below, filter_neg_above, remneg_fil, time_addgreennofil, timesimlandjoin) %>%
  #   gather("step", "value") %>%
  #   group_by(step) %>%
  #   summarise(sum = sum(value, na.rm = T))

  loop %>%
    select(RedFilter, GreenFilter, RxODE) %>%
    gather("step", "value") %>%
    group_by(step) %>%
    summarise(sum = sum(value, na.rm = T)) %>%
    mutate(sum = as.double(sum)) -> TimeSimplified

  TimeSimplified %>%
    add_row(step = "Other", sum = as.numeric(total, units = "days") * (24 * 3600) - sum(TimeSimplified$sum))

}


files <- list.files()


toread <- files[!grepl("^Ref", files)]
toread <- toread[! toread %in%  c("timeReference.RDS", "full_analysis.RDS")]


str_split(toread, pattern = "_", simplify = T) %>%
  as_tibble() %>%
  mutate(V3 = gsub("\\.RDS","", V3 )) %>%
  rename(nparam = V1, pct = V2, iteration = V3, meth= V4) %>%
  mutate(pct = as.double(pct), nparam = as.double(nparam)) %>%
  mutate(file = toread) %>%
  mutate(results = map(file, function(x){

    # x <- "1_0.5_3.RDS"
# print(x)
    obj <- readRDS(x)

    timetable(obj) %>%
      spread(step, sum) %>%
      mutate(VPfound= obj$poolVP %>% nrow, nsimul = obj$timeTrack$poolVP_compteur %>% pull(nsimul) %>% sum,
             timeTotal = obj$timeTrack$tTOTAL %>% as.double(), Tgreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalGreenFilter)%>% as.double(),
             Tred1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalRedFilter)%>% as.double(),
             nremRed1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(NremovedRedFilter)%>% as.double(),
             ndoneGreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(nextrapoGreen)%>% as.double())


  })) %>%
  unnest() %>%
saveRDS("full_analysis.RDS")




manual <- function(target, cohort){

  cohort <- crossing(cohort,  proto  = unique(target$protocol))


  self$protocols[] %>%
    bind_rows() %>%
    mutate(proto = names( self$protocols)) -> protocols

  eventsadmin  <- crossing(id = 1:2000, proto  = unique(target$protocol)) %>%
    left_join(protocols, by = "proto" ) %>%
    mutate(evid = 1)

  eventsadmin  <-   eventsadmin %>%
    bind_rows(

      eventsadmin  %>%
        mutate(evid = 0, amt = 0) %>%
        select(-time) %>%
        crossing(time = self$times)
    ) %>%
    arrange(id, time)


  demo <- cohort %>%
    rowid_to_column("ids") %>%
    mutate(group = floor(ids/2000)+1)


  idtorems <- double()

  resultsap <- list()


  for(a in unique(demo$group)){


    x <- demo %>% filter(group == a) %>%
      rowid_to_column("id")

    events <- eventsadmin %>%
      left_join(x%>% select(id, ids, group, proto) , by = c("id", "proto")) %>%
      filter(!is.na(ids))



    res <- self$model$solve(x %>% select(-proto), events , c(X2 = 0)) %>%
      as_tibble

    if(max(x$id) == 1) res <- res %>% mutate(id = 1)

    res <- res %>%
      left_join(x %>% select(id, proto, ids), by = "id")

    res %>%
      filter(time %in% target$time) %>%
      rename(protocol = proto) %>%
      gather("cmt", "value", unique(target$cmt)) %>%
      left_join(target, by = c("time", "protocol", "cmt")) %>%
      filter(value > max | value < min) %>% pull(ids ) -> idtorem

    idtorems <- c(idtorems, idtorem)


    resultsap[[a]]<-  res %>%   filter( ! (ids %in% idtorem))



  }

  resultsap <- bind_rows(resultsap)

  demo <- demo %>%
    filter(! (ids %in% idtorems)) %>%
    left_join( resultsap, by = c("ids", "proto") )


  # difftime(Sys.time(),t00, units = "s")

}



target <- tibble(protocol = "dose50", cmt = "tumVol", time = 48, min = -Inf, max = Inf)

cohort <- cohort_creator(nmodif = 2)


# Compute the reference

mbref <- microbenchmark(ref <- manual(target, cohort), times = 5)

saveRDS(mbref, "timeReference.RDS")

#  plot A - making: Simeoni analysis -----------------------------------------------
setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Simeoni")

allTimes <- readRDS("full_analysis.RDS") #%>%
  # mutate(meth = if_else(meth == "", "", "alt")) %>%
  # mutate(nparam = paste0(nparam, meth))

readRDS("timeReference.RDS") %>%
  as.data.frame() %>%
  mutate(time = time * 10E-10) %>%
  summarise(min = min(time), max = max(time), median = median(time)) -> ref

# Verifying the percentage

allTimes %>%
  mutate(pcteffe = VPfound / 200000) %>%
  filter(pct != pcteffe) # note: 1 or 2 of differences can be explained by internal rounding
                         # After check the additional/missing VP are from 1E-10 mm3 or similar


# Main plot time benefice

plotA <- allTimes %>%

  mutate(pct = as.double(pct)) %>%
  mutate(nparam = as.character(nparam)) %>%
  mutate(meth = if_else(meth == "", "Centered", "Lowest")) %>%
  group_by(nparam, pct, meth) %>%
  summarise(timeTotal = median(timeTotal)) %>%
  ungroup() %>%
  ggplot()+
  geom_rect(aes(xmin = 48, xmax = 52, ymin = 39, ymax = 41.5), col = "red",alpha = 0)+
  geom_segment(aes(x = 56, xend = 52, y = 45, yend = 42), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+
  # geom_text(aes(x = 62, y = 46, label = "next plots"), col = "red", size = 3)+ # Note: last 3 brut way, not reproducible
  geom_point(aes(x = (1-pct) * 100, y = timeTotal, col = nparam))+
  geom_line(aes(x = (1-pct)* 100, y = timeTotal, col = nparam, lty = meth))+
  geom_hline(data = ref, aes(yintercept = min), lty = 1)+
  geom_hline(data = ref, aes(yintercept = max), lty = 1)+
  geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max, fill = "Time of\nreference"), lty = 1,alpha = 0.2)+
  geom_hline(data = ref, aes(yintercept = median), lty = 2)+
  theme_bw()+
  scale_linetype_manual(values = c(1,2))+
  scale_fill_manual(values = "grey")+
  labs(x = "Percentage of rejected VP", y = "Total Time analysis (sec)", col = "Number\nvarying\nparameters", fill = "",
       lty = "Target\nmethod"); plotA

# wallTimes %>%
#   ggplot()+
#   geom_point(aes(x = (1-pct), y = median, col = factor(id)))+
#   geom_line(aes(x = (1-pct), y = median, col = factor(id)))+
#   geom_hline(data = ref, aes(yintercept = min), lty = 1)+
#   geom_hline(data = ref, aes(yintercept = max), lty = 1)+
#   geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.2)+
#   geom_hline(data = ref, aes(yintercept = median), lty = 2)+
#   geom_hline(yintercept =  16.7)+
#   # geom_errorbar(data = ref, aes(x = 0, ymin = min, ymax = max))+
#   theme_bw()
#
#
# allTimes %>%
#   mutate(pct = as.double(pct)) %>%
#   ggplot()+
#   geom_point(aes(x = (1-pct), y = nsimul, col = nparam))+
#   geom_line(aes(x = (1-pct), y = nsimul, col = nparam))+
#   geom_hline(data = ref, aes(yintercept = min), lty = 1)+
#   geom_hline(data = ref, aes(yintercept = max), lty = 1)+
#   geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.2)+
#   geom_hline(data = ref, aes(yintercept = median), lty = 2)+
#   # geom_hline(yintercept =  16.7)+
#   # geom_errorbar(data = ref, aes(x = 0, ymin = min, ymax = max))+
#   theme_bw()
#
#
# allTimes %>%
#   mutate(pct = as.double(pct)) %>%
#   ggplot()+
#   geom_boxplot(aes(x = factor((1-pct)), y = as.double(timeTotal)))



# Plot B ------------------------------------------------------------------
# First, take only the median time total

allTimes %>%
  mutate(timeTotal= as.double(timeTotal)) %>%
  left_join(
allTimes %>%
  group_by(nparam, pct,meth) %>%
  summarise(  timeTotal = median(as.double(timeTotal)) ) %>%
  mutate(test = T)
  ) %>%
  filter(test) %>%
  select(nparam, pct, meth,iteration,test) -> iteration_to_keep







plotB <- allTimes %>%
  filter(nparam == 5 & pct == 0.5 & meth == "") %>%
  left_join(iteration_to_keep) %>%
  filter(test) %>%
  gather("time", "value",GreenFilter, Other, RedFilter, RxODE) %>%
  mutate(pct = "New") %>%
  bind_rows( tibble(pct = "Old", time = c("Other", "RxODE"), value = c(7.1,47.9)) %>% crossing(nparam = 5)) %>%
  group_by(pct) %>%
  mutate(total = sum(value)) %>%
  ungroup() %>%
  ggplot()+
  geom_segment(aes(x = "New", xend = "New", y = 53, yend = 50), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+
  geom_col(aes(pct, value, fill = fct_reorder(time, value, .desc = F)), alpha = 0.5)+
  # geom_col(data =,aes(pct, value, fill = time), alpha = 0.4)+
  theme_bw()+
  # facet_wrap(~nparam)+
  scale_fill_manual(values = c("grey", "red", "darkgreen", "blue"))+
  # geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.2)+
  # geom_hline(data = ref, aes(yintercept = median), lty = 2)+
  labs(x = "Method", y = "Time of analysis (sec)", fill = "Step")+
  geom_text(aes(pct, value, fill = fct_reorder(time, value, .desc = F), label = as.double(value) %>% round(1)), position = position_stack(vjust = .5))+
  geom_label(aes(pct, total+ 5 ,  label = paste0("Total:",as.double(total) %>% round(1), "s"))); plotB

cowplot::plot_grid(plotA, plotB, nrow =   1)



# plot S A ----------------------------------------------------------------


suplementA <- allTimes %>%
  # filter(nparam !=1) %>%
  left_join(iteration_to_keep) %>%
  filter(test) %>%
  gather("time", "value",GreenFilter, Other, RedFilter, RxODE) %>%
  mutate(pct = as.factor((1-pct) * 100 )) %>%
  mutate(label = if_else(value < 4, "", as.character(value))) %>%
  mutate(nparam = paste0(nparam, meth) ) %>%
  mutate(nparam = if_else(nparam == "5alt.RDS", "5 alt", nparam)) %>%
  # bind_rows( tibble(pct = "Ref", time = c("Other", "RxODE"), value = c(7.1,47.9)) %>% crossing(nparam = 2:5)) %>%
  ggplot()+
  geom_col(aes(pct, value, fill = fct_reorder(time, value, .desc = F)), alpha = 0.5)+
  # geom_col(data =,aes(pct, value, fill = time), alpha = 0.4)+
  theme_bw()+
  facet_wrap(~nparam)+
  scale_fill_manual(values = c("grey", "blue", "red", "darkgreen"))+
  geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.2)+
  geom_hline(data = ref, aes(yintercept = median), lty = 2)+
  labs(x = "Percentage of rejected VP", y = "Time of analysis (sec)", fill = "")+
  geom_text(aes(pct, value, fill = fct_reorder(time, value, .desc = F), label = as.double(label) %>% round(1)), position = position_stack(vjust = .5)); suplementA


# plot C - data generation ----------------------------------------------------------------

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/time_impact")


times_to_try <- list(seq(0,48,48), seq(0,48,4),seq(0,48,2),seq(0,48,1),seq(0,48,0.5),seq(0,48,0.1))


cohort <- cohort_creator(nmodif = 5)

alltumVol <- readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Simeoni/Ref_5.RDS") %>%
  pull(tumVol)


# determine new targets
prototemp <- tibble(protocol = "dose50", cmt = "tumVol", time = 48, min =  quantile(alltumVol, probs = 0.25),
                    max = quantile(alltumVol, probs = 0.75))# and apply the new min and max


for(a in 1:length(times_to_try)){


  for(b in 1:5){


namefile <- paste0(a, "_", b,".RDS" )

if(!file.exists(namefile)){

  self <- VP_proj_creator$new()
  self$set_targets(manual = prototemp)
  self$times <- times_to_try[[a]]

  self$add_VP(VP_df = cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0, saveVPRej = F, time_compteur = T)
  self$timeTrack$ttotal

  saveRDS(object = self, namefile)
  }



  } # end for b

} # end for a

files <- list.files()


toread <- files[!grepl("^Ref", files)]



str_split(toread, pattern = "_", simplify = T) %>%
  as_tibble() %>%
  mutate(V2 = gsub("\\.RDS","", V2 )) %>%
  rename(Step = V1 , iteration = V2) %>%
  # mutate(pct = as.double(pct), nparam = as.double(nparam)) %>%
  mutate(file = toread) %>%
  mutate(results = map(file, function(x){

    # x <- "1_1.RDS"
    # print(x)
    obj <- readRDS(x)

    timetable(obj) %>%
      spread(step, sum) %>%
      mutate(VPfound= obj$poolVP %>% nrow, nsimul = obj$timeTrack$poolVP_compteur %>% pull(nsimul) %>% sum,
             timeTotal = obj$timeTrack$tTOTAL)


  })) %>%
  unnest() %>%
  saveRDS("full_analysis.RDS")

self <- VP_proj_creator$new() # create a new VP project
self$set_targets(manual = tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min = -1E99, max = 1E99))
target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min = -1E99, max = 1E99)

time_Impact <- function(times = 1:40){

  cohort <- cohort_creator(nmodif = 2)



  self$protocols[] %>%
    bind_rows() %>%
    mutate(proto = names( self$protocols)) -> protocols

  eventsadmin  <- crossing(id = 1:2000, proto  = unique(target$protocol)) %>%
    left_join(protocols, by = "proto" ) %>%
    mutate(evid = 1)

  eventsadmin  <-   eventsadmin %>%
    bind_rows(

      eventsadmin  %>%
        mutate(evid = 0, amt = 0) %>%
        select(-time) %>%
        crossing(time = times)
    ) %>%
    arrange(id, time)


  demo <- cohort %>%
    rowid_to_column("ids") %>%
    mutate(group = floor(ids/2000)+1)


  idtorems <- double()

  resultsap <- list()


  # for(a in unique(demo$group)){
  a <- 1

  x <- demo %>% filter(group == a) %>%
    rowid_to_column("id")

  events <- eventsadmin%>% filter(id <2000)

  t0 <- Sys.time()
  res <- self$model$solve(x, events  , c(X2 = 0)) %>%
    as_tibble

  difftime(Sys.time(), t0)

}



allTimesRxODE <-  map(times_to_try, function(x){

  time_Impact(x) * 100

})
# plot C - plot  -------

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/time_impact")



times_to_try <- list(seq(0,48,48), seq(0,48,4),seq(0,48,2),seq(0,48,1),seq(0,48,0.5),seq(0,48,0.1))

names(allTimesRxODE) <- map(times_to_try, ~ length(.x) %>% as.character)

temp <- readRDS("full_analysis.RDS") %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  group_by(Step) %>%
  mutate(median = median (timeTotal)) %>%
  filter(median == timeTotal) %>%
  ungroup() %>%
  mutate( n = map_dbl(times_to_try, ~ length(.x)) ) %>%
  mutate(timecomput =  map_dbl(allTimesRxODE, ~ as.double(.x)/ 100) ) %>%
  mutate(Old = 10 + allTimesRxODE %>% reduce(c) %>% as.double()) %>%
  rename(New = timeTotal) %>%
  gather("method", "value", Old, New) %>%
  mutate(value= as.double(value))

plotC <- temp %>%
  ggplot()+
  geom_vline(aes(xintercept = 0.375, lty= "Used in\nmain\nanalysis"))+
  geom_line(aes(timecomput, value, col = method), size = 2)+
  geom_point(aes(timecomput, value, col = method), size = 3)+
  geom_segment(aes(x = 0.45, xend = 0.39, y = 35, yend = 39), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+

  # geom_line(aes(timecomput, value, col = method), size = 2)+
  geom_ribbon(data = temp %>% spread(method, value), aes(x = timecomput, ymin = New, ymax =  Old), alpha = 0.3)+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Time to compute 2000 VPs with RxODE (sec)",y = "Time for performing 200.000 VPs (sec)", lty = "", col = "Method")+
  scale_linetype_manual(values = 2);plotC


plot_grid(plotA, plotB, plotC, nrow = 1)
temp <- tibble(n = map_dbl(times_to_try, ~ length(.x)), brut = allTimesRxODE %>% reduce(c),
               new  = allTimes %>% reduce(c)  ) %>%
  mutate(n2 = n , n = as.double(brut) / 2000 ) %>%
  mutate(brut = brut  * 0.7+ 7.1) %>%
  mutate(brut= as.double(brut), new = as.double(new))



# plot D ------------------------------------------------------------------
temp <- allTimes %>%
  # mutate(meth = if_else(meth == "", "Centered", "Lowest")) %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  group_by(nparam, pct, meth) %>%
  mutate(median = median(timeTotal)) %>%
  ungroup() %>%
  filter(timeTotal == median) %>%
  mutate(nparam = paste0(nparam, if_else(meth == "", "", "alt")))


plotD <- temp %>%
  ggplot()+
  geom_point(aes(nparam,nsimul, col = factor(1-pct)))+

  geom_line(data = temp %>% filter(nparam != "5alt"), aes(nparam, nsimul, col =factor(1-pct), group = pct))+
  geom_line(data = temp %>% filter(nparam %in% c("5", "5alt")), aes(nparam, nsimul, col =factor(1-pct), group = pct), lty = 2)+
   theme_bw()+
  geom_segment(aes(x = 5.5, xend = 5.1, y = 49000, yend = 47532), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+

  scale_y_log10()+
  labs(x = "Number of varying parameter", y = "Number of simulations performed",col = "Percentage\nRejection");plotD


# plot E ------------------------------------------------------------------
temp <- allTimes %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  group_by(nparam, pct, meth) %>%
  mutate(median = median(timeTotal)) %>%
  mutate(filtreTime = (Tgreen1 + Tgreen1)) %>%
  ungroup() %>%
  filter(timeTotal == median) %>%
  mutate(nparam = paste0(nparam, if_else(meth == "", "", "alt")))

plotE<- temp %>%
  ggplot()+

  geom_point(aes(nparam,filtreTime, col = factor(1-pct)))+
  geom_segment(aes(5.4, xend = 5.1, y = 16.5, yend = 14.7), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+

  geom_line(data = temp %>% filter(nparam != "5alt"), aes(nparam, filtreTime, col =factor(1-pct), group = pct))+
  geom_line(data = temp %>% filter(nparam %in% c("5", "5alt")), aes(nparam, filtreTime, col =factor(1-pct), group = pct), lty = 2)+
  theme_bw()+
  scale_y_log10()+
  labs(x = "Number of varying parameter", y = "Filtering time first iteration (sec)",col = "Percentage\nRejection");plotE


# plot F ------------------------------------------------------------------


plotF <- allTimes %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  group_by(nparam, pct) %>%
  mutate(median = median(timeTotal)) %>%
  mutate(filtreTime = (nremRed1)) %>%
  ungroup() %>%
  filter(timeTotal == median) %>%
  ggplot()+
  geom_point(aes(nparam,filtreTime, col = factor(1-pct)))+
  geom_line(aes(nparam, filtreTime, col =factor(1-pct)))+
  theme_bw()+
  # scale_y_log10()+
  labs(x = "Number of varying parameter", y = "Number VP removed first iteration",col = "Percentage\nRejection");plotF




# Plot G data generation------------------------------------------------------------------

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/ImpactSizeCohort")

# Step 1: create a 1 million database

base <- crossing( k2 = 0.5,
                  lambda0 = 0.1,
                  Vd =  40,
                  lambda1 = 12,
                  ke = 1
)


nperparam <- ceiling(1E6^(1/5)) # count how many parameter values per paramer, superior round


list <-   map(1:5, function(x){ # for each parameter, create the sequence function

  min <- base[[x]]/5
  max <- base[[x]]*5
  step <- (max - min)/(nperparam-1)

  expr(seq( !!min,  !!max, !!step  ))
})

names(list) <- names(base[1:5])

cohort <- crossing(!!!list) %>%
  slice(1:1E6) %>%
  mutate(psi = 20, w0 = 50, k1 = 0.5)

self <- VP_proj_creator$new() # create a new VP project
self2 <- VP_proj_creator$new() # create a new VP project
target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min = -1E99, max = 1E99)
self$set_targets(manual =target ) # Note: can't use (-) Inf
self2$set_targets(manual =target )
self$add_VP(cohort %>% slice(1:500000), use_green_filter = F)
self2$add_VP(cohort %>% slice(500001:1E6), use_green_filter = F,keep = "tumVol" )

saveRDS(self, "self_1sthalf")
saveRDS(self2, "self_2sthalf")

self <- readRDS("self_1sthalf")


tofill <- tibble()
for(a in 0:4){
print(a)
  tofill <-  bind_rows(tofill,


self$poolVP %>%
slice((a * 1E5+ 1):((a + 1)*1E5)) %>%
  unnest() %>%
  filter(time == 48)
  )


}

rm(self)

saveRDS(tofill, "alltumVol")

self2 <- readRDS("self_2sthalf")

for(a in 0:4){
  print(a)
  tofill <-  bind_rows(tofill,


                       self2$poolVP %>%
                         slice((a * 1E5+ 1):((a + 1)*1E5)) %>%
                         unnest() %>%
                         filter(time == 48)
  )


}


tofill %>%
  select(id, k2, lambda0, Vd, lambda1, ke, tumVol, time) %>%
  arrange(tumVol) %>%
  select(-id) %>%
  rowid_to_column("id") %>%
saveRDS( "alltumVol")
rm(list = ls())

# Step 2: sample different tumvol

  allTumVol <- readRDS("alltumVol")


  allSizeCohort <- c(10, 50 , 100,250,500,750,1000) * 1000

  for(a in allSizeCohort){
print(a)

    for(b in 1){
    name <- paste0(a,"_", b,".RDS")

    if(!file.exists(name)){

   cohort <- allTumVol %>%
      slice(1:a) %>%
      select(-id, -time) %>%
      mutate(k1 = 0.5, psi = 20, w0 = 50)


    target <- tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min =  quantile(cohort$tumVol  , probs = 0.25),
                        max =  quantile(   cohort$tumVol, probs = 0.75)) # and apply the new min and max


    self <- VP_proj_creator$new() # create a new VP project
    self$set_targets(manual =target) # No
    self$add_VP(cohort, use_green_filter = T, time_compteur = T, npersalve = 2000, pctActivGreen = 0, saveVPRej = T)


    saveRDS(self,name )
    }
    }# end for b
  }

  # Step 3: compute the stats

  toread <- list.files()


  toread <- toread[! toread %in%  c("self_2sthalf", "self_1sthalf","alltumVol", "full_analysis.RDS")]


  str_split(toread, pattern = "_", simplify = T) %>%
    as_tibble() %>%
    mutate(V2 = gsub("\\.RDS","", V2)) %>%
    rename(nCohort = V1, iteration = V2) %>%
    mutate(file = toread) %>%
    mutate(nCohort = as.double(nCohort)) %>%
    mutate(results = map(file, function(x){

      # x <- "1_0.5_3.RDS"
      # print(x)
      obj <- readRDS(x)

      timetable(obj) %>%
        spread(step, sum) %>%
        mutate(VPfound= obj$poolVP %>% nrow, nsimul = obj$timeTrack$poolVP_compteur %>% pull(nsimul) %>% sum,
               timeTotal = obj$timeTrack$tTOTAL %>% as.double(), Tgreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalGreenFilter)%>% as.double(),
               Tred1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalRedFilter)%>% as.double(),
               nremRed1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(NremovedRedFilter)%>% as.double(),
               ndoneGreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(nextrapoGreen)%>% as.double())


    })) %>%
    unnest() %>%
    saveRDS("full_analysis.RDS")



 # Plot G -plot------------------------------------------------------------------

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/ImpactSizeCohort")


  datas <- readRDS("full_analysis.RDS") %>%
    arrange(nCohort)


 plotG <-  datas %>%
    mutate(Old = 55.1 * nCohort / 2E5) %>%
    mutate(timeTotal = if_else(timeTotal < 4, timeTotal * 60, timeTotal)) %>%  #transform into minutes whats in sec (carefull reproductibiliy...)
    rename(New = timeTotal) %>%
    {pregather <<- .} %>%
    gather("Method", "value", Old, New) %>%
   mutate(pctExtra = round((nCohort-nsimul)/ nCohort * 100 )) %>%
   mutate(label = if_else(Method == "Old", "", paste0(pctExtra, "%"))) %>%
    ggplot()+
    geom_vline(aes(xintercept = 2E5, lty = "Used in\nmain\nanalysis"))+
    geom_segment(aes(x = 2.6E5, xend = 2.2E5, y = 35, yend = 42), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+
   # geom_text(aes(nCohort, value, label = label), nudge_y = -20)+
   geom_line(aes(nCohort, value, col = Method), size = 2) +
    geom_ribbon(data = pregather, aes(x = nCohort, ymin = Old, ymax = New), alpha = 0.3)+
    geom_point(aes(nCohort, value, col = Method), size = 3)+
     theme_bw() +
   scale_linetype_manual(values = 2)+
  labs(x = "Number of VPs in original cohort", y = "Time of computation (sec)");plotG



# Final merge for fig 4 ---------------------------------------------------

library(cowplot)

plot_grid(plotA, plotB, plotD, plotE, plotC, plotG, labels = LETTERS, nrow = 2)
