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
## GitHub: https://github.com/Thibaudpmx/PaSM
## ---------------------------
##
## For plot A: first a script to generate the data, stored into various RDS files
## Then a second script to read all those RDS files and produce the plot.
##
##
## ---------------------------

library(PaSM)



#  plot A - data generation: Simeoni analysis -----------------------------------------------





# Goal: vary percentage of acceptance/rejection and number of varying parameter
# Method: generate 200.000 VPs, compute them manually, then use the quantiles function
# to adjust the target, determining the pct of acceptance needed


# First, put where you want to save all the data (to heavy to put in GitHub, you will have to regenerate them yourself ! )
setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_PaSM/data/Simeoni")

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

  # Compute alternative methods
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

  namesCol <- names(loop)

  if(grepl("Treduc_filter_green",  namesCol) %>% sum > 1){

    naemesReduc <- namesCol[grepl("Treduc_filter_green", namesCol)]
    naemesReduc2 <- namesCol[grepl("TapplyGreenFilter", namesCol)]



    loop[ , naemesReduc] %>% map_df(~ as.double(.x)) %>% sum(na.rm = T)

    loop <- loop %>%
      mutate( Treduc_filter_pos_both  =  !!parse_expr(paste0(naemesReduc, collapse = "+"))) %>%
      mutate( TapplyGreenFilter  = !!parse_expr(paste0(naemesReduc2, collapse = "+")),
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
setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_PaSM/data/Simeoni_no_green")

allTimes <- readRDS("full_analysis.RDS") #%>%
# mutate(meth = if_else(meth == "", "", "alt")) %>%
# mutate(nparam = paste0(nparam, meth))


ref <- readRDS("../Simeoni_ref/timeReference.RDS")

# Verifying the percentage

# After check the additional/missing VP are from 1E-10 mm3 or similar



allTimesWith_green <- readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_PaSM/data/Simeoni/full_analysis.RDS") %>%
  mutate(pct = as.double(pct)) %>%
  mutate(nparam = as.character(nparam)) %>%
  mutate(meth = if_else(meth == "", "Centered", "Lowest")) %>%
  group_by(nparam, pct, meth) %>%
  summarise(timeTotal = median(timeTotal)) %>%
  ungroup() %>%
  rename(Yes = timeTotal)
# Main plot time benefice

plotA <- allTimes %>%

  mutate(pct = as.double(pct)) %>%
  mutate(nparam = as.character(nparam)) %>%
  mutate(meth = if_else(meth == "", "Centered", "Lowest")) %>%
  group_by(nparam, pct, meth) %>%
  summarise(timeTotal = median(timeTotal)) %>%
  ungroup() %>%
  ggplot()+
  geom_line(data = ref, aes(x= (1- pct)*100, y = medianTotal ), lty = 2)+
  geom_ribbon(data = ref, aes(x= (1- pct)*100, ymin = minTotal, ymax = maxTotal, fill = "Time of\nreference" ),alpha = 0.2)+
  geom_line(data = ref, aes(x= (1- pct)*100, y = minTotal ))+
  geom_line(data = ref, aes(x= (1- pct)*100, y = maxTotal ))+
  geom_point(aes(x = (1-pct) * 100, y = timeTotal, col = nparam))+
  geom_line(aes(x = (1-pct)* 100, y = timeTotal, col = nparam, lty = meth))+
  geom_line(data = allTimesWith_green , aes(x = (1-pct)* 100, y =  Yes, col = nparam, lty = meth), alpha = 0.5)+

  theme_bw()+
  scale_linetype_manual(values = c(1,2))+
  scale_fill_manual(values = "grey")+
  labs(x = "Percentage of rejected VP", y = "Total Time analysis (sec)", col = "Number\nvarying\nparameters", fill = "",
       lty = "Target\nmethod"); plotA




# plotB <- allTimes %>%
#   mutate(pct = as.double(pct)) %>%
#   mutate(nparam = as.character(nparam)) %>%
#   mutate(meth = if_else(meth == "", "Centered", "Lowest")) %>%
#   group_by(nparam, pct, meth) %>%
#   summarise(timeTotal = median(timeTotal)) %>%
#   ungroup() %>%
#   rename(No = timeTotal) %>%
#   left_join(allTimesWith_green  ) %>%
#   mutate(dif = No - Yes) %>%
#   ggplot()+
#   # geom_rect(aes(xmin = 48, xmax = 52, ymin = 39, ymax = 41.5), col = "red",alpha = 0)+
#   # geom_segment(aes(x = 56, xend = 52, y = 45, yend = 42), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+
#   # geom_text(aes(x = 62, y = 46, label = "next plots"), col = "red", size = 3)+ # Note: last 3 brut way, not reproducible
#
#   geom_point(aes(x = (1-pct) * 100, y = dif, col = nparam))+
#   geom_line(aes(x = (1-pct)* 100, y = dif, col = nparam, lty = meth))+
#   theme_bw()+
#   scale_linetype_manual(values = c(1,2))+
#   scale_fill_manual(values = "grey")+
#   labs(x = "Percentage of rejected VP", y = "Total Time analysis (sec)", col = "Number\nvarying\nparameters", fill = "",
#        lty = "Target\nmethod"); plotB



# Plot B DATA generation -------------------------------------------------------

# Changing the number of  dimension
setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_PaSM/data/Simeoni_less_usable_param")

alltumVol <- readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_PaSM/data/Simeoni_no_green/Ref_5.RDS") %>%
   pull(tumVol)


target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48,
                  min = quantile(alltumVol, probs = 0.25), max =  quantile(alltumVol, probs = 0.75))

cohort <- cohort_creator(nmodif = 5) # create the cohort of 200.000 VPs


pcttargets <- 0.5

a <- 5; b <- 1

for(a in 0:5){ # For each number of varying parameter to try

  cat(paste0("a = ", a))


  self <- VP_proj_creator$new() # create a new VP project
  self$set_targets(manual =target) # Note: can't use (-) Inf

  namepct <- paste0("Ref_",a, ".RDS") # Name of all simulations file


for(b in 0:1){ #activatio greenfilter

    for(d in 1:1){ # because each analyse repeated five time

      newnames <- paste0(a,"_", b, "_", d, ".RDS") # compute the name of the file (nparamvarying_pcttarget_iteration)

      if(!file.exists(newnames)){ # If the file does not exist yet
        # determine new targets



        self2 <- self$clone(deep = T) # create a copy of self


       if(a > 0) self2$param_increase$tumVol <- self2$param_increase$tumVol[-c(1:a)]
       if(a == 5) self2$param_reduce$tumVol <- character()

        # and do the computation ! With time_compteur activated for doing further analyses
        self2$add_VP(time_compteur = T, cohort, use_green_filter = b, npersalve = 2000, pctActivGreen = 0, saveVPRej = F)


        saveRDS(self2, file = newnames)
      }

    }

  }
}



files <- list.files()


toread <- files[!grepl("^Ref", files)]
toread <- toread[! toread %in%  c("timeReference.RDS", "full_analysis.RDS")]


str_split(toread, pattern = "_", simplify = T) %>%
  as_tibble() %>%
  mutate(V3 = gsub("\\.RDS","", V3 )) %>%
  rename(nparam = V1, pct = V2, iteration = V3) %>%
  mutate(pct = as.double(pct), nparam = as.double(nparam)) %>%
  mutate(file = toread) %>%
  mutate(results = map(file, function(x){

    # x <- "0_0_1.RDS"
    print(x)
    obj <- readRDS(x)

    timetable(obj) %>%
      spread(step, sum) %>%
      mutate(VPfound= obj$poolVP %>% nrow, nsimul = obj$timeTrack$poolVP_compteur %>% pull(nsimul) %>% sum,
             timeTotal = obj$timeTrack$tTOTAL %>% as.double(),
             Tred1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalRedFilter)%>% as.double(),
             nremRed1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(NremovedRedFilter)%>% as.double())


  })) %>%
  unnest() %>%
  saveRDS("full_analysis.RDS")



# Plot B  and C generation -------------------------------------------------------

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_PaSM/data/Simeoni_less_usable_param")


allTimes <- readRDS("full_analysis.RDS") %>%
  mutate(labelx = paste0( 5- nparam, "/5")) %>%
  mutate(GreenFilterlab = if_else(pct == 0, "Without", "With"))



plotB <- allTimes %>%
  ggplot()+
  geom_line(aes(x = 5-nparam, y = timeTotal, col =   GreenFilterlab))+
  geom_point(aes(x = 5-nparam, y = timeTotal, col =   GreenFilterlab))+
   geom_rect(data = mtcars, aes(xmin = -Inf, xmax = Inf, ymin = ref$minTotal[ref$pct == 0.5], ymax =  ref$maxTotal[ref$pct == 0.5], fill = "Time of\nreference"), lty = 1,alpha = 0.15)+
  geom_hline(data = mtcars, aes(yintercept = ref$medianTotal[ref$pct == 0.5]), lty = 2)+
  geom_hline(data = mtcars, aes(yintercept = ref$minTotal[ref$pct == 0.5]), lty = 1)+
  geom_hline(data = mtcars, aes(yintercept = ref$maxTotal[ref$pct == 0.5]), lty = 1)+

  scale_fill_manual(values = "grey")+
  # scale_y_log10()
  theme_bw()+
  scale_x_continuous(labels = c(paste0(0:5, "/5")))+
  labs( x = "Number of parameters with monotonicity ",y = "Time of analysis (s)", fill = "", col = "Green Filter");plotB


plotC <- allTimes %>%
  filter(pct == 1) %>%
  gather("param", "value", GreenFilter, RedFilter, RxODE,Other) %>%
  ggplot(aes(x = labelx, y = value, fill = param))+
  geom_col(aes(x = labelx, y = value, fill = param), alpha = 0.5)+
  geom_text(aes(x = labelx, y = value, label = paste0(round(value),"s"), group = param),       position = position_stack(vjust = .5))+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "grey", "red", "blue"))+
  labs(x = "Number of parameters with monotonicity", y = "Cumulative time for each step (sec)", fill = "Step"); plotC


cowplot::plot_grid(plotA, plotB, plotC, nrow = 1, labels = LETTERS)
