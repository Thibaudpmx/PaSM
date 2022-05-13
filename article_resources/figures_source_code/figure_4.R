## ---------------------------
## Script name: Figure R article
##
## Purpose of script: production of figure 4
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
## Plots require first to generate data, which in total
## take long time (especially du to the X5 reproduciton of analyses)
## Therefore, the script contains first data generation store in a desired folder
## Then the recuperation of these data to perform th eplots
##
## ---------------------------

library(QSPVP)
library(cowplot)

# Please provide the folder where you want the generated data to be stored
root <- "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data"

#  plot A - data generation: Simeoni analysis -----------------------------------------------

# Goal: vary percentage of acceptance/rejection and number of varying parameter
# Method: generate 200.000 VPs, compute them manually, then use the quantiles function
# to adjust the target, determining the pct of acceptance needed


# First, put where you want to save all the data (to heavy to put in GitHub, you will have to regenerate them yourself ! )
new_wd <- file.path(root, "Simeoni")
if(!file.exists(new_wd))  dir.create(new_wd)

setwd(new_wd)

# We first need a function to automatically create the 200.000 VPs according the number of dimension, so here it is
# Of note, this function contains a internal seed to be fully reproducible

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

# First, here are the percentage of acceptation desired
pcttargets <- c(1,0.99, 0.875, 0.75,0.5,0.25,0.125,0.01,0)

a <- 5; b <- 1

for(a in 1:5){ # For each number of varying parameter to try

  cat(paste0("a = ", a))

  cohort <- cohort_creator(nmodif = a) # create the cohort of 200.000 VPs with according dimensoin

  self <- VP_proj_creator$new() # create a new VP project
  target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min = -1E99, max = 1E99)
  self$set_targets(manual =target) # Note: can't use (-) Inf or it will set the information to TRUE

  # name of the file containing all simulations
  namepct <- paste0("Ref_",a, ".RDS") # Name of all simulations file

  # if it does not exists, create it
  if(!file.exists(namepct)){ # If it does not exist yet, compute it

    self2 <- self$clone(deep = T) # Create a clone
    self2$add_VP(cohort, use_green_filter = F, npersalve = 2000) # Compute, not using the green filter
    self2$poolVP %>%
      unnest(simul) %>%
      filter(time %in% target$time) %>% # take only the tumVol at the time of interest
      saveRDS(file = namepct )

    rm(self2)
  }

  # Read the fil containing all the simulations
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

  # If a == 5, we also want the alternative methods, so here it is
  if(a == 5){

    for(b in pcttargets){ # for every percentage

      for(d in 1:5){ # because each analyse repeated five time

        newnames <- paste0(a,"_", b, "_", d, "_alt.RDS") # compute the name of the file (nparamvarying_pcttarget_iteration)

        if(!file.exists(newnames)){ # If the file does not exist yet
          # determine new targets (with alternative methods)
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
# To do so, we create a funciton that analyse the time and other statistics for each analysis (object)
timetable <- function(self){

  tt <- self$timeTrack

  total <-tt$tTOTAL


  loop <- tt$poolVP_compteur %>%
    mutate( Treduc_filter_neg_both  =  Treduc_filter_neg_below +  Treduc_filter_neg_above) %>%
    mutate( RedFilter = TapplyRedFilterBoth + Treduc_filter_neg_both, RxODE = timemodel) %>%
    mutate(patientRemovedperSec = as.double(NremovedRedFilter) / as.double((RedFilter + timemodel )),
           refpatientRemovedperSec =  as.double(nsimul) / as.double((timemodel )))

  namesCol <- names(loop)

  # Different wether the green filter was activated or not
  if(grepl("Treduc_filter_green",  namesCol) %>% sum > 1){

    naemesReduc <- namesCol[grepl("Treduc_filter_green", namesCol)]
    naemesReduc2 <- namesCol[grepl("TapplyGreenFilter", namesCol)]



    loop[ , naemesReduc] %>% map_df(~ as.double(.x)) %>% sum(na.rm = T)

    loop <- loop %>%
      mutate( Treduc_filter_pos_both  =  !!parse_expr(paste0(naemesReduc, collapse = "+"))) %>%
      mutate( TapplyGreenFilter  = !!parse_expr(paste0(naemesReduc2, collapse = "+")),
              GreenFilter =  Treduc_filter_pos_both + TapplyGreenFilter)



  }else{ # if not activated, then set the time of green filter to 0

    loop <- loop %>%
      mutate( Treduc_filter_pos_both  =  0) %>%
      mutate( TapplyGreenFilter  =  0,
              GreenFilter =  0)

  }

  # Gather the time for red flter, green filter and RxODE
  loop %>%
    select(RedFilter, GreenFilter, RxODE) %>%
    gather("step", "value") %>%
    group_by(step) %>%
    summarise(sum = sum(value, na.rm = T)) %>%
    mutate(sum = as.double(sum)) -> TimeSimplified

  # And finally add the remaining time as "other"
  TimeSimplified %>%
    add_row(step = "Other", sum = as.numeric(total, units = "days") * (24 * 3600) - sum(TimeSimplified$sum))

}


# Now that we have the time analysing function, lets read all the files previously computed

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

    # we also add different metrics
    # (of note, not all of them are used in th eplots)
    timetable(obj) %>%
      spread(step, sum) %>%
      mutate(VPfound= obj$poolVP %>% nrow, nsimul = obj$timeTrack$poolVP_compteur %>% pull(nsimul) %>% sum,
             timeTotal = obj$timeTrack$tTOTAL %>% as.double(), Tgreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalGreenFilter)%>% as.double(),
             Tred1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalRedFilter)%>% as.double(),
             nremRed1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(NremovedRedFilter)%>% as.double(),
             ndoneGreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(nextrapoGreen)%>% as.double())


  })) %>%
  unnest() %>%
  saveRDS("full_analysis.RDS") # and finally save the results in one giant files



# Time reference ----------------------------------------------------------
# We also need to compute the reference time
# First, put where you want to save all the data (to heavy to put in GitHub, you will have to regenerate them yourself ! )
new_wd <- file.path(root, "Simeoni_ref")
if(!file.exists(new_wd))  dir.create(new_wd)

setwd(new_wd)


# here is the function used to perform manually the analyses
manual <- function(cohort, self){



  t0 <- Sys.time()
  cohort <- crossing(cohort,  proto  = unique(self$targets$protocol))


  self$protocols[] %>%
    bind_rows() %>%
    mutate(proto = names( self$protocols)) -> protocols

  eventsadmin  <- crossing(id = 1:2000, proto  = unique(self$targets$protocol)) %>%
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

timeModel <- double()


  for(a in unique(demo$group)){


    x <- demo %>% filter(group == a) %>%
      rowid_to_column("id")

    events <- eventsadmin %>%
      left_join(x%>% select(id, ids, group, proto) , by = c("id", "proto")) %>%
      filter(!is.na(ids))


    tmodel <- Sys.time()
    res <- self$model$solve(x %>% select(-proto), events) %>%
      as_tibble
    timeModel <- c(timeModel, difftime(Sys.time(),tmodel, units = 's' ) %>% as.double)

    if(max(x$id) == 1) res <- res %>% mutate(id = 1)

    res <- res %>%
      left_join(x %>% select(id, proto, ids), by = "id")

    res %>%
      filter(time %in% self$targets$time) %>%
      rename(protocol = proto) %>%
      gather("cmt", "value", unique(self$targets$cmt)) %>%
      left_join(self$targets, by = c("time", "protocol", "cmt")) %>%
      filter(value > max | value < min) %>% pull(ids ) -> idtorem

    idtorems <- c(idtorems, idtorem)


    resultsap[[a]]<-  res %>%   filter( ! (ids %in% idtorem))



  }

  resultsap <- bind_rows(resultsap)

  demo <- demo %>%
    filter(! (ids %in% idtorems)) %>%
    left_join( resultsap, by = c("ids", "proto") )


  TTotal <- difftime(Sys.time(),t0, units = "s") %>% as.double


  return(tibble(timeModel = sum(timeModel), timeTotal = TTotal))

}

# Here are the target to try
targets <- c(0,0.5,1)
# And the cohort (number of dimension not impactfull)
cohort <- cohort_creator(nmodif = 5) # doesn't matter

namepct <- file.path(root, "Simeoni", "Ref_5.RDS")

alltumVol <- readRDS(namepct)  %>% pull(tumVol)

# Perfom the same similar loop as before to compute the reference times
for(a in targets){

  for(b in 1:5){


# Name of all simulations file

    newnames <- paste0(a,"_", b, ".RDS") # compute the name of the file (nparamvarying_pcttarget_iteration)

    if(!file.exists(newnames)){ # If the file does not exist yet

      self <- VP_proj_creator$new()
      quant <- ( 1-a)/2
      target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48,
                        min = quantile(alltumVol, probs = quant), max =  quantile(alltumVol, probs = 1 - quant))
      self$set_targets(manual =target) # Note: can't use (-) Inf

      res <- manual( cohort, self)
      saveRDS(res, file = newnames)
    }

}

}


# And similar data manipulation to have the time of references
toread <- list.files()

tibble(file = toread[!grepl("timeReference", toread)]) %>%
  mutate(a = map(file, ~readRDS(.x))) %>%
  unnest() %>%
  mutate(pct = gsub("_.+", "", file)) %>%
  group_by(pct) %>%
  summarise(meanModel = mean(timeModel), medianModel = median(timeModel),minModel = min(timeModel),
            maxModel = max(timeModel),
            meanTotal = mean(timeTotal), medianTotal = median(timeTotal),minTotal = min(timeTotal),
            maxTotal = max(timeTotal)) %>%
  mutate(pct = as.double(pct)) %>%
 saveRDS("timeReference.RDS")


#  plot A - making: Simeoni analysis -----------------------------------------------

allTimes <- readRDS(file.path(root, "Simeoni", "full_analysis.RDS" ))

ref <-  readRDS(file.path(root, "Simeoni_ref", "timeReference.RDS" ))

# Verifying the percentage

allTimes %>%
  mutate(pcteffe = VPfound / 200000) %>%
  mutate(dif = pct  -  pcteffe) %>%
  group_by(dif) %>%
  tally
  # note: all the small of differences can be explained by internal rounding
  # After check the additional/missing VP are from 1E-10 mm3 or similar

# Main plot time benefice

plotA <- allTimes %>%

  mutate(pct = as.double(pct)) %>%
  mutate(nparam = as.character(nparam)) %>%
  mutate(meth = if_else(meth == "", "Centered", "Lowest")) %>%
  group_by(nparam, pct, meth) %>% # compute the median times of computatoin  (from x5 analyses)
  summarise(timeTotal = median(timeTotal)) %>%
  ungroup() %>%
  ggplot()+
  geom_rect(aes(xmin = 48, xmax = 52, ymin = 39, ymax = 41.5), col = "red",alpha = 0)+
  geom_segment(aes(x = 56, xend = 52, y = 45, yend = 42), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+
  # geom_text(aes(x = 62, y = 46, label = "next plots"), col = "red", size = 3)+ # Note: last 3 brut way, not reproducible
  geom_point(aes(x = (1-pct) * 100, y = timeTotal, col = nparam))+
  geom_line(aes(x = (1-pct)* 100, y = timeTotal, col = nparam, lty = meth))+
  geom_line(data = ref, aes(x= (1- pct)*100, y = medianTotal ), lty = 2)+
  geom_ribbon(data = ref, aes(x= (1- pct)*100, ymin = minTotal, ymax = maxTotal, fill = "Time of\nreference" ),alpha = 0.2)+
  geom_line(data = ref, aes(x= (1- pct)*100, y = minTotal ))+
  geom_line(data = ref, aes(x= (1- pct)*100, y = maxTotal ))+
  theme_bw()+
  scale_linetype_manual(values = c(1,2))+
  scale_fill_manual(values = "grey")+
  labs(x = "Percentage of rejected VP", y = "Total Time analysis (sec)", col = "Number\nvarying\nparameters", fill = "",
       lty = "Target\nmethod"); plotA



# Sinuosity ---------------------------------------------------------------
# Number of sinuosity (filters) written in the article

normal <-readRDS(file.path(root, "Simeoni", "5_0_5.RDS" ))
normal$n_filter_reduc()
normal # see number of filter after reduction

alt <- readRDS(file.path(root, "Simeoni", "5_0_5_alt.RDS" ))
alt$n_filter_reduc()
alt# see number of filter after reduction

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



# Take as reference the 50% acceptation
ref0.5 <- ref %>% filter(pct == 0.5)
ref0.5 <- tibble(pct = "Old", time = c("Other", "RxODE"), value = c(ref0.5$medianTotal - ref0.5$medianModel,ref0.5$medianModel)) %>% crossing(nparam = 5)

plotB <- allTimes %>%
  filter(nparam == 5 & pct == 0.5 & meth == "") %>%
  left_join(iteration_to_keep) %>%
  filter(test) %>%
  gather("time", "value",GreenFilter, Other, RedFilter, RxODE) %>%
  mutate(pct = "New") %>%
  bind_rows(ref0.5 ) %>%
  group_by(pct) %>%
  mutate(total = sum(value)) %>%
  ungroup() %>%
  ggplot()+
  geom_segment(aes(x = "New", xend = "New", y = 53, yend = 50), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+
  geom_col(aes(pct, value, fill = fct_reorder(time, value, .desc = F)), alpha = 0.5)+
  theme_bw()+
  scale_fill_manual(values = c("grey", "red", "darkgreen", "blue"))+
  labs(x = "Method", y = "Time of analysis (sec)", fill = "Step")+
  geom_text(aes(pct, value, fill = fct_reorder(time, value, .desc = F), label = as.double(value) %>% round(1)), position = position_stack(vjust = .5))+
  geom_label(aes(pct, total+ 5 ,  label = paste0("Total:",as.double(total) %>% round(1), "s"))); plotB

# cowplot::plot_grid(plotA, plotB, nrow =   1)



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
  ggplot()+
  geom_col(aes(pct, value, fill = fct_reorder(time, value, .desc = F)), alpha = 0.5)+
  theme_bw()+
  facet_wrap(~nparam)+
  scale_fill_manual(values = c("grey", "blue", "red", "darkgreen"))+
  labs( x= "Percentage of rejected VP", y = "Time of analysis (sec)", fill = "Step")+
  geom_text(aes(pct, value, fill = fct_reorder(time, value, .desc = F), label = as.double(label) %>% round(1)), position = position_stack(vjust = .5)); suplementA



# plot C ------------------------------------------------------------------
temp <- allTimes %>%
  # mutate(meth = if_else(meth == "", "Centered", "Lowest")) %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  group_by(nparam, pct, meth) %>%
  mutate(median = median(timeTotal)) %>%
  ungroup() %>%
  filter(timeTotal == median) %>%
  mutate(nparam = paste0(nparam, if_else(meth == "", "", "alt")))


plotC <- temp %>%
  ggplot()+
  geom_point(aes(nparam,nsimul, col = factor(1-pct)))+

  geom_line(data = temp %>% filter(nparam != "5alt"), aes(nparam, nsimul, col =factor(1-pct), group = pct))+
  geom_line(data = temp %>% filter(nparam %in% c("5", "5alt")), aes(nparam, nsimul, col =factor(1-pct), group = pct), lty = 2)+
  theme_bw()+
  geom_segment(aes(x = 5.5, xend = 5.1, y = 49000, yend = 47532), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+

  scale_y_log10()+
  labs(x = "Number of varying parameter", y = "Number of simulations performed",col = "Percentage\nRejection");plotC


# plot D ------------------------------------------------------------------
temp <- allTimes %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  group_by(nparam, pct, meth) %>%
  mutate(median = median(timeTotal)) %>%
  mutate(filtreTime = (Tgreen1 + Tgreen1)) %>%
  ungroup() %>%
  filter(timeTotal == median) %>%
  mutate(nparam = paste0(nparam, if_else(meth == "", "", "alt")))

plotD <- temp %>%
  ggplot()+

  geom_point(aes(nparam,filtreTime, col = factor(1-pct)))+
  geom_segment(aes(5.4, xend = 5.1, y = 16.5, yend = 14.7), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+

  geom_line(data = temp %>% filter(nparam != "5alt"), aes(nparam, filtreTime, col =factor(1-pct), group = pct))+
  geom_line(data = temp %>% filter(nparam %in% c("5", "5alt")), aes(nparam, filtreTime, col =factor(1-pct), group = pct), lty = 2)+
  theme_bw()+
  scale_y_log10()+
  labs(x = "Number of varying parameter", y = "Filtering time first iteration (sec)",col = "Percentage\nRejection");plotD



# plot E - data generation ----------------------------------------------------------------

setwd( file.path(root, "time_impact"))

# Times to try to print
times_to_try <- list(seq(0,48,48), seq(0,48,4),seq(0,48,2),seq(0,48,1),seq(0,48,0.5),seq(0,48,0.1))


# Same set-up as reference: 5 dimension
cohort <- cohort_creator(nmodif = 5)
alltumVol <- readRDS(file.path(root, "Simeoni/Ref_5.RDS"))%>%
  pull(tumVol)


# determine new targets using the center 50% acceptaiton set-up (quantiles 25 to 75)
prototemp <- tibble(protocol = "dose50", cmt = "tumVol", time = 48, min =  quantile(alltumVol, probs = 0.25),
                    max = quantile(alltumVol, probs = 0.75))# and apply the new min and max

# Compute all files
for(a in 1:length(times_to_try)){ # for each time to try


  for(b in 1:5){


    namefile <- paste0(a, "_", b,".RDS" )

    if(!file.exists(namefile)){

      self <- VP_proj_creator$new()
      self$set_targets(manual = prototemp)
      self$times <- times_to_try[[a]] # use this time instead of the default one

      self$add_VP(VP_df = cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0, saveVPRej = F, time_compteur = T)
      self$timeTrack$ttotal

      saveRDS(object = self, namefile)
    }



  } # end for b

} # end for a

# Then perform the analysis (timetable funciton needed)
# very similar to plot A so no need to reexplain

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


# plot E - plot  -------
setwd( file.path(root, "time_impact"))

# same list
times_to_try <- list(seq(0,48,48), seq(0,48,4),seq(0,48,2),seq(0,48,1),seq(0,48,0.5),seq(0,48,0.1))
target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min = -1E99, max = 1E99)

# Compute the time to of RxODE solving for each set-up, first by creating a function
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

# then mapping this function of all times-to try: this is our reference time !
allTimesRxODE <-  map(times_to_try, function(x){

  time_Impact(x) * 100

})

names(allTimesRxODE) <- map(times_to_try, ~ length(.x) %>% as.character)

# gather all times
temp <- readRDS("full_analysis.RDS") %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  group_by(Step) %>%
  mutate(median = median (timeTotal)) %>%
  filter(median == timeTotal) %>%
  ungroup() %>%
  mutate( n = map_dbl(times_to_try, ~ length(.x)) ) %>%
  mutate(timecomput =  map_dbl(allTimesRxODE, ~ as.double(.x)/ 100) ) %>%
  mutate(Old = 7.4 + allTimesRxODE %>% reduce(c) %>% as.double()) %>% # for reference time, we add fixed 7.4 sec
  # corresponding to the "other" time seens in  plot B (assumed the same)
  rename(New = timeTotal) %>%
  gather("method", "value", Old, New) %>%
  mutate(value= as.double(value))

# Perform the plot
plotE <- temp %>%
  ggplot()+
  geom_vline(aes(xintercept = (ref$medianModel/100) %>% median, lty= "Used in\nmain\nanalysis"))+
  geom_line(aes(timecomput, value, col = method), size = 2)+
  geom_point(aes(timecomput, value, col = method), size = 3)+
  geom_segment(aes(x = 0.52, xend = 0.45, y = 36, yend = 39), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+

  # geom_line(aes(timecomput, value, col = method), size = 2)+
  geom_ribbon(data = temp %>% spread(method, value), aes(x = timecomput, ymin = New, ymax =  Old), alpha = 0.3)+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Time to compute 2000 VPs with RxODE (sec)",y = "Time for performing 200.000 VPs (sec)", lty = "", col = "Method")+
  scale_linetype_manual(values = 2);plotE


plot_grid(plotA, plotB, plotC, nrow = 1)



# Plot F data generation------------------------------------------------------------------
setwd( file.path(root, "ImpactSizeCohort"))

# Step 1: create a 1 million database
# was complex in terms of overlaod (R crashing,...), so we divided in several sections

base <- crossing( k2 = 0.2,
                  lambda0 = 0.1,
                  Vd =  40,
                  lambda1 = 12,
                  ke = 1
)


sizeTotal <- 4E6
nperparam <- ceiling(sizeTotal^(1/5)) # count how many parameter values per paramer, superior round


list <-   map(1:5, function(x){ # for each parameter, create the sequence function

  min <- base[[x]]/10
  max <- base[[x]]*10
  step <- (max - min)/(nperparam-1)

  expr(seq( !!min,  !!max, !!step  ))
})

names(list) <- names(base[1:5])

cohort <- crossing(!!!list) %>%
  slice(1:sizeTotal) %>%
  mutate(psi = 20, w0 = 50, k1 = 0.5)

target <-  tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min = -1E99, max = 1E99)

# Here is the division we made
for(a in 0:(sizeTotal/5E5 - 1) ){


  if(!file.exists(paste0("self_", a))){
    self <- VP_proj_creator$new() # create a
    self$set_targets(manual =target )
    self$add_VP(cohort %>% slice((a * 500000+1):((a+1) * 500000)), use_green_filter = F, npersalve = 2000,timeSave = 48, keep = "tumVol")

    self$poolVP %>%
      unnest() %>%
      saveRDS(paste0("self_", a))
  }

}

files <- list.files()


toread <- files[grepl("^self_", files)]


map(toread, ~ readRDS(.x)) %>%
  bind_rows() %>%
  select(-rowid, - id,- tumVol_BU,- tumVol_AL) %>%
  saveRDS( "alltumVol")

rm(list = ls())

# Step 2: sample different tumvol (same process than before)
# this time changing the size of the cohort

allTumVol <- readRDS("alltumVol")


allSizeCohort <- c(10, 50 , 100,200,500,750,1000,2000,3000,4000) * 1000

for(a in allSizeCohort){
  print(a)

  for(b in 1:5){
    name <- paste0(a,"_", b,".RDS")

    if(!file.exists(name)){

      set.seed(b)

      cohort <- allTumVol %>%
        sample_n(a) %>%
        select(-time) %>%
        mutate(k1 = 0.5, psi = 20, w0 = 50)


      target <- tibble(protocol = "dose50", cmt  = "tumVol", time = 48, min =  quantile(cohort$tumVol  , probs = if_else(b == 2, 0.49, 0.25)),
                       max =  quantile(   cohort$tumVol, probs =  if_else(b == 2, 0.51, 0.75) )) # and apply the new min and max


      self <- VP_proj_creator$new() # create a new VP project
      self$set_targets(manual =target) # No
      self$add_VP(VP_df = cohort, use_green_filter = T, time_compteur = T, npersalve = 2000, pctActivGreen = 0, saveVPRej = F)

      self$n_filter_reduc()
      saveRDS(self,name )
    }
  }# end for b
}


# Step 3: compute the stats

toread <- list.files()


toread <- toread[! grepl("(self)|(alltumVol)|(full_analysis)", toread)]


str_split(toread, pattern = "_", simplify = T) %>%
  as_tibble() %>%
  mutate(V2 = gsub("\\.RDS","", V2)) %>%
  rename(nCohort = V1, iteration = V2) %>%
  mutate(file = toread) %>%
  mutate(nCohort = as.double(nCohort)) %>%
  mutate(results = map(file, function(x){

    # x <- "10000_1.RDS"
    print(x)
    obj <- readRDS(x)

    obj$n_filter_reduc()
    # saveRDS(obj, x)

    timetable(obj) %>%
      spread(step, sum) %>%
      mutate(VPfound= obj$poolVP %>% nrow, nsimul = obj$timeTrack$poolVP_compteur %>% pull(nsimul) %>% sum,
             timeTotal = obj$timeTrack$tTOTAL %>% as.double(), Tgreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalGreenFilter)%>% as.double(),
             Tred1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(TimeTotalRedFilter)%>% as.double(),
             nremRed1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(NremovedRedFilter)%>% as.double(),
             ndoneGreen1 = obj$timeTrack$poolVP_compteur %>%  slice(1) %>% pull(nextrapoGreen)%>% as.double(),
             nfilterabove=nrow(obj$filters_neg_above), nfilterbelow = nrow(obj$filters_neg_below),
             sinuosity = nfilterabove + nfilterbelow)


  })) %>%
  unnest() %>%
  saveRDS("full_analysis.RDS")



# Plot F -plot------------------------------------------------------------------
setwd( file.path(root, "ImpactSizeCohort"))


datas <- readRDS("full_analysis.RDS") %>%
  arrange(nCohort) %>%
  mutate(timeTotal = case_when(nCohort %in% (c(10,50,100,200,2000,3000,1000,4000) * 1000 )& iteration == 1 ~ timeTotal /60,
                               T ~timeTotal)) %>%
  group_by(nCohort) %>%
  summarise(timeTotal = median(timeTotal), nsimul = median(nsimul))


plotF <-  datas %>%
  mutate(Old = 50.1 * nCohort / 2E5 / 60) %>% # considering a strict proportionality
  rename(New = timeTotal) %>%
  mutate(time =  round(New / Old,1)) %>%
  {pregather <<- .} %>%
  gather("Method", "value", Old, New) %>%
  mutate(pctExtra = round((nCohort-nsimul)/ nCohort * 100 )) %>%

  mutate(label = if_else(Method == "Old", "",as.character(time))) %>%
  ggplot()+
  scale_y_log10()+
  scale_x_log10()+
  # stat_function(fun = function(x)55.1 * x / 2E5, aes(col = "Blue") , size = 2)+
  geom_vline(aes(xintercept = 2E5, lty = "Used in\nmain\nanalysis"))+
  # geom_segment(aes(x = 2.6E5, xend = 2.2E5, y = 35, yend = 42), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+
  # geom_text(aes(nCohort, value* 1E19, label = label), nudge_y = -20)+
  geom_line(aes(nCohort, value, col = Method), size = 2) +
  geom_ribbon(data = pregather, aes(x = nCohort, ymin = Old, ymax = New), alpha = 0.3)+
  geom_point(aes(nCohort, value, col = Method), size = 3)+
  theme_bw() +
  # geom_segment(aes(x = 2.6E5, xend = 2.2E5, y = 0.36, yend = 0.44), col = "red",  arrow = arrow(length = unit(0.2, "cm")))+

  # facet_wrap(~iteration)+
  scale_linetype_manual(values = 2)+
  labs(x = "Number of VPs in original cohort", y = "Time of computation (minute)");plotF



# datas %>%
#   mutate(Old = 55.1 * nCohort / 2E5 / 60) %>%
#   mutate(timeTotal = case_when(nCohort <= 200000 & iteration == 1 ~ timeTotal /60,
#                                nCohort <= 1000000 & iteration == 2 ~ timeTotal /60,
#                                T ~ timeTotal)) %>%
#   rename(New = timeTotal) %>%
#   mutate(New/Old) %>%
#   mutate(TperID = RxODE / nsimul ) %>%
#   select(TperID, everything()) %>%
#   arrange(iteration, nCohort)


# Final merge for fig 4 ---------------------------------------------------

library(cowplot)

plot_grid(plotA, plotB, plotC, plotD, plotE, plotF, labels = LETTERS, nrow = 2)



# Supplemental Graph ------------------------------------------------------



setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/ImpactSizeCohort")


datas <- readRDS("full_analysis.RDS") %>%
  arrange(nCohort) %>%
  mutate(timeTotal = case_when(nCohort %in% (c(10,50,100,200,2000,3000,1000,4000) * 1000 )& iteration == 1 ~ timeTotal /60,
                               T ~timeTotal)) %>%
  group_by(nCohort) %>%
  mutate(mediant = median(timeTotal)) %>%
  filter(timeTotal == mediant)

plotSupA <- datas %>%
  mutate( Old = 55.1 * nCohort / 2E5 / 60) %>%
  mutate(New = timeTotal) %>%
  mutate(timeTotal=timeTotal*60) %>%
  mutate(time =  round(New / Old,1)) %>%
  # filter(iteration == 1) %>%
  select(-iteration) %>%
  mutate(delta = nCohort - nsimul, Rxodeid =  RxODE /nsimul ) %>%

  mutate(TodeSaved = delta * Rxodeid, deltaratio = delta/nCohort ) %>%
  select(TodeSaved, Rxodeid,delta, deltaratio, everything()) %>%
  gather("step", "value", TodeSaved, GreenFilter, Other, RedFilter,RxODE, timeTotal ) %>%
  # filter(step == "t2") %>%
  mutate(value = value / nCohort) %>%
  # filter(nCohort >= 2E5) %>%
  ggplot(aes(x = nCohort, y = value, col = step))+
  geom_point()+
  geom_line()+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  labs(x= "Size cohort (VPs)", y = "Time per size cohort (sec / VP)"); plotSupA


plot_grid(
  plotSupA,
  datas %>%
    ggplot()+
    scale_x_log10()+
    geom_point(aes(nCohort, sinuosity))+
    geom_line(aes(nCohort, sinuosity))+
    theme_bw()+
    labs(x = "Size cohort (VPs)", y = "Sinuosity (number of filters)"),


  datas %>%
    mutate(ratioRxode =  nsimul/nCohort) %>%
    ggplot()+
    scale_x_log10()+
    geom_point(aes(nCohort, ratioRxode * 100))+
    geom_line(aes(nCohort, ratioRxode*100))+
    theme_bw()+
    labs(x = "Size cohort (VPs)", y=  "Percentage of VPs with ODE solved"),

  nrow = 1, labels = LETTERS, rel_widths = c(1.4,1,1)

)
