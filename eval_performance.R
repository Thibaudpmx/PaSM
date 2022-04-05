library(peccary)
library(QSPVP)
library(RxODE)
library(progress)
library(R6)
library(crayon)
library(profvis)
library(microbenchmark)



# Compute our methodology

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
# cohorts definitions --------------------------------

cohorts <- list()

cohort_creator <- function(nmodif){

  base <- crossing( k2 = 1,
                    lambda0 = 0.1,
                    Vd =  40,
                    lambda1 = c(12),
                    ke = 1 ,#*  seq(0.6,1.4,0.2),
                    w0 = 50
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
   mutate(k1 = 0.5, psi = 20)


}


cohort_creator(6)


# Targets -----------------------------------------------------------------

Targets <- function(proto = 1, cmt = c("tumVol"), time = c(45) ){

  proto <- c("dose50", "dose100", "dose0")[proto]
  cmt <- c("tumVol", "Conc")[cmt]

  crossing(protocol = proto, cmt = cmt , time = time) %>%
    mutate(min = -Inf, max = Inf)


}

Targets(proto = 1:3,cmt =  1:2,time =  c(25,30,50))



# Manual do it ------------------------------------------------------------

# target <- Targets(proto = 1:2,cmt =  1,time =  c(48))
# cohort <- cohort_creator(nmodif = 2)
self <- VP_proj_creator$new()

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
  mutate(group = ceiling(ids/2000)) %>%
  # filter(group == 2) -> x
  group_split(group) %>%
  map(function(x){



    x <-  x %>%
      rowid_to_column("id")

     events <- eventsadmin %>%
       left_join(x %>% select(id, ids, group, proto), by = c("id", "proto")) %>%
       filter(!is.na(ids))



    res <- self$model$solve(x %>% select(-proto), events , c(X2 = 0)) %>%
      as_tibble

    res %>%
      filter(time %in% target$time) %>%
      left_join(x %>% select(id, proto), by = "id") %>%
      rename(protocol = proto) %>%
      gather("cmt", "value", unique(target$cmt)) %>%
      left_join(target, by = c("time", "protocol", "cmt")) %>%
      filter(value > max | value < min) %>% pull(id) -> idtorem

    x %>%
      filter(! (id %in% idtorem)) %>%
      left_join(res %>%   filter( ! (id %in% idtorem))  %>% group_by(id) %>% nest(), by = "id")

  }) %>%
  invoke(.fn = bind_rows)
# difftime(Sys.time(),t00, units = "s")

}



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

manual(target, cohort)

# Main algorithm to compute -----------------------------------------------

target <- Targets(proto = 1,cmt =  1,time =  c(48))

target$min <- -1E99
target$max <- 1E99

cohort <- cohort_creator(nmodif = 2)
self$set_targets(manual = target)

self$add_VP(VP_df = cohort)



# Compute the reference
t0 <- Sys.time()
ref <- manual(target, cohort)
t0 <- difftime(Sys.time(),t0)

mbref <- microbenchmark(ref <- manual(target, cohort), times = 5)

# saveRDS(mbref, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS")
# readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS")

# Compute our methodology

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Simeoni")

pcttargets <- c(1, 0.75,0.5,0.25,0.125,0.01,0)

for(a in 1:5){
print(a)
cohort <- cohort_creator(nmodif = a)
self <- VP_proj_creator$new()
self$set_targets(manual = target)

namepct <- paste0("simulAlln_",a, "_proto_1_cmt_1_time1.RDS")
if(!file.exists(namepct)){

  self2 <- self$clone(deep = T)
  self2$add_VP(cohort, use_green_filter = F, npersalve = 2000)
  self2$poolVP %>%
    unnest(simul) %>%
    filter(time %in% target$time) %>%
    saveRDS(file = namepct )

  rm(self2)
}

allpct <- readRDS(namepct)

alltumVol <- allpct %>% pull(tumVol)

for(b in pcttargets){

  newnames <- paste0("tRG_n",a,"_", b, "pct_proto_1_cmt_1_time1.RDS")
  newnamessingle <-  paste0("SingletRG_n",a,"_", b, "pct_proto_1_cmt_1_time1.RDS")

  if(!file.exists(newnamessingle)){
    # determine new targets
    quant <- ( 1-b)/2
    prototemp <- target
    prototemp$min <- quantile(alltumVol, probs = quant)
    prototemp$max <- quantile(alltumVol, probs = 1 - quant)


    self2 <- self$clone(deep = T)
    self2$set_targets(manual = prototemp)


      self3 <- self2$clone(deep = T)
      set.seed(123)
      self3$add_VP(time_compteur = T, cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0.6)


      saveRDS(self3, file = newnamessingle)
}

#   if(!file.exists(newnames)){
#   # determine new targets
#   quant <- ( 1-b)/2
#   prototemp <- target
#   prototemp$min <- quantile(alltumVol, probs = quant)
#   prototemp$max <- quantile(alltumVol, probs = 1 - quant)
#
#
# self2 <- self$clone(deep = T)
# self2$set_targets(manual = prototemp)
#
# test <- function(){
#   self3 <- self2$clone(deep = T)
#   self3$add_VP(cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0.1)
# }
#
# mbref <- microbenchmark(test(), times = 5)
#
# saveRDS(mbref, newnames)
#   }

}# end for pcttargets

}




# Exploring the results ---------------------------------------------------





files <- list.files("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data")
setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data")

toread <- files[grep("^tRG_n", files)]

times <- map(toread, ~ readRDS(.x))
names(times) <- toread

readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS") %>%
  as.data.frame() %>%
  mutate(time = time * 10E-10) %>%
 summarise(min = min(time), max = max(time), median = median(time)) -> ref

allTimes <-

  map(times, function(x) x %>% as.data.frame) %>%
  bind_rows(.id = "id") %>%
  mutate(time = time * 10E-10) %>%
  group_by(id) %>%
  summarise(min = min(time), max = max(time), median = median(time)) %>%
  mutate(pct = gsub("pct.+", "", id)) %>% mutate(pct = gsub(".+_", "", pct) %>% as.double()) %>%
  mutate(id = gsub("tRG_n", "", id)) %>%
  mutate(id = gsub("_.+", "", id))




allTimes %>%
  ggplot()+
  geom_point(aes(x = (1-pct), y = median, col = factor(id)))+
  geom_line(aes(x = (1-pct), y = median, col = factor(id)))+
  geom_hline(data = ref, aes(yintercept = min), lty = 1)+
  geom_hline(data = ref, aes(yintercept = max), lty = 1)+
  geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.2)+
  geom_hline(data = ref, aes(yintercept = median), lty = 2)+
  geom_hline(yintercept =  16.7)+
  # geom_errorbar(data = ref, aes(x = 0, ymin = min, ymax = max))+
  theme_bw()

# stats for PAGE abstract

allTimes %>%
  summarise(min = min(median), max = max(median), mean = mean(median), q10 = quantile(median, 0.1), q90 = quantile(median, 0.9), median = median(median))

# allTimes %>%
  # filter(median > 16.7)

# Time of Analysis --------------------------------------------------------



files <- list.files("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Simeoni")

toread <- files[grep("^SingletRG_n", files)]

x <- "SingletRG_n1_0.01pct_proto_1_cmt_1_time1.RDS"

allTimes <- temp <- tibble(a = toread) %>%
  mutate(b = gsub("SingletRG_n", "", a )) %>%
  mutate(nparam = gsub("_.+", "", b)) %>%
  mutate(b = gsub("^._", "", b)) %>%
  mutate(pct  = gsub("pct.+", "", b)) %>%
  select(-b) %>%
  # slice(1:10) %>%
  mutate(nsimul = map_dbl(a, function(x){

    obj <- readRDS(x)

    as.double(obj$timeTrack$tTOTAL)


  }))



allTimes %>%
  mutate(pct = as.double(pct)) %>%
  ggplot()+
  geom_point(aes(x = (1-pct), y = nsimul, col = nparam))+
  geom_line(aes(x = (1-pct), y = nsimul, col = nparam))+
  geom_hline(data = ref, aes(yintercept = min), lty = 1)+
  geom_hline(data = ref, aes(yintercept = max), lty = 1)+
  geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.2)+
  geom_hline(data = ref, aes(yintercept = median), lty = 2)+
  # geom_hline(yintercept =  16.7)+
  # geom_errorbar(data = ref, aes(x = 0, ymin = min, ymax = max))+
  theme_bw()


temp <- tibble(a = toread) %>%
  mutate(b = gsub("SingletRG_n", "", a )) %>%
  mutate(nparam = gsub("_.+", "", b)) %>%
  mutate(b = gsub("^._", "", b)) %>%
  mutate(pct  = gsub("pct.+", "", b)) %>%
  select(-b) %>%
  # slice(1:10) %>%
  mutate(nsimul = map_dbl(a, function(x){

  obj <- readRDS(x)

  obj$timeTrack$poolVP_compteur %>%
    pull(nsimul) %>%
    sum

  }))




temp %>%
  filter(pct != 1) %>%
  mutate(pct_done =  nsimul/200000 * 100) %>%
  summarise(min = min(pct_done), max = max(pct_done), mean = mean(pct_done), median = median(pct_done))



temp %>%
  ggplot()+
  geom_point(aes(1-as.double(pct), nsimul/200000 * 100))+
  geom_line(aes(1-as.double(pct), nsimul/200000 * 100,col = nparam))

times <- map(toread, ~ readRDS(.x))
names(times) <- toread



readRDS("tRef_proto_1_cmt_1_time1.RDS")


testSingle <- readRDS("SingletRG_n6_0.5pct_proto_1_cmt_1_time1.RDS")


testSingle$timeTrack$poolVP_compteur %>% pull(nsimul) %>% sum

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



temp2 <- tibble(a = toread) %>%
  mutate(b = gsub("SingletRG_n", "", a )) %>%
  mutate(nparam = gsub("_.+", "", b)) %>%
  mutate(b = gsub("^._", "", b)) %>%
  mutate(pct  = gsub("pct.+", "", b)) %>%
  select(-b) %>%
  # slice(1:10) %>%
  mutate(nsimul = map(a, function(x){

    obj <- readRDS(x)
    timetable(obj)

  }))


temp2 %>%
  # slice(1) %>%
  unnest() %>%
  spread(step, sum) %>%
  left_join(temp) %>%
  mutate(timeExtrapol = GreenFilter + RedFilter) %>%
  mutate(timeNORxode = timeExtrapol + Other) %>%
  mutate(timeTotal = timeNORxode + Other) %>%
  filter(pct != 1) %>%
  mutate(pctExtra = (1 - nsimul/200000) *100) %>%
  saveRDS("stattable.RDS")


temp3 <- readRDS("stattable.RDS")


temp3 %>%
  arrange(desc(timeNORxode)) # recompute 0.5 pct with 6 param !

temp3 %>%
  summarise(min(timeExtrapol),min(timeNORxode), max(timeExtrapol), max(timeNORxode), median(timeExtrapol))


library(ggpubr)

# ggboxplot(data = temp3, x = "nparam", y = "nsimul", )

# Pct Extra vs naparam

temp3 %>%
  ggplot()+
  geom_boxplot(aes(nparam, pctExtra))

temp3 %>%
  ggplot()+
  geom_point(aes(nparam, pctExtra))+
  scale_y_log10()

cor.test(as.double(temp3$nparam), temp3$pctExtra, method=c("pearson", "kendall", "spearman"))


# Pct Extra vs pct rejection


temp3 %>%
  ggplot()+
  geom_point(aes(1 - as.double(temp3$pct), pctExtra))+
  geom_line(aes(1 - as.double(temp3$pct), pctExtra))+
  scale_y_log10()+
  facet_wrap(~nparam, scales = "free")


temp3 %>%
  ggplot()+
  geom_point(aes(1 - as.double(temp3$pct), timeNORxode))+
  scale_y_log10()



# timeNORxode vs pct rejection

temp3 %>%
  ggplot()+
  geom_point(aes(1 - as.double(temp3$pct), timeNORxode))+
  geom_line(aes(1 - as.double(temp3$pct), timeNORxode))+
  scale_y_log10()+
  facet_wrap(~nparam, scales = "free")


map(1:6, ~cor.test(as.double(temp3$pct[temp3$nparam == .x]), temp3$timeNORxode[temp3$nparam == .x], method=c("pearson", "kendall", "spearman")))


# timeNORxode vs pct naparam

temp3 %>%
  ggplot()+
  geom_point(aes(1 - as.double(temp3$nparam), timeNORxode))+
  geom_line(aes(1 - as.double(temp3$nparam), timeNORxode))+
  scale_y_log10()+
  facet_wrap(~pct, scales = "free")


map(unique(temp3$pct), ~cor.test(as.double(temp3$nparam[temp3$pct == .x]), temp3$timeNORxode[temp3$pct == .x], method=c("pearson", "kendall", "spearman")))


# Impact of the time of QSP ------------------------------------------
source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/time_impact")

    target <- Targets(proto = 1,cmt =  1,time =  c(48))
target$min <- 200
target$max <- 2794

times_to_try <- list(seq(0,48,48), seq(0,48,4),seq(0,48,2),seq(0,48,1),seq(0,48,0.5),seq(0,48,0.1))
x <- times_to_try[[1]]
cohort <- cohort_creator(nmodif = 4)

allTimes <- map(times_to_try, function(x){

  self <- VP_proj_creator$new()
  self$set_targets(manual = target)
  self$times <- x

  self$add_VP(VP_df = cohort, use_green_filter = T)
  self$timeTrack$ttotal

})
names(allTimes) <- map(times_to_try, ~ length(.x) %>% as.character)
saveRDS(allTimes, file = "allTimes")

allTimes <- readRDS( "allTimes")




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



names(allTimesRxODE) <- map(times_to_try, ~ length(.x) %>% as.character)

temp <- tibble(n = map_dbl(times_to_try, ~ length(.x)), brut = allTimesRxODE %>% reduce(c),
       new  = allTimes %>% reduce(c)  ) %>%
       mutate(n2 = n , n = as.double(brut) / 2000 ) %>%
  mutate(brut = brut  * 0.7+ 7.1) %>%
  mutate(brut= as.double(brut), new = as.double(new))


readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS") %>%
  as.data.frame() %>%
  mutate(time = time * 10E-10) %>%
  summarise(min = min(time), max = max(time), median = median(time)) -> ref2


temp %>%
  rename(Old = brut, New = new) %>%
  gather("method", "value", Old, New) %>%
  mutate(value= as.double(value)) %>%
  ggplot()+
  geom_line(aes(n, value, col = method), size = 2)+
  geom_vline(xintercept = 0.035)+
  geom_line(aes(n, value, col = method), size = 2)+
  geom_ribbon(data = temp, aes(x = n, ymin = brut, ymax = new), alpha = 0.3)+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  geom_hline(data = ref2, aes(yintercept = max), lty = 1)+
  geom_rect(data = ref2, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.3)+
  geom_hline(data = ref2, aes(yintercept = median), lty = 2)+
  geom_hline(data = ref2, aes(yintercept = min), lty = 1)+
  labs(x = "Time to compute 2000 VPs with RxODE (sec)",y = "Time for performing 200.000 VPs (sec)")


temp %>%
  mutate(benefice = (brut - new)/brut * 100, b2 = brut/new)

# Time analyse ------------------------------------------------------------

# ref


timeTrack2 <- list()

t00 <- t0 <- Sys.time()

    cohort <- crossing(VP_df,  proto  = unique(self$targets$protocol)) %>%
      mutate(psi = 20)


    target <- self$targets

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


    timeTrack2$prefor <- difftime(Sys.time(), t0, units = "s")

    loop2 <-  tibble()

    tpreloop <- Sys.time()

    for(a in unique(demo$group)){


      t0 <- Sys.time()

      x <- demo %>% filter(group == a) %>%
        rowid_to_column("id")

      events <- eventsadmin %>%
        left_join(x %>% select(id, ids, group, proto) , by = c("id", "proto")) %>%
        filter(!is.na(ids))


      perrun  <- tibble(pre_simul = difftime(Sys.time(), t0, units = "s"));  t0 <- Sys.time()

      res <- self$model$solve(x %>% select(-proto), events, c(X2 = 0)) %>%
        as_tibble

      if(max(x$id) == 1) res <- res %>% mutate(id = 1)

      perrun$simul <- difftime(Sys.time(), t0, units = "s") ;  t0 <- Sys.time()


      res <- res %>%
        left_join(x %>% select(id, proto, ids), by = "id")

      perrun$post_simul_join1 <- difftime(Sys.time(), t0, units = "s") ;  t0 <- Sys.time()




      res %>%
        filter(time %in% self$targets$time) %>%
        rename(protocol = proto) %>%
        gather("cmt", "value", unique(self$targets$cmt)) %>%
        left_join(target, by = c("time", "protocol", "cmt")) %>%
        filter(value > max | value < min) %>% pull(ids ) -> idtorem

      perrun$post_simul_join2 <- difftime(Sys.time(), t0, units = "s") ;  t0 <- Sys.time()


      idtorems <- c(idtorems, idtorem)


      resultsap[[a]]<-  res %>%   filter( ! (ids %in% idtorem))


      perrun$post_simul_join3 <- difftime(Sys.time(), t0, units = "s") ;  t0 <- Sys.time()


      loop2 <-  bind_rows(loop2,perrun )
    }

    timeTrack2$loop <- loop2 ;  t0 <- Sys.time()
    timeTrack2$tloop <- difftime(Sys.time(), tpreloop, units = "s")

    resultsap <- bind_rows(resultsap)



    demo <- demo %>%
      filter(! (ids %in% idtorems)) %>%
      left_join( resultsap, by = c("ids", "proto") )


    timeTrack2$finalmerge <- difftime(Sys.time(), t0, units = "s") ;  t0 <- Sys.time()
    timeTrack2$ttotal <- difftime(Sys.time(), t00, units = "s") ;  t0 <- Sys.time()


self <- readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/SingletRG_n1_0.5pct_proto_1_cmt_1_time1.RDS")


#protocoldf <- sample next lines and compute protocol


# Main times
loop %>%
  select(timesampleline, protocoldf, timemodel, Treduc_filter_neg_both ,TapplyRedFilterBoth, TremoveNegDirect,
         Treduc_filter_pos_both, TapplyGreenFilter, timesimlandjoin,res2) %>%
# loop %>%
  # select(timesampleline, protocoldf, timemodel, res2, filter_neg_below, filter_neg_above, remneg_fil, time_addgreennofil, timesimlandjoin) %>%
  gather("step", "value") %>%
  group_by(step) %>%
  summarise(sum = sum(value, na.rm = T)) -> temp; print(temp); temp %>%
  summarise(sum(sum))



loop %>%
  select(RedFilter, GreenFilter, RxODE) %>%
  gather("step", "value") %>%
  group_by(step) %>%
  summarise(sum = sum(value, na.rm = T)) -> TimeSimplified


dataredfilter <-  temp %>% filter(step %in% c("TapplyRedFilterBoth", "Treduc_filter_neg_both")) %>%
  mutate(name = if_else(step == "TapplyRedFilterBoth", "Apply", "Reduce")) %>%
  mutate(sum = as.double(sum))

datagreenfilter <-  temp %>% filter(step %in% c("TapplyGreenFilter", "Treduc_filter_pos_both")) %>%
  mutate(name = if_else(step == "TapplyGreenFilter", "Apply", "Reduce")) %>%
  mutate(sum = as.double(sum))


fig1 <- TimeSimplified %>%
  add_row(step = "other", sum = total - sum(TimeSimplified$sum)) %>%
  ggplot()+
  geom_col(aes(fct_reorder(step, sum, .desc = T), sum, fill = step))+
  geom_text(aes(fct_reorder(step, sum, .desc = T), sum, label = round(sum,1)), nudge_y = 0.3)+
  scale_fill_manual(values = c("darkgreen", "grey", "red", "blue"))+
  theme_bw()+
  geom_col(data = dataredfilter, aes("RedFilter", sum), col = "black", alpha = 0)+
  geom_text(data = dataredfilter %>% mutate(sum2 = sum, sum = dataredfilter %>% slice(1) %>% pull(sum), lag = c(0, TimeSimplified$sum[TimeSimplified$step == "RedFilter"] )),
            aes("RedFilter", (sum + lag)/2, label = paste0(name, "\n(",  round(sum2,1),")")), col = "black")+
  geom_col(data = datagreenfilter, aes("GreenFilter", sum), col = "black", alpha = 0)+
  geom_text(data = datagreenfilter %>% mutate(sum2 = sum, sum = datagreenfilter %>% slice(1) %>% pull(sum), lag = c(0, TimeSimplified$sum[TimeSimplified$step == "GreenFilter"] )),
            aes("GreenFilter", (sum + lag)/2, label = paste0(name, "\n(",  round(sum2,1),")")), col = "black"); fig1



ref <- timeTrack2$loop %>%
  summarise(RxODE = sum(simul))

ref <- ref %>%
  mutate(other = timeTrack2$ttotal - RxODE) %>%
  gather("step", "value") %>%
  group_by(step)


TimeSimplified %>%
  add_row(step = "other", sum = total - sum(TimeSimplified$sum)) %>%
  mutate(method = "New") %>%
  bind_rows(ref %>% rename(sum = value)%>% mutate(method = "Old")) %>%
  ggplot()+
  geom_col(aes(method, sum, fill = fct_reorder(step, sum, .desc = F)), alpha = 0.4)+
  geom_text(aes(method, sum, fill = fct_reorder(step, sum, .desc = F), label = as.double(sum) %>% round(1)), position = position_stack(vjust = .5))+
  scale_fill_manual(values = c("grey", "darkgreen", "red", "blue"))+
  # geom_col(data= ref , aes(x = "Old", y = value, fill =step))+
  # geom_text(data= ref, aes("New", sum, fill = step, label = as.double(sum) %>% round(1)), position = position_stack(vjust = .5))+
  geom_label(data = mtcars, aes("New", as.double(total) %>% round(1), label = paste0("Total:", as.double(total) %>% round(1), "s")), nudge_y = 5)+
  geom_label(data = mtcars, aes("Old", as.double(ref$value %>% sum) %>% round(1), label = paste0("Total:", as.double(ref$value %>% sum) %>% round(1), "s")), nudge_y = 5)+
  theme_bw()+
  labs(fill = "step")+
  labs(caption = paste0(loop %>% pull(nsimul) %>% sum, " VPs done instaed vs 200.000 (" , loop %>% pull(nsimul) %>% sum / 200000, "%)"))

loop %>% pull(nsimul) %>% sum / 200000

# Time comptation total

tt$nVP * loop %>%
  filter(nsimul == max(loop$nsimul)) %>%
  summarise(mean = mean(timemodel)) %>%
  pull(mean) / max(loop$nsimul)


loop %>%
  filter(!is.na(NremovedRedFilter)) %>%
  ggplot()+
  geom_point(aes(n, NremovedRedFilter))+
  geom_line(aes(n, NremovedRedFilter))

# Patients removed per seconds
loop %>%
  slice(1:2) %>%
  # select(n, NremovedRedFilter, RedFilter, timesampleline)
  # filter(!is.na(NremovedRedFilter)) %>%
  ggplot()+
  geom_line(aes(n, patientRemovedperSec))+
  geom_line(aes(n, refpatientRemovedperSec))



loop %>%
  slice(1:2) %>%
  ggplot()+
  geom_line(aes(n, TSavedRedFilter))+
  geom_line(aes(n, TimeTotalRedFilter))



  loop %>%
  mutate(patientRemovedperSec = as.double(NremovedRedFilter) / as.double((RedFilter + timesampleline )),
         refpatientRemovedperSec =  as.double(nsimul) / as.double((timesampleline )))




# Impact ncohort ----------------------------------------------------------


# Lindner Main algorithm to compute -----------------------------------------------




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





#
#   cohort <- cohort_creator(nmodif = 2)
#   self$set_targets(manual = target)
#
#   self$add_VP(VP_df = cohort)
#
#
#
#   # Compute the reference
#   t0 <- Sys.time()
#   ref <- manual(target, cohort)
#   t0 <- difftime(Sys.time(),t0)
#
#   mbref <- microbenchmark(ref <- manual(target, cohort), times = 5)

  # saveRDS(mbref, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS")
  # readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS")

  # Compute our methodology

  source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

  setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/Lindner")

  pcttargets <- c(1, 0.75,0.5,0.25,0.125,0.01,0)

  for(a in 1:6){
    print(a)
    cohort <- cohort_creator_Lindner(nmodif = a)
    self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")


    self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                      "unique",90,"Pore", 0, 200

    ))

    self$times <- c(0, 50:90)

    namepct <- paste0("simulAlln_",a, "_proto_1_cmt_1_time1.RDS")
    if(!file.exists(namepct)){

      self2 <- self$clone(deep = T)
      self2$add_VP(VP_df = cohort, use_green_filter = F, npersalve = 2000, keep = "Pore",time_compteur = T )

      self2$poolVP %>%
        unnest(simul) %>%
        filter(time %in% self2$target$time) %>%
        saveRDS(file = namepct )

      saveRDS(self2, "timeRefLindner")

      rm(self2)
    }

    allpct <- readRDS(namepct)

    alltumVol <- allpct %>% pull(tumVol)

    for(b in pcttargets){

      newnames <- paste0("tRG_n",a,"_", b, "pct_proto_1_cmt_1_time1.RDS")
      newnamessingle <-  paste0("SingletRG_n",a,"_", b, "pct_proto_1_cmt_1_time1.RDS")

      if(!file.exists(newnamessingle)){
        # determine new targets
        quant <- ( 1-b)/2
        prototemp <- target
        prototemp$min <- quantile(alltumVol, probs = quant)
        prototemp$max <- quantile(alltumVol, probs = 1 - quant)


        self2 <- self$clone(deep = T)
        self2$set_targets(manual = prototemp)


        self3 <- self2$clone(deep = T)
        self3$add_VP(time_compteur = T, cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0.75)


        saveRDS(self3, file = newnamessingle)
      }

      if(!file.exists(newnames)){
        # determine new targets
        quant <- ( 1-b)/2
        prototemp <- target
        prototemp$min <- quantile(alltumVol, probs = quant)
        prototemp$max <- quantile(alltumVol, probs = 1 - quant)


        self2 <- self$clone(deep = T)
        self2$set_targets(manual = prototemp)

        test <- function(){
          self3 <- self2$clone(deep = T)
          self3$add_VP(cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0.75)
        }

        mbref <- microbenchmark(test(), times = 5)

        saveRDS(mbref, newnames)
      }
    }# end for pcttargets

  }




self2$timeTrack$tTOTAL

self2$timeTrack$poolVP_compteur %>%
  pull(timemodel) %>%
  sum()/60

View(self2$timeTrack$poolVP_compteur )

self2$timeTrack$poolVP_compteur %>%
  pull(timesimlandjoin) %>%
  sum()

