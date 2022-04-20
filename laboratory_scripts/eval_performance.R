library(peccary)
library(QSPVP)
library(RxODE)
library(progress)
library(R6)
library(crayon)
library(profvis)
library(microbenchmark)
source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")




# Main algorithm to
analyzediscrepency <- function(input = "SingletRG_n4_0.5pct_proto_1_cmt_1_time1.RDS"){


  test <- readRDS(input)

  n <- gsub("SingletRG_n", "", input) %>%
    gsub(pattern = "_.+", replacement = "")

  b <- gsub("pct.+", "", input) %>%
    gsub(pattern = ".+_", replacement = "") %>%
    as.double

  namepct <- paste0("simulAlln_",n, "_proto_1_cmt_1_time1.RDS")

  allpct <- readRDS(namepct)


  quant <- ( 1-b)/2

  tmin <- quantile(allpct$tumVol, probs = quant)
  tmax <- quantile(allpct$tumVol, probs = 1 - quant)


  allTheoRaws <- allpct %>%
    filter(tumVol > tmin & tumVol < tmax)


  # Missing
  allTheoRaws %>%
    left_join(test$poolVP %>% mutate(test = T)) %>%
    filter(is.na(test)) -> missingrows

  VP_df <- cohort_creator(as.double(n))

  prototemp <- target
  prototemp$min <- tmin
  prototemp$max <- tmax

  self <- VP_proj_creator$new()
  self$set_targets(manual = prototemp)

  missingrows$rowid ->idmi

}




# Time of Analysis --------------------------------------------------------




allTimes %>%
  filter(nparam >1) %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  mutate(pct = (200000-nsimul) / 200000) %>%
  mutate(timeExtra = GreenFilter + RedFilter) %>%
  summarise(minT = min(timeTotal), maxT = max(timeTotal),medianT =  median(timeTotal),
            minpct = min(pct), maxpct = max(pct),medianpct =  median(pct),
            mintimeExtra= min(timeExtra), maxtimeExtra = max(timeExtra),mediantimeExtra=  median(timeExtra))

# Correlation nparam and time
cor.test(as.double(allTimes$nparam), as.double(allTimes$nparam %>% as.double ), method=c("pearson", "kendall", "spearman"))

# Correlation nparam and time extrapolation
cor.test(as.double(allTimes$nparam), as.double((allTimes$GreenFilter + allTimes$RedFilter) %>% as.double ), method=c("pearson", "kendall", "spearman"))

# Correlation nparam and nsimul
cor.test(as.double(allTimes$nparam), allTimes$nsimul , method=c("pearson", "kendall", "spearman"))




# Correlation pct rejetcion and time
allTimes %>%
  filter(pct %in% c(0.01,0.5)) %>%
  select(pct, nparam, timeTotal) %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
 spread(pct , timeTotal)  -> forwilcoxowvsmed

names(forwilcoxowvsmed) <- c("nparam", "lowhigh", "medium")

wilcox.test(forwilcoxowvsmed$lowhigh, forwilcoxowvsmed$medium, paired = TRUE)

allTimes %>%
  filter(pct %in% c(0.5,1)) %>%
  select(pct, nparam, timeTotal) %>%
  mutate(timeTotal = as.double(timeTotal)) %>%
  spread(pct , timeTotal)  -> forwilcoxmedvshigh

names(forwilcoxmedvshigh) <- c("nparam",  "medium","lowhigh")

wilcox.test(forwilcoxmedvshigh$lowhigh, forwilcoxmedvshigh$medium, paired = TRUE)

combined <- bind_rows(forwilcoxmedvshigh, forwilcoxowvsmed)
wilcox.test(combined$lowhigh, combined$medium, paired = TRUE)




allTimes %>%
  filter(pct %in% c(0.01,0.5)) %>%
  ggplot()+
  geom_boxplot(aes(factor(pct), as.double(timeTotal)))+
  geom_line(aes(factor(pct), as.double(timeTotal), group = nparam))


# Correlation pct rejetcion and nsimul
allTimes %>%
  filter(pct %in% c(0.01,0.5)) %>%
  select(pct, nparam, nsimul) %>%
  mutate(v = as.double(nsimul)) %>%
  spread(pct , nsimul)  -> forwilcoxowvsmed

names(forwilcoxowvsmed) <- c("nparam", "lowhigh", "medium")

wilcox.test(forwilcoxowvsmed$lowhigh, forwilcoxowvsmed$medium, paired = TRUE)

allTimes %>%
  filter(pct %in% c(0.5,1)) %>%
  select(pct, nparam, nsimul) %>%
  mutate(nsimul = as.double(nsimul)) %>%
  spread(pct , nsimul)  -> forwilcoxmedvshigh

names(forwilcoxmedvshigh) <- c("nparam",  "medium","lowhigh")

wilcox.test(forwilcoxmedvshigh$lowhigh, forwilcoxmedvshigh$medium, paired = TRUE)

combined <- bind_rows(forwilcoxmedvshigh, forwilcoxowvsmed)
wilcox.test(combined$lowhigh, combined$medium, paired = TRUE)



# Correlation pct rejetcion and nextrapolation time

allTimes %>%
  filter(pct %in% c(0.01,0.5)) %>%
  select(pct, nparam,  GreenFilter, RedFilter) %>%
  mutate(nsimul = as.double( GreenFilter + RedFilter)) %>%
  select(-GreenFilter, - RedFilter) %>%
  spread(pct , nsimul)  -> forwilcoxowvsmed

names(forwilcoxowvsmed) <- c("nparam", "lowhigh", "medium")

wilcox.test(forwilcoxowvsmed$lowhigh, forwilcoxowvsmed$medium, paired = TRUE)

allTimes %>%
  filter(pct %in% c(0.5,1)) %>%
  select(pct, nparam, nsimul) %>%
  mutate(nsimul = as.double(nsimul)) %>%
  spread(pct , nsimul)  -> forwilcoxmedvshigh

names(forwilcoxmedvshigh) <- c("nparam",  "medium","lowhigh")

wilcox.test(forwilcoxmedvshigh$lowhigh, forwilcoxmedvshigh$medium, paired = TRUE)

combined <- bind_rows(forwilcoxmedvshigh, forwilcoxowvsmed)
wilcox.test(combined$lowhigh, combined$medium, paired = TRUE)


# Correlation nparam and nsimul
cor.test(as.double(allTimes$nparam), allTimes$nsimul , method=c("pearson", "kendall", "spearman"))






t.test()


allTimes %>%
  filter(! VPfound %in% c(2000,25000,50000,100000,150000,0,200000))


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




# Lindner base ------------------------------------------------------------
  source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
  self <-   VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

  self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                    "unique",90,"Pore", 99, 100

  ))

  self$times <- c(0, 50:90)
cohort <-  cohort_creator_Lindner(6)

self$add_VP(VP_df = cohort, use_green_filter = T, npersalve = 2000, time_compteur = F)

self$timeTrack
# Lindner Main algorithm to compute -----------------------------------------------


self$poolVP %>%
  slice(1:100) %>%
  unnest() %>%
  filter(time >50) %>%
  ggplot()+
  geom_line(aes(time, Pore))

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



# Find Noptimal -----------------------------------------------------------


self <-   VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")


self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",90,"Pore", 10, 12

))

ntotest <- 100



cohort <- cohort_creator_Lindner(nmodif = 6) %>%
  rowid_to_column("id")


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


test_Lindner_speed <- function(ntotest){
t0 <- Sys.time()

res <- self$model$solve(cohort %>% slice(1:ntotest), events %>% filter(id <= ntotest), c(X2 = 0)) %>%
  as_tibble

if(max(x$id) == 1) res <- res %>% mutate(id = 1)
as.double(difftime(Sys.time(), t0, units = "s"))
}

test <- tibble(n = c(10,100,200,500,1000,2000)) %>%
  mutate(time = map_dbl(n, ~ test_Lindner_speed(.x)))


test %>%
 mutate(pern = time/n)

# Lindner Comparison to brut force ----------------------------------------

timeTrack2 <- list()

t00 <- t0 <- Sys.time()

cohort <- cohort_creator_Lindner(nmodif = 6)

self <-   VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner_origin.r")

self$set_targets(manual = tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                                  "unique",90,"Pore", 10, 12

))

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
    rowid_to_column("id") %>%
    mutate(proto = "unique")

  events <- eventsadmin %>%
    left_join(x %>% select(id, ids, group, proto) , by = c("id", "proto")) %>%
    filter(!is.na(ids))


  perrun  <- tibble(pre_simul = difftime(Sys.time(), t0, units = "s"));  t0 <- Sys.time()
  t0 <- Sys.time()
  res <- self$model$solve(x %>% select(-proto) %>% slice(1:100), events %>% filter(id <= 100), c(X2 = 0)) %>%
    as_tibble

  if(max(x$id) == 1) res <- res %>% mutate(id = 1)

  perrun$simul <- difftime(Sys.time(), t0, units = "s") ;  t0 <- Sys.time() ; perrun$simul


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
