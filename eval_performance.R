library(peccary)
library(QSPVP)
library(RxODE)
library(progress)
library(R6)
library(crayon)
library(profvis)
library(microbenchmark)



# cohorts definitions --------------------------------

cohorts <- list()

cohort_creator <- function(nmodif){

  base <- crossing( k2 = 1,
                    lambda0 = 0.1,
                    w0 = 50,
                    Vd =  40,
                    lambda1 = c(12),
                    ke = 1 #*  seq(0.6,1.4,0.2),
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

target <- Targets(proto = 1:2,cmt =  1,time =  c(48))
cohort <- cohort_creator(nmodif = 2)
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
       left_join(x , by = c("id", "proto")) %>%
       filter(!is.na(ids))



    res <- self$model$solve(x %>% select(-proto), events, c(X2 = 0)) %>%
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
      left_join(x , by = c("id", "proto")) %>%
      filter(!is.na(ids))



    res <- self$model$solve(x %>% select(-proto), events, c(X2 = 0)) %>%
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



cohort <- cohort_creator(nmodif = 2)


# Compute the reference
t0 <- Sys.time()
ref <- manual(target, cohort)
t0 <- difftime(Sys.time(),t0)

mbref <- microbenchmark(ref <- manual(target, cohort), times = 5)

# saveRDS(mbref, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS")
# readRDS("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/tRef_proto_1_cmt_1_time1.RDS")

# Compute our methodology

source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")

setwd("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data/")

pcttargets <- c(0.75,0.5,0.25,0.125,0.01,0)

for(a in 1:6){
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
  self3$add_VP(cohort, use_green_filter = T, npersalve = 2000, pctActivGreen = 0.1)
}

mbref <- microbenchmark(test(), times = 5)

saveRDS(mbref, newnames)
}
}# end for pcttargets

}




# Exploring the results ---------------------------------------------------

files <- list.files("D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/data")

toread <- files[grep("tRG_n", files)]

times <- map(toread, ~ readRDS(.x))
names(times) <- toread

readRDS("tRef_proto_1_cmt_1_time1.RDS") %>%
  as.data.frame() %>%
  mutate(time = time * 10E-10) %>%
 summarise(min = min(time), max = max(time), median = median(time)) -> ref


map(times, function(x) x %>% as.data.frame) %>%
  bind_rows(.id = "id") %>%
  mutate(time = time * 10E-10) %>%
  group_by(id) %>%
  summarise(min = min(time), max = max(time), median = median(time)) %>%
  mutate(pct = gsub("pct.+", "", id)) %>% mutate(pct = gsub(".+_", "", pct) %>% as.double()) %>%
  mutate(id = gsub("tRG_n", "", id)) %>%
  mutate(id = gsub("_.+", "", id)) %>%
  ggplot()+
  geom_point(aes(x = (1-pct), y = median, col = factor(id)))+
  geom_line(aes(x = (1-pct), y = median, col = factor(id)))+
  geom_hline(data = ref, aes(yintercept = min), lty = 1)+
  geom_hline(data = ref, aes(yintercept = max), lty = 1)+
  geom_rect(data = ref, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax =  max), lty = 1,alpha = 0.2)+
  geom_hline(data = ref, aes(yintercept = median), lty = 2)+
  # geom_errorbar(data = ref, aes(x = 0, ymin = min, ymax = max))+
  theme_bw()


  cohort_creator(nmodif = 6) %>%
    distinct(k1)
    distinct(w0 )

