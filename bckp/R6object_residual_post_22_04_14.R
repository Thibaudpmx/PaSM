
# zone maybe --------------------------------------------------------------



# Algorithm 2 -------------------------------------------------------------

# seq(0,3,0.0015) %>% length() *
  # + seq(0,1.4,0.000820) %>% length()
# source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
#
prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100, 431),
                    max = c(100.05, 431.05))

self <- VP_proj_creator$new()

self$set_targets(manual = prototiny)

#
# domain <- tribble(~param, ~from, ~to, ~digits,
#                   "k2", 0, 3, 3 ,
#                   "lambda0", 0, 1.4, 4
#                   )
# fix <-c(k1 = 0.5, ke = 1, Vd = 40, lambda1 = c(12))


domain <- tribble(~param, ~from, ~to, ~digits,
                  "k2", 0, 3, 2 ,
                  "lambda0", 0, 1.4, 2,
                  "ke", 0, 2,1,
                   "Vd", 0,40,0,
                    "lambda1", 0,24,1
)

fix <-c(k1 = 0.5)

ndomain <- function(domain){


  domain %>%
    mutate(how = pmap_dbl(list(from, to, digits), function(from, to , digits){

      length(seq(from, to, 10^(-digits)))
    })) %>%
    pull(how) %>%
    reduce(`*`)


}

# ndomain(domain) / 2000 * 0.5 / 3600 / 24

# blocs <- zone_maybe
ndomain2 <- function(blocs){


 temp <-  blocs %>%
   slice(1:2) %>%
    rowid_to_column("id") %>%
    gather("key", "value", -id , -contains("blocsPool")) %>%
    mutate(param = gsub("(min$)|(max$)", "", key)) %>%
    mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
    select(-key) %>%
    spread(a, value) %>%
    rename(from = min, to = max) %>%
    left_join(domain %>% distinct(param, digits), by = "param")


 temp %>%
   filter(!is.na(digits)) %>%
  mutate(how = pmap_dbl(list(from, to, digits), function(from, to , digits){

    max(length(seq(from, to, 10^(-digits))) - 2,1)


  })) %>%
    group_split(id) -> temp2


temp3 <-  temp2 %>%
    map_dbl(~  .x %>% pull(how) %>%
              reduce(`*`))
# calcul bord ?

# temp2 %>% #[[1]] -> x
#   map(function(x){
#
#     x %>%
#   mutate(values = map2(to, from, ~c(.x, .y))) %>%
#   pull(values) -> valu
# names(valu )<- x$param
#
# invoke(crossing, valu)
#   })


sum(temp3)
}


# (ndomain(domain)- sum(temp2)) / 2000 * 0.5 / 3600 / 24

# sum(temp2)/ 2000 * 0.5 / 3600 / 24

VP_proj_creator$set("public", "algo2", function(domain, fix = NULL){


nVPs0  <- ndomain(domain) - 2^nrow(domain)

## First division
t0 <- Sys.time()


nperparam <- floor(200000^(1/nrow(domain))) # + 1

namesparam <- c(paste0(domain$param,  "max"),paste0(domain$param,  "min") )

DFfix <- as.data.frame(fix) %>%
  rownames_to_column() %>%
  spread(rowname , fix)

blocs <- domain %>%
  mutate(id = "")
firstbloc <- T



blocs %>%
  mutate(sampl = pmap(list(from, to, digits), function(from, to, digits){

    # from = 0; to = 3; digits = 4
    # seq(from, to, (to-from)/(nperparam-1)) %>% round(digits) %>% unique()
   temp <- seq(from, to, (to-from)/(nperparam+1)) %>% round(digits) %>% unique()

   temp[-c(1, length(temp))]

  })) -> VPsparam

# Cross parameter (per bloc) and add fixed values
newVPs <- VPsparam %>%
  group_split(id) %>%
  map( function(x){
    temp <- invoke(.fn = crossing, x$sampl )
    names(temp) <- x$param
    temp
  }) %>%
  bind_rows() %>%
  add_column(DFfix)



# Compute the new VPs
taddvp <- Sys.time()

self$add_VP(newVPs, keepRedFiltaftDis = T, reducefilteratend = F)
taddvp <- difftime(Sys.time() ,  taddvp, units = "s")

self$n_filter_reduc()
# Compute zone_maybe
self$compute_zone_maybe()

# Which become the new blocs
maybe <- self$zone_maybe

message("Reduce maybe")
maybe <- reduce_maybe2(maybe)
message("End reduce maybe")


gatherblocsPool <- 100

blocs <- maybe %>%
  rowid_to_column("blocsPool") %>%
  mutate(blocsPool = floor(blocsPool/gatherblocsPool) + 1)



nVPs <- ndomain2(blocs)

timeperun <- 0.5/2000

self$algo2list[["tree"]] <- tibble(Name = "first", size = length(unique(blocs$blocsPool)), todo = 1,
                                   before = nVPs0, after = nVPs, ratio = nVPs/nVPs0,
                                   time = difftime(Sys.time(), t0, units = "s"), what = "zoom")

self$algo2list[["domain"]] <- domain

self$algo2list[["first"]] <- blocs

nextstep <-  self$algo2list[["tree"]] %>%
               filter(!is.na(todo)) %>%
                slice(1)



nVPs/nVPs0 * 100
nVPs/2000 * 0.5 / 3600 / 24

nVPs0/2000 * 0.5 / 3600 / 24
# Now the deep dive
while( nrow(nextstep) == 1){


  t0 <- Sys.time()

  todo <- nextstep$todo

  print(nextstep)

  maybe <-  self$algo2list[[nextstep$Name]]

  if(nextstep$Name !="first") maybe <- maybe$algo2list[[1]]

  if(nextstep$what == "all"){

    # Compute all patient


    blocs <- maybe[namesparam] %>%
      rowid_to_column("id") %>%
      gather("key", "value", -id) %>%
      mutate(param = gsub("(min$)|(max$)", "", key)) %>%
      mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
      select(-key) %>%
      spread(a, value) %>%
      rename(from = min, to = max) %>%
      left_join(domain %>% distinct(param, digits), by = "param")



    blocs %>%
      mutate(sampl = pmap(list(from, to, digits), function(from, to, digits){

        # from = 0; to = 3; digits = 4
        seq(from, to, 10^(-digits) )  %>% unique()


      })) -> VPsparam

    # VPsparam

    # Cross parameter (per bloc) and add fixed values
    newVPs <- VPsparam %>%
      group_split(id) %>%
      map( function(x){
        temp <- invoke(.fn = crossing, x$sampl )
        names(temp) <- x$param
        temp
      }) %>%
      bind_rows() %>%
      add_column(DFfix)


    newVPs <- newVPs %>% distinct()


    newVPs <-   newVPs %>%
      rowid_to_column("group") %>%
      mutate(group = floor(group / 400000) + 1)




    tempVP <- self$clone(deep = T)

    tempVP$filters_neg_above <- tibble()
    tempVP$filters_neg_below <- tibble()
    tempVP$algo2list  <- list()

for(a in newVPs$group %>% unique()){

  print(paste0(a, "/", length(newVPs$group %>% unique())))

  newVPs2 <- newVPs %>%
    filter(group == a)

    tempVP$add_VP(newVPs2)
    #
}

    print(tempVP)


  # tempVP$n_filter_reduc()

  newname <- paste0(nextstep$Name,"_final")

  self$algo2list[[newname]] <- tempVP

  self$algo2list$tree  <- self$algo2list$tree %>%
    add_row(Name = newname, size =  0, todo = 0, before = nrow(newVPs), after = 0, ratio = 0, what = "Done",
            time =   difftime(Sys.time(), t0, units = "s")
    )

  self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] <- 0


  saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")

  nextstep <-  self$algo2list[["tree"]] %>%
    filter(todo != 0 ) %>%
    filter(what == "all") %>%
    filter(!is.na(todo)) %>%
    slice(1)


  # nextstep <-  self$algo2list[["tree"]] %>%
  #   filter(todo != 0 ) %>%
  #   # filter(Name == "first") %>%
  #   filter(!is.na(todo)) %>%
  #   slice(1)
  #
  # next

  }else{


  maybe <-  maybe %>%
    # self$algo2list[[nextstep$Name]]
    filter(blocsPool == todo) %>%
    select(-blocsPool)

  nVP0s <- ndomain2(maybe)



 blocs <- maybe[namesparam] %>%
  rowid_to_column("id") %>%
  gather("key", "value", -id) %>%
  mutate(param = gsub("(min$)|(max$)", "", key)) %>%
  mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
  select(-key) %>%
  spread(a, value) %>%
  rename(from = min, to = max) %>%
  left_join(domain %>% distinct(param, digits), by = "param")

 nperparam <- ceiling((200000/length(unique(blocs$id)))^(1/nrow(domain)))


blocs %>%
  mutate(sampl = pmap(list(from, to, digits), function(from, to, digits){

# from = 0; to = 3; digits = 4
    seq(from, to, (to-from)/(nperparam-1)) %>% round(digits) %>% unique()


  })) -> VPsparam

# Cross parameter (per bloc) and add fixed values
newVPs <- VPsparam %>%
  group_split(id) %>%
  map( function(x){
   temp <- invoke(.fn = crossing, x$sampl )
   names(temp) <- x$param
   temp
    }) %>%
  bind_rows() %>%
  add_column(DFfix)


tempVP <- self$clone(deep = T)

tempVP$filters_neg_above <- tibble()
tempVP$filters_neg_below <- tibble()
tempVP$algo2list  <- list()
#
# if(nrow(self$poolVP)>0 ){
#
# newVPs <- newVPs %>%
#   left_join(self$poolVP %>% select(!!!parse_exprs(names(newVPs)),  id)) %>%
#   filter(is.na(id)) %>%
#   select(-id)
# }
#
# if(! firstbloc) newVPs <- newVPs %>%
#   left_join(prevVPs %>% mutate(test = 1)) %>%
#   filter(is.na(test)) %>%
#   select(-test)
#
# prevVPs <-newVPs



# newVPs <- invoke(.fn = crossing, .args = VPsparam)%>%
#   add_column(DFfix)

# Compute the new VPs
tempVP$add_VP(newVPs, keepRedFiltaftDis = T, reducefilteratend = T)




# self$n_filter_reduc()
# Compute zone_maybe
sizetemp <- NA
try({
tempVP$compute_zone_maybe()

# Which become the new blocs
# tempVP$zone_maybe
 zone_maybe <- tempVP$zone_maybe
 zone_maybe <- reduce_maybe2(zone_maybe)

tempVP$algo2list[["blocs"]] <- zone_maybe %>%
rowid_to_column("blocsPool") %>%
  mutate(blocsPool = floor(blocsPool/gatherblocsPool) + 1)

sizetemp <- max(tempVP$algo2list[["blocs"]]$blocsPool)

# ndomain2(tempVP$zone_maybe)

})


  # maybe

# tempVP$zone_maybe
nVPs <- ndomain2(zone_maybe)

tloop <- difftime(Sys.time(), t0, units = "s")

# Deactivate zoom system if needed
timesaved <- (nVP0s - nVPs) * as.double(taddvp) / 200000

nextsteptype <- "zoom"
if(tloop < 2 * timesaved ) nextsteptype <- "all"

# nVPs * as.double(taddvp) / 200000
#   as.double(tloop)



newname <- paste0(c(nextstep$Name, nextstep$todo), collapse = "_")
self$algo2list[[newname]] <- tempVP

self$algo2list$tree  <- self$algo2list$tree %>%
  add_row(Name = newname, size =  sizetemp, todo = 1, before = nVP0s, after = nVPs, ratio = nVPs/nVP0s,
          time = difftime(Sys.time(), t0, units = "s"), what = nextsteptype
          )

self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] <- self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] + 1
if(self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] > self$algo2list$tree$size[self$algo2list$tree$Name == nextstep$Name]){

  self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name]  <- 0
}

nextstep <-  self$algo2list[["tree"]] %>%
  filter(todo != 0 ) %>%
  filter(what == "all") %>%
  filter(!is.na(todo)) %>%
  slice(1)

saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")
#
# nextstep <- self$algo2list[["tree"]] %>%
#   filter(Name == newname)
# self$plot_2D(k2, lambda0)
# print("a")
  }

}

})

# self$algo2list$tree  <- self$algo2list$tree %>%
#   mutate(todo = if_else(todo >1,1, todo)) %>%
#   mutate(what = if_else(Name != "first", "all", "zoom"))
#   filter(what == "zoom" & todo != 0)
#   mutate(what = if_else())
# # saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")


# Reduce maybe ------------------------------------------------------------

reduce_maybe <- function(maybe){


  maybe2 <- maybe %>%
    # filter(blocsPool %in% pool) %>%
    rowid_to_column("group") %>%
    mutate(group = floor(group/ 2000)+1) %>%
    rowid_to_column("bloc")
t0 <- Sys.time()
  for(a in unique(maybe2$group)){

print(a)
    ids <-

      maybe2 %>%
      filter(group == a) %>%
      # filter(blocsPool %in% pool) %>%
      rowid_to_column() %>%
      group_split(rowid) %>%
      map(~  tibble(k2 = c(.x$k2min,.x$k2max), lambda0 = c(.x$lambda0max,.x$lambda0min),
                    ke = c(.x$kemax,.x$kemin), Vd = c(.x$Vdmax,.x$Vdmin), lambda1 = c(.x$lambda1max,.x$lambda1min)) %>%
            mutate(psi = 20, k1 = 0.5, bloc = .x$bloc)
      ) %>%
      bind_rows() %>%
      rowid_to_column("id")

    proto <- self$protocols$dose50 %>%
      mutate(evid = 1) %>%
      bind_rows(

        self$protocols$dose50 %>%
          mutate(evid = 0) %>%
          select(-time) %>%
          crossing(time = self$times) %>%
          mutate(amt = 0)

      ) %>%
      crossing(id = 1:nrow(ids)) %>%
      arrange(id, time)

  res <-  self$model$solve(ids, proto, c(X1 = 50)) %>%
      as.tibble() %>%
      left_join(ids %>% distinct(id, bloc))

  res %>%
    filter(time %in% self$targets$time) %>%
    gather("cmt", "value", !!!parse_exprs( self$targets$cmt)) %>%
    left_join(self$targets) %>%
    mutate(belowmin = value  < min, abovemax = value > max) %>%
    # filter(!testmin |  !testmax) %>%
    group_by(bloc, time) %>%
    summarise(min = sum(belowmin), max = sum(abovemax)) %>%
    filter(min ==2 | max == 2) %>%
    pull(bloc) %>%
    unique() -> bloctorem


  maybe2 <- maybe2 %>%
    filter(! bloc %in% bloctorem)
  # res %>%
  #   filter(bloc %in% bloctorem[1:20]) %>%
  #   # filter(rowid <20) %>%
  #   ggplot()+
  #   geom_line(aes(time, tumVol, group = id))+
  #   scale_y_log10()+
  #   geom_point(data = self$targets, aes(time, min), col = "red")+
  #   # geom_line(data=ref %>% select(-rowid), aes(time, tumVol, group = id), col = "red")+
  #   facet_wrap(~bloc)

  }

difftime(Sys.time(), t0, "s")
return(maybe2)
}




reduce_maybe2 <- function(maybe, obj = self){


  testabove <- obj$clone(deep = T)

  testbelow <- obj$clone(deep = T)
  #Which are above

  # tempreduce$targets$min <-  tempreduce$targets$max
  testabove$targets$max <- Inf
  testbelow$targets$min <-- Inf

  maybeabove <- maybe
  maybebelow <- maybe
  for(a in obj$param){

if(a %in% obj$param_increase[[1]]){

  names(maybeabove)[names(maybeabove) == paste0(a, "max")] <- a
  names(maybebelow)[names(maybebelow) == paste0(a, "min")] <- a
}else if(a %in% obj$param_reduce[[1]] ){

  names(maybeabove)[names(maybeabove) == paste0(a, "min")] <- a
  names(maybebelow)[names(maybebelow) == paste0(a, "max")] <- a
  }else{

    #to handle the rest...
    names(maybeabove)[names(maybeabove) == paste0(a, "min")] <- a
    names(maybebelow)[names(maybebelow) == paste0(a, "max")] <- a
  }
}

  testabove$add_VP(maybeabove %>% select(!!!parse_exprs(obj$param)), use_green_filter =   T, npersalve = 2000)
  testbelow$add_VP(maybebelow %>% select(!!!parse_exprs(obj$param)), use_green_filter = T, npersalve = 2000)


  testabove$poolVP %>% select(id) %>%
    left_join(   testbelow$poolVP %>% select(id, rowid)   ) %>%
    filter(!is.na(rowid)) %>%
    pull(id) -> idtokeep


  maybe2 <- maybe %>%
    slice(idtokeep)


  # difftime(Sys.time(), t0, "s")
  return(maybe2)
}

# Massive screening -------------------------------------------------------



# 1000^6
# 10^18



# Test_reduce bloc --------------------------------------------------------


# ndomain2 <- function(blocs){
#
#
#
#
#   blocs <- tempVP$zone_maybe %>%
#     rowid_to_column("id") %>%
#     gather("key", "value", -id) %>%
#     mutate(param = gsub("(min$)|(max$)", "", key)) %>%
#     mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
#     select(-key) %>%
#     spread(a, value) %>%
#     rename(from = min, to = max) %>%
#     left_join(domain %>% distinct(param, digits), by = "param")
#
#   nperparam <- ceiling((200000/length(unique(blocs$id)))^(1/nrow(domain)))
#   nperparam <- 100
#
#   blocs %>%
#     mutate(sampl = pmap(list(from, to, digits), function(from, to, digits){
#
#       # from = 0; to = 3; digits = 4
#       seq(from, to, (to-from)/(nperparam-1)) %>% round(digits) %>% unique()
#
#
#     })) -> VPsparam
#
#   # Cross parameter (per bloc) and add fixed values
#   newVPs <- VPsparam %>%
#     group_split(id) %>%
#     map( function(x){
#       temp <- invoke(.fn = crossing, x$sampl )
#       names(temp) <- x$param
#       temp
#     }) %>%
#     bind_rows()
#
#
#   newVPs %>%
#     distinct()
#
#
#   temp <-  blocs %>%
#     rowid_to_column("id") %>%
#     gather("key", "value", -id , -contains("blocsPool")) %>%
#     mutate(param = gsub("(min$)|(max$)", "", key)) %>%
#     mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
#     select(-key) %>%
#     spread(a, value) %>%
#     rename(from = min, to = max) %>%
#     left_join(domain %>% distinct(param, digits), by = "param")
#
#   blocs
#   temp%>%
#     filter(!is.na(digits))
#
#   temp %>%
#     filter(!is.na(digits)) %>%
#     mutate(how = pmap_dbl(list(from, to, digits), function(from, to , digits){
#
#       length(seq(from, to, 10^(-digits)))
#     })) %>%
#     group_split(id) -> temp
#
#
#   temp2 <-  temp %>%
#     map_dbl(~  .x %>% pull(how) %>%
#               reduce(`*`))
#
#   sum(temp2)
# }


reduce_bloc <- function(){


blocsid <- blocs %>%
  rowid_to_column("id")


param <- self$param

line <-  blocsid %>%
    slice(1)

temp <- "lambda1"

maxt  <- line[[paste0(temp, "max")]]


temp2 <- blocsid %>%
  left_join(line %>% select(-lambda1min, - lambda1max) %>% mutate(test=1)) %>%
  filter(test ==1)

temp2 %>%
  filter(lambda1min == maxt) -> to_merge

# linenew <- line
blocsid$lambda1max[blocsid$id == line$id]<- to_merge$lambda1max

blocsid <- blocsid %>% filter(id != to_merge$id)

}

