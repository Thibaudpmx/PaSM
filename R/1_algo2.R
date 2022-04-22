# Algorithm 2 -------------------------------------------------------------
# library(QSPVP)
# seq(0,3,0.0015) %>% length() *
# + seq(0,1.4,0.000820) %>% length()
# source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")




reduce_maybe2 <- function(maybe, obj = self){


  testabove <- obj$clone(deep = T)
  testabove$poolVP <- tibble()

  testbelow <- obj$clone(deep = T)
  testbelow$poolVP <- tibble()
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

  testabove$add_VP(maybeabove %>% select(!!!parse_exprs(obj$param)) , use_green_filter =   T, npersalve = 2000)

  maybeabove %>% select(!!!parse_exprs(obj$param)) %>%
    rowid_to_column("tokeep") %>%
    left_join(testabove$poolVP %>% distinct(!!!parse_exprs(obj$param)) %>% mutate(test = T)) %>%
    filter(test) %>%
    pull(tokeep) -> idpostabove


  maybebelow2 <- maybebelow %>%  rowid_to_column("tokeep") %>% slice(idpostabove) %>% select(!!!parse_exprs(obj$param), tokeep)

  testbelow$add_VP(maybebelow %>% slice(idpostabove) %>% select(!!!parse_exprs(obj$param)), use_green_filter = T, npersalve = 2000)



  testbelow$poolVP  %>%
    select(!!!parse_exprs(obj$param)) %>%
    distinct() %>%
    left_join(maybebelow2) %>%  pull(tokeep) -> idtokeep

  # testabove$poolVP %>%
  #   select(-rowid) %>%
  #   left_join(   testbelow$poolVP %>% select(id, rowid)   ) %>%
  #   filter(!is.na(rowid)) %>%
  #   pull(id) %>%
  #   unique()




  maybe2 <- maybe %>%
    slice(idtokeep)
    # rowid_to_column("id") %>%
    # filter(id %in% idtokeep)


  # difftime(Sys.time(), t0, "s")
  return(maybe2)
}
#

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

fix <-c(k1 = 0.5, w0 = 50)

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
add_nvp_bloc <- function(blocs){


  temp <-  blocs %>%
    # slice(1:2) %>%
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

      max(length(seq(from, to, 10^(-digits))) - 2,1) # ici -2 to avoid border?


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

  blocs %>%
    mutate(temp3)
  # sum(temp3)
}


# (ndomain(domain)- sum(temp2)) / 2000 * 0.5 / 3600 / 24

# sum(temp2)/ 2000 * 0.5 / 3600 / 24

prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100, 431),
                    max = c(100.05, 431.05))

prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(50, 200),
                    max = c(60, 210))

prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(50, 200),
                    max = c(70, 230)) %>%
  bind_rows(tibble(
    protocol = "dose100", cmt = "tumVol", time = c(12,40), min = c(10, 125),
    max = c(55, 150))
  ) %>%
  bind_rows(tibble(
    protocol = "dose0", cmt = "tumVol", time = c(12,40), min = c(200, 600),
    max = c(250, 6000))
  )


self <- VP_proj_creator$new()

self$set_targets(manual = prototiny)


 npersalve = 2E5
 npersalveFinal = 1E6
 fix <-c(k1 = 0.5, w0 = 50)

file <- "D:/these/Second_project/QSP/modeling_work/VT_simeoni/testtwodose.RDS"
self <- readRDS(file)
file <- ""
save_every = 2
# Main functin ------------------------------------------------------------


VP_proj_creator$set("public", "algo2", function(domain, fix = NULL, npersalve = 2E5, npersalveFinal = 1E6, file = "", save_every = 5){


  nVPs0  <- ndomain(domain) - 2^nrow(domain) # Compute the number of VPs

  ## First division
  t0 <- Sys.time()




  namesparam <- c(paste0(domain$param,  "max"),paste0(domain$param,  "min") ) # names of parameter with min and max

  # The parameter that will not vary
  DFfix <- as.data.frame(fix) %>%
    rownames_to_column() %>%
    spread(rowname , fix)

  blocs <- domain %>% # See after?
    mutate(id = "")

  firstbloc <- T # Used in the loop


  # Compute the first the first vectors of parameters

  nperparam <- floor(npersalve^(1/nrow(domain))) # how many division per parameter per iteration

  blocs %>%
    mutate(sampl = pmap(list(from, to, digits), function(from, to, digits){

      # from = 0; to = 3; digits = 4
      # seq(from, to, (to-from)/(nperparam-1)) %>% round(digits) %>% unique()
      temp <- seq(from, to, (to-from)/(nperparam+1)) %>% round(digits) %>% unique()

      temp[-c(1, length(temp))]

    })) -> VPsparam

  # Cross parameter (per bloc) and add fixed values: perfom the first cohort
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

  self$add_VP(newVPs, keepRedFiltaftDis = T, reducefilteratend = F, use_green_filter = T, pctActivGreen = 0.1)
  taddvp <- difftime(Sys.time() ,  taddvp, units = "s")

  self$n_filter_reduc()
  # Compute zone_maybe
  self$compute_zone_maybe()

  # self$compute_zone_sure()
  # Which become the new blocs
  maybe <- self$zone_maybe


  for(a in names(maybe)[grepl("max", names(maybe))]){

  if(length(   maybe[[a]][ maybe[[a]] ==Inf]) >0)  maybe[[a]][ maybe[[a]] ==Inf]  <- domain$to[domain$param == gsub("max", "", a)]
  }

  for(a in names(maybe)[grepl("min", names(maybe))]){


    if(length(   maybe[[a]][ maybe[[a]] ==0]) >0)  maybe[[a]][ maybe[[a]] ==0]  <- domain$from[domain$param == gsub("min", "", a)]
  }

  self$zone_maybe   <- maybe
#
#   maybe %>%
#     filter(k2min >1.6 & k2max < 2.5 & lambda0min >0.06 & lambda0max < 1.2 & kemin >1 & kemax < 1.5&
#              Vdmin >33 & Vdmax < 38)
#
#
#   maybe %>%
#     filter(k2min <= 2 & k2max >= 2 & lambda0min <= 0.09 & lambda0max >= 0.09 & kemin <= 1.3 & kemax >= 1.3 &
#              Vdmin <= 35 & Vdmax >= 35 & lambda1min <= 16 & lambda1max >= 16)
#
#
#   self$poolVP %>% filter(k2 == 2)

  cat("Reduce maybe")
  maybe <- reduce_maybe2(maybe)
  cat("End reduce maybe")

  maybe <- add_nvp_bloc(maybe)

  self$zone_maybe   <- maybe


to_do_final <-   maybe %>%
    filter(temp3 < npersalveFinal)

gatherblocsPool <- (npersalveFinal / median(maybe$temp3)) %>% floor()

to_do_final <- to_do_final %>%
    rowid_to_column("blocsPool") %>%
    mutate(blocsPool = floor(blocsPool/gatherblocsPool) + 1) %>%
    mutate(todo = "final")

to_dive <-   maybe %>%
  filter(temp3 >= npersalveFinal) %>%
  rowid_to_column("blocsPool") %>%
  mutate(blocsPool = max(to_do_final$blocsPool)+blocsPool) %>%
  mutate(todo = "dive")


maybeFinal <- bind_rows(to_do_final, to_dive)

  nVPs <- sum(maybeFinal$temp3)

  timeperun <- 0.5/2000

  self$algo2list[["tree"]] <- tibble(Name = "first", size = length(unique(maybeFinal$blocsPool)), todo = 1,
                                     before = nVPs0, after = nVPs, ratio = nVPs/nVPs0,
                                     time = difftime(Sys.time(), t0, units = "s"), what = "zoom")

  self$algo2list[["domain"]] <- domain

  self$algo2list[["first"]] <- maybeFinal


  # maybeFinal[namesparam] %>% distinct()

  if(file != "") saveRDS(self, file)

  tree <- self$algo2list[["tree"]]
  # nextstep <-  self$algo2list[["tree"]] %>%
  #   filter(!is.na(todo)) %>%
  #   slice(1)
  #
  #
  #
  # nVPs/nVPs0 * 100
  # nVPs/2000 * 0.5 / 3600 / 24
  #
  # nVPs0/2000 * 0.5 / 3600 / 24
  # Now the deep dive
  nsave <- 0
  while(sum( self$algo2list[["tree"]]$todo) > 0 ){


    t0 <- Sys.time()

    nextstep <- self$algo2list[["tree"]] %>% filter(todo > 0) %>% slice(1)

    print(nextstep)

    maybe <-  self$algo2list[[nextstep$Name]] %>%
      filter(blocsPool == nextstep$todo)

    # if(nextstep$Name !="first") maybe <- maybe$algo2list[[1]]

    nVP0s <- sum(maybe$temp3)

    if(unique(maybe$todo) == "final"){

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
          temp <- seq(from, to, 10^(-digits) )  %>% unique()
          temp[2:max(length(temp)-1,2)]

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


      # newVPs <- newVPs %>% distinct()

#
#       newVPs <-   newVPs %>%
#         rowid_to_column("group") %>%
#         mutate(group = floor(group / 400000) + 1)




      tempVP <- self$clone(deep = T)

      tempVP$filters_neg_above <- tibble()
      tempVP$filters_neg_below <- tibble()
      tempVP$algo2list  <- list()
      tempVP$poolVP  <- tibble()

      # for(a in newVPs$group %>% unique()){

        # print(paste0(a, "/", length(newVPs$group %>% unique())))

        # newVPs2 <- newVPs %>%
          # filter(group == a)

        tempVP$add_VP(newVPs, keepRedFiltaftDis = T, reducefilteratend = F, use_green_filter = T, pctActivGreen = 0.1,npersalve = 2000)
        #
      # }

      # print(tempVP)


      # tempVP$n_filter_reduc()

      newname <- paste0(nextstep$Name,"_", nextstep$todo)


      tempVPstorage <- list("neg_above" = tempVP$filters_neg_above, "neg_below" = tempVP$filters_neg_below, "poolVP" = tempVP$poolVP)

      self$algo2list[[newname]] <- tempVPstorage

      self$algo2list$tree  <- self$algo2list$tree %>%
        add_row(Name = newname, size =  0, todo = 0, before = nrow(newVPs), after = nrow( tempVP$poolVP), ratio = 0, what = "Done",
                time =   difftime(Sys.time(), t0, units = "s")
        )

      self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] <- if_else(nextstep$todo == nextstep$size, 0, nextstep$todo + 1)



      # nextstep <-  self$algo2list[["tree"]] %>%
      #   filter(todo != 0 ) %>%
      #   filter(what == "all") %>%
      #   filter(!is.na(todo)) %>%
      #   slice(1)


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

      nVP0s <- sum(maybe$temp3)



      blocs <- maybe[namesparam] %>%
        rowid_to_column("id") %>%
        gather("key", "value", -id) %>%
        mutate(param = gsub("(min$)|(max$)", "", key)) %>%
        mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
        select(-key) %>%
        spread(a, value) %>%
        rename(from = min, to = max) %>%
        left_join(domain %>% distinct(param, digits), by = "param")

      nperparam <- ceiling((npersalve/length(unique(blocs$id)))^(1/nrow(domain)))


      blocs %>%
        mutate(sampl = pmap(list(from, to, digits), function(from, to, digits){

          # from = 0; to = 3; digits = 4
          # seq(from, to, (to-from)/(nperparam-1)) %>% round(digits) %>% unique()
          withbord <- seq(from, to, 10^(-digits)) %>% round(digits) %>% unique()
          withbord[2:max(length(withbord)-1,2)]

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
      tempVP$add_VP(newVPs, keepRedFiltaftDis = T, reducefilteratend = F)




      # self$n_filter_reduc()
      # Compute zone_maybe

   test_zone_maybe <-    try({   aezatempVP$compute_zone_maybe()})

        # Which become the new blocs
        # tempVP$zone_maybe
   if(class(test_zone_maybe) != "try-error"){
        zone_maybe <- tempVP$zone_maybe
        zone_maybe <- reduce_maybe2(zone_maybe)

        tempVP$algo2list[["blocs"]] <- zone_maybe %>%
          rowid_to_column("blocsPool") %>%
          mutate(blocsPool = floor(blocsPool/gatherblocsPool) + 1)

        sizetemp <- max(tempVP$algo2list[["blocs"]]$blocsPool)

        # ndomain2(tempVP$zone_maybe)
   }else{

     zone_maybe <- tibble()

     tempVP$algo2list[["blocs"]] <- "too big"
     sizetemp <- NA
     nVPs <- NA
   }



      # maybe

      # tempVP$zone_maybe
      if(nrow(zone_maybe)> 0 ){
      zone_maybe <- add_nvp_bloc(zone_maybe)
      nVPs <- sum(zone_maybe$temp3)
      }else if(!is.na(sizetemp)){

        nVPs <- 0
        sizetemp <- 0
      }
      tloop <- difftime(Sys.time(), t0, units = "s")

      # Deactivate zoom system if needed
      timesaved <- (nVP0s - nVPs) * as.double(taddvp) / 200000

      if(is.na(timesaved)){

        nextsteptype <- "not_computable"
      }else if(tloop < 2 * timesaved ) {

        nextsteptype <- "all"
      }else{
      nextsteptype <- "zoom"
      }

      # nVPs * as.double(taddvp) / 200000
      #   as.double(tloop)



      newname <- paste0(c(nextstep$Name, nextstep$todo), collapse = "_")
      self$algo2list[[newname]] <- tempVP

      self$algo2list$tree  <- self$algo2list$tree %>%
        add_row(Name = newname, size =  sizetemp, todo = min(nVPs,1), before = nVP0s, after = nVPs, ratio = nVPs/nVP0s,
                time = difftime(Sys.time(), t0, units = "s"), what = nextsteptype
        )

      self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] <- self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] + 1
      if(self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] > self$algo2list$tree$size[self$algo2list$tree$Name == nextstep$Name]){

        self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name]  <- 0
      }

      nextstep <-  self$algo2list[["tree"]] %>%
        filter(todo != 0 ) %>%
        # filter(what == "all") %>%
        filter(!is.na(todo)) %>%
        slice(1)

      # saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")
      #
      # nextstep <- self$algo2list[["tree"]] %>%
      #   filter(Name == newname)
      # self$plot_2D(k2, lambda0)
      # print("a")
    } # end else


    nsave <- nsave + 1
    if(file != "" &   nsave == save_every){

      cat(silver("Saving..."))
      saveRDS(self, file)
      nsave <- 0
    }

    # saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")


  } # end while loop



})

# self$algo2list$tree  <- self$algo2list$tree %>%
#   mutate(todo = if_else(todo >1,1, todo)) %>%
#   mutate(what = if_else(Name != "first", "all", "zoom"))
#   filter(what == "zoom" & todo != 0)
#   mutate(what = if_else())
# # saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")



# Tries -------------------------------------------------------------------

# try to compute new patients

  nnew <- (600000 / nrow(blocs))%>% ceiling

  temp <-  blocs %>%
    # slice(1:2) %>%
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
    # slice(1) %>%
    mutate(how = pmap(list(from, to, digits), function(from, to , digits){


      tempp <- seq(from, to, (to-from)/(3+1)) %>% round(digits) %>% unique()

      # print(tempp)
      tibble(n = list(tempp[-c(1, length(tempp))]), ntotal =     max(length(seq(from, to, 10^(-digits))) - 2,1))






    })) -> temp2

  temp2 %>%
  # temp2 %>%
    unnest() %>%
    filter(id == 1) -> x
    group_split(id) %>%
    map(function(x){

      print(unique(unique(x$id)))
      x %>%
        mutate(exp = paste0(n)) %>%
        pull(exp) -> exp

      names(exp ) <- x$param

      crossing(!!!parse_exprs(exp)) %>%
        mutate(id = unique(x$id))
    }) -> temp3




    temp3 %>%
    bind_rows() %>%
    group_by(id) %>%
    tally %>%
    group_by(n) %>%
    tally

    newcohort <-
      temp3 %>%
      bind_rows() %>%
      select(-id) %>%
      mutate(k1 = 0.5, w0 = 50)

    self2 <- VP_proj_creator$new()

    self2$set_targets(manual = prototiny)

    self2$add_VP(newcohort, keepRedFiltaftDis = T, reducefilteratend = F)

    self2$n_filter_reduc()

    self2

    # Which become the new blocs
    self2$compute_zone_maybe()
    maybe <- self$zone_maybe


    self3 <- self2$clone(deep = T)

    self2$filters_neg_above %>%
      distinct(lambda0) %>%
      pull(lambda0)

    self3$filters_neg_above  <- self3$filters_neg_above  %>% filter(lambda0 < 0.52)
    self3$filters_neg_below  <- self3$filters_neg_below  %>% filter(lambda0 < 0.52)
    self3$compute_zone_maybe()

    self4$filters_neg_above  <- self3$filters_neg_above  %>% filter(lambda0 >=  0.52)
    self4$filters_neg_below  <- self3$filters_neg_below  %>% filter(lambda0 >=  0.52)
    self4$compute_zone_maybe()


    self2$compute_zone_maybe()

    message("Reduce maybe")
    maybe <- reduce_maybe2(maybe)
    message("End reduce maybe")


    gatherblocsPool <- 100

    blocs <- maybe %>%
      rowid_to_column("blocsPool") %>%
      mutate(blocsPool = floor(blocsPool/gatherblocsPool) + 1)


    blocs <- add_nvp_bloc(blocs)