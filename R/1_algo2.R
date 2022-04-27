


# Helping functions -------------------------------------------------------
#' Title
#'
#' @param ind_param
#' @param add_events
#'
#' @return
#' @export
#'
#' @examples
#'
cleanObj <- function(obj){

  obj$filters_neg_above <- tibble()
  obj$filters_neg_below <- tibble()
  obj$poolVP <- tibble()
obj
}

#' Title
#'
#' @param ind_param
#' @param add_events
#'
#' @return
#' @export
#'
#' @examples
#'
reduce_maybe2 <- function(maybe, obj = self){





  for(OoI in obj$targets$cmt %>% unique){



    testabove <- obj$clone(deep = T) %>% cleanObj()

    testbelow <- obj$clone(deep = T) %>% cleanObj()

    #Which are above

    # tempreduce$targets$min <-  tempreduce$targets$max
    testabove$targets$max <- Inf
    testbelow$targets$min <-- Inf

    # kept only the target of desired OoI
    testabove$targets <- testabove$targets %>% filter(cmt == OoI)
    testbelow$targets <- testbelow$targets %>% filter(cmt == OoI)

    maybeabove <- maybe
    maybebelow <- maybe

  if(nrow(maybe) > 0){

  for(a in obj$param){

    if(a %in% obj$param_increase[[OoI]]){

      names(maybeabove)[names(maybeabove) == paste0(a, "max")] <- a
      names(maybebelow)[names(maybebelow) == paste0(a, "min")] <- a
    }else if(a %in% obj$param_reduce[[OoI]] ){

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

  if(nrow(maybebelow %>% slice(idpostabove)) > 0){
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




  maybe <- maybe %>%
    slice(idtokeep)
  # rowid_to_column("id") %>%
  # filter(id %in% idtokeep)

  }else{

    maybe <- maybe %>% slice(0)
  }

  }

  }# end if nrow(maybe > 0)
  # difftime(Sys.time(), t0, "s")
  return(maybe)
}
#



#' Title
#'
#' @param ind_param
#' @param add_events
#'
#' @return
#' @export
#'
#' @examples
#'
ndomain <- function(domain){


  domain %>%
    mutate(how = pmap_dbl(list(from, to, by), function(from, to , by){

      length(seq(from, to, by))
    })) %>%
    pull(how) %>%
    reduce(`*`)


}




#' Title
#'
#' @param ind_param
#' @param add_events
#'
#' @return
#' @export
#'
#' @examples
#'
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
    left_join(domain %>% distinct(param, by), by = "param")

  temp %>%
    filter(!is.na(by)) %>%
    mutate(how = pmap_dbl(list(from, to, by), function(from, to , by){

      max(length(seq(from, to, by)) - 2,1) # ici -2 to avoid border?


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

# Algorithm 2 -------------------------------------------------------------
# library(QSPVP)
# seq(0,3,0.0015) %>% length() *
# + seq(0,1.4,0.000820) %>% length()
# source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")



#
# domain <- tribble(~param, ~from, ~to, ~digits,
#                   "k2", 0, 3, 3 ,
#                   "lambda0", 0, 1.4, 4
#                   )
# fix <-c(k1 = 0.5, ke = 1, Vd = 40, lambda1 = c(12))


# domain <- tribble(~param, ~from, ~to, ~digits,
#                   "k2", 0, 3, 2 ,
#                   "lambda0", 0, 1.4, 2,
#                   "ke", 0, 2,1,
#                   "Vd", 0,40,0,
#                   "lambda1", 0,24,1
# )
#
# fix <-c(k1 = 0.5, w0 = 50)
#
# # ndomain(domain) / 2000 * 0.5 / 3600 / 24
#
# # blocs <- zone_maybe
#
#
#
# # (ndomain(domain)- sum(temp2)) / 2000 * 0.5 / 3600 / 24
#
# # sum(temp2)/ 2000 * 0.5 / 3600 / 24
#
# prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(100, 431),
#                     max = c(100.05, 431.05))
#
# prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(50, 200),
#                     max = c(60, 210))
#
# prototiny <- tibble(protocol = "dose50", cmt = "tumVol", time = c(12,40), min = c(50, 200),
#                     max = c(70, 230)) %>%
#   bind_rows(tibble(
#     protocol = "dose100", cmt = "tumVol", time = c(12,40), min = c(10, 125),
#     max = c(55, 150))
#   ) %>%
#   bind_rows(tibble(
#     protocol = "dose0", cmt = "tumVol", time = c(12,40), min = c(200, 600),
#     max = c(250, 6000))
#   )
#
#
# self <- VP_proj_creator$new()
#
# self$set_targets(manual = prototiny)
#
#
#  npersalve = 2E5
#  npersalveFinal = 1E6
#  fix <-c(k1 = 0.5, w0 = 50)
#
# # file <- "D:/these/Second_project/QSP/modeling_work/VT_simeoni/testtwodose.RDS"
# self <- readRDS(file)
# file <- ""
# save_every = 2
# Main functin ------------------------------------------------------------

#' Title
#'
#' @param ind_param
#' @param add_events
#'
#' @return
#' @export
#'
#' @examples
#'

VP_proj_creator$set("public", "algo2", function(domain, fix = NULL, npersalve = 2E5, npersalveFinal = 1E6, file = "", save_every = 5, method = 1){


  nVP0s  <- ndomain(domain) - 2^nrow(domain) # Compute the number of VPs


  timeTrak <- list()
  ## First division
  t0 <- Sys.time()

  timeTrak$Start <- t0


  DFfix <- as.data.frame(fix) %>%
    rownames_to_column() %>%
    spread(rowname , fix)

  DFfix2 <- DFfix %>%
    gather("param", "value") %>%
    crossing(a = c("min", "max")) %>%
    mutate(param = paste0(param, a)) %>%
    select(-a) %>%
    spread(param, value)

  namesparam <- c(paste0(domain$param,  "max"),paste0(domain$param,  "min") ) # names of parameter with min and max



  blocs <- domain %>% # See after?
    mutate(id = "")

  firstbloc <- T # Used in the loop



  # Compute the first vectors of parameters

  nperparam <- floor(npersalve^(1/nrow(domain))) # how many division per parameter per iteration




algo1 <- function(obj = self, VPsparam ){




  newVPs <-  rlang::invoke(.fn = crossing, VPsparam$sampl )
  names(newVPs) <- VPsparam$param
  newVPs <- newVPs %>%
    add_column(DFfix)

  timeTrakTemp <- list()

  # Compute the VPs using algo 1
  taddvp <- Sys.time()

  obj$add_VP(newVPs, keepRedFiltaftDis = T, reducefilteratend = F, use_green_filter = T, pctActivGreen = 0.1)

  taddvp <- difftime(Sys.time() ,  taddvp, units = "s")

  timeTrakTemp$FirstAlgo1 <- taddvp

  # Filter reduction
  tFilterReduc <- Sys.time()

  obj$n_filter_reduc()

  timeTrakTemp$FilterReduc <- difftime(Sys.time(),  tFilterReduc, units = "s")


  # Compute zone_maybe
  tMaybe <- Sys.time()
  obj$compute_zone_maybe(limits = blocs)
  timeTrakTemp$ZoneMaybe <- difftime(Sys.time(),  tMaybe, units = "s")

  maybe <- obj$zone_maybe


  for(a in names(maybe)[grepl("max", names(maybe))]){

    if(length(   maybe[[a]][ maybe[[a]] ==Inf]) >0)  maybe[[a]][ maybe[[a]] ==Inf]  <- domain$to[domain$param == gsub("max", "", a)]
  }

  for(a in names(maybe)[grepl("min", names(maybe))]){


    if(length(   maybe[[a]][ maybe[[a]] ==0]) >0)  maybe[[a]][ maybe[[a]] ==0]  <- domain$from[domain$param == gsub("min", "", a)]
  }


  obj$zone_maybe   <- maybe

  return(list(obj = obj, timeTrack = timeTrakTemp))
}




  # maybeFinal[namesparam] %>% distinct()

  if(file != "") saveRDS(self, file)


  # Now the deep dive
  nsave <- 0
  while(sum( self$algo2list[["tree"]]$todo) > 0 | is.null( self$algo2list[["tree"]]) ){


    # Determine what to do
    t0 <- Sys.time()

    if(is.null( self$algo2list[["tree"]])){
      maybe <- list()
      maybe$todo <- "First0"
    }else{
    nextstep <- self$algo2list[["tree"]] %>% filter(todo > 0) %>% slice(1)

    print(nextstep)

    maybe <-  self$algo2list[[nextstep$Name]] %>%
      filter(blocsPool == nextstep$todo)

    # if(nextstep$Name !="first") maybe <- maybe$algo2list[[1]]

    nVP0s <- sum(maybe$temp3)
    }



    tempVP <- self$clone(deep = T)

    tempVP$filters_neg_above <- tibble()
    tempVP$filters_neg_below <- tibble()
    tempVP$algo2list  <- list()
    ### Here do differently wether we should cut a domain or completely explore it

# Zoom --------------------------------------------------------------------


    if(unique(maybe$todo) != "final"){ # In case we need to cut a domain




      if( unique(maybe$todo) == "First0"){
        blocs <- domain
      }else{

        blocs <- maybe %>%
          select(-todo) %>%
          gather("param", "value") %>%
          filter(grepl("(min$)|(max$)", param)) %>%
          mutate(minmax = if_else(grepl("min$", param), "min", "max")) %>%
          mutate(param = gsub("(min$)|(max$)", "", param)) %>%
          filter(param %in% domain$param) %>%
          spread(minmax, value) %>%
          rename(from = min, to = max) %>%
          left_join(domain %>% select(param, by), by = "param") %>%
          select(param, from, to, by)

      }


      nperparam <- floor(npersalve^(1/nrow(domain)))


      blocs %>%
        mutate(sampl = pmap(list(from, to, by), function(from, to, by){

          # from = 0; to = 75; digits = 5
          temp <- seq(from, to, (to-from)/(nperparam+1))

          temp <- temp - temp %% by

          temp[-c(1, length(temp))]

        })) -> VPsparam






      map2(VPsparam$param, VPsparam$sampl, function(x, y){

        namemax <- parse_expr(paste0(x,"max"))

        tibble(!!namemax := c(0, y, domain$to[domain$param == x])) %>%
          mutate(!!parse_expr(paste0(x,"min")) := lag(!!namemax)) %>%
          slice(-1)

      }) %>%
        rlang::invoke(.fn = tidyr::crossing) %>%
        mutate(w0min = 50, w0max = 50, k1min = 0.5, k1max = 0.5)-> allBlocs


      if(method == 1){

        usealgo1 <- algo1(tempVP, VPsparam)

        tempVP <- usealgo1$obj

        maybe <- tempVP$zone_maybe

        timeTrak$round1  <- usealgo1$timeTrack
      }else{


        maybe <- map2(VPsparam$param, VPsparam$sampl, function(x, y){

          namemax <- parse_expr(paste0(x,"max"))

          tibble(!!namemax := c( blocs$from[blocs$param == x], y, blocs$to[blocs$param == x])) %>%
            mutate(!!parse_expr(paste0(x,"min")) := lag(!!namemax)) %>%
            slice(-1)

        }) %>%
          rlang::invoke(.fn = tidyr::crossing) %>%
          mutate(DFfix2)

      }

      timeTrak$sizeMaybeBeforeReduc <- nrow(maybe)
      cat("Reduce maybe")
      tMaybeReduce <- Sys.time()
      maybe <- reduce_maybe2(maybe, obj = tempVP)
      timeTrak$tMaybeReduce <- difftime(Sys.time(),  tMaybeReduce, units = "s")
      timeTrak$sizeMaybeAfterReduc <- nrow(maybe)
      cat("End reduce maybe")


      if(nrow(maybe) > 0 ){
      tSizeBloc<- Sys.time()
      maybe <- add_nvp_bloc(maybe)
      timeTrak$tSizeBlock <- difftime(Sys.time(),  tSizeBloc, units = "s")

      tempVP$zone_maybe   <- maybe


      to_do_final <-   maybe %>%
        filter(temp3 < npersalveFinal)

      gatherblocsPool <- (npersalveFinal / median(to_do_final$temp3)) %>% floor()

      to_do_final <- to_do_final %>%
        rowid_to_column("blocsPool") %>%
        mutate(blocsPool = floor(blocsPool/gatherblocsPool) ) %>%
        mutate(todo = "final")

      if(nrow(to_do_final) >0) if(min(to_do_final$blocsPool) == 0) to_do_final$blocsPool <- to_do_final$blocsPool + 1


      to_dive <-   maybe %>%
        filter(temp3 >= npersalveFinal) %>%
        rowid_to_column("blocsPool") %>%
        mutate(blocsPool = max(to_do_final$blocsPool,0)+blocsPool) %>%
        mutate(todo = "dive")


      maybeFinal <- bind_rows(to_do_final, to_dive)

      nVPs <- sum(maybeFinal$temp3)
      ntodo <- 1
      }else{ #if every blocs deleated

        maybeFinal <- maybe

        nVPs <- 0
    ntodo <- 0
      }

if(unique(maybe$todo) == "First0"){

      self <- tempVP

      self$algo2list[["tree"]] <- tibble(Name = "first", size = length(unique(maybeFinal$blocsPool)), todo = 1,
                                         before = nVP0s, after = nVPs, ratio = nVPs/nVP0s,
                                         time = difftime(Sys.time(), t0, units = "s"), what = "zoom")

      self$algo2list[["domain"]] <- domain

      self$algo2list[["first"]] <- maybeFinal

}else{


  newname <- paste0(nextstep$Name,"_", nextstep$todo)


  if(nrow(tempVP$poolVP) > 0  ) self$poolVP <- bind_rows(self$poolVP, tempVP$poolVP %>% mutate(from = newname))

  self$algo2list[[newname]] <- maybeFinal

  self$algo2list$tree  <- self$algo2list$tree %>%
    add_row(Name = newname, size =  max(c(maybeFinal$blocsPool,0)), todo = ntodo, before = nVP0s, after = nVPs, ratio = nVPs/nVP0s, what = "Done",
            time =   difftime(Sys.time(), t0, units = "s")
    )

  self$algo2list$tree$todo[self$algo2list$tree$Name == nextstep$Name] <- if_else(nextstep$todo == nextstep$size, 0, nextstep$todo + 1)



}



# dig ---------------------------------------------------------------------


    }else{




      # Compute all patient


      blocs <- maybe[namesparam] %>%

        rowid_to_column("id") %>%
        gather("key", "value", -id) %>%
        mutate(param = gsub("(min$)|(max$)", "", key)) %>%
        mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
        select(-key) %>%
        spread(a, value) %>%
        rename(from = min, to = max) %>%
        left_join(domain %>% distinct(param, by), by = "param")

      blocs %>%
        mutate(sampl = pmap(list(from, to, by), function(from, to, by){

          # from = 0; to = 75; digits = 5
          temp <- seq(from, to, (to-from)/(nperparam+1))

          temp <- temp - temp %% by

          temp[-c(1, length(temp))] %>% unique()

        })) -> VPsparam


      # VPsparam

      # Cross parameter (per bloc) and add fixed values
      newVPs <- VPsparam %>%
        group_split(id) %>%
        map( function(x){
          temp <- rlang::invoke(.fn = crossing, x$sampl )
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

      if(nrow(tempVP$poolVP) > 0  ) self$poolVP <- bind_rows(self$poolVP, tempVP$poolVP %>% mutate(from = newname))

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



      # saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")
      #
      # nextstep <- self$algo2list[["tree"]] %>%
      #   filter(Name == newname)
      # self$plot_2D(k2, lambda0)
      # print("a")
    } # end else


    nsave <- nsave + 1

    cat(paste0("Nsave = ", nsave))
    if(file != "" &   nsave == save_every){

      cat(silver("Saving..."))
      saveRDS(self, file)
      nsave <- 0
    }

    # saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")


  } # end while loop


  saveRDS(self, file)
})

# self$algo2list$tree  <- self$algo2list$tree %>%
#   mutate(todo = if_else(todo >1,1, todo)) %>%
#   mutate(what = if_else(Name != "first", "all", "zoom"))
#   filter(what == "zoom" & todo != 0)
#   mutate(what = if_else())
# # saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/algo2.RDS")



# Tries -------------------------------------------------------------------
#
# maybe
#
# VPsparam
#
#
# domain
# # newVPs <-
# map2(VPsparam$param, VPsparam$sampl, function(x, y){
#
#   namemax <- parse_expr(paste0(x,"max"))
#
#   tibble(!!namemax := c(0, y, domain$to[domain$param == x])) %>%
#     mutate(!!parse_expr(paste0(x,"min")) := lag(!!namemax)) %>%
#     slice(-1)
#
# }) %>%
#   invoke(.fn = tidyr::crossing) %>%
#   mutate(w0min = 50, w0max = 50, k1min = 0.5, k1max = 0.5)-> allBlocs
#
# allBlocs <- allBlocs
#
# T00 <- Sys.time()
# allBlocs <- reduce_maybe2(allBlocs)
# talt <-  difftime( Sys.time(),T00)
#
#
# allBlocs <- add_nvp_bloc(allBlocs)
#
#
# allBlocs %>%
#   rename(nnew= temp3) %>%
#   left_join(maybeFinal %>% select(-blocsPool)) %>%
#   filter(is.na(temp3))
#
# # try to compute new patients
# #
# #   nnew <- (600000 / nrow(blocs))%>% ceiling
# #
# #   temp <-  blocs %>%
# #     # slice(1:2) %>%
# #     rowid_to_column("id") %>%
# #     gather("key", "value", -id , -contains("blocsPool")) %>%
# #     mutate(param = gsub("(min$)|(max$)", "", key)) %>%
# #     mutate(a = map2_chr(param, key, ~ gsub(.x, "", .y))) %>%
# #     select(-key) %>%
# #     spread(a, value) %>%
# #     rename(from = min, to = max) %>%
# #     left_join(domain %>% distinct(param, digits), by = "param")
# #
# #
# #   temp %>%
# #     filter(!is.na(digits)) %>%
# #     # slice(1) %>%
# #     mutate(how = pmap(list(from, to, digits), function(from, to , digits){
# #
# #
# #       tempp <- seq(from, to, (to-from)/(3+1)) %>% round(digits) %>% unique()
# #
# #       # print(tempp)
# #       tibble(n = list(tempp[-c(1, length(tempp))]), ntotal =     max(length(seq(from, to, 10^(-digits))) - 2,1))
# #
# #
# #
# #
# #
# #
# #     })) -> temp2
# #
# #   temp2 %>%
# #   # temp2 %>%
# #     unnest() %>%
# #     filter(id == 1) -> x
# #     group_split(id) %>%
# #     map(function(x){
# #
# #       print(unique(unique(x$id)))
# #       x %>%
# #         mutate(exp = paste0(n)) %>%
# #         pull(exp) -> exp
# #
# #       names(exp ) <- x$param
# #
# #       crossing(!!!parse_exprs(exp)) %>%
# #         mutate(id = unique(x$id))
# #     }) -> temp3
# #
# #
# #
# #
# #     temp3 %>%
# #     bind_rows() %>%
# #     group_by(id) %>%
# #     tally %>%
# #     group_by(n) %>%
# #     tally
# #
# #     newcohort <-
# #       temp3 %>%
# #       bind_rows() %>%
# #       select(-id) %>%
# #       mutate(k1 = 0.5, w0 = 50)
# #
# #     self2 <- VP_proj_creator$new()
# #
# #     self2$set_targets(manual = prototiny)
# #
# #     self2$add_VP(newcohort, keepRedFiltaftDis = T, reducefilteratend = F)
# #
# #     self2$n_filter_reduc()
# #
# #     self2
# #
# #     # Which become the new blocs
# #     self2$compute_zone_maybe()
# #     maybe <- self$zone_maybe
# #
# #
# #     self3 <- self2$clone(deep = T)
# #
# #     self2$filters_neg_above %>%
# #       distinct(lambda0) %>%
# #       pull(lambda0)
# #
# #     self3$filters_neg_above  <- self3$filters_neg_above  %>% filter(lambda0 < 0.52)
# #     self3$filters_neg_below  <- self3$filters_neg_below  %>% filter(lambda0 < 0.52)
# #     self3$compute_zone_maybe()
# #
# #     self4$filters_neg_above  <- self3$filters_neg_above  %>% filter(lambda0 >=  0.52)
# #     self4$filters_neg_below  <- self3$filters_neg_below  %>% filter(lambda0 >=  0.52)
# #     self4$compute_zone_maybe()
# #
# #
# #     self2$compute_zone_maybe()
# #
# #     message("Reduce maybe")
# #     maybe <- reduce_maybe2(maybe)
# #     message("End reduce maybe")
# #
# #
# #     gatherblocsPool <- 100
# #
# #     blocs <- maybe %>%
# #       rowid_to_column("blocsPool") %>%
# #       mutate(blocsPool = floor(blocsPool/gatherblocsPool) + 1)
# #
# #
# #     blocs <- add_nvp_bloc(blocs)
