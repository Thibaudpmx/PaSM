# library(RxODE)
# library(R6)


# https://adv-r.hadley.nz/r6.html


# VP_proj_creator ---------------------------------------------------------

#'  Create poolVP
#' @export
#'
#'
#'

VP_proj_creator <- R6Class("VT",

  public = list( model = NULL,

  param = NULL,
  filters_neg_above = tibble(),
  filters_neg_below  = tibble(),
  filters_pos_above = NULL,
  filters_pos_below = NULL,
  data = NULL,
  parameters_default_values = NULL,# c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
  initial_cmt_values = NULL, #c(X1 = 50) # initial compartment values. At least one, and every missing cmt name would be set to 0
  times = NULL, #seq(0,52, 1)
  poolVP = tibble(),
  protocols = NULL,
  param_reduce = NULL,
  param_increase = NULL,
  param_no_impact = NULL,
  targets = NULL,
  timesaver = NULL,

  initialize = function(sourcefile= "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config.r"){
    myEnv <- new.env()
    source(sourcefile, local=myEnv)
    # attach(myEnv, name="sourced_scripts")

    self$model <- myEnv$model_RxODE

    self$data <- myEnv$data_VT
    self$parameters_default_values <- myEnv$parameters_default_values
    self$initial_cmt_values <- myEnv$initial_cmt_values
    self$times <- myEnv$times
    self$protocols <- myEnv$protocols
    self$param_reduce <- myEnv$param_reduce
    self$param_increase <- myEnv$param_increase
    self$param_no_impact <- myEnv$param_no_impact
    self$filters_neg_above <- tibble()
    self$filters_neg_below <- tibble()

    # get param

    param <- myEnv$model_RxODE$params
    param  <-  param[!param%in%  names(myEnv$parameters_default_values)]

    lines_model <- str_split(deparse(myEnv$model_RxODE$model), pattern = "\\\\n")[[1]] %>%
      gsub(pattern = "==", replacement = "nop")

    already_computed <- lines_model[grepl("=", lines_model)] %>%
      gsub(pattern = "=.+", replacement = "")
    param <- param[!param %in% already_computed]
    self$param <- param
  },
  print = function(...) {
    cat("Person: \n")
    cat("  Name: ", self$times, "\n", sep = "")
    cat("  Age: ")
    invisible(self)
  }

  )
)




# Set targets -------------------------------------------------------------

#
#
VP_proj_creator$set("public", "set_targets", function(..., filter = NULL, ntime = 3, manual = NULL, timeforce = NULL){

filter <- enexpr(filter)

if(!is.null(self$self$poolVP))  return("Virtual Patients have already been generated - It is thus no possible to modify the targets.")
if(!is.null(self$filters_neg_below)) if(map_dbl(self$filters_neg_below, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")
if(!is.null(self$filters_neg_above)) if(map_dbl(self$filters_neg_above, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")
if(!is.null(self$filters_pos_above)) if(map_dbl(self$filters_neg_below, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")
if(!is.null(self$filters_pos_below)) if(map_dbl(self$filters_pos_below, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")



if(!is.null(manual)){

  targets <- manuel


}else{

 targets <- data_segment(data = self$data, protocol, cmt, filter = !!filter, ntime = ntime, timeforce = timeforce) %>%
    filter(max != "0")

}


print( as.data.frame(targets))


    # update targets
    self$targets <- targets

    # initialise filters

    filters_list= list()

    for(a in unique(targets$cmt))  filters_list[[a]] <- tibble()

    self$filters_pos_above  <- filters_list
    self$filters_pos_below  <- filters_list



    return(  message("Once any VP's or filters have been created, targets will be locked."))

})


# VP_production -----------------------------------------------------------
VP_proj_creator$set("public", "add_VP", function(VP_df,  saven = 50, drug = NULL, update_at_end = T, time_compteur = F,  fillatend = F, reducefilteratend = F, npersalve = 1000, use_green_filter = F, pctActivGreen = 0.75){

  # protocols = self$protocols

  toadd <- VP_df

  poolVP <-     VP_df %>%
    rowid_to_column("cellid") %>%
    crossing(protocol = unique(self$targets$protocol)) %>%
    rowid_to_column("rowid")


  if(time_compteur == T)  timesaver <- list()


  #
  #   poolVP %>%
  #     filter(round(k2,1) == 1.2 & round(lambda0,2) == 0.07)

  # Add the columns for each output
  self$targets %>%
    filter(protocol %in% unique(poolVP$protocol)) %>%
    pull(cmt) %>%
    unique -> cmts

  col_to_add <- c(paste0(cmts, "_BU"),
                  paste0(cmts, "_AL"))

  for(a in col_to_add) poolVP[a] <- NA


  ## Apply the filters

  if(time_compteur == T){
    t0 <- Sys.time()
    nbef <- nrow(poolVP)
  }

  # remove filter neg

  if(nrow(self$filters_neg_above) > 0){

    # remove neg_above
    # temp <-  self$filters_neg_above[[cmt]]
    # filter_temp <- paste("!(", filter_template[[1]], ")")
    # if(nrow(temp) > 0){
    #   for(a in 1:nrow(temp)){
    #
    #     ref <-temp %>% slice(a)
    #
    #     poolVP <- poolVP %>%
    #       filter(!!parse_expr(filter_temp))
    #
    #   }
    # }
    # # remove neg_below
    # temp <-  self$filters_neg_below[[cmt]]
    # filter_temp <- paste("!(", filter_template[[2]], ")")
    # if(nrow(temp) > 0){
    #   for(a in 1:nrow(temp)){
    #
    #     ref <-temp %>% slice(a)
    #
    #     poolVP <- poolVP %>%
    #       filter(!!parse_expr(filter_temp))
    #
    #   }
    # }
    # poolVP <- poolVP %>%
    #   select(-rowid) %>%
    #   rowid_to_column()

  }

  # rm filter pos
  for(cmt in cmts){

    filter_template <- self$make_filters(cmt)%>%
      gsub(pattern = "line\\$", replacement = "")



    # set pos abov
    temp <-  self$filters_pos_above[[cmt]]
    if(is.null(temp)) temp <- tibble()
    filter_temp <- paste("(", filter_template[[1]], ")")
    if(nrow(temp) > 0){
      for(a in 1:nrow(temp)){

        ref <-temp %>% slice(a)

        poolVP %>%
          filter(!!parse_expr(filter_temp)) %>%
          pull(rowid) ->rowidstemp

        poolVP[[paste0(cmt, "_AL")]][rowidstemp] <- rep(T, length(rowidstemp))

      }
    }
    # set pos below
    temp <-  self$filters_pos_below[[cmt]]
    filter_temp <- paste("(", filter_template[[2]], ")")
    if(is.null(temp)) temp <- tibble()
    if(nrow(temp) > 0){
      for(a in 1:nrow(temp)){

        ref <-temp %>% slice(a)

        poolVP %>%
          filter(!!parse_expr(filter_temp)) %>%
          pull(rowid) ->rowidstemp

        poolVP[[paste0(cmt, "_BU")]][rowidstemp] <- T

      }
    }
  }
  if(time_compteur == T){
    timesaver$filtre_prev <- tibble(time = difftime(Sys.time(), t0, units = "s"), nrem =nbef -  nrow(poolVP) )

  }


  if(nrow(self$poolVP) > 0 ) poolVP <- poolVP %>%
    mutate(cellid = cellid + max(self$poolVP$cellid),
           rowid = rowid + max(self$poolVP$rowid))



  ### eviter les loupes infinis si un protocole n'as pas d'observion (ex no PK for control PD group)..

  if(time_compteur == T) t0 <- Sys.time()
  crossing(protocols = unique(toadd$protocol), cmt =  self$targets %>%
             filter(protocol %in% unique(poolVP$protocol)) %>%
             pull(cmt) %>%
             unique) %>%
    full_join(self$targets) %>%
    filter(is.na(time)) -> torem

  if(time_compteur == T) timesaver$toremfilter <- difftime(Sys.time(),t0, units = "s")

  if(nrow(torem)>0){

    for(a in 1:nrow(torem)){

      cmt_to_rm <- torem$cmt[[a]]
      pro <- torem$protocol[[a]]

      poolVP[[paste0(cmt_to_rm,"_BU")]][poolVP$protocol == pro] <- FALSE
      poolVP[[paste0(cmt_to_rm,"_AL")]][poolVP$protocol == pro] <- FALSE
    }

  }

  if(time_compteur == T) timesaver$toremfilter2 <- difftime(Sys.time(),t0, units = "s")


  # Handling death and survival agents, to greatly accelerate the process
  all_param <- self$param


  filters <-  list()

  for(a in unique(self$targets$cmt)){

    filters[[a]] <- self$make_filters(cmt = a)   %>%
      gsub(pattern = "line\\$", replacement = "" )

    }


  # line_compar <- paste0("which(", line_compar,")")
  # saveRDS(object = poolVP, file = gsub("\\.RDS", "_todetermine.RDS", file))


  # Time compteur

  if(time_compteur == T){

    timesaver$poolVP_compteur <- tibble(n = NA, time  = NA, nelim =NA, ninfo =NA_real_, computmodel = NA)
    n_compteur <- 0
  }



  siml <- tibble(cellid = integer(), protocol = character())
  # just in case we never enter into the loop (if already filled, almost always useless)

  ntotal <- nrow(poolVP)
  t00 <- Sys.time()


  maxinfo <- is.na(poolVP[, col_to_add]) %>% sum

  pb <- progress_bar$new(
    format = "  VP creation [:bar] :current/:total (:percent) in :elapsed",
    total = maxinfo, clear = FALSE, width= 60)
  # pb <- progress_bar$new(total = )



  slice0 <- poolVP   %>% slice(0) %>% select(!!!parse_exprs(all_param[all_param != "protocol"])) # filre_neg_above

  neg_below <- neg_above <-  pos_below <- pos_above <- list() # filtre_neg_below


  for(a in self$targets$cmt){

    pos_below[[a]] <- slice0
    pos_above[[a]] <- slice0

  }
  use_red_filter <- T

  message(green("Start main loop"))
  # begining while lopp----------------------------------------------------------
  while(is.na(poolVP[, col_to_add]) %>% sum > 0){

    newratio <- is.na(poolVP[, col_to_add]) %>% sum
    pb$update(ratio = (maxinfo -newratio )/maxinfo)

    if(time_compteur == T){

      n_compteur <- n_compteur + 1

      poolVP_compteur_new <- timesaver$poolVP_compteur %>%
        slice(1) %>%
        mutate(n = n_compteur)
    }

    # Just compute some stat...
    t0 <- Sys.time()


    # Sample one rows among the not done yet

    # filter =  # which(is.na(poolVP[, col_to_add]) %>% apply(1, sum) != 0 )

    filter_to_use <-  paste0("is.na(", col_to_add,")") %>%
      paste0(collapse = "|")

    line <- poolVP %>%
      filter(!!parse_expr(filter_to_use))

    line <- line %>%
      filter(protocol == sample(line$protocol, 1))

    line <- line %>%
      sample_n(min(npersalve, nrow(line))) %>%
      rowid_to_column("id")

    if(time_compteur == T) poolVP_compteur_new$timesampleline <- difftime(Sys.time(), t0, units = "s")
    # Now we need to handle the administrations
    # by making a temporar copy
    # protocol  <- self$protocols[[line$protocol[[1]]]]
    if(time_compteur == T) t02 <- Sys.time()

    protocol <-  line %>%
      mutate(protocol2 = map(protocol, ~ self$protocols[[.x]])) %>%
      select(id, protocol2) %>%
      unnest(protocol2)


    if(time_compteur == T) poolVP_compteur_new$protocoldf <- difftime(Sys.time(), t02, units = "s")
    # add_events_line$amt[is.na(add_events_line$amt )] <- 0

    # And now we can make the simulation and extract the result !
    time_simulations <- Sys.time()

    b <- Sys.time()
    res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                       icmt = self$initial_cmt_values, time_vec =self$times,
                       pardf = self$parameters_default_values, model = self$model)#;res

    res <- as.tibble(res)

    time_simulations <-  difftime(Sys.time(), b, units = "s")
    n_simulations <-  nrow(line)
    if(time_compteur == T){
      poolVP_compteur_new$timemodel <- time_simulations
      poolVP_compteur_new$nsimul <- n_simulations
    }
    # if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
      # mutate(computmodel = as.double(difftime(Sys.time(),b, units = "sec")))



      # filter(protocol %in% line$protocol)


    remv <- F
    cellidtorem <- double()


      # targets_temp2 <- targets_temp %>%
        # filter(cmt == cmtt, protocol == line$protocol[[1]])
      # below upper?
      # res %>%
      #   filter(time %in% targets_temp2$time) %>%
      #   pull(!!parse_expr(cmtt)) -> values

      if(time_compteur == T) t02 <- Sys.time()



    # # get targets for this patients
    targets_temp <- self$targets %>%
      filter(protocol == unique(line$protocol))#%>%

      res2 <- res %>%
        gather("cmt", "value", unique(targets_temp$cmt)) %>%
        filter(time %in% targets_temp$time) %>%
        left_join( targets_temp, by = c("time", "cmt")) %>%
        filter(!is.na(min)) %>%
        mutate(be_up = value <= max) %>%
        mutate(ab_low = value >= min)

# Red filter --------------------------------------------------------------
      # If no red filter
      if(use_red_filter == F){


        t02 <- Sys.time()
        res2 %>%
          left_join(line %>% distinct(id, cellid), by  = "id") %>%
          filter(be_up == F | ab_low == F) %>%
          pull(cellid) %>%
          unique -> cellidstorem

        poolVP <- poolVP %>%
          filter(! cellid %in% cellidstorem)

        difftime(Sys.time(), t02, units = "s")


      }else
        {

      # if red filter
      if(time_compteur == T) poolVP_compteur_new[paste0("res2_", cmtt)] <- difftime(Sys.time(), t02, units = "s")


##### Remove patients neg  and compute red filters if activated ########
      # lines output neg_ below
      if(time_compteur == T) t02 <- Sys.time()

      t0filternegbelow <- Sys.time()

      res2 %>%
        filter(ab_low == 0) %>%
        group_by(id, cmt) %>%
        slice(1) %>% ungroup() -> idsbelow

      filters_neg_below <- idsbelow %>%
        left_join(line,by = c("id", "protocol")) %>%
        select(!!!parse_exprs(all_param), "cmt","cellid")

      if(nrow(filters_neg_below) > 0){
      filters_neg_below_up_reduc <- filter_reduc(filters_neg_below ,obj = self, direction = "below")
      #%>%
      # arrange(k2, desc(lambda0))

      }else{

        filters_neg_below_up_reduc <- filters_neg_below
      }

      if(time_compteur == T){
        poolVP_compteur_new$filter_neg_below <- difftime(Sys.time(), t02, units = "s")
        poolVP_compteur_new$nfilter_negbel_bef <- nrow(filters_neg_below)
        poolVP_compteur_new$nfilters_negbel_af <- nrow(filters_neg_below_up_reduc)
      }
      #

      timefilternegbelowmake <- difftime(Sys.time(),t0filternegbelow , units = "s")

      # lines output neg_ above
      if(time_compteur == T) t02 <- Sys.time()

      t0filternegabovemake <- Sys.time()

      res2 %>%
        filter(be_up == 0) %>%
        group_by(id, cmt) %>%
        slice(1) %>% ungroup() ->idsabove


      filters_neg_above <-  idsabove %>%
        left_join(line,by = c("id", "protocol")) %>%
        select(!!!parse_exprs(all_param), "cmt" , "cellid")

      if(nrow(filters_neg_above) > 0){
      filters_neg_above_reduc <- filter_reduc(filters_neg_above ,obj = self, direction =  "above")
# %>%
  # arrange(desc(k2))
      }else{

        filters_neg_above_reduc <- filters_neg_above
      }

       if(time_compteur == T){
        poolVP_compteur_new$filter_neg_above <- difftime(Sys.time(), t02, units = "s")
        poolVP_compteur_new$nfilter_negab_bef <- nrow(filters_neg_above)
        poolVP_compteur_new$nfilters_negab_af <- nrow(filters_neg_above_reduc)
       }

      timefilternegabovemake <- difftime(Sys.time(),t0filternegabovemake , units = "s")

      if(time_compteur == T) t02 <- Sys.time()
      ## We remove all the one not accepted
      poolVP <- poolVP %>%
        filter(! cellid %in% unique(filters_neg_above$cellid) & !  cellid %in% unique(filters_neg_below$cellid))





#### Apply red filter if activated

      t02 <- Sys.time()
      nref <- nrow(poolVP)
      if(time_compteur == T){
        t02 <- t02
        nref <- nref
      }
      if(nrow(filters_neg_above_reduc) > 0){
        for(a in 1:nrow(filters_neg_above_reduc)){

          ref <- filters_neg_above_reduc %>% slice(a)

          poolVP <- poolVP %>%
            mutate(test = !!parse_expr(filters[[ref$cmt]][["above"]])) %>%
            filter(test == F)
        }
      }

      if(time_compteur == T){
        poolVP_compteur_new$remneg_above_fil <- difftime(Sys.time(), t02, units = "s")
        poolVP_compteur_new$nremoved_above_fil <-   nref - nrow(poolVP)
      }


      if(nrow(filters_neg_below_up_reduc)> 0){
        for(a in 1:nrow(filters_neg_below_up_reduc)){

          ref <- filters_neg_below_up_reduc %>% slice(a)

          poolVP <- poolVP %>%
            mutate(test = !!parse_expr(filters[[ref$cmt]][["below"]])) %>%
            filter(test == F)
        }
      }
      if(time_compteur == T){
        poolVP_compteur_new$remneg_below_fil <- difftime(Sys.time(), t02, units = "s")
        poolVP_compteur_new$nremoved_below_fil <-   nref - nrow(poolVP)
      }

      nrem <- nref -  nrow(poolVP)


      time_filter_neg_apply <- difftime(Sys.time(), t02, units = "s")
      nremoved_below_fil <-  nrow(poolVP)

       if(time_compteur == T){
                 poolVP_compteur_new$remneg_below_fil <- time_filter_neg_above_apply
                 poolVP_compteur_new$nremoved_below_fil <- nremoved_below_fil
              }
##### Compute the rendement of red filter and disable if negative

      totaltimeredfilter <- time_filter_neg_apply + timefilternegabovemake + timefilternegbelowmake
      totalsave <-  nrem * time_simulations /  n_simulations

      if(totalsave < totaltimeredfilter){

      message(red("\nRed filter system disabled."))
      use_red_filter <- F
}
      }


      # gree filter --------------------------------------------------------------



      if(use_green_filter == F | nrow(res2 %>% filter(be_up == T & ab_low== T)) < pctActivGreen * nrow(res2)  ){


        if(time_compteur == T) t02 <- Sys.time()

       targets_temp %>%
          # filter(protocol == unique(line$protocol), cmt == cmtt) %>%
          pull(time) %>% length -> nmax


        res2 %>%
          filter(be_up == 1 & ab_low == 1) %>%
          group_by(id) %>%
          tally -> nabove

        res2 %>%
          left_join(nabove, by = "id") %>%
          filter(n == nmax) %>%
          group_by(id) %>%
          slice(1) %>% pull(id) -> idsgood

        line %>% filter(id %in% idsgood) %>% pull(cellid) -> idsgood

        poolVP[poolVP$cellid %in% idsgood & poolVP$protocol == unique(line$protocol), col_to_add     ] <- T


        if(time_compteur == T)  poolVP_compteur_new$time_addgreennofil <- difftime(Sys.time(), t02, units = "s")




        if(time_compteur == T) t02 <- Sys.time()

        res %>%
          left_join(line %>% distinct(id, cellid), by = "id") %>%
          filter(cellid %in% idsgood) %>%
          mutate(protocol = unique(line$protocol)) %>%
          group_by(cellid,protocol) %>%

          nest() -> forjoin



        siml <- siml %>%
          bind_rows(forjoin)

        if(time_compteur == T) poolVP_compteur_new$timesimlandjoin <- difftime(Sys.time(), t02, units = "s")





      }else {

    for(cmtt in unique(self$targets$cmt)){
      if(time_compteur == T){

        befusegreen <-   poolVP %>%
            group_by(tumVol_BU, tumVol_AL) %>%
            tally
       t0green   <- Sys.time()
      }
        # lines output below_up

        if(time_compteur == T) t02 <- Sys.time()
        self$targets %>%
          filter(protocol == unique(line$protocol), cmt == cmtt) %>%
          pull(time) %>% length -> nmax


        res2 %>%
          filter(be_up == 1 & ab_low == 1) %>%
          group_by(id) %>%
          tally -> nabove

        res2 %>%
          left_join(nabove, by = "id") %>%
          filter(n == nmax) %>%
          group_by(id) %>%
          slice(1) %>% pull(id) -> idsbelowpos

        line %>% filter(id %in% idsbelowpos) -> filters_pos_below_up
        filters_pos_below_up$cmt = cmtt
        filters_pos_below_up_reduc <- filter_reduc(filters_pos_below_up , filtre =  "below")

        # %>%
        #   arrange(k2, desc(lambda0))

        pos_below[[cmtt]] <- bind_rows(pos_below[[cmtt]], filters_pos_below_up_reduc[names(pos_below[[1]])])


        if(time_compteur == T){
          poolVP_compteur_new$filter_pos_above <- difftime(Sys.time(), t02, units = "s")
          poolVP_compteur_new$nfilter_pobe_bef <- nrow(filters_pos_below_up)
          poolVP_compteur_new$nfilters_posbe_af <- nrow(filters_pos_below_up_reduc)
        }

        # lines output bpostif above lo
        if(time_compteur == T) t02 <- Sys.time()

        res2 %>%
          filter(be_up == 1 & ab_low == 1) %>%
          group_by(id) %>%
          tally -> nbelow

        res2 %>%
          left_join(nbelow, by = "id") %>%
          filter(n == nmax) %>%
          group_by(id) %>%
          slice(1) %>% pull(id) -> idsabovepos


        line %>% filter(id %in% idsabovepos) -> filters_ps_above_lo

        filters_ps_above_lo$cmt = cmtt

        filters_ps_above_lo_reduc <- filter_reduc(filters_ps_above_lo, filtre = "above")
        # %>%
        #   arrange(desc(k2), lambda0)

        pos_above[[cmtt]] <- bind_rows(pos_above[[cmtt]], filters_ps_above_lo_reduc[names(pos_above[[1]])])

        if(time_compteur == T){
          poolVP_compteur_new$filter_pos_above <- difftime(Sys.time(), t02, units = "s")
          poolVP_compteur_new$nfilter_posab_bef <- nrow(filters_ps_above_lo)
          poolVP_compteur_new$nfilters_posab_af <- nrow(filters_ps_above_lo_reduc)
        }

        ###### Use the filters

        if(time_compteur == T){
          t02 <- Sys.time()
          reff <- sum(is.na(poolVP$tumVol_AL))

        }
        if(nrow(filters_ps_above_lo_reduc) >0){

          for(a in 1:nrow(filters_ps_above_lo_reduc)){

            ref <- filters_ps_above_lo_reduc %>% slice(a)

            poolVP %>%
              mutate(test = !!parse_expr(filters[[ref$cmt]][["above"]])) %>%
              filter(test == T) %>% pull(cellid) -> id_temp

            poolVP$tumVol_AL[poolVP$cellid %in% id_temp & poolVP$protocol == unique(line$protocol)] <- T
          }

        }
        if(time_compteur == T){
          poolVP_compteur_new$time_add_above_fil <- difftime(Sys.time(), t02, units = "s")
          poolVP_compteur_new$n_add_above_fil <-   sum(is.na(poolVP$tumVol_AL)) - reff
        }


        if(time_compteur == T){
          t02 <- Sys.time()
          reff <- sum(is.na(poolVP$tumVol_BU))

        }
        if(nrow(filters_ps_above_lo_reduc) >0){
          for(a in 1:nrow(filters_pos_below_up_reduc)){

            ref <- filters_pos_below_up_reduc %>% slice(a)

            poolVP %>%
              mutate(test = !!parse_expr(filters[[ref$cmt]][["below"]])) %>%
              filter(test == T) %>% pull(cellid) -> id_temp

            poolVP$tumVol_BU[poolVP$cellid %in% id_temp & poolVP$protocol == unique(line$protocol)] <- T
          }

          if(time_compteur == T){
            poolVP_compteur_new$time_add_below_fil <- difftime(Sys.time(), t02, units = "s")
            poolVP_compteur_new$n_add_below_fil <-   sum(is.na(poolVP$tumVol_BU)) - reff
          }
      }



        if(time_compteur == T){

          afusegreen <-   poolVP %>%
            group_by(tumVol_BU, tumVol_AL) %>%
            tally  %>%
            rename(n2 = n)

          left_join(afusegreen, befusegreen) %>%
            mutate(Diff = n2 - n) %>%
            mutate(colname = paste0(tumVol_BU, "-", tumVol_AL)) -> temp3

          for(a in 1:nrow(temp3)){

            poolVP_compteur_new[temp3$colname[[a]]] <-  temp3$Diff[[a]]

          }
          poolVP_compteur_new$wholegreenfilter <- difftime(Sys.time(), t0green, units = "s")
        }

        if(time_compteur == T) t02 <- Sys.time()

        res %>%
          left_join(line %>% distinct(id, cellid), by = "id") %>%
          filter(cellid %in% unique(poolVP$cellid)) %>%
          mutate(protocol = unique(line$protocol)) %>%
          group_by(cellid,protocol) %>%

          nest() -> forjoin



        siml <- siml %>%
          bind_rows(forjoin)

        if(time_compteur == T) poolVP_compteur_new$timesimlandjoin <- difftime(Sys.time(), t02, units = "s")


  } # end use green filter


    }# end for each compartment



    # if() siml <- tibble(cellid = integer(), protocoles = character(), simul= list())

    if(time_compteur == T)  timesaver$poolVP_compteur <- bind_rows(timesaver$poolVP_compteur,poolVP_compteur_new)



    # print(nn)
  }# fin while 1


  ## Need to remove filter pos for which some have been delated (because of multiple protocol,one ok, one bad)

  # pos_below <- map(pos_below, function(x) x %>% filter(cellid %in% poolVP$cellid) %>% select(-cellid))
  # pos_above <- map(pos_above, function(x) x %>% filter(cellid %in% poolVP$cellid) %>% select(-cellid))

  ## Here happen after nsave iteration
  # print("########################### SAVNG RDS #############################")

  # if(time_compteur == T)  poolVP_compteur <<- poolVP_compteur
  print(Sys.time() - t00)

  if(time_compteur == T)  t02 <- Sys.time()
  poolVP <- poolVP %>%
    select(-starts_with("simul")) %>%
    left_join(siml  %>%  rename(simul = data))

  if(time_compteur == T)  timesaver$joinsimul <- difftime(Sys.time(), t02, units = "s")

  # Update filters
  if(time_compteur == T)  t02 <- Sys.time()
  for(a in unique(self$targets$cmt)){

    if(is.null(self$filters_pos_above[[a]])) self$filters_pos_above[[a]] = tibble()
  self$filters_pos_above[[a]] <- bind_rows(self$filters_pos_above[[a]] %>% mutate(PrimFilter = T) , pos_above[[a]] )

  if(is.null(self$filters_pos_below[[a]])) self$filters_pos_below[[a]] = tibble()
   self$filters_pos_below[[a]] <- bind_rows(self$filters_pos_below[[a]] %>% mutate(PrimFilter = T) , pos_below[[a]] )

  }


  self$filters_neg_above <- bind_rows(self$filters_neg_above %>% mutate(PrimFilter = T) , filters_neg_above_reduc )
  self$filters_neg_below <- bind_rows(self$filters_neg_below %>% mutate(PrimFilter = T) , filters_neg_below_up_reduc )


  if(time_compteur == T)  timesaver$updatefilters <- difftime(Sys.time(), t02, units = "s")


  # bind_rows(self$poolVP, poolVP )


  # save pooVP
  if(is.null(self$poolVP)){

    self$poolVP <- poolVP

  }else{
    self$poolVP <- bind_rows(self$poolVP, poolVP )
  }


  if(time_compteur == T)   self$timesaver <- timesaver
  # Fill missing profile if required
  if(fillatend){

    print("filling")
    self$fill_simul()
  }

  if(reducefilteratend){

    print("filter reduction")
    self$n_filter_reduc()
  }
  # # Recompute the whole poolVP

})


# Complete VP simul -------------------------------------------------------

VP_proj_creator$set("public", "fill_simul", function(){



self$poolVP %>%
    mutate(missing = map_lgl(simul, ~if_else(is.null(.x), TRUE, FALSE ))) %>%
    filter(missing) %>%
    rowid_to_column("id")-> line

  if(nrow(line) == 0) return(message("All VP's have simulations already available"))

  protocol <-  line %>%
    mutate(protocol2 = map(protocol, ~ self$protocols[[.x]])) %>%
    select(id, protocol2) %>%
    unnest()
  # add_events_line$amt[is.na(add_events_line$amt )] <- 0

  # And now we can make the simulation and extract the result !

  b <- Sys.time()
  res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                      icmt = self$initial_cmt_values, time_vec =self$times,
                      pardf = self$parameters_default_values, model = self$model);res


  res <- as.data.frame(res)

  res %>%
    left_join(line %>% distinct(id, cellid,protocol)) %>%
    group_by(cellid, protocol) %>%
    nest() -> newsimuls


  temp <- self$poolVP %>%
    mutate(test= map_lgl(simul, ~ is.null(.x)))

  bind_rows(

    temp %>% filter(test == F),

    temp %>% filter(test == T) %>%
      select(-simul) %>%
      left_join(newsimuls %>% rename(simul = data))

  ) %>%
    arrange(cellid) %>%
    select(-test) -> output

  self$poolVP <- output

})


# VP plot -----------------------------------------------------------------

VP_proj_creator$set("public", "plot_VP", function(){


  self$poolVP %>%
    # slice(1:100) %>%
    unnest(simul) %>%
    gather("cmt", "value", !!!parse_exprs(unique(self$targets$cmt))) %>%
    filter( ! (cmt == "Conc" & value == 0)) %>%
    filter( ! (cmt == "Conc" & time > 15)) %>%
    ggplot()+
    geom_line(aes(time, value, group = cellid))+
    facet_wrap(protocol~cmt, scales = "free")+
    scale_y_log10()+
    geom_segment(data = self$targets ,
                 aes(x = time, xend = time, y = min, yend = max), col ="red", size = 2)


})


# make_filters ------------------------------------------------------------


VP_proj_creator$set("public", "make_filters", function(cmt = "tumVol"){

  all_param <- self$param

  line_compar <-   paste0("line$", all_param, " == ref$", all_param) %>%
    paste0(collapse = " & ")


  # above
  line_above <- line_compar
  line_below <- line_compar


  # self$filters_neg_above

  for(a in self$param_increase[[cmt]]){

    line_above <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above)
    line_below <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below)

  }


  for(a in self$param_reduce[[cmt]]){

    line_above <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above)
    line_below <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below)
  }

  for(a in  self$param_no_impact[[cmt]]){

    line_above <-  gsub(paste0("&? * line\\$", a, " *== *ref\\$",a),"", line_above)
    line_below <- gsub(paste0("&? * line\\$", a, " *== *ref\\$",a),"", line_below)
  }

  return(c(above = line_above, below = line_below))

})


# plot 2D -----------------------------------------------------------------
# x = expr(k2)
# y = expr(lambda0)
# y = expr(ke)
# toaddneg = VP_df
VP_proj_creator$set("public", "plot_2D", function(x, y , cmt_green = "tumVol", toaddneg = NULL, plotMain = F, add_point =F , IDref = NULL, plotoreturn = 3){

  x <- enexpr(x)
  y <-  enexpr(y)



  # Apply filter to keep only on plan

  # if( is.null(IDref)) IDref <- sample(  self$poolVP$cellid, size = 1)
  #
  # line <- self$poolVP %>%
  #   filter(cellid == IDref)
  #
  # namesparam <- names(line)
  # namesparam <- namesparam[! namesparam %in% c("rowid", "cellid", "protocol", "simul", deparse(x), deparse(y))]
  # namesparam <- namesparam[!grepl("(_BU$)|(_AL$)", namesparam)]
  #
  # filtre <- line[, namesparam] %>%
  #   unlist() %>%
  #   imap_chr(~ paste0(.y, " == " ,.x)) %>%
  #   paste0(collapse = " & ") %>%
  #   parse_expr
  # end apply filter

 ytype <- names(self$filters_pos_above)[[1]]#Conc

  neg_above <- invoke(self$filters_neg_above, .fn = bind_rows) #%>% filter(!!filtre)
  neg_below <- invoke(self$filters_neg_below, .fn = bind_rows) #%>% filter(!!filtre)
  pos_above <- invoke(self$filters_pos_above, .fn = bind_rows) #%>% filter(!!filtre)
  pos_below <- invoke(self$filters_pos_below, .fn = bind_rows) #%>% filter(!!filtre)




    plot_dot <-

      ggplot()+
      geom_point(data = pos_above, aes(!!x, !!y), col = "darkgreen") +
      geom_point(data = pos_below, aes(!!x, !!y), col = "darkgreen") +
      geom_point(data = neg_below, aes(!!x, !!y, shape = cmt), col = "chocolate") +
      geom_point(data = neg_above , aes(!!x, !!y, shape = cmt), col = "red", alpha = 1)+
      # geom_point(data = self$filters_neg_below, aes(k2, lambda0), col = "red", alpha = 1)+
      theme_bw()+
      ggtitle( "VP tested")+
      theme(plot.title = element_text(hjust = 0.5)); plot_dot

if(length(unique(self$targets$cmt) == 1))  plot_dot <- plot_dot + guides(shape = NULL)

if(!is.null(toaddneg)){

  toaddneg %>%
    left_join(self$poolVP) %>%
    filter(is.na(rowid)) -> addneg

  plot_dot <- plot_dot +
    geom_point(data = addneg, aes(!!x, !!y), col = "red")+
    geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen")

}


    if(add_point == T) plot_dot <-   plot_dot+
      geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen")

    rectangles_above <- tibble()
    rectangles_below <- tibble()

  for(a in unique(self$targets$cmt) ){


    xana <- case_when(deparse(x) %in% self$param_increase[[a]] ~ "inc",
                      deparse(x) %in% self$param_reduce[[a]] ~ "dec",
                      deparse(x) %in% self$param_no_impact[[a]] ~ "no_impact",
                      T ~ "None")


    yana <- case_when(deparse(y) %in% self$param_increase[[a]] ~ "inc",
                      deparse(y) %in% self$param_reduce[[a]] ~ "dec",
                      deparse(y) %in% self$param_no_impact[[a]] ~ "no_impact",
                      T ~ "None")


    rectangles_above  <- neg_above %>%
      filter(cmt == a) %>%
      mutate(xmin = case_when(xana %in% c("dec", "no_impact") ~ 0, T ~ !!x),
             xmax = case_when(xana %in% c("dec","None") ~ !!x, T ~Inf),
             ymin = case_when(yana %in% c("dec", "no_impact") ~ 0, T ~!!y),
             ymax = case_when(yana  %in% c("dec","None") ~ !!y, T ~Inf)) %>%
            bind_rows(rectangles_above)

    rectangles_below  <- neg_below %>%
      filter(cmt == a) %>%
      mutate(xmin = case_when(xana %in% c("inc","no_impact") ~ 0, T ~ !!x),
             xmax = case_when(xana %in% c("inc","None") ~ !!x, T ~Inf),
             ymin = case_when(yana %in% c("inc","no_impact") ~ 0, T ~!!y),
             ymax = case_when(yana %in% c("inc","None") ~ !!y, T ~Inf))%>%
      bind_rows(rectangles_below)




  }

rectangles <- bind_rows(rectangles_above %>% mutate(where = "above"), rectangles_below %>% mutate(where = "below"))


  plot1 <-  plot_dot +
    geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,  fill = where, col = where), alpha = 0.3)+
    ggtitle( "zone rejection")+
    scale_color_manual(values = c("red", "chocolate"))+
    scale_fill_manual(values = c("red", "chocolate"));plot1

  if(plotoreturn == 1) return(plot1)
  # plot1 +
  #   geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen", alpha = 0.2)
  #



  rectangles_above_lower <-  pos_above %>%
          mutate(xmin = case_when(xana == "dec" ~ 0, T ~ !!x),
         xmax = case_when(xana == "dec" ~ !!x, T ~Inf),
         ymin = case_when(yana == "dec" ~ 0, T ~!!y),
         ymax = case_when(yana == "dec" ~ !!y, T ~Inf))

  rectangles_below_upper  <- pos_below%>%
    mutate(xmin = case_when(xana == "inc" ~ 0, T ~ !!x),
           xmax = case_when(xana == "inc" ~ !!x, T ~Inf),
           ymin = case_when(yana == "inc" ~ 0, T ~!!y),
           ymax = case_when(yana == "inc" ~ !!y, T ~Inf))





 plot2 <-   plot_dot +
    geom_rect(data = rectangles_below_upper, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "darkgreen",  col = "darkgreen")+
   ggtitle( "zone below upper limit")

 plot3 <-  plot_dot +
    geom_rect(data = rectangles_above_lower, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "darkgreen",  col = "darkgreen")+
   ggtitle( "zone above lower limit")

 # plot_dot +
 #   geom_rect(data = rectangles_above_lower, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "red",  col = "red")+
 #   geom_rect(data = rectangles_below_upper, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "blue",  col = "blue")
 #   ggtitle( "zone above lower limit")

if(plotMain == F) return(plot_grid(plot_dot, plot1, plot2, plot3))


 self$poolVP %>%
   select(!!x, !!y) -> temp


 crossing(a = unique(temp[[deparse(x)]]), b = unique(temp[[deparse(x)]])) %>%
   filter(b >= a) %>%
   mutate(height = map2(a, b, function(a,b){

     temp %>%
       filter(!!x>=a & !!x<=b) %>%
       group_by(!!x) %>%
       summarise(max = max(!!y), min = min(!!y)) -> temp2

     tibble( floor = max(temp2$min), ceiling = min(temp2$max))
     # temp2

   })) %>%
   unnest() %>%
   filter(ceiling >= floor) %>%
   arrange(desc(b)) %>%
   group_by(a, floor, ceiling) %>%
   slice(1) %>%
   rename(xmin = a, xmax = b, ymin = floor, ymax = ceiling)-> col_green
#
#
 plot4<-
   ggplot()+
     geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,  fill = where, col = where), alpha = 1)+
   scale_color_manual(values = c("red", "chocolate"))+
   scale_fill_manual(values = c("red", "chocolate"))+
   geom_rect(data = col_green, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 1,  fill = "darkgreen",  col = "darkgreen")+
   theme_bw()+
   ggtitle( "Total")+
   theme(plot.title = element_text(hjust = 0.5))


 if(!is.null(toaddneg)){


     plot4 <- plot4 +
     geom_point(data = addneg, aes(!!x, !!y), col = "red", alpha = 0)+
     geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen", alpha = 0)

 }

plot_grid(plot4,plot_grid(plot_dot, plot1, plot2, plot3))

})


# plot 3 D ----------------------------------------------------------------
# x = expr(k2)
# y = expr(lambda0)
# z = expr(Vd)
VP_proj_creator$set("public", "plot_3D", function(x, y, z,  toaddneg = NULL,add_point =F ){

  x <- enexpr(x)
  y <- enexpr(y)
  z <- enexpr(z)

  neg_above <- invoke(self$filters_neg_above, .fn = bind_rows) #%>% filter(!!filtre)
  neg_below <- invoke(self$filters_neg_below, .fn = bind_rows) #%>% filter(!!filtre)
  pos_above <- invoke(self$filters_pos_above, .fn = bind_rows) #%>% filter(!!filtre)
  pos_below <- invoke(self$filters_pos_below, .fn = bind_rows) #%>% filter(!!filtre)


  allpoints <- bind_rows(self$poolVP %>% mutate(test = "Accepted"),
                         neg_above %>% mutate(test = "Rejected_Above"),
                         neg_below %>% mutate(test = "Rejected_Below")

                         )


 pltly <-  plot_ly()%>%
    add_markers(type = "scatter3d",
                mode = "markers",
                data = allpoints,
                x = ~lambda0,
                y = ~k2,
                z = ~Vd,
                color = ~test,
                opacity = 1,
                colors = c('darkgreen', "orange", 'red'))





  # Add above

  # neg_above %>%
  #   slice(1)

 neg_above %>%
   mutate(a = pmap(list(k2, lambda0, Vd), function(k2, lambda0, Vd){

     filtrecube <- crossing(k2 = c(0, k2),
                            lambda0 = c(lambda0, 0.18),
                            Vd  = c(Vd, 50)) %>%
                    slice(1,2,6,5,3,4,8,7)

     pltly <<- pltly %>%
       add_trace(type = 'mesh3d',
                 data = filtrecube,
                 x = ~lambda0,
                 y = ~k2,
                 z = ~Vd,
                 i = c(7, 0, 0, 0, 4, 4, 6, 1, 4, 0, 3, 6),
                 j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3),
                 k = c(0, 7, 2, 3, 6, 7, 1, 6, 5, 5, 7, 2),
                 facecolor = rep("orange", 12),
                 opacity = 0.4
       )

     "r"

   }))

 # Add below


 neg_below %>%
   mutate(a = pmap(list(k2, lambda0, Vd), function(k2, lambda0, Vd){

     filtrecube <- crossing(k2 = c(k2, 3),
                            lambda0 = c(0,  lambda0),
                            Vd  = c(0, Vd)) %>%
       slice(1,2,6,5,3,4,8,7)


     pltly <<- pltly %>%
 add_trace(type = 'mesh3d',
           data = filtrecube,
           x = ~lambda0,
           y = ~k2,
           z = ~Vd,
           i = c(7, 0, 0, 0, 4, 4, 6, 1, 4, 0, 3, 6),
           j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3),
           k = c(0, 7, 2, 3, 6, 7, 1, 6, 5, 5, 7, 2),
           facecolor = rep("red", 12),
           opacity = 0.4
 )
 "r"

   }))



print(pltly)
"Done"

})
# filter reduc -----------------------------------------------------------------

#'  Create poolVP
#' @export
#'
#'
filter_reduc <- function(df = filters_neg_up, filtre = NULL , obj = NULL, direction = "below") { #filters[[1]]

cmts <- unique(df$cmt)


if(!is.null(obj)){


  filters_per_cmt <- list()

  ## look first for cmt with ni ==o impacting parameters
cmt_w_indp <- self$param_no_impact[map_lgl(obj$param_no_impact, ~length(.x) > 0) ]
cmt_w_indp <- cmt_w_indp[names(cmt_w_indp) %in% cmts]

for(a in  names(cmt_w_indp)){


  filters_per_cmt[[a]] <-  filter_reduc(df = df %>% filter(cmt == a), direction = direction, obj = NULL, filtre = obj$make_filters(cmt = a)[[direction]]) %>%
    mutate(cmt = a)
}

### remaining cmt

if(length(cmt_w_indp) >0) cmts <- cmts[!cmts %in% names(cmt_w_indp)]


for(a in cmts){

  dftemp <- df %>%
    filter(cmt == a)

### module to apply first the previous filters
  if(length(cmt_w_indp) >0){
for(b in names(filters_per_cmt)){



  filter_temp <- obj$make_filters(b)[[direction]] %>%
    gsub(pattern = "line\\$", replacement = "")

  filter_temp <- paste0("!(", filter_temp, ")")

for(d in 1:nrow(filters_per_cmt[[b]])){

ref <-filters_per_cmt[[b]] %>% slice(d)

dftemp <- dftemp %>%
  filter(!!parse_expr(filter_temp))

} # end for d

} # end for b
  }
#### end module reduce datafilter

  filters_per_cmt[[a]] <-  filter_reduc(df = dftemp, direction = direction, obj = NULL, filtre = obj$make_filters(a)[[direction]])%>%
    mutate(cmt = a)
}

return(filters_per_cmt %>% reduce(bind_rows))
}


  df2 <- df %>%
    select(-starts_with("iddummy")) %>%
    rowid_to_column("iddummy")

  filtre <- gsub("line\\$", "", filtre)

  for(a in df2$iddummy){
    # print(a)

    if(a %in% df2$iddummy){
      ref <- df2 %>%
        filter(iddummy == a)

      df2 %>%
        mutate(test = !!parse_expr(filtre))  -> df2

      df2$test[df2$iddummy == a] <- F
      df2 <- df2 %>%
        filter(test == F)
    }
  }
  df2
}

VP_proj_creator$set("public", "n_filter_reduc", function(){


  filters <- self$make_filters() %>%
    gsub(pattern = "line\\$", replacement = "" )



  # Handle negative above


  self$filters_neg_above %>%
    select(-starts_with("rowid")) %>%
    # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
    rowid_to_column() -> temp

  message("Reducing negative above filters")

  if(nrow(temp) > 0 ){
  for(a in 1:nrow(temp)){
    # print(a)

    if(a %in% temp$rowid){
      ref <- temp %>%
        filter(rowid == a)

      temp %>%
        mutate(test = !!parse_expr(filters[[1]]))  -> temp

      temp$test[temp$rowid == a] <- F
      temp <- temp %>%
        filter(test == F)
    }
    self$filters_neg_above <- temp %>% mutate(PrimFilter = T) %>% select(-test)
  }
  }




  # Handle negative below

  self$filters_neg_below %>%
    select(-starts_with("rowid")) %>%
    # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
    rowid_to_column() -> temp

  message("Reducing negative below filters")

  if(nrow(temp) > 0 ){
  for(a in 1:nrow(temp)){
    # print(a)

    if(a %in% temp$rowid){
      ref <- temp %>%
        filter(rowid == a)

      temp %>%
        mutate(test = !!parse_expr(filters[[2]]))  -> temp

      temp$test[temp$rowid == a] <- F
      temp <- temp %>%
        filter(test == F)
    }

    self$filters_neg_below <- temp %>% mutate(PrimFilter = T) %>% select(-test)
  }
}





  # for each cmtt
for(cmtt in unique((self$targets$cmt))){

# filter pos above <low



self$filters_pos_above[[cmtt]] %>%
  select(-starts_with("rowid")) %>%
  # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
  rowid_to_column() -> temp


if(nrow(temp) > 0){

  message("Reducing positive above filters")
for(a in 1:nrow(temp)){
  # print(a)

  if(a %in% temp$rowid){
    ref <- temp %>%
      filter(rowid == a)

    temp %>%
      mutate(test = !!parse_expr(filters[[1]]))  -> temp

    temp$test[temp$rowid == a] <- F
    temp <- temp %>%
      filter(test == F)
  }
}

  self$filters_pos_above[[cmtt]] <- temp %>% mutate(PrimFilter = T) %>% select(-test)

}


# handle pos belw

self$filters_pos_below[[cmtt]] %>%
  select(-starts_with("rowid")) %>%
  # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
  rowid_to_column() -> temp

message("Reducing positive above filters")

if(nrow(temp) > 0){
  message("Reducing positive above filters")

for(a in 1:nrow(temp)){
  # print(a)

  if(a %in% temp$rowid){
    ref <- temp %>%
      filter(rowid == a)

    temp %>%
      mutate(test = !!parse_expr(filters[[2]]))  -> temp

    temp$test[temp$rowid == a] <- F
    temp <- temp %>%
      filter(test == F)
  }
}
  self$filters_pos_below[[cmtt]] <- temp %>% mutate(PrimFilter = T) %>% select(-test)
}



}


})


# time analyser -----------------------------------------------------------


VP_proj_creator$set("public", "time_analyserr", function(){

  timetoana <- self$timesaver

  timetoanaloop <- timetoana$poolVP_compteur %>%
    filter(!is.na(n))%>%
    select(-time, -nelim, -ninfo, -computmodel)

  # green filter

  timetoanaloop %>%
    distinct(n, wholegreenfilter , `TRUE-NA`,  `TRUE-TRUE`, `NA-TRUE`, time_add_above_fil, time_add_below_fil, time_addgreennofil)




  timetoanaloop
names(timetoanaloop)


# Impact filter neg

timetoanaloop %>%
  select(timemodel, nremoved_above_fil, filter_neg_above, remneg_above_fil, nfilter_negab_bef,nfilters_negab_af ) %>%
  mutate(timesaved = nremoved_above_fil * timemodel / 2000)

timetoanaloop %>%
  select(timemodel, nremoved_below_fil, filter_neg_below, remneg_below_fil, nfilter_negbel_bef,nfilters_negbel_af ) %>%
  mutate(timesaved = nremoved_below_fil * timemodel / 2000)

timetoanaloop %>%
  gather("key", "value", filter_neg_above, remneg_above_fil)

sum(timetoanaloop$remneg_above_fil)







})




# Massive screening -------------------------------------------------------



# 1000^6
# 10^18
VP_proj_creator$set("public", "big_screening", function(domain){


  if(!"param" %in% names(domain) & nrow(domain) == 2){

  domain <-   imap(domain, function(x,y){

      values <- domain[[y]]
      values <- values[order(values)]
      tibble(param = y, from = values[[1]], to = values[[2]]) %>%
       mutate(step = if_else(from == to, 0,1))

    }) %>%
      invoke(.fn = bind_rows)
  }
# nsize = 10^18

param_fluct <- domain %>%
  filter(step != 0)

# How many param for having 20000
nperparam <- 200000^(1/nrow(param_fluct)) %>% ceiling()



param_fluct %>%
  mutate(values = pmap(list(from, to , step), function(from, to , step){

    seqq <- seq(from, to, step)
  indic <- length(seqq) * 0:(nperparam-1)/(nperparam-1)
  indic[indic == 0 ] <- 1
  round(seqq[indic])

  })) %>%
  unnest() -> allparamcut

allparamcut %>%
  group_split(param) %>%
  map(function(x){

    list()

    x %>%
      select(values) -> temp
    names(temp) <- unique(x$param)
    temp
  }) %>%
  invoke(.fn = crossing) -> toadd

# fixed one
domain %>%
  filter(step == 0) -> fixed

if(nrow(fixed)>0){
for(a in 1:nrow(fixed)){

  toadd[[fixed$param[[a]]]] <- fixed$from[[a]]
}
}

nsize/nrow(toadd)

VP_temp <- self#VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner.r")
VP_temp$filters_neg_above <- tibble()
VP_temp$filters_neg_below <- tibble()

VP_temp$add_VP(toadd)

VP_temp$n_filter_reduc()
if(nrow(VP_temp$poolVP) >0) return(VP_temp$poolVP)
# VP_temp$plot_VP()
group_split(VP_temp$filters_neg_above,rowid)[[1]] -> x
nrowss <- nrow(VP_temp$filters_neg_above)
whereVP <- group_split(VP_temp$filters_neg_above,rowid)[1:min(20,nrowss)] %>%
  map(function(x){

    y <- x
    for(a in unique(allparamcut$param)){

      temp <- allparamcut %>% filter(param == a)
      currentloc <- which(temp$values == y[[a]])

      if(currentloc == 1){
        value_if_increase <- 0
        value_if_decrease <- temp$values[[currentloc+1]]
      }else if (currentloc == nrow(temp)){

        value_if_increase <- temp$values[[currentloc-1]]
        value_if_decrease <-  max(temp$values)

      }else{
        value_if_increase <- temp$values[[currentloc-1]]
        value_if_decrease <- temp$values[[currentloc+1]]
      }


      y[[a]] <- case_when(a %in% self$param_increase$Pore ~ value_if_increase ,
                          a %in% self$param_reduce$Pore ~ value_if_decrease)

    }

    bind_rows(x, y)[, domain$param]

  })

return(whereVP)

})

