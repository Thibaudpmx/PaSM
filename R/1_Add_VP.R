## ---------------------------
## Script name: add_VP
##
## Purpose of script: First main algorithm
##
## Author: Thibaud Derippe
##
## Date Created: 2022-04-04 (day of code repartition from R6object.R)
##
## Under GPL-3 License
## Email: thibaud.derippe@gmail.com
## GitHub: https://github.com/Thibaudpmx/QSPVP
## ---------------------------
##
## Notes:
##
##
## ---------------------------

 fix_df = NULL; saven = 50; time_compteur = T;  fillatend = F; reducefilteratend = F; npersalve = 2000; use_green_filter = T;
pctActivGreen = 0.1; keepRedFiltaftDis = T; methodFilter = 2; use_red_filter = T; cmtalwaysbelow = NULL; keep = NULL; saveVPRej = T; timeSave = NULL;PreviousResults = tibble()
RedFilterDisAllProt = F; GreenFilterDisAllProt = F

VP_proj_creator$set("public", "add_VP", function(VP_df, fix_df = NULL, saven = 50, time_compteur = F,  fillatend = F, reducefilteratend = F, npersalve = 1000, use_green_filter = F,
                                                 pctActivGreen = 0.1, keepRedFiltaftDis = T, methodFilter = 2, use_red_filter = T, cmtalwaysbelow = NULL, keep = NULL, saveVPRej = T, timeSave = NULL, PreviousResults = tibble(),
                                                 RedFilterDisAllProt = F, GreenFilterDisAllProt = F){

  tTOTAL <- Sys.time()

  if(time_compteur == T){

    timesaver <- list()
    timesaver$nVP0 <- nrow(VP_df)

  }


  ### Handle the plannification if there are fixed parameter
  if(!is.null(fix_df)){

    self$action_programmed$VP_df <- VP_df
    self$action_programmed$fix_df <- fix_df

    self$launch_programmed()

    return(cat(red("Done")))
  }


  # Create the pool of virtual patient
  # each VP should have a unique id and one row per protocol
  # Two possibilities: use crossing or map protocole followed by bind_rows
  # See test microbrenchmark 1 -> second option clear winner
  # (crossing is a powerfull yet slow function)



  poolVP <-     VP_df %>%
    rowid_to_column("id") # add row equal to id

  poolVP <-  map( unique(self$targets$protocol), ~ poolVP %>% mutate(protocol = .x)) %>% # one df per protocole
    bind_rows() %>% # bind all df
    rowid_to_column("rowid") # add rowid



  # Simply modify the id if already stored in the R6 object
  # Such as former and new VP don't have the same IDs
  if(nrow(self$poolVP) > 0 ) poolVP <- poolVP %>%
    mutate(id = id + max(self$poolVP$id),
           rowid = rowid + max(self$poolVP$rowid))





  # Extract the cmts or interest (ytype / observation type)
  self$targets %>%
    filter(protocol %in% unique(poolVP$protocol)) %>%
    pull(cmt) %>%
    unique -> cmts

  # For each cmts, add a column "Below upper" or "above lower"
  # In fact useless if not use green filter but fast anyway
  col_to_add <- c(paste0(cmts, "_BU"),
                  paste0(cmts, "_AL"))

  for(a in col_to_add) poolVP[a] <- NA



  ### Some protocol can have less obsevation than other
  ### Ex:  PK and PD for dose group, but only PD for control
  ### We need to tell them in such case to note perform PK for control
  ### Use crossing but speed anyway

  if(time_compteur == T) t0 <- Sys.time()


  crossing(protocols = unique(poolVP$protocol), cmt =  self$targets %>%
             filter(protocol %in% unique(poolVP$protocol)) %>%
             pull(cmt) %>%
             unique) %>%
    full_join(self$targets) %>%
    filter(is.na(time)) -> torem


  if(nrow(torem)>0){

    for(a in 1:nrow(torem)){

      cmt_to_rm <- torem$cmt[[a]]
      pro <- torem$protocol[[a]]

      poolVP[[paste0(cmt_to_rm,"_BU")]][poolVP$protocol == pro] <- -1
      poolVP[[paste0(cmt_to_rm,"_AL")]][poolVP$protocol == pro] <- -1
    }

  }

  if(time_compteur == T) timesaver$cleanProto <- difftime(Sys.time(),t0, units = "s")


  # Handling death and survival agents, to greatly accelerate the process
  all_param <- self$param


  # Compute the filters now and once for all and not during the loop
  # Do it for each
  filters <-  list()

  filters <- map(unique(self$targets$cmt), ~   self$make_filters(cmt = .x)   %>%
                   gsub(pattern = "line\\$", replacement = "" ) )
  names(filters) <- unique(self$targets$cmt)

  for(a in unique(self$targets$cmt)){

    filters[[a]] <- self$make_filters(cmt = a)   %>%
      gsub(pattern = "line\\$", replacement = "" )

  }


  # line_compar <- paste0("which(", line_compar,")")
  # saveRDS(object = poolVP, file = gsub("\\.RDS", "_todetermine.RDS", file))


  # Time compteur

  if(time_compteur == T){

    timesaver$poolVP_compteur <- tibble(n = NA) %>% slice(0)

  }

  n_compteur <- 0
  # Create a table saving the siml, preparing a final join
  siml <-  list() #tibble(id = integer(), protocol = character())
  # just in case we never enter into the loop (if already filled, almost always useless)

  ntotal <- nrow(poolVP)

  maxinfo <- nrow(poolVP)

  pb <- progress_bar$new(
    format = "  VP creation [:bar] :current/:total (:percent) in :elapsed",
    total = maxinfo, clear = FALSE, width= 60)



  # Preparing some object to store elements (filters etc)

  slice0 <- poolVP   %>% slice(0) %>% select(!!!parse_exprs(all_param[all_param != "protocol"])) # filre_neg_above

  neg_below <- neg_above <-  pos_below <- pos_above <- list() # filtre_neg_below


  for(a in self$targets$cmt){

    pos_below[[a]] <- slice0
    pos_above[[a]] <- slice0

  }

  # Add Inf handling
  nrowcolInf <- 0

  maxinf <- self$targets %>% filter(max == Inf)
  if(nrow(maxinf)>0){

    for(a in 1:nrow(maxinf)){

      line <- maxinf %>%
        slice(a)

      idsInf <- which(poolVP$protocol == line$protocol)
      poolVP[[paste0(line$cmt, "_BU")]][idsInf] <- Inf

      nrowcolInf <- nrowcolInf + length(idsInf)

    }


  }

  maxminf <- self$targets %>% filter(min == -Inf)
  if(nrow(maxminf)>0){

    for(a in 1:nrow(maxminf)){

      line <- maxminf %>%
        slice(a)

      idsInf <- which(poolVP$protocol == line$protocol)
      poolVP[[paste0(line$cmt, "_AL")]][idsInf] <- -Inf

      nrowcolInf <- nrowcolInf + length(idsInf)

    }


  }




  message(green("Start main loop"))
  newratio <- is.na(poolVP[, col_to_add]) %>% sum


  if(time_compteur == T){

    timesaver$tprewhile <- difftime(Sys.time(), tTOTAL, units = "s")
  }

  VP_rej <- tibble()

  nextrapoGreen0 <- 0


  if(nrow( PreviousResults) > 0){

    PreviousResults <- PreviousResults %>%
      select(-starts_with("id")) %>%
      left_join(poolVP)

  }
  PrevLoop <- F


  use_red_filter <- tibble(protocol = unique(self$targets$protocol), use = use_red_filter)
  use_green_filter <- tibble(protocol = unique(self$targets$protocol), use = use_green_filter)
  # begining while lopp----------------------------------------------------------
  while(newratio %>% sum > 0){

    # Just compute some stat...
    t0 <- Sys.time()

    if(time_compteur == T) poolVP_compteur_new <- tibble(n = n_compteur)

    n_compteur <- n_compteur + 1
    pb$update(ratio = (maxinfo -newratio / length(col_to_add))/maxinfo)


    if(nrow(PreviousResults) > 0){ # First use the previous results if provided

      PrevLoop <- T

     proto   <- sample(PreviousResults$protocol %>% unique, 1)

      res <-  PreviousResults %>%
        filter(protocol == proto)

      res <- res %>%
        slice(1:min(2000, nrow(PreviousResults)))


      line <- res %>%
        mutate(id_origin = id)


      if(time_compteur == T){
        poolVP_compteur_new$timemodel <- NA
        poolVP_compteur_new$nsimul <- NA
      }

    }else{ # then complete the rest

      if( PrevLoop == T) cat(red("\nEnd of using previous results\n"))

      PrevLoop <- F



    # set.seed(89651)
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
      rename(id_origin = id) %>%
      rowid_to_column("id")



    if(time_compteur == T) poolVP_compteur_new$timesampleline <- difftime(Sys.time(), t0, units = "s")
    # Now we need to handle the administrations
    # by making a temporar copy
    # protocol  <- self$protocols[[line$protocol[[1]]]]

    protocol <-  line %>%
      mutate(protocol2 = map(protocol, ~ self$protocols[[.x]])) %>%
      select(id, protocol2) %>%
      unnest(protocol2)


    if(time_compteur == T) poolVP_compteur_new$protocoldf <- difftime(Sys.time(), t0, units = "s") - poolVP_compteur_new$timesampleline
    # add_events_line$amt[is.na(add_events_line$amt )] <- 0

    # And now we can make the simulation and extract the result !
    time_simulations <- Sys.time()

    b <- Sys.time()

    res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                        icmt = self$initial_cmt_values, time_vec =self$times,
                        pardf = self$parameters_default_values, model = self$model) %>%
      as_tibble() # %>%
    # get back the real id
    # left_join(line %>% distinct(id, id_origin), by = "id") %>%
    # select(-id) %>%
    # rename(id = id_origin)


    if(!"id" %in% names(res)) res <- res %>% mutate(id = 1)

    if(!is.null(keep)){

      tokeep <- c(self$targets$cmt, keep) %>% unique()
      res <- res %>%
        select(id, time,!!!parse_exprs(tokeep))
    }

    time_simulations <-  difftime(Sys.time(), b, units = "s")


    n_simulations <-  nrow(line)
    if(time_compteur == T){
      poolVP_compteur_new$timemodel <- time_simulations
      poolVP_compteur_new$nsimul <- n_simulations
    }



    }  # end if_else previous results vs new ode system


    if(time_compteur == T) t02 <- Sys.time()
    # # get targets for this patients
    targets_temp <- self$targets %>%
      filter(protocol == unique(line$protocol))#%>%



    # # Compute the tests
    res2 <- res %>%
      select(id, time,  unique(targets_temp$cmt)) %>%
      gather("cmt", "value", unique(targets_temp$cmt)) %>%
      filter(time %in% targets_temp$time) %>%
      left_join( targets_temp, by = c("time", "cmt")) %>%
      filter(!is.na(min)) %>%
      mutate(be_up = value <= max) %>%
      mutate(ab_low = value >= min) %>%
      left_join(line %>% distinct(id, id_origin), by = c("id"))


    if(time_compteur == T){

      poolVP_compteur_new$res2 <- difftime(Sys.time(), t02)
      t02 <- Sys.time()


    }

    ## Addition 04/05/22

    list_error_computation <- list()

    analytical_problem <-  res2 %>% filter(value < 0)

    if(nrow(analytical_problem) >0){

      list_error_computation[[n_compteur]] <- analytical_problem

      analytical_problem %>%
        pull(id_origin) -> toLabelErrorAna

      poolVP <- poolVP %>%
        filter(! id %in% toLabelErrorAna)

      # poolVP_id <- poolVP_id %>%
        # filter(! id %in% toLabelErrorAna)

    }






    # Red filter --------------------------------------------------------------
    # If no red filter


    if(use_red_filter$use[use_red_filter$protocol == unique(line$protocol)] == F & PrevLoop == F){


      t02 <- Sys.time()
      res2 %>%
        # left_join(line %>% distinct(id, id), by  = "id") %>%
        filter(be_up == F | ab_low == F) %>%
        pull(id_origin) %>%
        unique -> idstorem

      if(keepRedFiltaftDis == T){


        newfilters <- res2 %>%
          filter(be_up == F) %>%
          distinct(id, cmt) %>%
          left_join(line, by  = "id") %>%
          distinct(!!!parse_exprs(all_param), cmt,id)


        self$filters_neg_above <- bind_rows(self$filters_neg_above,
                                            newfilters )

        newfilters <- res2 %>%
          filter( ab_low == F) %>%
          left_join(line, by  = "id") %>%
          # filter( id %in% idstorem) %>%
          distinct(!!!parse_exprs(all_param), cmt,id)

        self$filters_neg_below <- bind_rows(self$filters_neg_below, newfilters )

      }


      if(saveVPRej == T){


        VP_rej <- bind_rows(VP_rej, res2 %>%
                              filter(id_origin %in% idstorem) %>%
                              mutate(n = n_compteur)
        )
      }

      poolVP <- poolVP %>%
        filter(! id %in% idstorem)

      # poolVP_id <- poolVP_id %>%
        # filter(! id %in% idstorem)

      if(time_compteur == T){

        poolVP_compteur_new$ElimRedUnac <-  difftime(Sys.time(), t02, units = "s")

      }


    }else{

      # if red filter
      # if(time_compteur == T) poolVP_compteur_new[paste0("res2_", cmtt)] <- difftime(Sys.time(), t02, units = "s")


      ##### Remove patients neg  and compute red filters if activated ########
      # lines output neg_ below


      t0filternegbelow <- Sys.time()

      res2 %>%
        filter(ab_low == 0) %>%
        group_by(id_origin, cmt) %>%
        slice(1) %>% ungroup() -> idsbelow

      filters_neg_below <- idsbelow %>%
        left_join(line,by = c("id", "protocol", "id_origin")) %>%
        select(!!!parse_exprs(all_param), "cmt","id", "id_origin")

      if(nrow(filters_neg_below) > 0){

        filters_neg_below_up_reduc <- filter_reduc(df = filters_neg_below ,obj = self, direction = "below")


      }else{

        filters_neg_below_up_reduc <- filters_neg_below
      }

      timefilternegbelowmake <- difftime(Sys.time(),t0filternegbelow , units = "s")

      if(time_compteur == T){
        poolVP_compteur_new$Treduc_filter_neg_below <- timefilternegbelowmake
        poolVP_compteur_new$nfilter_negbel_bef <- nrow(filters_neg_below)
        poolVP_compteur_new$nfilters_negbel_af <- nrow(filters_neg_below_up_reduc)
      }
      #



      # lines output neg_ above


      t0filternegabovemake <- Sys.time()

      res2 %>%
        filter(be_up == 0) %>%
        group_by(id, cmt) %>%
        slice(1) %>% ungroup() ->idsabove


      filters_neg_above <-  idsabove %>%
        left_join(line,by = c("id", "protocol", "id_origin")) %>%
        select(!!!parse_exprs(all_param), "cmt" , "id", "id_origin")

      if(nrow(filters_neg_above) > 0){
        filters_neg_above_reduc <- filter_reduc(filters_neg_above ,obj = self, direction =  "above")

      }else{

        filters_neg_above_reduc <- filters_neg_above
      }

      timefilternegabovemake <- difftime(Sys.time(),t0filternegabovemake , units = "s")


      if(time_compteur == T){
        poolVP_compteur_new$Treduc_filter_neg_above <- timefilternegabovemake
        poolVP_compteur_new$nfilter_negab_bef <- nrow(filters_neg_above)
        poolVP_compteur_new$nfilters_negab_af <- nrow(filters_neg_above_reduc)
      }



      if(time_compteur == T) t02 <- Sys.time()

      ## SaveRej if needed

      if(saveVPRej == T){


        VP_rej <- bind_rows(VP_rej, res2 %>%
          filter(id_origin %in% c(unique(filters_neg_above$id_origin), unique(filters_neg_below$id_origin))) %>%
            mutate(n = n_compteur)
        )

      }

      ## We remove all the one not accepted
      tfilterODE <- Sys.time()

      idremoved <-  c(unique(filters_neg_above$id_origin), unique(filters_neg_below$id_origin))
      poolVP <- poolVP %>%
        filter(! id %in%idremoved)

if(PrevLoop == T) PreviousResults <- PreviousResults %>%
        filter(id %in% idremoved & protocol == line$protocol [[1]])

      tfilterODE <- difftime( Sys.time(), tfilterODE, units = "s")


      # which(c(unique(filters_neg_above$id_origin, unique(filters_neg_below$id_origin))) %in% poolVP$id[!is.na(poolVP$tumVol_BU)])

      # poolVP_id <- poolVP_id %>%
        # filter(! id %in% c(unique(filters_neg_above$id_origin), unique(filters_neg_below$id_origin)))


      if(time_compteur == T) poolVP_compteur_new$TremoveNegDirect <- difftime(Sys.time(),t02 , units = "s")

      #### Apply red filter if activated

      t02 <- Sys.time()
      nref <- nrow(poolVP)
      if(time_compteur == T){
        t02 <- t02
        nref <- nref
      }


      if(nrow(filters_neg_above_reduc) > 0){

        # method 1 with loop
        if(methodFilter == 1){
          for(a in 1:nrow(filters_neg_above_reduc)){

            ref <- filters_neg_above_reduc %>% slice(a)

            poolVP %>%
              mutate(test = !!parse_expr(filters[[ref$cmt]][["above"]])) %>%
              filter(test == T) %>%
              pull(id) ->idtorem

            if(saveVPRej == T)    VP_rej <- bind_rows(VP_rej, tibble(id_origin = idtorem) %>% mutate(n = n_compteur))

            poolVP <- poolVP %>%
              filter(!id %in% idtorem)

            if(PrevLoop == T) PreviousResults <- PreviousResults %>%
              filter(id %in% idtorem & protocol == line$protocol [[1]])
            # poolVP_id <- poolVP_id %>%
              # filter(!id %in% idtorem)
          }
        }else{ # method 2 with biiig filter


          temp <- filters_neg_above_reduc %>%
            mutate(filtre = map_chr(cmt,  ~ filters[[.x]][["above"]]))

          for(a in all_param){

            temp %>%
              mutate(filtre = map2_chr(filtre, !!parse_expr(a), function(filtre, x){

                gsub(paste0("ref\\$", a),x,  filtre)

              })) -> temp

          }

          temp %>%
            pull(filtre) -> filtres

          if(length(filtres) > 0 ){
            filtre_line <- paste0("(", filtres, ")") %>% paste0(collapse = "|")








            poolVP %>%
              mutate(test = !!parse_expr(filtre_line)) %>%
              filter(test == T) %>% pull(id) ->idtorem

            # idtorem[idtorem %in%idmi]
            # if(idtorem[idtorem %in%idsMissings] %>% sum > 0) stop("IDs missing ! ")

            if(saveVPRej == T)    VP_rej <- bind_rows(VP_rej, tibble(id_origin = idtorem) %>% mutate(n = n_compteur))


            poolVP <- poolVP %>%
              filter(! id %in% idtorem)


            if(PrevLoop == T) PreviousResults <- PreviousResults %>%
              filter(id %in% idtorem & protocol == line$protocol [[1]])

            # poolVP_id <- poolVP_id %>%
              # filter(! id %in% idtorem)
          }


        }
      }


      if(time_compteur == T) poolVP_compteur_new$TapplyRedFilterAbove <- difftime(Sys.time(),t02 , units = "s")


      if(nrow(filters_neg_below_up_reduc)> 0){

        # method 1 with loop
        if(methodFilter == 1){
          for(a in 1:nrow(filters_neg_below_up_reduc)){

            ref <- filters_neg_below_up_reduc %>% slice(a)


            poolVP %>%
              mutate(test = !!parse_expr(filters[[ref$cmt]][["below"]])) %>%
              filter(test == T) %>%
              pull(id)  ->idtorem

            if(saveVPRej == T)    VP_rej <- bind_rows(VP_rej, tibble(id_origin = idtorem) %>% mutate(n = n_compteur))

            poolVP <- poolVP %>%
              filter(!id %in% idtorem)


            if(PrevLoop == T) PreviousResults <- PreviousResults %>%
              filter(id %in% idtorem & protocol == line$protocol [[1]])

            # poolVP_id <- poolVP_id %>%
              # filter(!id %in% idtorem)
          }
        }else{ # method 2 with biiig filter


          temp <- filters_neg_below_up_reduc %>%
            mutate(filtre = map_chr(cmt,  ~ filters[[.x]][["below"]]))

          for(a in all_param){

            temp %>%
              mutate(filtre = map2_chr(filtre, !!parse_expr(a), function(filtre, x){

                gsub(paste0("ref\\$", a),x,  filtre)

              })) -> temp

          }

          temp %>%
            pull(filtre) -> filtres

          if(length(filtres) > 0 ){
            filtre_line <- paste0("(", filtres, ")") %>% paste0(collapse = "|")


            poolVP %>% #poolVP_id
              mutate(test = !!parse_expr(filtre_line)) %>%
              filter(test == T) %>% pull(id) ->idtorem





            if(saveVPRej == T)    VP_rej <- bind_rows(VP_rej, tibble(id_origin = idtorem) %>% mutate(n = n_compteur))


            poolVP <- poolVP %>%
              filter(! id %in% idtorem)


            if(PrevLoop == T) PreviousResults <- PreviousResults %>%
              filter(id %in% idtorem & protocol == line$protocol [[1]])

          }
        }

      }

      nrem <- nref -  nrow(poolVP)

      time_filter_neg_apply <- difftime(Sys.time(), t02, units = "s")
      nremoved_below_fil <-  nrow(poolVP)

      if(time_compteur == T){
        poolVP_compteur_new$TapplyRedFilterBoth <- time_filter_neg_apply
        poolVP_compteur_new$TapplyRedFilterBelow <-      poolVP_compteur_new$TapplyRedFilterBoth -  poolVP_compteur_new$TapplyRedFilterAbove
        poolVP_compteur_new$NremovedRedFilter <-   nrem
      }




      ### Update filters

      self$filters_neg_above <- bind_rows(self$filters_neg_above, filters_neg_above_reduc )
      self$filters_neg_below <- bind_rows(self$filters_neg_below, filters_neg_below_up_reduc )

      ##### Compute the rendement of red filter and disable if negative
if(PrevLoop == F){

      totaltimeredfilter <- time_filter_neg_apply + timefilternegabovemake + timefilternegbelowmake
      totalsave <-  nrem * time_simulations /  n_simulations + tfilterODE

      if(time_compteur == T){

        poolVP_compteur_new$TSavedRedFilter <- totalsave
        poolVP_compteur_new$TimeTotalRedFilter <- totaltimeredfilter
      }

      if(totalsave < totaltimeredfilter){

        if(RedFilterDisAllProt == T){

          message(red(paste0("\nRed filter system (all protocols) disabled.")))

          use_red_filter$use <- F
        }else{

          message(red(paste0("\nRed filter system for protocol ", unique(line$protocol) , "  disabled.")))

          use_red_filter$use[use_red_filter$protocol == unique(line$protocol)] <- F

        }


      }

}

    }# end use of red filter


    # gree filter --------------------------------------------------------------
    # dont use the freen filter if use_green_filter F
    # OR if percentage green is below  pctActivGreen

    # cat(green(paste0("use_green_filter is equa to :" , use_green_filter)))
    # cat(green(paste0(nrow(res2 %>% filter(be_up == T & ab_low== T)), "<", pctActivGreen * nrow(res2) )))

    if( (use_green_filter$use[use_red_filter$protocol == unique(line$protocol)] == F | nrow(res2 %>% filter(be_up == T & ab_low== T)) < pctActivGreen * nrow(res2)) & PrevLoop == F  ){


      if(time_compteur == T) t02 <- Sys.time()

      targets_temp %>%
        # filter(protocol == unique(line$protocol), cmt == cmtt) %>%
        pull(time) %>% length -> nmax


      res2 %>%
        filter(be_up == 1 & ab_low == 1) %>%
        group_by(id_origin) %>%
        tally -> nabove

      res2 %>%
        left_join(nabove, by = "id_origin") %>%
        filter(n == nmax) %>%
        group_by(id_origin) %>%
        slice(1) %>% pull(id_origin) -> idsgood

      # line %>% filter(id %in% idsgood) %>% pull(id) -> idsgood

      poolVP[poolVP$id %in% idsgood & poolVP$protocol == unique(line$protocol), col_to_add     ] <- 0


      if(time_compteur == T)  poolVP_compteur_new$time_addgreennofil <- difftime(Sys.time(), t02, units = "s")




      if(time_compteur == T) t02 <- Sys.time()

      if( !is.null(timeSave)) res <- res %>% filter(time %in% timeSave)

      res %>%
        left_join(line %>% distinct(id, id_origin), by = "id") %>%
        filter(id_origin %in% idsgood) %>%
        mutate(protocol = unique(line$protocol)) %>%
        select(-id) %>%
        rename(id = id_origin) %>%
        group_by(id,protocol) %>%

        nest() -> forjoin



      siml[[n_compteur]] <-  forjoin#siml %>%


      if(time_compteur == T) poolVP_compteur_new$timesimlandjoin <- difftime(Sys.time(), t02, units = "s")





    }else if (PrevLoop == F | nrow(PreviousResults) > 0){



      t0greenfilterFull <- Sys.time()

      cmttall <- unique(self$targets$cmt)

      while(length(cmttall) > 0){


        cmtt <- cmttall[[1]]

        t0greenfilter <- Sys.time()

        aboveExpr <- parse_expr(paste0(cmtt, "_AL"))
        belowExpr <- parse_expr(paste0(cmtt, "_BU"))


        if(time_compteur == T){

          befusegreen <-   poolVP %>%
            group_by(!!aboveExpr, !!belowExpr) %>%
            tally
          t0green   <- Sys.time()
        }
        # lines output below_up

        donebeforegreen <-  poolVP %>%
            filter(!is.na(!!belowExpr) & protocol %in% line$protocol[[1]] & !is.na(!!aboveExpr) ) %>%
            nrow()

        donebeforegreen <- donebeforegreen +  res2 %>%
          filter(cmt == cmtt) %>%
          filter(be_up == T &ab_low  == T ) %>%
          nrow()


        if(time_compteur == T) t02 <- Sys.time()

        self$targets %>%
          filter(protocol == unique(line$protocol), cmt == cmtt) %>%
          pull(time) %>% length -> nmax


        res2 %>%
          filter(be_up == 1 & ab_low == 1) %>%
          group_by(id_origin) %>%
          tally -> naccepted

        res2 %>%
          left_join(naccepted, by = "id_origin") %>%
          filter(n == nmax) %>%
          group_by(id_origin) %>%
          slice(1) %>% pull(id_origin) -> idsaccepted




        line %>% filter(id_origin %in% idsaccepted) -> filters_pos_below_up ->  filters_ps_above_lo
        filters_pos_below_up$cmt = cmtt
        filters_ps_above_lo$cmt = cmtt


        # automatic desactivation if Inf or - Inf
        below_green <- T
        above_green <- T
        if( nrow(self$targets %>%
                 filter(!(protocol == unique(res2$protocol) & cmt == cmtt & max == Inf))) == 0) below_green <- F
        if( nrow(self$targets %>%filter(!(protocol == unique(res2$protocol) & cmt == cmtt & min == -Inf))) == 0 )  above_green <- F

        if(nrow(filters_pos_below_up) == 0)  below_green  <- F
        if(nrow(filters_ps_above_lo) == 0)  above_green  <- F


        if(below_green == T){

          if(time_compteur == T) t03 <- Sys.time()


          filters_pos_below_up_reduc <- filter_reduc(filters_pos_below_up,obj = self , direction  =  "below")

          if(time_compteur == T) poolVP_compteur_new[paste0("Treduc_filter_green_below_", cmtt)] <- difftime(Sys.time(), t03, units = "s")

          pos_below[[cmtt]] <- bind_rows(pos_below[[cmtt]], filters_pos_below_up_reduc[names(pos_below[[1]])])

          if(time_compteur == T){
            t02 <- Sys.time()
            reff <- sum(is.na(poolVP[deparse(belowExpr)]))

          }

          if(nrow(filters_pos_below_up_reduc) >0){
            for(a in 1:nrow(filters_pos_below_up_reduc)){

              ref <- filters_pos_below_up_reduc %>% slice(a)

              poolVP %>%
                mutate(test = !!parse_expr(filters[[ref$cmt]][["below"]])) %>%
                filter(test == T) %>% pull(id) -> id_temp

              poolVP[[paste0(cmtt,"_BU")]][poolVP$id %in% id_temp & poolVP$protocol == unique(line$protocol)] <-  ref$rowid

          if(PrevLoop) PreviousResults[[paste0(cmtt,"_BU")]][PreviousResults$id %in% id_temp & PreviousResults$protocol == unique(line$protocol)] <-  ref$rowid
            }

            if(time_compteur == T){
              poolVP_compteur_new[paste0("TapplyGreenFilterBelow_", cmtt)] <- difftime(Sys.time(), t02, units = "s")
              poolVP_compteur_new[paste0("n_add_below_fil", cmtt)] <-   sum(is.na(poolVP[deparse(belowExpr)])) - reff
            }
          }



        } # end if below green

        if(above_green == T){


          if(time_compteur == T) t03 <- Sys.time()

          filters_ps_above_lo_reduc <- filter_reduc(filters_ps_above_lo,obj = self, direction = "above")

          if(time_compteur == T) poolVP_compteur_new[paste0("Treduc_filter_green_above_", cmtt)] <- difftime(Sys.time(), t03, units = "s")

          pos_above[[cmtt]] <- bind_rows(pos_above[[cmtt]], filters_ps_above_lo_reduc[names(pos_above[[1]])])


          if(time_compteur == T){
            t02 <- Sys.time()
            reff <- sum(is.na(poolVP[[deparse(aboveExpr)]]))


          }


          if(nrow(filters_ps_above_lo_reduc) >0){

            for(a in 1:nrow(filters_ps_above_lo_reduc)){

              ref <- filters_ps_above_lo_reduc %>% slice(a)

              poolVP %>%
                mutate(test = !!parse_expr(filters[[ref$cmt]][["above"]])) %>%
                filter(test == T) %>% pull(id) -> id_temp


              poolVP[[paste0(cmtt,"_AL")]][poolVP$id %in% id_temp & poolVP$protocol == unique(line$protocol)] <- ref$rowid

              if(PrevLoop) PreviousResults[[paste0(cmtt,"_AL")]][PreviousResults$id %in% id_temp & PreviousResults$protocol == unique(line$protocol)] <-  ref$rowid

            }

          }
          if(time_compteur == T){
            poolVP_compteur_new[[paste0("TapplyGreenFilterAbove_", cmtt)]] <- difftime(Sys.time(), t02, units = "s")
            poolVP_compteur_new[[paste0("n_add_above_fil_", cmtt)]] <-   sum(is.na(poolVP[[deparse(aboveExpr)]])) - reff
          }



        } # end if above green


        if(time_compteur == T){ # TODO adapt it with recent modification (counting effi)

          afusegreen <-   poolVP %>%
            group_by(!!aboveExpr, !!belowExpr) %>%
            tally  %>%
            rename(n2 = n)

          left_join(afusegreen, befusegreen) %>%
            mutate(Diff = n2 - n) %>%
            mutate(colname = paste0(deparse(belowExpr), "-", deparse(aboveExpr))) -> temp3

          if(nrow(temp3) >=1){
          for(a in 1:nrow(temp3)){

            poolVP_compteur_new[temp3$colname[[a]]] <-  temp3$Diff[[a]]

          }
          poolVP_compteur_new$wholegreenfilter <- difftime(Sys.time(), t0green, units = "s")
        }
        }


        if(time_compteur == T) t02 <- Sys.time()

        if( !is.null(timeSave)) res <- res %>% filter(time %in% timeSave)

        res %>%
          left_join(line %>% distinct(id, id_origin), by = "id") %>%
          filter(id_origin %in% unique(poolVP$id) & !id_origin %in% siml$id) %>%
          mutate(protocol = unique(line$protocol)) %>%
          select(-id) %>% rename(id = id_origin) %>%
          group_by(id,protocol) %>%

          nest() -> forjoin



        siml[[n_compteur]] <-  forjoin#sim


        if(time_compteur == T) poolVP_compteur_new$timesimlandjoin <- difftime(Sys.time(), t02, units = "s")


        # Count how many VPs are already good
     doneaftergreen <-  poolVP %>%
          filter(!is.na(!!belowExpr) & protocol %in% line$protocol  & !is.na(!!aboveExpr) ) %>%
          nrow()

     nextrapoGreen <- doneaftergreen - donebeforegreen

     timegreen <- difftime(Sys.time(), t0greenfilter, units = "s")




     if(PrevLoop == F){

       totalsaveGreen <-  nextrapoGreen * time_simulations /  n_simulations

       # cat(green(paste0("Ratio Green filter:",totalsaveGreen - timegreen , "\n")))
       if(totalsaveGreen < timegreen){


         if(GreenFilterDisAllProt == T){

           message(green(paste0("\nGreen filter system (all protocols) disabled.")))

           use_green_filter$use <- F
         }else{

           message(green(paste0("\nGreen filter system for protocol ", unique(line$protocol) , "  disabled.")))

           use_green_filter$use[use_green_filter$protocol == unique(line$protocol)] <- F

         }

         cmttall <- character()

       }else{

         cmttall <- cmttall[-1]

       }


       if(time_compteur == T){

         # poolVP_compteur_new$TSavedGreenFilter <- totalsaveGreen
         # poolVP_compteur_new$TimeTotalGreenFilter <- timegreen
         # poolVP_compteur_new$nextrapoGreen <- thisTurn
       }
     }

     if(PrevLoop == T) PreviousResults <- PreviousResults %>%   filter(!!parse_expr(filter_to_use))


        # end for each OoI
      }




      if(time_compteur == T){
        poolVP_compteur_new$TimeTotalGreenFilter <- difftime(Sys.time(), t0greenfilterFull, units = "s")
      }



    }# end use green filter



    # if() siml <- tibble(id = integer(), protocoles = character(), simul= list())

    if(time_compteur == T)  timesaver$poolVP_compteur <- bind_rows(timesaver$poolVP_compteur,poolVP_compteur_new)

    # Compute newratio for knowing when to stop
    newratio <- is.na(poolVP[, col_to_add]) %>% sum
  }# fin while 1

  cat(red("End while loop"))
  ## Need to remove filter pos for which some have been delated (because of multiple protocol,one ok, one bad)

  # pos_below <- map(pos_below, function(x) x %>% filter(id %in% poolVP$id) %>% select(-id))
  # pos_above <- map(pos_above, function(x) x %>% filter(id %in% poolVP$id) %>% select(-id))


  if(time_compteur == T)  t02 <- Sys.time()

  siml <- bind_rows(siml)

  poolVP <- poolVP %>%
    select(-starts_with("simul")) %>%
    left_join(siml  %>%  rename(simul = data))

  if(time_compteur == T)  timesaver$joinsimul <- difftime(Sys.time(), t02, units = "s")

  if(length(list_error_computation) > 0) self$errorComputation <-  bind_rows(list_error_computation)



  # Update filters
  if(time_compteur == T)  t02 <- Sys.time()
  for(a in unique(self$targets$cmt)){

    if(is.null(self$filters_pos_above[[a]])) self$filters_pos_above[[a]] = tibble()
    self$filters_pos_above[[a]] <- bind_rows(self$filters_pos_above[[a]] %>% mutate(PrimFilter = T) , pos_above[[a]] )

    if(is.null(self$filters_pos_below[[a]])) self$filters_pos_below[[a]] = tibble()
    self$filters_pos_below[[a]] <- bind_rows(self$filters_pos_below[[a]] %>% mutate(PrimFilter = T) , pos_below[[a]] )

  }


  if(time_compteur == T)  timesaver$updatefilters <- difftime(Sys.time(), t02, units = "s")


  if(time_compteur == T)   self$timesaver <- timesaver
  # bind_rows(self$poolVP, poolVP )


  # save pooVP
  if(is.null(self$poolVP)){

    self$poolVP <- poolVP

  }else{
    self$poolVP <- bind_rows(self$poolVP, poolVP )
  }

  # Save Rej if needed
  if(saveVPRej == T)   self$VP_rejected <-  VP_rej %>% left_join(

    VP_df %>%
      rowid_to_column("id_origin"), by = "id_origin"

  )



# Final total time saving
  if(time_compteur == T){
    timesaver$tTOTAL <- difftime(Sys.time(), tTOTAL, units = 's')

    self$timeTrack <- timesaver
  }else{


    self$timeTrack <- list(ttotal = difftime(Sys.time(), tTOTAL))
  }





  # Fill missing profile if required
  if(fillatend){

    print("filling")
    self$fill_simul()
  }

  # Filter reduce if asked
  if(reducefilteratend){

    print("filter reduction")
    self$n_filter_reduc()
  }

  self
  # # Recompute the whole poolVP

})


