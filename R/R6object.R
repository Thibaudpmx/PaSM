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
  # targets = NULL,
  filters_neg_above = NULL,
  filters_neg_below = NULL,
  filters_pos_above = NULL,
  filters_pos_below = NULL,
  data = NULL,
  parameters_default_values = NULL,# c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
  initial_cmt_values = NULL, #c(X1 = 50) # initial compartment values. At least one, and every missing cmt name would be set to 0
  times = NULL, #seq(0,52, 1)
  poolVP = NULL,
  protocols = NULL,
  param_reduce = NULL,
  param_increase = NULL,
  param_no_impact = NULL,
  targets = NULL,
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
print(data_segment_plot(data = self$data, targets = targets))

    # update targets
    self$targets <- targets

    # initialise filters

    filters_list= list()

    for(a in unique(targets$cmt))  filters_list[[a]] <- tibble()


    self$filters_neg_above[[a]]  <- tibble()
    self$filters_neg_below[[a]]  <- tibble()
    self$filters_pos_above[[a]]  <- tibble()
    self$filters_pos_below[[a]]  <- tibble()



    return(  message("Once any VP's or filters have been created, targets will be locked."))

})


# VP_production -----------------------------------------------------------


VP_proj_creator$set("public", "add_VP", function(VP_df,  saven = 50, drug = NULL, update_at_end = T, time_compteur = F,  fillatend = F, reducefilteratend = F){

  # protocols = self$protocols

  toadd <- VP_df

  poolVP <-     VP_df %>%
    rowid_to_column("cellid") %>%
    crossing(protocol = unique(self$targets$protocol)) %>%
    rowid_to_column("rowid")





  #
#   poolVP %>%
#     filter(round(k2,1) == 1.2 & round(lambda0,2) == 0.07)

  # Add the columns for each output
  self$targets %>%
    filter(protocol %in% unique(poolVP$protocol)) %>%
    pull(cmt) %>%
    unique -> col_to_add

  col_to_add <- c(paste0(col_to_add, "_BU"),
                  paste0(col_to_add, "_AL"))

  for(a in col_to_add) poolVP[a] <- NA


  ## Apply the filters

  for(cmt in "tumVol"){

    filter_template <- self$make_filters(cmt)%>%
      gsub(pattern = "line\\$", replacement = "")

    # remove neg_above
    temp <-  self$filters_neg_above[[cmt]]
    filter_temp <- paste("!(", filter_template[[1]], ")")
    if(nrow(temp) > 0){
    for(a in 1:nrow(temp)){

      ref <-temp %>% slice(a)

      poolVP <- poolVP %>%
        filter(!!parse_expr(filter_temp))

    }
    }
    # remove neg_below
    temp <-  self$filters_neg_below[[cmt]]
    filter_temp <- paste("!(", filter_template[[2]], ")")
    if(nrow(temp) > 0){
    for(a in 1:nrow(temp)){

      ref <-temp %>% slice(a)

      poolVP <- poolVP %>%
        filter(!!parse_expr(filter_temp))

    }
}
    poolVP <- poolVP %>%
      select(-rowid) %>%
      rowid_to_column()

    # set pos abov
    temp <-  self$filters_pos_above[[cmt]]
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



 if(!is.null(self$poolVP$cellid)) poolVP <- poolVP %>%
    mutate(cellid = cellid + max(self$poolVP$cellid),
             rowid = rowid + max(self$poolVP$rowid))



  ### eviter les loupes infinis si un protocole n'as pas d'observion..
  crossing(protocols = unique(toadd$protocol), cmt =  self$targets %>%
             filter(protocol %in% unique(poolVP$protocol)) %>%
             pull(cmt) %>%
             unique) %>%
    full_join(self$targets) %>%
    filter(is.na(time)) -> torem

  if(nrow(torem)>0){

    for(a in 1:nrow(torem)){

      cmt_to_rm <- torem$cmt[[a]]
      pro <- torem$protocol[[a]]

      poolVP[[paste0(cmt_to_rm,"_BU")]][poolVP$protocol == pro] <- FALSE
      poolVP[[paste0(cmt_to_rm,"_AL")]][poolVP$protocol == pro] <- FALSE
    }

  }

  # test filterDf system
  #
  #   filterdf <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/saveFilterDf.RDS")
  #
  # poolVP$res <- NULL
  # for(a in 154:nrow(filterdf)){
  #   b <- Sys.time()
  #
  # filter <- filterdf$filter[[a]]
  # result <- filterdf$res[[a]]
  #
  # filter <- gsub("(line\\$BAK0)|(line\\$BAXc0)", "1000", filter) %>%
  #   gsub(pattern = " \\& poolVP\\$group== line\\$group", replacement = "")
  #
  # nrow <- eval(parse_expr(filter))
  # poolVP[nrow, "res"] <- result
  #
  # print(paste0(a, ":", length(nrow),":", difftime(Sys.time(), b)))
  #
  #
  # }


  #

  # Handling death and survival agents, to greatly accelerate the process
  all_param <- names(poolVP)
  all_param <- all_param[! all_param %in% c("cellid", "rowid",col_to_add)]


  line_compar <-   paste0("poolVP$", all_param, " == line$", all_param) %>%
    paste0(collapse = " & ")
  # line_compar <- paste0("which(", line_compar,")")
  # saveRDS(object = poolVP, file = gsub("\\.RDS", "_todetermine.RDS", file))


  # Time compteur

  if(time_compteur == T){

    poolVP_compteur <- tibble(n = NA, time  = NA, nelim =NA, ninfo =NA_real_, computmodel = NA)
    n_compteur <- 0
  }



  siml <- tibble(cellid = integer(), protocol = character(), simul= list())
  # just in case we never enter into the loop (if already filled, almost always useless)

  ntotal <- nrow(poolVP)
  t00 <- Sys.time()





  # poolVPprev <-poolVP
  # poolVP <- poolVPprev


  maxinfo <- is.na(poolVP[, col_to_add]) %>% sum

  pb <- progress_bar$new(
    format = "  VP creation [:bar] :current/:total (:percent) in :elapsed",
    total = maxinfo, clear = FALSE, width= 60)
  # pb <- progress_bar$new(total = )



  slice0 <- poolVP   %>% slice(0) %>% select(!!!parse_exprs(all_param[all_param != "protocol"])) # filre_neg_above

  neg_below <- neg_above <-  pos_below <- pos_above <- list() # filtre_neg_below


  for(a in targets$cmt){

    neg_below[[a]] <- slice0
    neg_above[[a]] <- slice0
    pos_below[[a]] <- slice0
    pos_above[[a]] <- slice0

    }

    # begining while lopp----------------------------------------------------------
    while(is.na(poolVP[, col_to_add]) %>% sum > 0){

      newratio <- is.na(poolVP[, col_to_add]) %>% sum
      pb$update(ratio = (maxinfo -newratio )/maxinfo)

      if(time_compteur == T){

        n_compteur <- n_compteur + 1

        poolVP_compteur_new <- poolVP_compteur %>%
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
        filter(!!parse_expr(filter_to_use)) %>%
        sample_n(1)


      # Now we need to handle the administrations
      # by making a temporar copy
      protocol  <- self$protocols[[line$protocol]]

      # add_events_line$amt[is.na(add_events_line$amt )] <- 0

      # And now we can make the simulation and extract the result !

      b <- Sys.time()
      res <- simulations(ind_param = line, add_events = protocol, returnSim = T,
                         icmt = self$initial_cmt_values, time_vec =self$times,
                         pardf = self$parameters_default_values, model = self$model);res

      if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
        mutate(computmodel = as.double(difftime(Sys.time(),b, units = "sec")))


      # get targets for this patients
      targets_temp <- self$targets %>%
        filter(protocol %in% line$protocol)

      # Now let's see if we can extrapolate some other results
      # cmtt <- "tumVol"
      # cmtt <- "Conc"
      cmt_to_update <- unique(targets_temp$cmt)

      line %>%
        gather("key", "value") %>%
        filter(is.na(value)) -> currentlyna

      cmt_to_update <- cmt_to_update[cmt_to_update %in% gsub("(_AL$)|(_BU$)", "", currentlyna$key)]

      remv <- F
      cellidtorem <- double()
      for(cmtt in cmt_to_update){

        targets_temp2 <- targets_temp %>%
          filter(cmt == cmtt)
        # below upper?
        res %>%
          filter(time %in% targets_temp2$time) %>%
          pull(!!parse_expr(cmtt)) -> values

        be_up <- if_else(min(values <= targets_temp2$max) == 0, F, T);be_up
        ab_low <-  if_else(min(values >= targets_temp2$min) == 0, F, T);ab_low

        pa_in_temp <- self$param_increase[[cmtt]]
        pa_in_temp <- pa_in_temp[pa_in_temp %in% all_param]

        pa_re_temp <-  self$param_reduce[[cmtt]]
        pa_re_temp <- pa_re_temp[pa_re_temp %in% all_param]

        pa_ni_temp <-  self$param_no_impact[[cmtt]]
        pa_ni_temp <- pa_ni_temp[pa_ni_temp %in% all_param]


        # if the line output are above
        if(be_up == F){

          # create a copy of the line_compar with everything "=="
          reject <- paste0("which(", line_compar,")")

          # Then replace "==" by "<=" for survival parameters
          for(a in pa_in_temp){

            reject <- gsub(paste0(a, " *=="), paste0(a, " >= "), reject)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


          }

          # Then replace "==" by ">=" for survival parameters
          for(a in pa_re_temp){

            reject <- gsub(paste0(a, " *=="), paste0(a, " <= "), reject)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


          }

          for(a in pa_ni_temp){

            reject <-
              gsub(paste0("&? * poolVP\\$", a, " *== *line\\$",a),"", reject)
          }

          reject <- gsub("line\\$protocol",  paste0("\"", line$protocol,"\""), reject)

          # Compute the test
          reject_eval <- eval(parse_expr(reject))
          cellidtorem_above <- poolVP[reject_eval, "cellid"]$cellid

          cellidtorem <- c(cellidtorem, cellidtorem_above)
          remv <- T
          # print(paste0(length(cellidtorem), " cells removed"))

          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate(nelim = length(cellidtorem))

          neg_above[[cmtt]] <- bind_rows(neg_above[[cmtt]], line %>% select(!!!parse_exprs(all_param[all_param != "protocol"])))
          ### if the line output is survival
        }



        if(ab_low == F){



          # create a copy of the line_compar with everything "=="
          # reject <- paste0(line_compar, "& is.na(poolVP$",    paste0(cmtt, "_AL"),")")
          reject <- paste0("which(", line_compar,")")

          # Then replace "==" by "<=" for death parameters
          for(a in pa_re_temp){


            reject <- gsub(paste0(a, " *=="), paste0(a, " >= "), reject)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


          }

          # Then replace "==" by ">=" for survival parameters
          for(a in pa_in_temp){

            reject <- gsub(paste0(a, " *=="), paste0(a, " <= "), reject)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


          }



          for(a in pa_ni_temp){

            reject <-
              gsub(paste0("&? * poolVP\\$", a, " *== *line\\$",a),"", reject)
          }

          reject <- gsub("line\\$protocol",  paste0("\"", line$protocol,"\""), reject)

          # Compute the test
          reject_eval <- eval(parse_expr(reject))
          cellidtorem_below <- poolVP[reject_eval, "cellid"]$cellid

          cellidtorem <- c(cellidtorem, cellidtorem_below)
          remv <- T
          # print(paste0(length(cellidtorem), " cells removed"))

          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate(nelim = length(cellidtorem))


          neg_below[[cmtt]] <- bind_rows(neg_below[[cmtt]] , line %>% select(!!!parse_exprs(all_param[all_param != "protocol"])))

        } # end if-else








        # Now update the lines
        if(ab_low == TRUE &   paste0(cmtt, "_AL") %in% currentlyna$key){

          # create a copy of the line_compar with everything "=="
          test_above_lower_lim <- paste0(line_compar, "& is.na(poolVP$",    paste0(cmtt, "_AL"),")")
          test_above_lower_lim <- paste0("which(", test_above_lower_lim,")")

          for(a in pa_re_temp){

            test_above_lower_lim <- gsub(paste0(a, " *=="), paste0(a, " <= "), test_above_lower_lim)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])
          }

          for(a in pa_in_temp){

            test_above_lower_lim <- gsub(paste0(a, " *=="), paste0(a, " >= "), test_above_lower_lim)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])

          }

          for(a in pa_ni_temp){

            test_above_lower_lim  <-
              gsub(paste0("&? * poolVP\\$", a, " *== *line\\$",a),"", test_above_lower_lim)
          }



          whichaboveloweer <- eval(parse_expr(test_above_lower_lim))

          poolVP[[paste0(cmtt, "_AL")]][whichaboveloweer]  <- TRUE


          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate( ninfo = if_else(is.na(poolVP_compteur_new$ninfo),0,poolVP_compteur_new$ninfo)+  length(whichaboveloweer))


          pos_above[[cmtt]] <- bind_rows(pos_above[[cmtt]], line %>% select(!!!parse_exprs(all_param[all_param != "protocol"])))
        } # end if ab_low == T



        if(be_up == TRUE &  paste0(cmtt, "_BU") %in% currentlyna$key){

          # create a copy of the line_compar with everything "=="
          test_below_upper_lim <- paste0(line_compar, "& is.na(poolVP$",    paste0(cmtt, "_BU"),")")
          test_below_upper_lim <- paste0("which(", test_below_upper_lim,")")

          # Then replace "==" by "<=" for death parameters
          for(a in pa_re_temp){


            test_below_upper_lim <- gsub(paste0(a, " *=="), paste0(a, " >= "), test_below_upper_lim)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


          }

          # Then replace "==" by ">=" for survival parameters
          for(a in pa_in_temp){

            test_below_upper_lim <- gsub(paste0(a, " *=="), paste0(a, " <= "), test_below_upper_lim)%>%
              gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


          }

          for(a in pa_ni_temp){

            test_below_upper_lim  <-
              gsub(paste0("&? * poolVP\\$", a, " *== *line\\$",a),"", test_below_upper_lim)
          }



          whichbelowupper <- eval(parse_expr(test_below_upper_lim))



          poolVP[[paste0(cmtt, "_BU")]][whichbelowupper]  <- TRUE

          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate( ninfo = if_else(is.na(poolVP_compteur_new$ninfo),0,poolVP_compteur_new$ninfo)+  length(whichbelowupper))

          pos_below[[cmtt]]  <- bind_rows(pos_below[[cmtt]], line %>% select(!!!parse_exprs(all_param[all_param != "protocol"])))

        } # end if be_up == T

      }# end for each compartment




      if(remv == F ){

        siml <- siml %>%
          add_row(cellid = line$cellid, protocol = line$protocol, simul = list(res))

      }else{

        poolVP <- poolVP %>%
          filter(!cellid %in%cellidtorem)

      }
      # if() siml <- tibble(cellid = integer(), protocoles = character(), simul= list())

      if(time_compteur == T)  poolVP_compteur <- bind_rows(poolVP_compteur,poolVP_compteur_new %>%
                                                             mutate(time = as.double(difftime(Sys.time(), t0, units = "sec"))))



      # print(nn)
    }# fin while 1

    # Just print some stuff
    # nnewlines <- sum(!is.na(poolVP$res)) - before
    # print( paste0(nnewlines,  " new lines proceeded"))
    # print(Sys.time() - t0)
    # print(Sys.time() - t00)



  ## Here happen after nsave iteration
  # print("########################### SAVNG RDS #############################")

  if(time_compteur == T)  poolVP_compteur <<- poolVP_compteur
  print(Sys.time() - t00)


   poolVP <- poolVP %>%
    left_join(siml)



# Update filters

for(a in unique(self$targets$cmt)){

  self$filters_neg_above[[a]] <- bind_rows(self$filters_neg_above[[a]] %>% mutate(PrimFilter = T) , neg_above[[a]] )
  self$filters_neg_below[[a]] <- bind_rows(self$filters_neg_below[[a]] %>% mutate(PrimFilter = T) , neg_below[[a]] )
  self$filters_pos_above[[a]] <- bind_rows(self$filters_pos_above[[a]] %>% mutate(PrimFilter = T) , pos_above[[a]] )
  self$filters_pos_below[[a]] <- bind_rows(self$filters_pos_below[[a]] %>% mutate(PrimFilter = T) , pos_below[[a]] )

}


   # bind_rows(self$poolVP, poolVP )


# save pooVP
if(is.null(self$poolVP)){

  self$poolVP <- poolVP

}else{
  self$poolVP <- bind_rows(self$poolVP, poolVP )
}



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
    mutate(fill = map(rowid, function(x){

      line <- self$poolVP %>% filter(rowid == x)
      proto <- self$protocols[[line$protocol]]

      sim <- simulations(ind_param = line %>% select(-simul), add_events = proto, icmt = self$initial_cmt_values,
                  time_vec = self$times, pardf = self$parameters_default_values,model = self$model)
      self$poolVP$simul[ self$poolVP$rowid == x] <- list(sim)
    }))

  print("Done")

})


# VP plot -----------------------------------------------------------------

VP_proj_creator$set("public", "plot_VP", function(){


  self$poolVP %>%
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

  all_param <- names(self$filters_neg_above[[1]])
  all_param <- all_param[! all_param %in% c("cellid", "rowid", "PrimFilter")]


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

  return(c(above = line_above, below = line_below))

})


# plot 2D -----------------------------------------------------------------
# x = expr(k2)
# y = expr(lambda0)
# toaddneg = VP_df
VP_proj_creator$set("public", "plot_2D", function(x, y , cmt_green = "tumVol", toaddneg = NULL, plotMain = F, add_point =F , IDref = NULL){

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


  neg_above <- invoke(self$filters_neg_above, .fn = bind_rows) #%>% filter(!!filtre)
  neg_below <- invoke(self$filters_neg_below, .fn = bind_rows) #%>% filter(!!filtre)
  pos_above <- invoke(self$filters_pos_above, .fn = bind_rows) #%>% filter(!!filtre)
  pos_below <- invoke(self$filters_pos_below, .fn = bind_rows) #%>% filter(!!filtre)




    plot_dot <-

      ggplot()+
      geom_point(data = pos_above, aes(!!x, !!y), col = "darkgreen") +
      geom_point(data = pos_below, aes(!!x, !!y), col = "darkgreen") +
      geom_point(data = neg_below, aes(!!x, !!y), col = "red") +
      geom_point(data = neg_above , aes(!!x, !!y), col = "red", alpha = 1)+
      # geom_point(data = self$filters_neg_below, aes(k2, lambda0), col = "red", alpha = 1)+
      theme_bw()+
      ggtitle( "VP tested")+
      theme(plot.title = element_text(hjust = 0.5)); plot_dot

if(!is.null(toaddneg)){

  toaddneg %>%
    left_join(self$poolVP) %>%
    filter(is.na(rowid)) -> addneg

  plot_dot <- plot_dot +
    geom_point(data = addneg, aes(!!x, !!y), col = "red")+
    geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen")

}


    if(add_point == T) plot_dot <-   plot_dot+
      geom_point(data = self$poolVP, aes(!!x, !!y))

 xana <- case_when(deparse(x) %in% self$param_increase$tumVol ~ "inc",
                   deparse(x) %in% self$param_reduce$tumVol ~ "dec",
                   T ~ "None")


 yana <- case_when(deparse(y) %in% self$param_increase$tumVol ~ "inc",
                   deparse(y) %in% self$param_reduce$tumVol ~ "dec",
                   T ~ "None")


  rectangles_above  <- neg_above %>%
    mutate(xmin = case_when(xana == "dec" ~ 0, T ~ !!x),
           xmax = case_when(xana == "dec" ~ !!x, T ~Inf),
           ymin = case_when(yana == "dec" ~ 0, T ~!!y),
           ymax = case_when(yana == "dec" ~ !!y, T ~Inf))

  rectangles_below  <- neg_below %>%
    mutate(xmin = case_when(xana == "inc" ~ 0, T ~ !!x),
           xmax = case_when(xana == "inc" ~ !!x, T ~Inf),
           ymin = case_when(yana == "inc" ~ 0, T ~!!y),
           ymax = case_when(yana == "inc" ~ !!y, T ~Inf))

  rectangles <- bind_rows(rectangles_above, rectangles_below)

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


 plot1 <-  plot_dot +
    geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "red", col = "red")+
    ggtitle( "zone rejection")



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
   geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 1,  fill = "red", col = "red")+

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


VP_proj_creator$set("public", "n_filter_reduc", function(){


  all_param <- names(self$filters_neg_above[[1]])
  all_param <- all_param[! all_param %in% c("cellid", "rowid", "PrimFilter")]


  line_compar <-   paste0("line$", all_param, " == ref$", all_param) %>%
    paste0(collapse = " & ")

# Handle negative above

for(cmtt in unique((self$targets$cmt))){

line_above_reject <- line_compar
line_above_rem_filtre_pre <- line_compar

for(a in self$param_increase[[cmtt]]){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_rem_filtre_pre)

}


for(a in self$param_reduce[[cmtt]]){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_rem_filtre_pre)
}

filters_neg_above_temp <- self$filters_neg_above[[cmtt]] %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column()

filters_neg_above_temp$PrimFilter[[1]] <- T

list_filter <- filters_neg_above_temp %>%
  filter(PrimFilter == T) %>%
  group_split(rowid)
# a <- a + 1
for(a in filters_neg_above_temp %>% filter(is.na(PrimFilter)) %>% pull(rowid)){
# print(a)

line <- filters_neg_above_temp %>%
        filter(rowid == a);line

# test if to add
 test <-  map_lgl(list_filter, function(ref){
    eval(parse_expr(line_above_reject))

  });test


 if(max(test) == 0){


   list_filter[length(list_filter) + 1]<-  list(line)

  # if he is to add, should we remove some?
  test2  <- map_lgl(list_filter, function(ref){
    eval(parse_expr(line_above_rem_filtre_pre))

  });test2

  test2[length(test2)] <- FALSE

  list_filter <- list_filter[!test2]
 }
# print("ici?")
}

self$filters_neg_above[[cmtt]] <- reduce(list_filter, bind_rows) %>% select(-rowid) %>% mutate(PrimFilter = T)


# Handle negative below



line_below_reject <- line_compar
line_below_rem_filtre_pre <- line_compar

for(a in self$param_increase[[cmtt]]){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_rem_filtre_pre)

}


for(a in self$param_reduce[[cmtt]]){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_rem_filtre_pre)
}

filters_neg_below_temp  <- self$filters_neg_below[[cmtt]] %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column()

filters_neg_below_temp$PrimFilter[[1]] <- T

list_filter <- filters_neg_below_temp %>%
  filter(PrimFilter == T) %>%
  group_split(rowid)

  # a <- a + 1
for(a in filters_neg_below_temp %>% filter(is.na(PrimFilter)) %>% pull(rowid)){
  # print(a)

  line <-  filters_neg_below_temp %>%
    filter(rowid == a);line

  # test if to add
  test <-  map_lgl(list_filter, function(ref){
    eval(parse_expr(line_below_reject))

  });test


  if(max(test) == 0){


    list_filter[length(list_filter) + 1]<-  list(line)

    # if he is to add, should we remove some?
    test2  <- map_lgl(list_filter, function(ref){
      eval(parse_expr(line_below_rem_filtre_pre))

    });test2

    test2[length(test2)] <- FALSE

    list_filter <- list_filter[!test2]
  }
  # print("ici?")
}

self$filters_neg_below[[cmtt]] <- reduce(list_filter, bind_rows) %>% select(-rowid) %>% mutate(PrimFilter = T)

# filter pos above <low



line_above_reject <- line_compar
line_above_rem_filtre_pre <- line_compar

for(a in self$param_increase[[cmtt]]){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_rem_filtre_pre)

}


for(a in self$param_reduce[[cmtt]]){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_rem_filtre_pre)
}

filters_pos_above_temp  <- self$filters_pos_above[[cmtt]] %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column()



filters_pos_above_temp$PrimFilter[[1]] <- T


list_filter <- filters_pos_above_temp %>%
                      filter(PrimFilter == T) %>%
                     group_split(rowid)
###
# a <- a + 1
for(a in filters_pos_above_temp %>% filter(is.na(PrimFilter)) %>% pull(rowid)){
  # print(a)

  line <-  filters_pos_above_temp %>%
    filter(rowid == a);line

  # test if to add
  test <-  map_lgl(list_filter, function(ref){
    eval(parse_expr(line_above_reject))

  });test


  if(max(test) == 0){


    list_filter[length(list_filter) + 1]<-  list(line)

    # if he is to add, should we remove some?
    test2  <- map_lgl(list_filter, function(ref){
      eval(parse_expr(line_above_rem_filtre_pre))

    });test2

    test2[length(test2)] <- FALSE

    list_filter <- list_filter[!test2]
  }
  # print("ici?")
}

self$filters_pos_above[[cmtt]] <-reduce(list_filter, bind_rows) %>% select(-rowid) %>% mutate(PrimFilter = T)

# handle pos belw

line_below_reject <- line_compar
line_below_rem_filtre_pre <- line_compar

for(a in self$param_increase[[cmtt]]){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_rem_filtre_pre)

}


for(a in self$param_reduce[[cmtt]]){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_rem_filtre_pre)
}

filters_pos_below_temp  <- self$filters_pos_below[[cmtt]] %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column()



filters_pos_below_temp$PrimFilter[[1]] <- T


list_filter <- filters_pos_below_temp %>%
                      filter(PrimFilter == T) %>%
                      group_split(rowid)
###
# a <- a + 1
for(a in filters_pos_below_temp %>% filter(is.na(PrimFilter)) %>% pull(rowid)){
  # print(a)

  line <- filters_pos_below_temp %>%
    filter(rowid == a);line

  # test if to add
  test <-  map_lgl(list_filter, function(ref){
    eval(parse_expr(line_below_reject))

  });test


  if(max(test) == 0){


    list_filter[length(list_filter) + 1]<-  list(line)

    # if he is to add, should we remove some?
    test2  <- map_lgl(list_filter, function(ref){
      eval(parse_expr(line_below_rem_filtre_pre))

    });test2

    test2[length(test2)] <- FALSE

    list_filter <- list_filter[!test2]
  }
  # print("ici?")
}

self$filters_pos_below[[cmtt]] <- reduce(list_filter, bind_rows) %>% select(-rowid) %>% mutate(PrimFilter = T)


}


})
