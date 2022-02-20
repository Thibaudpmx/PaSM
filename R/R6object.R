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
  filters_ab_lo = NULL,
  filters_be_up = NULL,
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
VP_proj_creator$set("public", "set_targets", function(..., filter = NULL, ntime = 3, manual = NULL, timeforce = NULL, force = F){

filter <- enexpr(filter)

if(!is.null(self$targets)){



  stop("Targets have already been set, it is not possible to modify them afterwards.
This is because filters used for determining parameter values spaces of VP
rejection / acceptation are strictly dependants on those targets.")

}


if(!is.null(manual)){

  print(manual)
  print(data_segment_plot(data = self$data, targets = manual))



}else{

 targets <- data_segment(data = self$data, protocol, cmt, filter = !!filter, ntime = ntime, timeforce = timeforce) %>%
    filter(max != "0")



 print( as.data.frame(targets))

}

if(force == T){

  self$targets <- targets

  return("Done ! ")

}

  message("See the above segmentation table plus the plot for checking")
  message("Defining targets is definitive.")
  message("Is it okay (y) or do you want to perform manual changes (n) ?")


  ask <- readline(prompt="")

  if(ask == "y"){

    self$targets <- targets



    return("Done !")
  }else{

    message ("Targets not set. Use manual arg for manual targets definition")


  }

})


# VP_production -----------------------------------------------------------


VP_proj_creator$set("public", "add_VP", function(VP_df,  saven = 50, drug = NULL, update_at_end = T, time_compteur = F,  fillatend = F){

  # protocols = self$protocols

  toadd <- VP_df

  poolVP <-   VP_df %>%
    crossing(protocol = unique(self$targets$protocol)) %>%
    group_by(!!!parse_exprs(names(toadd)[names(toadd) != "protocols"])) %>%
    nest() %>%
    rowid_to_column("cellid") %>%
    unnest() %>%
    rowid_to_column("rowid") %>%
    ungroup()
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


  poolVPbef <- poolVP

  above <- double() # filre_neg_above
  below <- double() # filtre_neg_below


    # begining while 2----------------------------------------------------------
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


      # Add the columns for each output
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


        # if the line output is death
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
          cellidtorem <- poolVP[reject_eval, "cellid"]$cellid

          poolVP <- poolVP %>%
            filter(!cellid %in%cellidtorem)

          remv <- T
          # print(paste0(length(cellidtorem), " cells removed"))

          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate(nelim = length(cellidtorem))

          if(cmtt == "tumVol") above <- c(above, cellidtorem)
          ### if the line output is survival
        }else if(ab_low == F){



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
          cellidtorem <- poolVP[reject_eval, "cellid"]$cellid

          poolVP <- poolVP %>%
            filter(!cellid %in%cellidtorem)


          remv <- T
          # print(paste0(length(cellidtorem), " cells removed"))

          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate(nelim = length(cellidtorem))


          if(cmtt == "tumVol") below <- c(below, cellidtorem)
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
          # poolVP %>%
          #   slice(whichaboveloweer) %>%
          #   filter(is.na(tumVol_AL)) %>% distinct(protocols)
          poolVP[[paste0(cmtt, "_AL")]][whichaboveloweer]  <- TRUE


          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate( ninfo = if_else(is.na(poolVP_compteur_new$ninfo),0,poolVP_compteur_new$ninfo)+  length(whichaboveloweer))
          # poolVP %>% slice(whichaboveloweer) %>% filter(is.na(tumVol_AL))
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

          # poolVP %>%
          #   slice(whichbelowupper) %>%
          #   filter(is.na(tumVol_BU))

          poolVP[[paste0(cmtt, "_BU")]][whichbelowupper]  <- TRUE

          if(time_compteur == T) poolVP_compteur_new <- poolVP_compteur_new %>%
            mutate( ninfo = if_else(is.na(poolVP_compteur_new$ninfo),0,poolVP_compteur_new$ninfo)+  length(whichbelowupper))

        } # end if be_up == T

      }# end for each compartment




      if(remv == F ){

        siml <- siml %>%
          add_row(cellid = line$cellid, protocol = line$protocol, simul = list(res))

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



# Handle filter neg above

print("icii")


newfilnegab <-  poolVPbef %>%
     filter(cellid %in% above) %>%
    select( !!!parse_exprs(all_param)) %>%
    select(-protocol) %>%
    distinct()

if(is.null(self$filters_neg_above )){


  self$filters_neg_above<- newfilnegab
}else{


  self$filters_neg_above <- bind_rows(  self$filters_neg_above, newfilnegab )
}


# Handle filter neg below

newfilnegbel <-  poolVPbef %>%
  filter(cellid %in% below) %>%
  select( !!!parse_exprs(all_param)) %>%
  select(-protocol) %>%
  distinct()

if(is.null(self$filters_neg_below )){
  self$filters_neg_below<- newfilnegbel
}else{


  self$filters_neg_below <- bind_rows(  self$filters_neg_below, newfilnegbel )
}

# save pooVP

self$poolVP <- poolVP


# Fill missing profile if required
if(fillatend){

print("filling")
  self$fill_simul()
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


# plot 2D -----------------------------------------------------------------
# x = expr(k2)
# y = expr(lambda0)
VP_proj_creator$set("public", "plot_2D", function(x, y){

  x <- enexpr(x)
  y <-  enexpr(y)

 plot_dot <- ggplot()+
    geom_point(data = self$filters_ab_lo, aes(k2, lambda0), col = "darkgreen") +
   geom_point(data = self$filters_be_up, aes(k2, lambda0), col = "darkgreen") +
    geom_point(data = self$filters_neg_above, aes(k2, lambda0), col = "red", alpha = 1)+
    geom_point(data = self$filters_neg_below, aes(k2, lambda0), col = "red", alpha = 1)+
   theme_bw()+
   ggtitle( "VP tested")+
   theme(plot.title = element_text(hjust = 0.5)); plot_dot

 xana <- case_when(deparse(x) %in% self$param_increase$tumVol ~ "inc",
                   deparse(x) %in% self$param_reduce$tumVol ~ "dec",
                   T ~ "None")


 yana <- case_when(deparse(y) %in% self$param_increase$tumVol ~ "inc",
                   deparse(y) %in% self$param_reduce$tumVol ~ "dec",
                   T ~ "None")


  rectangles_above  <- self$filters_neg_above %>%
    mutate(xmin = case_when(xana == "dec" ~ 0, T ~ !!x),
           xmax = case_when(xana == "dec" ~ !!x, T ~Inf),
           ymin = case_when(yana == "dec" ~ 0, T ~!!y),
           ymax = case_when(yana == "dec" ~ !!y, T ~Inf))

  rectangles_below  <- self$filters_neg_below %>%
    mutate(xmin = case_when(xana == "inc" ~ 0, T ~ !!x),
           xmax = case_when(xana == "inc" ~ !!x, T ~Inf),
           ymin = case_when(yana == "inc" ~ 0, T ~!!y),
           ymax = case_when(yana == "inc" ~ !!y, T ~Inf))

  rectangles <- bind_rows(rectangles_above, rectangles_below)

  rectangles_above_lower <-  self$filters_ab_lo %>%
          mutate(xmin = case_when(xana == "dec" ~ 0, T ~ !!x),
         xmax = case_when(xana == "dec" ~ !!x, T ~Inf),
         ymin = case_when(yana == "dec" ~ 0, T ~!!y),
         ymax = case_when(yana == "dec" ~ !!y, T ~Inf))

  rectangles_below_upper  <- self$filters_be_up %>%
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

 plot_dot +
   geom_rect(data = rectangles_above_lower, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "red",  col = "red")+
   geom_rect(data = rectangles_below_upper, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "blue",  col = "blue")
   ggtitle( "zone above lower limit")

#
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



plot_grid(plot4,plot_grid(plot_dot, plot1, plot2, plot3))

})

# filter reduc -----------------------------------------------------------------


VP_proj_creator$set("public", "n_filter_reduc", function(){


  all_param <- names(self$filters_neg_above)
  all_param <- all_param[! all_param %in% c("cellid", "rowid")]


  line_compar <-   paste0("line$", all_param, " == ref$", all_param) %>%
    paste0(collapse = " & ")

# Handle negative above



line_above_reject <- line_compar
line_above_rem_filtre_pre <- line_compar

for(a in self$param_increase$tumVol){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_rem_filtre_pre)

}


for(a in self$param_reduce$tumVol){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_rem_filtre_pre)
}

self$filters_neg_above <- self$filters_neg_above %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column()

list_filter <- list(self$filters_neg_above %>%
  slice(1))
# a <- a + 1
for(a in 2:nrow( self$filters_neg_above)){
# print(a)

line <-  self$filters_neg_above %>%
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

self$filters_neg_above <- invoke(bind_rows, list_filter)


# Handle negative below



line_below_reject <- line_compar
line_below_rem_filtre_pre <- line_compar

for(a in self$param_increase$tumVol){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_rem_filtre_pre)

}


for(a in self$param_reduce$tumVol){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_rem_filtre_pre)
}

self$filters_neg_below <- self$filters_neg_below %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column()

list_filter <- list(self$filters_neg_below %>%
                      slice(1))
# a <- a + 1
for(a in 2:nrow( self$filters_neg_below)){
  # print(a)

  line <-  self$filters_neg_below %>%
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

self$filters_neg_below <- invoke(bind_rows, list_filter)


# filter above low



line_above_reject <- line_compar
line_above_rem_filtre_pre <- line_compar

for(a in self$param_increase$tumVol){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_rem_filtre_pre)

}


for(a in self$param_reduce$tumVol){

  line_above_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_above_reject)
  line_above_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_above_rem_filtre_pre)
}

self$filters_ab_lo  <- self$poolVP %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column() %>%
  select(-simul)

list_filter <- list(self$filters_ab_lo %>%
                      slice(1))
# a <- a + 1
for(a in 2:nrow( self$filters_ab_lo)){
  # print(a)

  line <-  self$filters_ab_lo %>%
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

self$filters_ab_lo <- invoke(bind_rows, list_filter)

# handle below up

line_below_reject <- line_compar
line_below_rem_filtre_pre <- line_compar

for(a in self$param_increase$tumVol){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_rem_filtre_pre)

}


for(a in self$param_reduce$tumVol){

  line_below_reject <- gsub(paste0(a, " *=="), paste0(a, " >="), line_below_reject)
  line_below_rem_filtre_pre <- gsub(paste0(a, " *=="), paste0(a, " <="), line_below_rem_filtre_pre)
}

self$filters_be_up <- self$poolVP %>%
  select(-starts_with("rowid")) %>%
  rowid_to_column() %>%
  select(-simul)

list_filter <- list(self$filters_be_up %>%
                      slice(1))
# a <- a + 1
for(a in 2:nrow( self$filters_be_up)){
  # print(a)

  line <-  self$filters_be_up %>%
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

self$filters_be_up <- invoke(bind_rows, list_filter)





})
