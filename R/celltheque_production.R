#'
#'  Create celltheque
#' @export
#'
#'

# toadd <- crossing(ke_Venetoclax = c(0.1,0,1,10))
# drug <- c(1,4)
# add_events <- tibble(cmt = c("Venetoclax")) %>%
#   mutate(time = 0, amt = "conc1")

celltheque_produc  <- function(file.name = "first_try.RDS", toadd = NULL, saven = 50, drug = NULL, update_at_end = T, time_compteur = F, saveFilter  = F, saveSimul = T){



  # Compute the path of where the celltheque file will be
  file <- file.path(active_VT_project, "2_celltheques","celltheques", file.name)


    # Create the celltheque

  celltheque <-  toadd %>%
      group_by(!!!parse_exprs(names(toadd)[names(toadd) != "protocols"])) %>%
      nest() %>%
      rowid_to_column("cellid") %>%
      unnest() %>%
      rowid_to_column("rowid") %>%
      ungroup()


# Add the columns for each output
  targets %>%
    filter(protocols %in% unique(celltheque$protocols)) %>%
    pull(cmt) %>%
    unique -> col_to_add

  col_to_add <- c(paste0(col_to_add, "_BU"),
  paste0(col_to_add, "_AL"))

  for(a in col_to_add) celltheque[a] <- NA



### eviter les loupes infinis si un protocole n'as pas d'observion..
 crossing(protocols = unique(toadd$protocols), cmt =  targets %>%
            filter(protocols %in% unique(celltheque$protocols)) %>%
            pull(cmt) %>%
            unique) %>%
   full_join(targets) %>%
   filter(is.na(time)) -> torem

 if(nrow(torem)>0){

   for(a in 1:nrow(torem)){

cmt_to_rm <- torem$cmt[[a]]
pro <- torem$protocols[[a]]

     celltheque[[paste0(cmt_to_rm,"_BU")]][celltheque$protocols == pro] <- FALSE
     celltheque[[paste0(cmt_to_rm,"_AL")]][celltheque$protocols == pro] <- FALSE
   }

 }

# test filterDf system
#
#   filterdf <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/saveFilterDf.RDS")
#
# celltheque$res <- NULL
# for(a in 154:nrow(filterdf)){
#   b <- Sys.time()
#
# filter <- filterdf$filter[[a]]
# result <- filterdf$res[[a]]
#
# filter <- gsub("(line\\$BAK0)|(line\\$BAXc0)", "1000", filter) %>%
#   gsub(pattern = " \\& celltheque\\$group== line\\$group", replacement = "")
#
# nrow <- eval(parse_expr(filter))
# celltheque[nrow, "res"] <- result
#
# print(paste0(a, ":", length(nrow),":", difftime(Sys.time(), b)))
#
#
# }


#

# Handling death and survival agents, to greatly accelerate the process
all_param <- names(celltheque)
all_param <- all_param[! all_param %in% c("cellid", "rowid",col_to_add)]


line_compar <-   paste0("celltheque$", all_param, " == line$", all_param) %>%
  paste0(collapse = " & ")
# line_compar <- paste0("which(", line_compar,")")
  # saveRDS(object = celltheque, file = gsub("\\.RDS", "_todetermine.RDS", file))


  # Time compteur

if(time_compteur == T){

celltheque_compteur <- tibble(n = NA, time  = NA, nelim =NA, ninfo =NA_real_, computmodel = NA)
n_compteur <- 0
}

if(saveFilter == T) saveFilterDf <- character()

if(saveSimul == T) siml <- tibble(cellid = integer(), protocols = character(), simul= list())
  # just in case we never enter into the loop (if already filled, almost always useless)

  ntotal <- nrow(celltheque)
  t00 <- Sys.time()



  # begining while ----------------------------------------------------------

  # cellthequeprev <-celltheque
  # celltheque <- cellthequeprev
  cellthequeDone <- celltheque %>% slice(0)

  maxinfo <- is.na(celltheque[, col_to_add]) %>% sum

  pb <- progress_bar$new(
    format = "  VP creation [:bar] :current/:total (:percent) in :elapsed",
    total = maxinfo, clear = FALSE, width= 60)
  # pb <- progress_bar$new(total = )


  above <- double()
  below <- double()

  while(is.na(celltheque[, col_to_add]) %>% sum > 0 ){
    nn <- 0

    # ndone <- nrow(celltheque %>% filter(!is.na(res)))
    # print(ndone)


    # To gain time, we remove in the celltheque the line already done


    # celltheque <- celltheque  %>%slice(-indexdone)

    ## Allow to have intermediate save, usefull when we let computer run all night
    ## If the server crashes, we don't loose everything...

    # begining while 2----------------------------------------------------------
    while(is.na(celltheque[, col_to_add]) %>% sum > 0& nn < saven){

      newratio <- is.na(celltheque[, col_to_add]) %>% sum
     pb$update(ratio = (maxinfo -newratio )/maxinfo)

      if(time_compteur == T){

        n_compteur <- n_compteur + 1

        celltheque_compteur_new <- celltheque_compteur %>%
          slice(1) %>%
          mutate(n = n_compteur)
      }

      # Just compute some stat...
      t0 <- Sys.time()


      # Sample one rows among the not done yet

      # filter =  # which(is.na(celltheque[, col_to_add]) %>% apply(1, sum) != 0 )

       filter_to_use <-  paste0("is.na(", col_to_add,")") %>%
          paste0(collapse = "|")

       line <- celltheque %>%
         filter(!!parse_expr(filter_to_use)) %>%
         sample_n(1)


      # Now we need to handle the administrations
      # by making a temporar copy
      protocol  <- protocols[[line$protocols]]

      # add_events_line$amt[is.na(add_events_line$amt )] <- 0

      # And now we can make the simulation and extract the result !

      b <- Sys.time()
      res <- simulations(ind_param = line, add_events = protocol, returnSim = T);res

      if(time_compteur == T) celltheque_compteur_new <- celltheque_compteur_new %>%
        mutate(computmodel = as.double(difftime(Sys.time(),b, units = "sec")))


      # Add the columns for each output
      targets_temp <- targets %>%
        filter(protocols %in% line$protocols)

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

        pa_in_temp <- param_increase[[cmtt]]
        pa_in_temp <- pa_in_temp[pa_in_temp %in% all_param]

        pa_re_temp <- param_reduce[[cmtt]]
        pa_re_temp <- pa_re_temp[pa_re_temp %in% all_param]

        pa_ni_temp <- param_no_impact[[cmtt]]
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
              gsub(paste0("&? * celltheque\\$", a, " *== *line\\$",a),"", reject)
          }

          reject <- gsub("line\\$protocols",  paste0("\"", line$protocols,"\""), reject)

          # Compute the test
          reject_eval <- eval(parse_expr(reject))
          cellidtorem <- celltheque[reject_eval, "cellid"]$cellid

          celltheque <- celltheque %>%
            filter(!cellid %in%cellidtorem)

          if(saveFilter == T) saveFilterDf <- c(saveFilterDf, reject)
          remv <- T
          # print(paste0(length(cellidtorem), " cells removed"))

          if(time_compteur == T) celltheque_compteur_new <- celltheque_compteur_new %>%
            mutate(nelim = length(cellidtorem))

         if(cmtt == "tumVol") above <- c(above, cellidtorem)
          ### if the line output is survival
        }else if(ab_low == F){



          # create a copy of the line_compar with everything "=="
          # reject <- paste0(line_compar, "& is.na(celltheque$",    paste0(cmtt, "_AL"),")")
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
              gsub(paste0("&? * celltheque\\$", a, " *== *line\\$",a),"", reject)
          }

          reject <- gsub("line\\$protocols",  paste0("\"", line$protocols,"\""), reject)

          # Compute the test
          reject_eval <- eval(parse_expr(reject))
          cellidtorem <- celltheque[reject_eval, "cellid"]$cellid

          celltheque <- celltheque %>%
            filter(!cellid %in%cellidtorem)

          if(saveFilter == T) saveFilterDf <- c(saveFilterDf, reject)
          remv <- T
          # print(paste0(length(cellidtorem), " cells removed"))

          if(time_compteur == T) celltheque_compteur_new <- celltheque_compteur_new %>%
            mutate(nelim = length(cellidtorem))


          if(cmtt == "tumVol") below <- c(below, cellidtorem)
        } # end if-else








        # Now update the lines
        if(ab_low == TRUE &   paste0(cmtt, "_AL") %in% currentlyna$key){

          # create a copy of the line_compar with everything "=="
          test_above_lower_lim <- paste0(line_compar, "& is.na(celltheque$",    paste0(cmtt, "_AL"),")")
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
              gsub(paste0("&? * celltheque\\$", a, " *== *line\\$",a),"", test_above_lower_lim)
          }



          whichaboveloweer <- eval(parse_expr(test_above_lower_lim))
          # celltheque %>%
          #   slice(whichaboveloweer) %>%
          #   filter(is.na(tumVol_AL)) %>% distinct(protocols)
          celltheque[[paste0(cmtt, "_AL")]][whichaboveloweer]  <- TRUE


          if(time_compteur == T) celltheque_compteur_new <- celltheque_compteur_new %>%
            mutate( ninfo = if_else(is.na(celltheque_compteur_new$ninfo),0,celltheque_compteur_new$ninfo)+  length(whichaboveloweer))
          # celltheque %>% slice(whichaboveloweer) %>% filter(is.na(tumVol_AL))
        } # end if ab_low == T

        if(be_up == TRUE &  paste0(cmtt, "_BU") %in% currentlyna$key){

          # create a copy of the line_compar with everything "=="
          test_below_upper_lim <- paste0(line_compar, "& is.na(celltheque$",    paste0(cmtt, "_BU"),")")
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
              gsub(paste0("&? * celltheque\\$", a, " *== *line\\$",a),"", test_below_upper_lim)
          }



          whichbelowupper <- eval(parse_expr(test_below_upper_lim))

          # celltheque %>%
          #   slice(whichbelowupper) %>%
          #   filter(is.na(tumVol_BU))

          celltheque[[paste0(cmtt, "_BU")]][whichbelowupper]  <- TRUE

          if(time_compteur == T) celltheque_compteur_new <- celltheque_compteur_new %>%
            mutate( ninfo = if_else(is.na(celltheque_compteur_new$ninfo),0,celltheque_compteur_new$ninfo)+  length(whichbelowupper))

      } # end if be_up == T

      }# end for each compartment




if(remv == F & saveSimul == T){

   siml <- siml %>%
        add_row(cellid = line$cellid, protocols = line$protocols, simul = list(res))

}
  # if() siml <- tibble(cellid = integer(), protocoles = character(), simul= list())

      if(time_compteur == T)  celltheque_compteur <- bind_rows(celltheque_compteur,celltheque_compteur_new %>%
                                                                 mutate(time = as.double(difftime(Sys.time(), t0, units = "sec"))))


      nn <- nn +1
      # print(nn)
      }# fin while 1

      # Just print some stuff
      # nnewlines <- sum(!is.na(celltheque$res)) - before
      # print( paste0(nnewlines,  " new lines proceeded"))
      # print(Sys.time() - t0)
      # print(Sys.time() - t00)




print("#####################################################")
    indexdone <- which(is.na(celltheque[, col_to_add]) %>% apply(1, sum) == 0 )

    cellthequeDone <- bind_rows(cellthequeDone, celltheque %>% slice(indexdone))
    celltheque <- celltheque %>% slice(-indexdone)


    above <<- above
    below <<- below
    ## Allow to have intermediate save, usefull when we

      # pctdone <- (length(celltheque$rowid[!is.na(celltheque$res)]) +  nrow(cellthequeDone)) / ntotal
      # print(paste0("Percentage done: ", round(pctdone * 100, 3), "%"))
      # nn <- nn +1
      # print(nn)
    }#fin while 2

    ## Here happen after nsave iteration
    # print("########################### SAVNG RDS #############################")

  if(time_compteur == T)  celltheque_compteur <<- celltheque_compteur
  print(Sys.time() - t00)


  if(saveSimul == T) cellthequeDone <- cellthequeDone %>%
    left_join(siml)

  if(saveFilter == T) saveFilterDf <<- c(saveFilterDf, reject)
    return(cellthequeDone)
    # # Recompute the whole celltheque
    # celltheque <- bind_rows(celltheque, cellthequeDone)
    # cellthequeDone <- celltheque %>% slice(0)
    # # Save it
    # saveRDS(object = celltheque, file = file)

    # Start new session

  }


#   if(time_compteur == T){
#
#
#     celltheque_compteur <<- celltheque_compteur
#
#   }
#
#   if(saveFilter == T) saveFilterDf <<- saveFilterDf
#   # At the end of everything
#   # Recompute the whole celltheque
#   # celltheque <- bind_rows(celltheque, cellthequeDone)
# return(celltheque)
#
#   # And save the final and completely filled celltheque !
#   # saveRDS(object = celltheque, file = file)
#
#
#
#
#
# }
#'  Create celltheque
#' @export
#'
#'
extract_filter <- function(toadd, results = celltheque2){

# # above
# below

toadd2 <- toadd %>%
  group_by(!!!parse_exprs(names(toadd)[names(toadd) != "protocols"])) %>%
  nest() %>%
  rowid_to_column("cellid") %>%
  unnest() %>%
  rowid_to_column("rowid") %>%
  ungroup()

# Filter because above
abovedf <- toadd2 %>%
  filter(cellid %in% above)

belowdf <- toadd2 %>%
  filter(cellid %in% below)

for(a in param_increase$tumVol){

  abovedf[paste0(a, "sign")] <- ">="
  belowdf[paste0(a, "sign")] <- "<="
}

for(a in param_reduce$tumVol){

  abovedf[paste0(a, "sign")] <- "<="
  belowdf[paste0(a, "sign")] <- ">="
}

abovedf %>%
  bind_rows(belowdf) %>%
  select(-rowid, - cellid) -> temp

temp[ , order(names(temp))]

}





#'  Create celltheque
#' @export
#'
#'
reduce_filter <- function(filter = filtre_rouge){



  filter %>%
    select(k2, k2sign, lambda0, lambda0sign) %>%
    group_by(k2, k2sign, lambda0sign) %>%
    nest() %>%
    # pull(data)-> x; x <-x[[1]]
    mutate(lambda0 = map2(lambda0sign, data, function(sign, x){

       if(sign == ">=") return(min(x$lambda0))

      return(max(x$lambda0))
    })) %>%
    select(-data) %>%
    unnest() %>%
    group_by( lambda0, k2sign,lambda0sign ) %>%
    nest() %>%
    # pull(data)-> x; x <-x[[1]]
    mutate(k2 = map2(k2sign, data, function(sign, x){

      if(sign == ">=") return(min(x$k2))

      return(max(x$k2))
    })) %>%
    select(-data) %>%
    unnest()

}

#'  Create celltheque
#' @export
#'
#'
extra_filter_green <- function(results = celltheque2){

#là on parle de square
  results %>%
    select(k2, lambda0) -> temp

  crossing(a = unique(temp$k2), b = unique(temp$k2)) %>%
    filter(b >= a) %>%
    mutate(height = map2(a, b, function(a,b){

temp %>%
    filter(k2>=a & k2<=b) %>%
        group_by(k2) %>%
        summarise(max = max(lambda0), min = min(lambda0)) -> temp2

      tibble( floor = max(temp2$min), ceiling = min(temp2$max))
      # temp2

    })) %>%
    unnest() %>%
    filter(ceiling >= floor) %>%
    arrange(desc(b)) %>%
    group_by(a, floor, ceiling) %>%
    slice(1)

}

#' load_spread
#' @export
#'
#'

celltheque_produccomp_line_per_line <- function(toadd){


  toadd %>%
    slice(1:100) %>%
    rowid_to_column() %>%
    rename(protocol = protocols) %>%
    mutate(simul = map(rowid, function(x){
      print(x)
      toadd %>%
       slice(x) -> param
  # print(param)
  # print(param$protocols)
  # print(protocols[[param$protocols]])
      add_events <- protocols[[param$protocols]]

      print(add_events)
      simulations(param, add_events)

    })) %>%
    unnest()


}


# cellthequeLoad(toadd = toadd, fill = T, saven = 5)
# Load -------------------------------------------------------
#'
#' @author Thibaud Derippe
#' @export
#' Create celltheque
#'
