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

extract_filter <- function(toadd, results = celltheque2){

above
below

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

load_spread <- function(celltheque = NULL, path_celltheque, drug, returnOnIDperGroup = T, update = F){


  # if path_celltheque is just the name, compute the full path
  if(!grepl("/", path_celltheque)) path_celltheque <- file.path(active_VT_project, "2_celltheques", "celltheques",path_celltheque)

  # get the  name
  file.name <- gsub(".+/", "", path_celltheque)

  # and read the celltheque from RDS if not provided
  if(is.null(celltheque)) celltheque <- readRDS(path_celltheque)

  # need to handle drug as list  (ex list(c(1,4), c(2,4)))

  if(is.list(drug)){
    name <-   map_chr(drug, ~ paste0(.x, collapse = "_")) %>% paste0(collapse = "_AND_")
  }else{

    name <-  paste0(drug, collapse = "_")
  }

  # Now, where is supposed to be the containing folder?
  if(returnOnIDperGroup == T){
    folder_path <- file.path(active_VT_project, "2_celltheques",
                             paste0("/celltheque_one_per_cell_spread_drug_",name))

  }else{

    folder_path <- file.path(active_VT_project, "2_celltheques",
                             paste0("/celltheque_spread_drug_", name))


  }
  # create the folder if it does not exist yet
  if(!file.exists(folder_path)) dir.create(folder_path)


  # Compute the path of where the celltheque file will be
  file <- file.path(folder_path,
                    file.name)

  if(!file.exists(file) | update == T){

    temp_spread  <- celltheque_spread(celltheque, returnOnIDperGroup = returnOnIDperGroup, drug = drug)
    saveRDS(temp_spread, file)
  }else{

    temp_spread <- readRDS(file)

  }

  return(temp_spread)

}



cellthequeLoad <- function(drug = 1, update = F, return = 1){


    # First let's take the path where the celltheque are
    root <- file.path(active_VT_project, "2_celltheques")

    # If it does not exist, create a folder containing  shrinked celltheque per drug combination
   if(!file.exists(file.path(root, "celltheque_one_per_cell"))) dir.create(file.path(root, "celltheque_one_per_cell"))

    # Compute the path of this shrinked celltheque and its "spread" version
    path <- file.path(root, "celltheque_one_per_cell", paste0("Drug", paste0(drug, collapse = "_"), ".RDS"))
    # path_spread <- gsub(".RDS", "_spread.RDS", path)


    # See if this shrink version already exists
    celltheque <- try(readRDS(path) )


    # if yes, lets register which celltheque cells come from, and read the spread version
    if(class(celltheque)[[1]] != "try-error"){

    alreadytestest <- unique(celltheque$from)
    # spread_base <- readRDS(path_spread)
    spread_base <- celltheque_spread(celltheque, returnOnIDperGroup = T, drug = drug)
    toupdate <- F
    }else{

      alreadytestest <- ""
      spread_base <- NULL
      toupdate <- T
    }


  # Now make a list of all celltheques available


  # look if we need to update (not asked by user and spread_base already exists)
  if(update == F  & !is.null(spread_base)){
    listf <- character()
  }else{
    # else, make the likst of all celltheque
    listf <- list.files(file.path(root,"celltheques"))
  }

    cont <- 0
  # For each of this celltheque
  for(a in listf[!listf %in%alreadytestest]){


    print(a)
    celltheque_temp <- readRDS(file.path(root, "celltheques", a))

    # rm every non desired concentration

    for(b in (1:ndrug)[! (1:ndrug %in%drug)]){

      if(paste0("conc",b) %in% names(celltheque_temp)){

        expr_temp <- parse_expr(paste0("conc",b))

        celltheque_temp <-  celltheque_temp %>%
          filter(!! expr_temp == 0) %>%
          select(-!!expr_temp)


      }

    }

    # First we need to make sure this cell line is complete (no NA res)
    # AND all the asked drugs are there

    if(sum(is.na(celltheque_temp$res)) == 0 &
           sum(paste0("conc", drug) %in% names(celltheque_temp)) == length(drug)){


    # Get the spread
    temp_spread <- load_spread(path_celltheque = a, drug = drug, update = T, returnOnIDperGroup = T) # try(readRDS(path_spread))


    # Then, let's make the celltheques ! If it does not exist yet
    if(is.null(spread_base) & nrow(temp_spread) >0){

      # Compute the newIDS (remember we keep only one ID per result profiles)
      temp_spread %>% pull(cellid)-> newIDs

      # celltheque is initialized as this first celltheque file
      # keeping only one ID per result profile,
      # and keeping track on which file they come from and
      # waht was their ID in the original file
      celltheque <-  celltheque_temp %>%
        filter(cellid %in% newIDs) %>%
        mutate(cellid0 = cellid, from = a)

      # spread base is directly temp_spread
      spread_base <- temp_spread

    }else if(nrow(temp_spread) >0){

      # Now imagine we already had a celltheque,
      # We don't want to add profile we already had, right?
      # So let's first compute profile we already have

      spread_base %>%
        select(starts_with("Conc")) %>%
        mutate(torem = T) -> already_have

      # and join the new profile to see those we already have (tagged by "torem")
    temp_spread %>% # taking the new spread
      select(cellid, starts_with("Conc")) %>% # keeping only cell id and all conc result columns
      left_join(already_have) %>% # making the left_join
      filter(is.na(torem)) %>%  # removing profiles we already have
      pull(cellid) -> newIDs # extracting only the new ids !


     # Now update the cellthque
      celltheque <-  celltheque %>% # take the previous one and add
        bind_rows(
          celltheque_temp %>% # new rows
                    filter(cellid %in% newIDs) %>% # but only if new profiles
                    mutate(cellid0 = cellid, from = a,  # keep track of provenance
                           cellid = cellid+max(celltheque$cellid)) # and create a new cell id
          )

     # And do the same thigs with spread_base
      spread_base <- spread_base %>%
        bind_rows(temp_spread %>%
                    filter(cellid %in% newIDs) %>%
                    mutate(from = a, cellid = cellid+max(spread_base$cellid)))
    }
    }else{

      if(sum(is.na(celltheque_temp$res)) >  0)  print("Incomplete celltheque, can not be included")
      if(sum(paste0("conc", drug) %in% names(celltheque_temp)) != length(drug)){


        missi <- paste0("conc", drug)[!paste0("conc", drug) %in% names(celltheque_temp)]
        print(paste0("Concentration needed in celltheque missing: ", paste0(missi, collapse = ", ")))

      }

      }

  }

  # Now take

  path_spread <- gsub(".RDS", "_spread.RDS", path)



# Just arrange

  drugs <- parse_exprs(paste0("conc", drug))


  celltheque <- celltheque %>%
    arrange(cellid, !!!drugs)


  if( ! (update == F  & toupdate == F)){


    saveRDS(celltheque, path)
    saveRDS(spread_base, path_spread)

  }


  # name <- paste0(gsub("(.+/)|(\\.RDS)", "", path),"_theque_minimal")
  # if(! exists(name)) eval(expr(!!parse_expr(name) <<- try(celltheque, silent = T)))
  #
  # name <- paste0(gsub("(.+/)|(\\.RDS)", "", path),"_theque_minimal_spread")
  # if(! exists(name)) eval(expr(!!parse_expr(name) <<- try(spread_base, silent = T)))

  if(length(return) == 1 & return[[1]] == 1) return(celltheque)
  if(length(return) == 1 & return[[1]] == 2) return(spread_base)
  return(list(celltheque, spread_base))

}

# cellthequeLoad(drug = 2, equilibrium = F, return = 3)
# toadd <- crossing(Bak = c(1000), Bax = c(1000), Bcl2 = seq(20,1020,120), Bclxl = seq(20,1020,120),
#          Mcl1 = seq(10,130,40), BIM = seq(0,150,50), PUMA = seq(0,150,50), NOXA =  seq(0,150,50)); toadd
# file = "D:/these/Second_project/QSP/modeling_work/celltheque_equilibrium_full"
# with equilibrium
# cellthequeFill_equilibrium(file, toadd = toadd, saven = 200)


