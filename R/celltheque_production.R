#'
#'  Create celltheque
#' @export
#'
#'

# toadd <- crossing(ke_Venetoclax = c(0.1,0,1,10))
# drug <- c(1,4)
# add_events <- tibble(cmt = c("Venetoclax")) %>%
#   mutate(time = 0, amt = "conc1")

celltheque_produc  <- function(file.name = "first_try.RDS", toadd = NULL, saven = 50, drug = NULL,
                               add_events = NULL, update_at_end = T,
                               add_events_new_drug = NULL, new_drug = NULL, time_compteur = F, saveFilter  = F){

  if(is.null(drug) & is.null(new_drug )) stop("Need drug or new_drug")
  if(is.null(drug)) drug <- new_drug
if(!is.list(drug)) drug <- list(drug)
  if(is.null(add_events) & is.null(add_events_new_drug )) stop("Need add_events or add_events_new_drug")
  if(is.null(add_events)) add_events <- add_events_new_drug


  # Expression conc: usefull to handle columns with unknown conc numbers
  name <-   drug %>% reduce(c) %>% unique()
  name <- name[order(name)]

  conc_expr <- paste0("conc", drug) %>% parse_exprs()
  conc_expr <- paste0("conc", name) %>% parse_exprs()

  # Compute the path of where the celltheque file will be
  file <- file.path(active_VT_project, "2_celltheques","celltheques", file.name)


  # The celltheque file can already exist or not, so let's try !
  celltheque <-  try(readRDS(file), silent = T)

  # If the celltheque does not exist yet, we have to create it
  if(class(celltheque)[[1]] == "try-error"){

    # Create an empty data frame that will contains the cellthque


    celltheque <- tibble(cellid = double(), rowid =double(),
                         above_lower = logical(), below_upper = logical(), source = double())

    # Add all desired concentrations umber
    for(a in conc_expr){

      celltheque <- celltheque %>%
        mutate(!!a := double())
    }


  }

  # then handle what we want to add to the cellthque
  if(!is.null(toadd)){


    # compute the new lines to add to the cellthque
    toaddtest <- try({

      map(drug, function(x){

        conc_expr_temp <- paste0("conc", x) %>% parse_exprs()

        crossing( !!! conc_expr_temp) %>%
          mutate(group = paste0("Drug", paste0( x, collapse = "_")))

      }) %>%
        bind_rows() %>%
        map_df(function(x){

          x[is.na(x)] <- 0
          x

        }) -> conc_df

      toadd %>%
        distinct() %>%
        crossing(conc_df) %>%
        # crossing( !!! conc_expr) %>%
        # next three lines to remove those already tested
        left_join(celltheque %>% mutate(Test = T)) %>%
        filter(is.na(Test)) %>%
        select(-Test, -cellid, - rowid) %>%
        mutate(res = NA, source = NA)

    }, silent = T)

    # Display a message if it faile
    if(class(toaddtest)[[1]] == "try-error"){

      # if addition left, just let a warning to say you
      warning(paste0("Failed to add the new line, please use the format : crossing(Bak = 1000, Bax = 1000, Bcl2 = seq(20,650,100), Bclxl = seq(100,1000,100),
                     Mcl1 = seq(5,80,10), eta = seq(100,1200,300))"))

    }else{

      # Otherwise lets add those new cells

      # Firstwe need to know from which value cellid and rowid starts,
      # Depending on if the celltheque already existed
      if(nrow(celltheque) == 0){

        cellidstart <- 0
        rowidstart <- 0

      }else{

        cellidstart <-  max(celltheque$cellid)
        rowidstart <- max(celltheque$rowid)
      }



      # Now we need to compute the new cell id by biological featuring (independant of concentrations)
      # these are all the columns of to add minus res, source, and all concX

      biol_feat <- names(toaddtest)[!names(toaddtest) %in% c("res", "source","group" )]
      biol_feat <- biol_feat[!grepl("^conc[[:digit:]]*$", biol_feat)]
      biol_feat <- parse_exprs(biol_feat)


      cell_id_for_join <-  toaddtest  %>%
        group_by(!!! biol_feat) %>%
        nest() %>%
        select(-data) %>%
        distinct() %>%
        rowid_to_column("cellid")


      # And that's it, now we add the new rows
      toaddtest <- left_join(

        # Joining the new dataset o add
        toaddtest%>%
          rowid_to_column(),

        # With the cell id reference table
        cell_id_for_join
      ) %>%
        select(cellid, rowid,  res, everything()) %>%
        # increment nrw rowid and cell id if there was previous value
        mutate(rowid = rowid + rowidstart, cellid = cellid + cellidstart)


      # And finally, merge the old (empty or not) and new celltheque !
      celltheque <- celltheque %>%
        bind_rows(toaddtest)

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
all_param <- all_param[! all_param %in% c("cellid", "rowid", "res", "source", "above_lower", "below_upper")]
line_compar <-   paste0("celltheque$", all_param, "== line$", all_param) %>%
  paste0(collapse = " & ")
line_compar <- paste0("which(", line_compar," & is.na(celltheque$res))")
  # saveRDS(object = celltheque, file = gsub("\\.RDS", "_todetermine.RDS", file))


  # Time compteur

if(time_compteur == T){

celltheque_compteur <- tibble(n = integer(), nline = double(), time  = double())
n_compteur <- 0
}

if(saveFilter == T) saveFilterDf <- tibble(filter = character(), res = logical())

  # just in case we never enter into the loop (if already filled, almost always useless)
  cellthequeDone <- celltheque %>%  slice(0)

  ntotal <- nrow(celltheque)
  t00 <- Sys.time()


# begining while ----------------------------------------------------------

  # cellthequeprev <-celltheque
  # celltheque <- cellthequeprev

  while(sum(is.na(celltheque$above_lower) |  is.na(celltheque$below_upper)) != 0 ){
    nn <- 0

    ndone <- nrow(celltheque %>% filter(!is.na(res)))
    print(ndone)


    # To gain time, we remove in the celltheque the line already done
    # cellthequeDone <- celltheque %>% filter(!is.na(res))
    # celltheque <- celltheque  %>% filter(is.na(res))

    ## Allow to have intermediate save, usefull when we let computer run all night
    ## If the server crashes, we don't loose everything...
    while(sum(is.na(celltheque$above_lower) |  is.na(celltheque$below_upper)) != 0 & nn < saven){



      # Just compute some stat...
      before <- sum(!is.na(celltheque$res))
      t0 <- Sys.time()


      # Sample one rows among the not done yet
      filter = which(is.na(celltheque$above_lower) | is.na(celltheque$below_upper))
      n <- sample(x =c(filter,filter), size = 1);n
      # and extract the line to be tested !
      line <- celltheque %>%
        slice(n)



      # Now we need to handle the administrations
      # by making a temporar copy
      add_events_line <- add_events

      # and replace all "concX" by corresponding value found in the line
      for(a in drug %>% reduce(c) %>% unique){

        add_events_line$amt[ add_events_line$amt == paste0("conc", a)] <- as.character(line[[paste0("conc", a)]])

      }

      add_events_line$amt <- as.double(add_events_line$amt)

      add_events_line$amt[is.na(add_events_line$amt )] <- 0

      # And now we can make the simulation and extract the result !
      dose <<- line$conc1
      res <- simulations(ind_param = line, add_events = add_events_line, returnSim = F);res

      # Finally, we update the celltheque results
      # celltheque$res[[n]] <- res
      # celltheque$source[[n]] <- line$rowid


      # Now let's see if we can extrapolate some other results

      # if the line output is death
      if(res$be_up == F){

        # create a copy of the line_compar with everything "=="
        also_dead_line <- line_compar

        # Then replace "==" by "<=" for survival parameters
        for(a in param_survive[param_survive %in% all_param]){

          also_dead_line <- gsub(paste0(a, " *=="), paste0(a, " >= "), also_dead_line)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


        }

        # Then replace "==" by ">=" for survival parameters
        for(a in param_death[param_death %in% all_param]){

          also_dead_line <- gsub(paste0(a, " *=="), paste0(a, " <= "), also_dead_line)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


        }


        # Compute the test
        also_dead <- eval(parse_expr(also_dead_line))
        cellidtorem <- celltheque[also_dead, "cellid"]$cellid


        # and modify accordingly the celltheque

         celltheque <- celltheque %>% filter(! cellid %in% cellidtorem)

         print(paste0(length(cellidtorem), " cells removed"))


        if(saveFilter == T) saveFilterDf <-saveFilterDf  %>%
          add_row(filter = also_dead_line)

       ### if the line output is survival
      }else if(res$ab_low == F){


        # create a copy of the line_compar with everything "=="
        also_survive_line <- line_compar

        # Then replace "==" by "<=" for death parameters
        for(a in param_death[param_death %in% all_param]){


          also_survive_line <- gsub(paste0(a, " *=="), paste0(a, " >= "), also_survive_line)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


        }

        # Then replace "==" by ">=" for survival parameters
        for(a in param_survive[param_survive %in% all_param]){

          also_survive_line <- gsub(paste0(a, " *=="), paste0(a, " <= "), also_survive_line)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


        }


        # Compute the test
        also_survive <- eval(parse_expr(also_survive_line))
        cellidtorem <- celltheque[also_survive, "cellid"]$cellid


        # and modify accordingly the celltheque

        celltheque <- celltheque %>% filter(! cellid %in% cellidtorem)

        print(paste0(length(cellidtorem), "cells removed"))


        if(saveFilter == T) saveFilterDf <-saveFilterDf  %>%
          add_row(filter = also_dead_line)


      }

      # Now update the lines
      if(res$ab_low == TRUE){

        # create a copy of the line_compar with everything "=="
        test_above_lower_lim <- line_compar

        for(a in param_death[param_death %in% all_param]){

          test_above_lower_lim <- gsub(paste0(a, " *=="), paste0(a, " <= "), test_above_lower_lim)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])
        }

        for(a in param_survive[param_survive %in% all_param]){

            test_above_lower_lim <- gsub(paste0(a, " *=="), paste0(a, " >= "), test_above_lower_lim)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])

        }


        whichaboveloweer <- eval(parse_expr(test_above_lower_lim))
        celltheque$above_lower[whichaboveloweer]  <- TRUE
      }

      if(res$be_up == TRUE){
        test_below_upper_lim <- line_compar

        # Then replace "==" by "<=" for death parameters
        for(a in param_death[param_death %in% all_param]){


          test_below_upper_lim <- gsub(paste0(a, " *=="), paste0(a, " >= "), test_below_upper_lim)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


        }

        # Then replace "==" by ">=" for survival parameters
        for(a in param_survive[param_survive %in% all_param]){

          test_below_upper_lim <- gsub(paste0(a, " *=="), paste0(a, " <= "), test_below_upper_lim)%>%
            gsub(pattern = paste0("line\\$", a), replacement = line[[a]])


        }


        whichbelowupper <- eval(parse_expr(test_below_upper_lim))


        celltheque$below_upper[whichbelowupper]  <- TRUE


      }

      print(Sys.time() - t00)
      }

      # Just print some stuff
      # nnewlines <- sum(!is.na(celltheque$res)) - before
      # print( paste0(nnewlines,  " new lines proceeded"))
      print(Sys.time() - t0)
      print(Sys.time() - t00)

      if(time_compteur == T){

        n_compteur <- n_compteur + 1
        celltheque_compteur <- celltheque_compteur %>%
          add_row(n = n_compteur, nline = nnewlines, time  =  as.double(difftime(Sys.time() ,t0 , units='hours')))

      }

      # pctdone <- (length(celltheque$rowid[!is.na(celltheque$res)]) +  nrow(cellthequeDone)) / ntotal
      # print(paste0("Percentage done: ", round(pctdone * 100, 3), "%"))
      # nn <- nn +1
      # print(nn)
    }

    ## Here happen after nsave iteration
    print("########################### SAVNG RDS #############################")

    return(celltheque)
    # # Recompute the whole celltheque
    # celltheque <- bind_rows(celltheque, cellthequeDone)
    # cellthequeDone <- celltheque %>% slice(0)
    # # Save it
    # saveRDS(object = celltheque, file = file)

    # Start new session

  }


  if(time_compteur == T){


    celltheque_compteur <<- celltheque_compteur

  }

  if(saveFilter == T) saveFilterDf <<- saveFilterDf
  # At the end of everything
  # Recompute the whole celltheque
  # celltheque <- bind_rows(celltheque, cellthequeDone)
return(celltheque)

  # And save the final and completely filled celltheque !
  # saveRDS(object = celltheque, file = file)





}


#' load_spread
#' @export
#'
#'

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


