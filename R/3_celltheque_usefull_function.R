# Lindner equations -------------------------------------------------------
#'
#' Create celltheque
#' @author Thibaud Derippe
#' @export
plot_iv <- function(data = iv_data, ..., shortcut = 0, conc2x = F){

    filtre <- enexprs(...)

    if(shortcut == 1) filtre$additional <- expr(Drug==2 & Cell_line_bin == 1)
    if(shortcut == 2) filtre$additional <- expr( Cell_line_bin == 1)


  if(!exists("iv_data")) data("iv_data")

    if(conc2x == F){

      conc1 <- expr(conc1)
      conc2 <- expr(conc2)
    }else{


      conc1 <- expr(conc2)
      conc2 <- expr(conc1)
    }

    data %>%
    filter(!!!unname(filtre)) %>%
    ggplot(aes(!!conc1, Value,group =!!conc2, col= as.factor(!!conc2)))+
    geom_line()+
    geom_point()+
    facet_grid(Cell_line~Drug1)+
    scale_x_log10()+
    theme_bw()+
    labs(col = "Dose")


}


# verification output -----------------------------------------------------

#'
#' Create celltheque
#' @author Thibaud Derippe
#' @export
celltheque_verification <- function( celltheque, add_events, row_id = NULL, drug = NULL){

  # celltheque <- readRDS("D:/these/Second_project/QSP/modeling_work/Lind_eq_VTpckg/2_celltheques/celltheques/test_1.RDS")

  if(is.null(row_id)){

    # if no row_id are provided, sample every rows where soruce are different from cellid
    # (meaning the result has been "guessed" based on provided death and survival agents)
    # in that case, the idea of the function is to let it run and it will stop  if he
    # detects a discrepencies....

      row_id <- celltheque$cellid[celltheque$source != celltheque$cellid]
      row_id <- sample(row_id)
  }


  if(length(row_id) > 1){

    return(map(row_id, ~ celltheque_verification(.x, celltheque = celltheque,add_events = add_events, drug = drug)) %>%
      bind_rows())

  }

  if(is.null(drug)) drug <- names(celltheque)[grepl("^conc", names(celltheque))] %>%
      gsub(pattern = "conc", replacement = "") %>%
      as.double


  line <- celltheque %>%
    filter(rowid == row_id)


  add_events_line <- add_events

  # and replace all "concX" by corresponding value found in the line
  for(a in drug %>% reduce(c) %>% unique){

    add_events_line$amt[ add_events_line$amt == paste0("conc", a)] <- as.character(line[[paste0("conc", a)]])

  }

  add_events_line$amt <- as.double(add_events_line$amt)

  result <- simulations(ind_param = line, add_events = add_events_line, returnSim = F)


  if(result == line$res){

    print(paste0("Rowid ", row_id, " okay"))

  }else{

    error(paste0("Rowid ", row_id, " Problem !!!"))

  }

#
#   line %>% mutate(pct = result) %>%
#     select(res, pct, everything())
}




# matrix_cell(3600)
#
# cellthequeLoad() %>%
#   filter(cellid == 3600) %>%
#   filter(conc1 == 0)

# test <- celltheque_verification(1:100)
# test %>%
#   mutate(test2 = if_else(pct < 10, F, T)) %>%
#   filter(res != test2)


# Matrix cells -------------------------------------------------------
#'
#' Create celltheque
#' @author Thibaud Derippe
#' @export
matrix_cell <- function( celltheque = cellthequeLoad(), cellID, drug = NULL, plot = T, title = T){


  if(is.list(drug)){

    return(
      drug %>%
        map(~ matrix_cell(celltheque = celltheque, cellID = cellID, drug = .x, plot = T)+
              guides(fill = F)) %>%
        invoke(.fn = plot_grid)

    )

  }



  titlee <- paste0("cell ", cellID)
  if(class(celltheque)[[1]] == "VirtualTumor"){

    if(is.null(drug)) drug <- celltheque@drugs

    titlee <- paste0(titlee, " (x", sum(celltheque@harvest == cellID), ")")
    celltheque <- celltheque@celltheque
  }



  concn1 <- parse_expr(paste0("conc", drug[[1]]))

  if(length(drug) == 1){

    concn2 <- expr(concn2)
    temp2 <- temp2 %>% mutate(concn2 = 0)

  }else{

    concn2 <- parse_expr(paste0("conc", drug[[2]]))

  }


  if(plot == T){

    if(!"group" %in% names(celltheque)) celltheque$group <- paste0("Drug", paste0(drug, collapse = "_"))

    celltheque %>%
      filter(group == paste0("Drug", paste0(drug, collapse = "_"))) %>%
      filter( cellid == cellID )  %>% #, Drug %in% drug
      mutate(res = if_else(res == T, "death", "survive"))  -> temp2




    # print(conc1)
    # print(conc2)
    temp2 %>%
      ggplot()+
      geom_tile(aes(factor(!!concn1), factor(!!concn2), fill = res), col  = "black")+
      theme_bw()+
      labs(x = paste(deparse(concn1)," (uM)"), y = "A-1210477 (uM)", fill = "")+
      # facet_wrap(~Drug)+
      theme(
            line = element_blank(),rect =  element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.x=element_text(angle = 90, size = 8)
      ) -> temp

    if(title == T) temp <- temp + ggtitle(titlee)
    if(sum(temp2$res == "death") == 0) temp <- temp + scale_fill_manual(values = "#00BFC4")

    return(temp)

  }

  celltheque %>%
    filter( cellid == cellID) %>%
    select(!!concn1, !!concn2, res) %>%
    # group_split(Drug) %>%
    # map(~ .x %>%
    spread(key = !!concn2, value = res) %>%
    arrange(desc(!!concn1))
            # )

}


#'

cellsneeded_pdf <- function(harvest, ..., celltheque = cellthequeLoad(), drug = 2, plot = F){


  spreadd <- celltheque_spread(celltheque_theo_perfect, harvest = harvest)
  temp <-
    spreadd  %>%
    left_join(
    celltheque_spread(cellthequeLoad(), returnOnIDperGroup = T) %>% select(-cellid)
    )



librar

}

# Critical conc1  -------------------------------------------------

#' Critical concentration
#'
#' @author Thibaud Derippe
#' @export
crit_conc <- function(...,  conc2_vs_1 = F, justCount = F, drug = 2){

  # filtre <- exprs(Drug == 2)
  filtre <- enexprs(... )
  filtre$Drug <- expr(Drug ==drug)
  filtre <- unname(filtre)
  # print(filtre)
  if(conc2_vs_1 == F){

    conc1expr <- expr(conc1)
    conc2expr <- expr(conc2)
  }else{

    conc1expr <- expr(conc2)
    conc2expr <- expr(conc1)
  }
  # print(conc1)
  # print(conc2)
  if(justCount == T){

    cellthequeLoad() %>%
      filter(!!!filtre) %>%
      bind_rows({

        cellthequeLoad() %>%
          filter(!!!filtre) %>%
          group_by(cellid, !!conc2expr) %>%
          slice(1) %>%
          mutate(res = TRUE) -> temp

        temp[deparse(conc1expr)] <- Inf
        temp
      }) %>%
      filter(res == TRUE) %>%
      group_by(cellid, !!conc2expr) %>%
      slice(1) %>%
      mutate(firstconc = !!conc1expr) %>%
      group_by(firstconc, !!conc2expr) %>%
      tally() %>%
      spread(key = !!conc2expr, value = n) %>%
      arrange(desc(firstconc)) -> to_return
  }else{

    cellthequeLoad() %>%
      filter(!!!filtre) %>%
      bind_rows({

        cellthequeLoad() %>%
          filter(!!!filtre) %>%
          group_by(cellid, !!conc2expr) %>%
          slice(1) %>%
          mutate(res = TRUE) -> temp

        temp[deparse(conc1expr)] <- Inf
        temp
      }) %>%
      # distinct(cellid)
      # filter(cellid == 1) %>%
      filter(res == TRUE) %>%
      # distinct(cellid)
      group_by(cellid, !!conc2expr) %>%
      slice(1) %>%
      mutate(firstconc = !!conc1expr) %>%
      group_by(firstconc, !!conc2expr) %>%
      nest() %>%
      mutate(data = map(data, ~ .x$cellid)) %>%
      # tally() %>%
      spread(key = !!conc2expr, value = data) %>%
      arrange(desc(firstconc)) -> to_return


  }

  names(to_return)  <- c(deparse(conc1expr), paste0(deparse(conc2expr),"_",  names(to_return)[-1]))
  to_return %>%
    ungroup



}




# Sampling cells -------------------------------------------------------
#'
#' Create celltheque
#' @author Thibaud Derippe
#' @export
#'
# filter <- expr(ID == 1)
# filter_data <- expr(Drug == 1 & conc2 == 0)
# filter_data <- expr(Drug == 1  & Cell_line_bin == 1)
# celltheque =celltheques[[1]]
# drug <-c(1,4)
curve_sampling <- function(data = data_VT, filter_data , uniqueIDperOutcome = T, celltheque = cellthequeLoad(), drug = 1:2,  drugspread = NULL){
#
  # allconc1 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
  # allconc2 <- c(0,5,10,15)


  filter_data <- enexpr(filter_data)


  if(is.null(drugspread)) drugspread <- drug


  # Just to know if we want to keep all id or just the minmal
  # Unless doing thing manual, users already works with minimal celltheque anyway
  if(uniqueIDperOutcome == T){
  uniqueID <- celltheque_spread(returnOnIDperGroup = T, celltheque = celltheque, drug = drug)$cellid
  }else{

    uniqueID <- unique(celltheque$cellid)
  }



# Ok let's start

  # Take the observation data
  data %>%
    # apply the filter (cell line, drug,....)
    filter(!!filter_data) %>%
    # reduce and nest
    group_by(ID) %>%
    nest() %>%
    # pull(data) -> x; x <- x[[1]]
    mutate(data = map(data, function(x){


      # we need to guess the main drug (one with on different conc per line)
      x %>%
        select(starts_with("conc")) %>%
        gather("key","value") %>%
        distinct() %>%
        group_by(key) %>%
        tally %>%
        filter(n == nrow(x)) %>%
        pull(key) -> maindrug



      coltemp <- parse_expr(maindrug)

      # get all the concentrations, pus (for surviving)
      allconc1 <-  x[[maindrug]] %>% unique()
      allconc1 <- c(allconc1, Inf)


      # Ok so x is one curve profile
      x %>%
        # first replace value above 100 by 100
        mutate(Value = if_else(Value > 100, 100, Value)) -> x

      # Then prevent value to reincrease by comparing a value and the previous one
      # If the value is higher, replace by the previous one
      for(a in 2:nrow(x)){

        if(x$Value[[a]] >  x$Value[[a - 1]]) x$Value[[a]] <-  x$Value[[a - 1]]

      }

      # take the data
      x %>%
      # add a new row to handle From last conc to Inf
        bind_rows({

          x %>%
            slice(1) %>%
            mutate(Value = 0) -> tempp

          tempp[maindrug] <- Inf


          tempp

        }) %>%
        # make a lag to compute the difference (providing n, pct of cell that died during the step)
        mutate(Value2 = lag(Value)) %>%
        mutate(Value2 = if_else(is.na(Value2), 100, Value2)) %>%
        select(Value, Value2, everything()) %>%
        mutate(n = Value2 - Value) %>%
        # slice(-1) %>%
        mutate(n = round(n)) -> temp



      # now constitute the bags, hardest part to code!

      temp %>%
        rowid_to_column() %>%
        group_by(rowid) %>%
        nest() %>%
        # pull(data) -> y; y <- y[[2]]
        mutate(sample = map(data, function(y){

          # for each row


            # take the concentration on the line
            nconc <- which(allconc1 == y[[maindrug]])

            # and the concentration before
            nconcc <- allconc1[c(nconc-1, nconc)]

            # take our celltheque
            celltheque_temp <- celltheque

            # see all concentration present in our celltheque and data
            allconc <- names(y)[grepl("^conc", names(y))]
            allconc <- allconc[allconc %in% names(celltheque_temp)]
            # for all other concentration that is not the main,
            # filter the celltheque to keep the right secondary concentration
            # the same as in our line

            for(a in allconc[allconc != maindrug]){

            celltheque_temp <- celltheque_temp[celltheque_temp[[a]] ==  y[[a]] , ]

            }

          # so let's capture all the cells dying at the conc but surviving the conc before:

            # target
            # trick here is if there is  celltheque does not have  conc Inf, so n= 1
            # means only one row was taken (at max conc) and was false
            # while for all others rows, it means 2 rows were taken and one of the two
            # is false, one is true
            # For conc = 0 (first row), we look at cell already dead, else 1 TRUE

            test_to_do <- if_else(length(nconcc) == 1 & nconcc[[1]] == 0, TRUE, FALSE)

              celltheque_temp %>%
                filter(cellid %in% uniqueID) %>%
                ungroup() -> celltheque_temp
                #  take only the conc and the one before

              celltheque_temp[celltheque_temp[[maindrug]] %in% nconcc ,] %>%
                # now take only the cell with res:
                # we want to keep only the cell ID with exactly one remaining rows
                # see explanation above
                filter(res == test_to_do) %>%
                group_by(cellid) %>%
                tally %>%
                filter(n == 1) %>%
                pull(cellid) -> temp



              temp


      }))



        })) %>%
    unnest() %>%
    unnest(data) %>%
    select(n , sample, everything())



  }
#
#
#

# cell_harvest(sample_df = testsampling2, ncol = n, samplecol = expr(samplefinal), returnPlot = T, wrap = T, plot_by_con2 = 2)$plot
# bag -------------------------------------------------------
#'
#' bag
#' @author Thibaud Derippe
#' @export
cell_harvest <- function( sample_df = curve_sampling(), ncol = n, samplecol = sample, returnPlot = F, xlog = T,  plot_by_con2 = F,wrap = T){


  # ncol <- expr(n); samplecol <- expr(sample)
  ncol<- enexpr(ncol)
  samplecol<- enexpr(samplecol)


  # sample_df %>%
  #   select(n, sample) %>%
  #   unnest()

  sample_df %>%
    # slice(6:7) %>%
    mutate(sample = map2(!!ncol, !!samplecol, function(nn, samplee){
      # print(nn)
      if(class(samplee) == "try-error" ) return(NA)
      if(length(samplee) == 0 ) return(NA)
      sample(c(samplee, samplee), size = max(nn,0), replace = T)

    } )) %>%
    unnest(sample) %>%
    filter(!is.na(sample)) %>%
    pull(sample) -> harvest
  # %>%
  # pull(sample)
  if(length(harvest) < 100){

    harvest <- c(harvest, sample(harvest, size = 100 - length(harvest), replace = T))

  }

  if(returnPlot == F) return(harvest)
  # Verification by reconstituting the plot



  return(list(bag = harvest, plot = reconstitution_plot(harvest, xlog = xlog,  plot_by_con2 = plot_by_con2, wrap = wrap)))
}





harvest_distrib_comparison <- function(celltheque1, harvest1, harvest2, celltheque2, name1 = "KARPAS-422", name2 = "SU-DHL4") {
  #
  # celltheque1 <- readRDS("D:/these/Second_project/QSP/modeling_work/harvest_to_beat_after_cell_by_cell_cellLine2_drug2_celltheque.RDS")
  # harvest1 <- readRDS("D:/these/Second_project/QSP/modeling_work/harvest_to_beat_after_cell_by_cell_cellLine2_drug2.RDS")
  #
  #
  # harvest2 <-    readRDS("D:/these/Second_project/QSP/modeling_work/onecellperbagharvest.RDS")
  # celltheque2 <-     readRDS("D:/these/Second_project/QSP/modeling_work/onecellperbagtheque.RDS")
  #



  tibble(cellid = harvest1) %>%
    left_join(celltheque1 %>%     distinct(cellid, Bak, Bax, Bcl2, Bclxl, Mcl1, eta) %>% mutate(tumor = name1) ) %>%
    bind_rows(


      tibble(cellid = harvest2) %>%
        left_join(celltheque2 %>%     distinct(cellid, Bak, Bax, Bcl2, Bclxl, Mcl1, eta) %>% mutate(tumor = name2) )

    ) %>%
    gather("key", "value", -cellid, -tumor) %>%
    ggplot() +
    geom_density(aes(value, fill = tumor), alpha = 0.3)+
    # geom_histogram(aes(value, y = ..density fill = tumor), alpha = 0.3)+
    facet_wrap(~key, scales = "free")+
    theme_bw()


}



#' Doublon
#'
#' @author Thibaud Derippe
#' @export
harvest_doublon <- function(harvest = harvest_to_beat, celltheque){

  if("bag" %in% names(harvest)) harvest <- harvest$bag

  ol <- list() # output list

  ol$all <- tibble(cellid = harvest) %>%
    group_by(cellid) %>%
    tally %>%
    arrange(desc(n))

  ol$sum <- ol$all %>% filter(n!=1)  %>% summarise(ncell = sum(n), pct = sum(n)/length(harvest));ol$sum

  return(ol)

  allids <- ol$all %>% pull(cellid)
  #
  # allids
  #
  # idtest <- 6000

  pdf("harvest_description_by_cell.pdf", width = 7,height = 6)
  ol$all %>%
    # slice(1:3) %>%
    mutate(plot = map2(cellid, n,  function(idtest,n){

      temp <- matrix_cell(cellID = idtest, drug = 2, celltheque = celltheque, plot = F) %>%
        gather(- conc1, key = "conc2", value = "value") %>%
        arrange(conc1, conc2) %>%
        filter(value) %>%
        group_by(conc2) %>%
        mutate(conc2 = as.double(conc2)) %>%
        slice(1) %>%
        left_join(
          reconstitution_plot(celltheque = celltheque,harvest = harvest, accesRes = T)
        )



      reconstitution_plot(celltheque = celltheque,harvest = harvest, wrap = F, comparisonData = F)+
        geom_point(data = temp, aes(conc1, res, group = res ), col = "red", size = 3)+
        facet_wrap(~paste0("cellid = ", idtest, ", n = ", n))+
        labs(col = "conc2") +
        annotation_custom(
          ggplotGrob(matrix_cell(cellID = idtest, plot = T, drug = 2, celltheque = celltheque)+ guides(fill = F)),
          xmin = log(1.7), xmax = log(4), ymin = 63, ymax = 103
        )



    })) %>%
    pull(plot)
  dev.off()


  # pdf("temp_removable2.pdf")
  # ggplot()
  # dev.off()

  shell.exec("harvest_description_by_cell.pdf")



}



#' celltheque theo
#'
#' @author Thibaud Derippe
#' @export
# All theoretical possibilities
celltheque_theo_full <- function(Drug = 2){

crossing(A = 0:10, B = 0:10, C = 0:10, D = 0:10) %>% # 14K: for each line/conc2, how many concentration survive
filter(B <= A) %>% # 7K
filter(C <= B) %>% # 3k2
filter(D <= C) %>% #991
  mutate(theores = pmap(list(A,B,C,D), function(A,B,C,D){

    proj_conc1 <- conc1

    proj_conc2 <- conc4
    # just to get the name (can be optimised)
    # temp <- celltheque %>%
    #   filter(cellid == 1) %>%
    #   mutate(conc_comb = paste0("C1=",conc1," C2=", conc2)) %>%
    #   select(-conc1, -conc2, -rowid, - source) %>%
    #   spread(conc_comb, res)
    #
    # names <-  paste0("`",colnames(allprofiles)[-c(1:8)],"`")
    #
    theoRes <- crossing(conc1 = proj_conc1, conc2 = proj_conc2) %>% mutate(res = NA)

    # replace
    # for(a in proj_conc2){}
    theoRes$res[theoRes$conc2 == 0] <- c(rep( F, A), rep(T, 10 - A))
    theoRes$res[theoRes$conc2 == 5] <- c(rep( F, B), rep(T, 10 - B))
    theoRes$res[theoRes$conc2 == 10] <- c(rep( F, C), rep(T, 10 - C))
    theoRes$res[theoRes$conc2 == 15] <- c(rep( F, D), rep(T, 10 - D))

    theoRes
    # theoRes %>%
    #   mutate(conc_comb = paste0("C1=",conc1," C2=", conc2)) %>%
    #   select(-conc1, -conc2) %>%
    #   spread(conc_comb, res)
  #   theoRes %>%
  #     mutate(Drug = 2) %>%
  #     mutate(cellid = 1) %>%
  #     {matrix_cell(plot = T, celltheque = ., cellID = 1)}
  #
  #  map() c(rep( F, A), rep(T, 10 - A))
  #
  # map(0:10, ~ c(rep( F, .x ), rep(T,  10 - .x)))
  #
  #   proj_conc1

  })) -> testallcombin

testallcombin %>%
  rowid_to_column("cellid") %>%
  # select(cellid, theores) %>%
  unnest -> theo_all_celltheque

theo_all_celltheque %>%
  crossing(Drug = Drug)
}



# print pdf with all
# testallcombin %>%
#   mutate(tounest= map(theores, function(theoRes){
#
#     theoRes %>%
#       mutate(conc_comb = paste0("C1=",conc1," C2=", conc2)) %>%
#       select(-conc1, -conc2) %>%
#       spread(conc_comb, res)
#
#
#
#   })) %>%
#   unnest(tounest) -> testallcombin2
#
#
# idobs <- all_profiles(returnOnIDperGroup = T)
#
# pdf("test_all.pdf")
#
# testallcombin2 %>%
#   select(-A, -B, -C, -D) %>%
#   left_join(idobs) %>%
#   mutate(test = is.na(n)) %>%
#   select(test, n, cellid, everything()) %>%
#   mutate(n = if_else(is.na(n), 0L , n)) %>%
#   # filter(test == T) %>%
#   # slice(1:100) %>%
#   # pull(theores)
#   mutate(plot = map2(theores, n, function(theoRes, n){
#
#
#       matrix_cell(celltheque =  theoRes %>% mutate(Drug = 2, cellid = 1), cellID = 1, plot = T) +
#       ggtitle(paste0("n = ", n))
#
#   })) %>%
#   pull(plot)
#
#   dev.off()
#   shell.exec("test_all.pdf")
#


# testsampling <- curve_sampling()
# For stochastic algorithm ------------------------------------------------





# filterOF <- expr(conc2 >= 0)
# harvest <- harvest_to_beat
#
# critical <- crit_conc(conc2_vs_1 = T)
# critical2 <- crit_conc()

moduleVertical <- function(harvest, filterOF =  conc2 >=0 ){

  filterOF <- enexpr(filterOF)




  if("bag" %in% names(harvest)) harvest <- harvest$bag

  conc1possible <- unique(iv_data2$conc1)
  conc1Sample <- sample(conc1possible, 1)
  conc1Sample <- 1.30

  iv_data2 %>%
    filter(conc1 == conc1Sample) %>%
    select(conc2, Value) %>%
    mutate(Valu2 = 100 - Value / max(Value) * 100) %>%
    mutate(lagg = lag(Valu2)) %>%
    mutate(pct = Valu2 - lagg) %>%
    slice(-1) %>%
    pull(pct) -> percentage # among the resistant at conc1, percentage that die at each remaining concentration

  percentage <- c(percentage, 100 - sum(percentage))

  cells_of_interest <-  critical[c("conc2", paste0("conc1_", conc1Sample))]
  names(cells_of_interest)[[2]] <- "conc1ofinter"


  testsampling %>%
    mutate(n = n * length(harvest)/100) %>%
    filter(conc2 >= 0 )  %>%
    mutate(from_harvest = map(sample, ~ harvest[harvest %in% .x])) %>%
    mutate(neff = map_dbl(sample, ~ sum(harvest %in% .x))) %>%
    mutate(ndif = neff-n) %>%
    filter(conc2 == 0)  %>%
    # slice(9) %>% pull(from_harvest) -> from_harvest ; from_harvest <- from_harvest[[1]]
    # slice(9) %>% pull(sample) -> sample ; sample <- sample[[1]]
    mutate(newharvest = pmap(list(sample, conc1, from_harvest), function(sample, conc1, from_harvest){

     if(conc1 < conc1Sample) return(from_harvest)
      try({
      # cells_of_interest
      ntoreplace <- length(from_harvest)
      nnew <- round(ntoreplace * percentage / 100)

      cells_of_interest %>%
        arrange(conc2) %>%
        mutate(nnew = c(0, nnew)) %>%
        mutate(conc1ofinter2 = map(conc1ofinter, ~ .x[.x %in% sample])) %>%
        mutate(infal = map2(conc1ofinter2, nnew, ~ sample(c(.x, .x), .y))) %>%
        pull(infal) %>% reduce(c)
      }, silent = T)
    })) %>%
    mutate(newharvest = map2(newharvest, from_harvest, function(x, y){

      if(class(x) == "try-error") return(y)
      x
    })) %>%
    pull(newharvest) %>% reduce(c) -> output

  print("done")

  output



}



# prev <- readRDS("D:/these/Second_project/QSP/modeling_work/onecellperbagtheque.RDS")
 # saveRDS(celltheque, "D:/these/Second_project/QSP/modeling_work/onecellperbagtheque.RDS")
#'
#' celltheque_merge
#' @author Thibaud Derippe
#' @export

celltheque_merge <- function(prev = NULL, drug = 1:2, equili = T){


  if(is.null(prev)){

  celltheque <- cellthequeLoad()
  celltheque_spread(celltheque, returnOnIDperGroup = T, drug = drug) -> spread_base

  celltheque <- celltheque %>%
    filter(cellid %in% spread_base$cellid) %>%
    mutate(from = "Base")

  alreadytestest <- "Base"
  }else{

    celltheque <- prev
    alreadytestest <- unique(celltheque$from)
    celltheque_spread(celltheque, returnOnIDperGroup = T, drug = drug) -> spread_base
    celltheque <- celltheque %>%
      filter(cellid %in% spread_base$cellid)
  }


  if(equili == T){

    root <- "D:/these/Second_project/QSP/modeling_work/Lindner_equilibrium"
  }else{

    root <- "D:/these/Second_project/QSP/modeling_work/Lindner_classique"
  }


  listf <- list.files(root)
  listf <- listf[grepl("Bax....?_", listf)]
  listf <- c( listf, "another_one.RDS") # 9 new profiles !



  for(a in listf[!listf %in%alreadytestest]){
    print(a)
    celltheque_temp <- readRDS(file.path(root, a))
    temp_spread  <- celltheque_spread(celltheque_temp, returnOnIDperGroup = T, drug = drug)

    temp_spread %>%
      select(-n,  -Bak, - Bax, - Bcl2, -Bclxl, - Mcl1, - eta) %>%
      left_join(spread_base %>%
                  select(-n, - cellid, -Bak, - Bax, - Bcl2, -Bclxl, - Mcl1, - eta) %>%
                  mutate(torem = T)) %>%
      filter(is.na(torem)) %>% pull(cellid) -> newIDs

    # update the celltheque
    celltheque <-  celltheque %>%
      bind_rows(celltheque_temp %>%
                  filter(cellid %in% newIDs) %>%
                  mutate(from = a, cellid = cellid+max(celltheque$cellid)))

    spread_base <- spread_base %>%
      bind_rows(temp_spread %>%
                  filter(cellid %in% newIDs) %>%
                  mutate(from = a, cellid = cellid+max(celltheque$cellid)))


  }

  print(celltheque %>%
    distinct(cellid, from) %>%
    group_by(from) %>%
    tally)

  return(celltheque)

}

celltheque_celltag <- function(celltheque, drug ){


  celltheque %>%
    filter(group == paste0("Drug", paste0(drug, collapse = "_"))) %>%
    group_by(cellid) %>%
    nest() %>%
    # pull(data) -> data; data <- data[[1]]
    mutate(type = map(data, function(data){

      #empty tibble to fill for futur unnest
      results <- tibble(a = NA)


      # look at first drug
      matrix_drug1 <- matrix_cell(celltheque = data %>% mutate(cellid = 1), cellID =  1, drug = drug, plot = F)
      drug1 <- sum(matrix_drug1$`0`)
      nconc1 <- nrow(matrix_drug1)

      results$drug1 <- case_when(drug1 == nconc1 ~ "Native_dead",
                                 drug1 >= nconc1 / 2 ~ "High_sensitiv",
                                 drug1 > 0 ~ "Sensitiv",
                                 drug1 == 0 ~ "Resistant"
      )
      results$drug1Level <- nconc1 - drug1 + 1
      if(results$drug1Level == nconc1 + 1) results$drug1Level <- Inf
      # look at second drug
      if(length(drug) == 2){

        matrix_drug2 <- matrix_cell(celltheque = data %>% mutate(cellid = 1), cellID =  1, drug = c(drug[[2]], drug[[1]]), plot = F)
        drug2 <- sum(matrix_drug2$`0`)
        nconc2 <- nrow(matrix_drug2)

        results$drug2 <- case_when(drug2 == nconc2 ~ "Native_dead",
                                   drug2 >= nconc2 / 2 ~ "High_sensitiv",
                                   drug2 > 0 ~ "Sensitiv",
                                   drug2 == 0 ~ "Resistant"
        )
        results$drug2Level <- nconc2 - drug2 + 1
        if(results$drug2Level == nconc2 + 1) results$drug2Level <- Inf
      }
      # Finally, the total

      results$total <- map_dbl(matrix_drug1[,-1 ], ~sum(.x )) %>% sum

      nconc <-  nrow(matrix_drug1) * (ncol(matrix_drug1) - 1)

      results %>%
        select(-a) %>%
        mutate(profile = case_when(drug1 %in% c("High_sensitiv", "Sensitiv") & drug2 == "Resistant" ~ "Drug1_dependant",
                                   drug2 %in% c("High_sensitiv", "Sensitiv") & drug1 == "Resistant" ~ "Drug2_dependant",
                                   drug1  %in% c("High_sensitiv", "Sensitiv") & drug2 %in% c("High_sensitiv", "Sensitiv") ~ "Both",
                                   # drug1 == "Resistant" & drug2 == "Resistant" & total < nconc/4 ~
                                   total == 0 ~ "Unkilable",
                                   total == nconc ~ "Die anyway",
                                   drug1 == "Resistant" &drug2 == "Resistant"~ "Synergism",
                                   T ~ "hm?"
        ))


    })) %>%
    select(-data) %>%
    unnest(type)



    }






# Simulations PK/PD style !  ----------------------------------------------
#
# celltheque
# harvest <- harvest_to_beat
#
# temp <- tibble(cellid = harvest_to_beat) %>%
#   left_join(
# celltheque %>%
#   filter(cellid %in% harvest) %>%
#   select(-rowid, - Drug, - conc1, - conc2, -res, - source) %>%
#   distinct()
#   ) %>%
#   gather("Prot", "Value", -cellid)
#
# stat <- temp %>%
#   group_by(Prot) %>%
#   summarise(mean = mean(Value),sd =  sd(Value), meanlog = mean(log(Value)), sdlog = sd(log(Value)), min = min(Value), max = max(Value))
#
#
#
# nsample <-  1000
#
# samples <- tribble(~Prot, ~Distrib,
#         "Bak", "const",
#         "Bax", "const",
#         "Bcl2", "normal",
#         "Bclxl", "uniform",
#         "Mcl1", "normal",
#         "eta", "lognormal") %>%
#   left_join(stat) %>%
#   group_by(Prot) %>%
#   nest %>%
#   mutate(Value = map(data, function(data){
#
#     if(data$Distrib == "const") return(rep(1000,nsample))
#     if(data$Distrib == "normal") return(rnorm(nsample, data$mean, data$sd))
#     if(data$Distrib == "lognormal") return(rnorm(nsample, data$meanlog, data$sdlog) %>% exp)
#     if(data$Distrib == "uniform") return(runif(n = nsample, min = data$min, max = data$max))
#   })) %>%
#   mutate(ncell = map(data, ~ 1:nsample)) %>%
#   unnest()
#
# samples <- samples %>%
#   mutate(Value = if_else(Value < 0 , 0, Value))
#
# temp %>%
#   ggplot()+
#   geom_density(aes(Value))+
#   geom_density(data=samples, aes(Value), col = "red")+
#   # geom_histogram(data=samples, aes(sampling), col = "red")+
#   facet_wrap(~Prot)
#
#
# cells <- samples %>%
#   select(ncell, Value, Prot) %>%
#   spread(Prot, Value)
#
#
# tofill <- cells %>%
#   crossing(conc1 = proj_conc1, conc2 = proj_conc2) %>%
#   mutate(res = NA) %>%
#   rowid_to_column("rowid") %>%
#   mutate(Drug = 2)
#
# # library(progress)
# #
# # # t0 <- Sys.time()
# # # nn <- 0
# # pb <- progress_bar$new(total = nrow(tofill %>% filter(is.na(res))))
# # ## while loop (very similar to celltheque production....)
# # while(sum(is.na(tofill$res)) != 0 ){
# #
# #   previous = sum(!is.na(tofill$res))
# #
# #     filter = tofill$rowid[which(is.na(tofill$res))]
# #     # filter <-which(is.na(celltheque$res) & celltheque$Drug == 2 & celltheque$conc2== 5 & celltheque$Bak == 1000 & celltheque$Bax == 1000)
# #     n <- sample(x =c(filter,filter), size = 1);n
# #     # print(n)
# #     line <- tofill %>%
# #       slice(n)
# #
# #     ind_param<- tibble(BAK0 = line$Bak, BAXc0 = line$Bax, Bcl20 = line$Bcl2, Bclxl0 = line$Bclxl, Mcl10 = line$Mcl1,  ke_BCl2_I = 0, k2_bcl2_I = 2,
# #                        ke_BClxl_I = 0, k2_Bclxl_I = 2,
# #                        ke_Mcl1_I = 0, k2_Mcl1_I = 2, nameparset = "")
# #
# #     add_events <- tibble(cmt = c("Bcl2_I", "Bclxl_I", "Mcl1_I")) %>%
# #       mutate(Proto = 1, time = 0, amt = c( (line$Drug == 1) * line$conc1, (line$Drug == 2) * line$conc1, line$conc2), method = "add")
# #
# #     res <- Lindner_simulations(ind_param, eta = line$eta, add_events = add_events, returnSim = F); res
# #     # print(res)
# #     # update the specific line
# #     tofill$res[[n]] <- res > 10
# #     # if it failes, every value with lower value will fail too
# #     if(res < 10){
# #
# #
# #       # those with higher Bcl2 will fail too
# #       tempwhich <- which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 > line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res))
# #
# #       # those with higher Bclxl will fail too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl > line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with higher Mcl1 will fail too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 > line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with lower eta will fail too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta < line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with lower conc1 will fail too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 < line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with lower conc2 will fail too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 < line$conc2 & is.na(tofill$res)))
# #
# #       tofill$res[tempwhich] <- FALSE
# #
# #
# #
# #
# #
# #     }else{
# #
# #       # those with lower Bcl2 will succeed too
# #       tempwhich <- which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 < line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res))
# #
# #       # those with lower Bclxl will succeed too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl < line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with lower Mcl1 will succeed too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 < line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with higher eta will succeed too
# #       tempwhich <- c(tempwhich, which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta > line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with higher conc1 will succeed too
# #       tempwhich <- c(tempwhich,  which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 > line$conc1 & tofill$conc2 == line$conc2 & is.na(tofill$res)))
# #
# #       # those with higher conc2 will succeed too
# #       tempwhich <- c(tempwhich,  which(tofill$Bak == line$Bak & tofill$Bax == line$Bax & tofill$Bcl2 == line$Bcl2 & tofill$Bclxl == line$Bclxl & tofill$Mcl1 == line$Mcl1 & tofill$eta == line$eta & tofill$Drug == line$Drug & tofill$conc1 == line$conc1 & tofill$conc2 > line$conc2 & is.na(tofill$res)))
# #
# #       tofill$res[tempwhich] <- TRUE
# #
# #
# #
# #     }
# #     new = sum(!is.na(tofill$res))
# #
# #     pb$tick(new-previous)
# #     # print( paste0(sum(!is.na(celltheque$res)) - before,  " new lines proceeded"))
# #     # print(Sys.time() - t0)
# #     # print(Sys.time() - t00)
# #     # print(paste0("Percentage done: ", round((1 - (length(celltheque$rowid[filter]) )/ (length(celltheque$rowid[filter2]) + ndone)) * 100, 3), "%"))
# #     # print(paste0())
# #     # print(nn)
# #     # nn <- nn +1
# #     # print(nn)
# #
# # }
# # #
# #
# #
# # reconstitution_plot(celltheque = tofill %>% rename(cellid = ncell), harvest = 1:100)
# #
# # temp <- tibble(cellid = harvest_to_beat) %>%
# #   left_join(celltheque %>% distinct(Bak, Bax, Bcl2, Bclxl, Mcl1, eta, cellid))
# #
# # cor(temp$Mcl1,temp$Bcl2, method="spearman")
