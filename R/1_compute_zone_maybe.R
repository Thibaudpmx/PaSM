## ---------------------------
## Script name: compute_zone_maybe
##
## Purpose of script: Compute zones we are sure the VPs are uncertains
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
#
#
#
# VP_proj_creator$set("public", "compute_zone_maybe", function(){
#
#
#   param_unique <- map(self$param, function(x){
#
#
#
#     temp <-  c(self$filters_neg_above[[x]], self$filters_neg_below[[x]] )  %>% unique
#
#
#
#     maxcol <- paste0(x, "max")
#     mincol <- paste0(x, "min")
#
#
#     if(length(temp) == 1){
#
#       df_temp <- tibble(max = temp, min = temp)
#
#     }else{
#
#       df_temp <-tibble(max =  temp) %>%
#         # add_row(lambda0max = 0) %>%
#         # add_row(lambda0max = 1) %>%
#         arrange(max) %>%
#         mutate(min = lag(max)) %>%
#         slice(-1)
#     }
#
#     df_temp <- df_temp %>%
#       select(min, max)
#
#     names(df_temp) <- c(mincol, maxcol )
#
#     df_temp
#   })
#
#   testnsquare <- map_dbl(param_unique, ~ nrow(.x)) %>% reduce(`*`)
#
#   if(testnsquare > 2000000) stop(paste0("Too many square (", testnsquare,")"))
#
#   allsquares <- rlang::invoke(.fn = crossing, .args = param_unique )
#
#   # if(direction == "below"){
#   # if(!is.null(param_increase)) list_arrange <- map(param_increase,~ expr(desc(!!parse_expr(.x))))
#   # if(!is.null(param_reduce)) list_arrange <- c(list_arrange, map(param_reduce,~ expr(!!parse_expr(.x))))
#   # }else{
#
#
#
#   for(x in names(self$param_increase)){
#
#     filters <- self$make_filters(x) %>%
#       map_chr(~ gsub("line\\$", "", .x))
#
#     filters <- paste0("!( ", filters, ")")
#
#     other <- self$param[! self$param %in% c(self$param_increase[[x]], self$param_reduce[[x]])]
#
#     for(a in c(self$param_increase[[x]], other)){
#
#       filters[1] <- gsub(paste0(" ",a, " "), paste0(" ", a, "min "), filters[1])
#       filters[2] <- gsub(paste0(" ",a, " "), paste0(" ", a, "max "), filters[2])
#
#     }
#
#     for(a in self$param_reduce[[x]]){
#
#       filters[1] <- gsub(paste0(" ",a, " "), paste0(" ", a, "max "), filters[1])
#       filters[2] <- gsub(paste0(" ",a, " "), paste0(" ", a, "min "), filters[2])
#
#     }
#
#
#
#     filter_above <- filters[1]
#
#     filters_neg_above <- self$filters_neg_above %>% filter(cmt == x)
#
#     # list_arrange <- list()
#     # if(!is.null(self$param_increase)) list_arrange <- map(self$param_increase[[1]],~ expr(!!parse_expr(.x)))
#     # if(!is.null(self$param_reduce))  list_arrange <- c(list_arrange, map(self$param_reduce[[1]],~ expr(desc(!!parse_expr(.x)))))
#
#     # if(!is.null(self$param_increase)) list_arrange <- map(self$param_increase[[1]],~ expr(desc(!!parse_expr(.x))))
#     # if(!is.null(self$param_reduce)) list_arrange <- c(list_arrange, map(self$param_reduce[[1]],~ expr(!!parse_expr(.x))))
#
#
#     # filters_neg_above <- filters_neg_above %>%
#     # arrange(!!!list_arrange)
#
#     # filters_neg_above %>%
#     #   arrange()
#
#     if(nrow(filters_neg_above) > 0){
#       for(a in 1:nrow(filters_neg_above)){
#
#         # print(a)
#         ref <- filters_neg_above %>% slice(a)
#
#         allsquares <- allsquares %>%
#           filter(!!parse_expr(filter_above))
#
#         # print(a)
#         # print(nrow(allsquares))
#       }
#     }
#
#     filter_below<- filters[2]
#     filters_neg_below <- self$filters_neg_below %>% filter(cmt == x)
#
#
#
#
#     if(nrow(filters_neg_below) > 0){
#       for(a in 1:nrow(filters_neg_below)){
#
#         ref <- filters_neg_below %>% slice(a)
#
#         allsquares <- allsquares %>%
#           filter(!!parse_expr(filter_below))
#
#       }
#     }
#   }
#
#   # can we reduce those square now?
#
#   self$zone_maybe <- allsquares
# })
#



# start functin -----------------------------------------------------------




VP_proj_creator$set("public", "compute_zone_maybe", function(zone_sure = F, keptSingleValue = T, limits = NULL){


  allcmt <- unique(self$targets$cmt) # to know if we add 0 and Inf

  param_unique <- map(self$param, function(x){



    temp <-  c(self$filters_neg_above[[x]], self$filters_neg_below[[x]])  %>% unique



    # if the paramter is included in the system
  if(map_lgl(allcmt, ~ x %in% c(self$param_no_impact[[.x]],self$param_increase[[.x]],self$param_reduce[[.x]] )) %>% max ){

    if((length(temp) == 1 &   keptSingleValue == F) |length(temp) > 1 ){ # if there is only one value of a param, do you
      # want to extrapolate to 0 et Inf, or keet it simple by saying you don't want any dimension on it
   if(is.null(limits)){
     temp <- c(temp, 0, Inf) %>% unique
   }else{

     temp <- c(temp, limits$from[limits$param == x], limits$to[limits$param == x]) %>% unique
   }

    }
  }


if(zone_sure == T ) temp <- c(temp, self$poolVP[[x]]) %>% unique

# if(max(temp %in% min(newVPs[[x]])) == 1) temp <- c(temp, 0)
# if(max(temp %in% max(newVPs[[x]])) == 1) temp <- c(temp, Inf)

# temp <- c(temp, 0, Inf)


    maxcol <- paste0(x, "max")
    mincol <- paste0(x, "min")


    if(length(temp) == 1){

      df_temp <- tibble(max = temp, min = temp)

    }else{

      df_temp <-tibble(max =  temp) %>%
        # add_row(lambda0max = 0) %>%
        # add_row(lambda0max = 1) %>%
        arrange(max) %>%
        mutate(min = lag(max)) %>%
        slice(-1)
    }

    df_temp <- df_temp %>%
      select(min, max)

    names(df_temp) <- c(mincol, maxcol )

    df_temp
  })



  testnsquare <- map_dbl(param_unique, ~ nrow(.x)) %>% reduce(`*`)

  if(testnsquare > 2000000) stop(paste0("Too many square (", testnsquare,")"))

  allsquares0 <- rlang::invoke(.fn = crossing, .args = param_unique )


  # if(direction == "below"){
  # if(!is.null(param_increase)) list_arrange <- map(param_increase,~ expr(desc(!!parse_expr(.x))))
  # if(!is.null(param_reduce)) list_arrange <- c(list_arrange, map(param_reduce,~ expr(!!parse_expr(.x))))
  # }else{
# perfor <- 1E6
# for(a in 0:(ceiling(nrow(allsquares0)/perfor)-1)){

  allsquares <-  allsquares0


  for(x in names(self$param_increase)){

    filters <- self$make_filters(x) %>%
      map_chr(~ gsub("line\\$", "", .x))

    filters <- paste0("!( ", filters, ")")

    other <- self$param[! self$param %in% c(self$param_increase[[x]], self$param_reduce[[x]])]

    for(a in c(self$param_increase[[x]], other)){

      filters[1] <- gsub(paste0(" ",a, " "), paste0(" ", a, "min "), filters[1])
      filters[2] <- gsub(paste0(" ",a, " "), paste0(" ", a, "max "), filters[2])

    }

    for(a in self$param_reduce[[x]]){

      filters[1] <- gsub(paste0(" ",a, " "), paste0(" ", a, "max "), filters[1])
      filters[2] <- gsub(paste0(" ",a, " "), paste0(" ", a, "min "), filters[2])

    }



    filter_above <- filters[1]

    filters_neg_above <- self$filters_neg_above %>% filter(cmt == x)

    # list_arrange <- list()
    # if(!is.null(self$param_increase)) list_arrange <- map(self$param_increase[[1]],~ expr(!!parse_expr(.x)))
    # if(!is.null(self$param_reduce))  list_arrange <- c(list_arrange, map(self$param_reduce[[1]],~ expr(desc(!!parse_expr(.x)))))

    # if(!is.null(self$param_increase)) list_arrange <- map(self$param_increase[[1]],~ expr(desc(!!parse_expr(.x))))
    # if(!is.null(self$param_reduce)) list_arrange <- c(list_arrange, map(self$param_reduce[[1]],~ expr(!!parse_expr(.x))))


    # filters_neg_above <- filters_neg_above %>%
    # arrange(!!!list_arrange)

    # filters_neg_above %>%
    #   arrange()






    if(nrow(filters_neg_above) > 0){

      pb <- progress_bar$new(total = nrow(filters_neg_above), format = paste0("Zone maybe - Neg above - ", x,
                                                                              "[:bar] :percent Remains: :eta"))

      for(a in 1:nrow(filters_neg_above)){

        pb$tick()
        # print(a)
        ref <- filters_neg_above %>% slice(a)

        allsquares <- allsquares %>%
          filter(!!parse_expr(filter_above))

        # print(a)
        # print(nrow(allsquares))
      }
    }

    filter_below<- filters[2]
    filters_neg_below <- self$filters_neg_below %>% filter(cmt == x)




    if(nrow(filters_neg_below) > 0){

      pb <- progress_bar$new(total = nrow(filters_neg_below), format = paste0("Zone maybe - Neg below - ", x,
                                                                              "[:bar] :percent Remains: :eta"))

      for(a in 1:nrow(filters_neg_below)){
        pb$tick()
        ref <- filters_neg_below %>% slice(a)

        allsquares <- allsquares %>%
          filter(!!parse_expr(filter_below))

      }
    }
  }

  # If zone sure

  if(zone_sure == T){

    self$poolVP %>%
      mutate(test= map_chr(simul, ~class(.x)[[1]])) %>%
      filter(test != "tbl_df") %>%
      distinct(tumVol_BU,tumVol_AL) -> extrapolated_ones


    for(a in 1:nrow(extrapolated_ones)){

      line <-  extrapolated_ones %>% slice(a)
      alref <- self$poolVP %>%
        filter(rowid %in% c(line$tumVol_AL, line$tumVol_BU))

     namepar <- self$param

     allvalues <-  map_chr(namepar, function(x){

        paste0( x, "min >= ", min(alref[[x]]), "& ", x, "max <= ", max(alref[[x]]))
      }) %>%
       paste0(collapse = " & ")


      allsquares %>%
        filter(!!parse_expr(allvalues))

    }

  }

# }
  # can we reduce those square now?

  self$zone_maybe <- allsquares


})
