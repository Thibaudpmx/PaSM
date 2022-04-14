## ---------------------------
## Script name: filter_reduc
##
## Purpose of script: Reduce the number of filter
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
## Notes:  hybrid use w/wo the R6 object
##
##
## ---------------------------

#'  filter_reduc
#' @export
#'
#'
filter_reduc <- function(df = filters_neg_up, filtre = NULL , obj = NULL, direction = "below",
                         param_increase = NULL, param_reduce = NULL, arrangeDf = T) { #filters[[1]]

  # take all cmt
  cmts <- unique(df$cmt)



  if(!is.null(obj)){# if a obj is inputted

    # see how many ytype there is
    filters_per_cmt <- list()

    ## look first for parameters with no impacte on cmt
    cmt_w_indp <- obj$param_no_impact[map_lgl(obj$param_no_impact, ~length(.x) > 0) ]
    cmt_w_indp <- cmt_w_indp[names(cmt_w_indp) %in% cmts]

    for(a in  names(cmt_w_indp)){


      filters_per_cmt[[a]] <-  filter_reduc(df = df %>% filter(cmt == a), direction = direction, obj = NULL, filtre = obj$make_filters(cmt = a)[[direction]]) %>%
        mutate(cmt = a)
    }

    ### remaining cmt

    if(length(cmt_w_indp) >0) cmts <- cmts[!cmts %in% names(cmt_w_indp)]

    ## For each cmts

    for(a in cmts){

      # keep only filter corresponding to the current cmt
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

      # Now make the filter
      filters_per_cmt[[a]] <-  filter_reduc(df = dftemp, direction = direction, obj = NULL, filtre = obj$make_filters(a)[[direction]],
                                            param_increase = obj$param_increase[[a]], param_reduce = obj$param_reduce[[a]], arrangeDf = arrangeDf)%>%
        mutate(cmt = a)
    }

    return(filters_per_cmt %>% reduce(bind_rows))
  }


  df2 <- df %>%
    select(-starts_with("iddummy")) %>%
    rowid_to_column("iddummy")

  filtre <- gsub("line\\$", "", filtre)


  ## Bloc to reorder the param to increase the speed (see test microbenchark number 2 - both theo and practical)
  # direction = "below", param_increase = NULL, param_no_impact = NULL, param_reduce = NULL

  if(arrangeDf == T){
    list_arrange <- list()

    if(direction == "below"){
      if(!is.null(param_increase)) list_arrange <- map(param_increase,~ expr(desc(!!parse_expr(.x))))
      if(!is.null(param_reduce)) list_arrange <- c(list_arrange, map(param_reduce,~ expr(!!parse_expr(.x))))
    }else{
      if(!is.null(param_increase)) list_arrange <- map(param_increase,~ expr(!!parse_expr(.x)))
      if(!is.null(param_reduce))  list_arrange <- c(list_arrange, map(param_reduce,~ expr(desc(!!parse_expr(.x)))))
    }

    df2 <- df2 %>%
      arrange(!!!list_arrange)


  }
  ## end bloc to reorder
  if(!"PrimFilter" %in% names(df2)) df2$PrimFilter <- NA

  for(a in df2$iddummy[is.na(df2$PrimFilter)]){
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
