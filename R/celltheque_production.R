
extract_filter <- function(toadd, results = poolVP2){

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





#'  Create poolVP
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

#'  Create poolVP
#' @export
#'
#'
extra_filter_green <- function(results = poolVP2){

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

poolVP_produccomp_line_per_line <- function(toadd){


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


# poolVPLoad(toadd = toadd, fill = T, saven = 5)
# Load -------------------------------------------------------
#'
#' @author Thibaud Derippe
#' @export
#' Create poolVP
#'
