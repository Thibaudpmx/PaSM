# celltheque <- celltheque_temp

#' celltheque_spread
#'
#' @param celltheque
#' @param harvest
#' @param returnOnIDperGroup
#'
#' @return
#' @export
#'
#' @examples
# celltheque <- celltheque_temp
# drug <- list(c(1,4), c(2,4))
celltheque_spread <- function(drug= 1, celltheque = cellthequeLoad(), harvest= NULL, returnOnIDperGroup = F){


  if(!is.null(harvest)) celltheque <-  celltheque %>% filter(cellid %in% harvest)


  # if we want only drug1, filter to not have the other drugs

  conc_here <- names(celltheque)[grepl("^conc",  names(celltheque))]


  conc_wanted <- paste0("conc", drug)

  extra_conc <- conc_here[!conc_here %in% conc_wanted]

  if(length(extra_conc) > 0){
  filtree <- paste0(extra_conc, " == 0 ") %>% paste0(collapse = "&")

  celltheque <- celltheque %>%
    filter(!!parse_expr(filtree))
  }

  # very ugly but works....
  paste0("\"", conc_wanted, "=\"  ,", conc_wanted) %>% paste0(collapse = ",") -> temp

  paste0("celltheque %>%
           mutate(conc_comb =paste0(", temp, "))") %>% parse_expr() %>% eval -> allprofiles

  allprofiles <- allprofiles %>%
    group_by(conc_comb, cellid) %>%
    slice(1) %>%  #& to remove redondant line with all conc == 0... %>%
    select(-starts_with(conc_here), - contains("rowid"), - contains("source"), - group) %>%
    spread(conc_comb, res)


  if(returnOnIDperGroup == F) return(allprofiles)

  names <-  paste0("`",colnames(allprofiles)[grepl("^conc.=", colnames(allprofiles))],"`")

  allprofiles %>%
    group_by(!!!parse_exprs(names)) %>%
    tally %>%
    select(n, everything()) -> ns

  allprofiles %>%
    group_by(!!!parse_exprs(names)) %>%
    slice(1) %>%
    left_join(ns) %>%
    select(n, cellid,  everything()) -> ns2

  if(returnOnIDperGroup == T) return(ns2)


}
