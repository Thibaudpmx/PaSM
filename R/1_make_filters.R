## ---------------------------
## Script name: make_filters
##
## Purpose of script: make the filter template for a model and an OoI
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

VP_proj_creator$set("public", "make_filters", function(cmt = "tumVol"){

  all_param <- self$param

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

  for(a in  self$param_no_impact[[cmt]]){

    line_above <-  gsub(paste0("&? * line\\$", a, " *== *ref\\$",a),"", line_above)
    line_below <- gsub(paste0("&? * line\\$", a, " *== *ref\\$",a),"", line_below)
  }

  # handle the IF structure for param increase
  if_to_handle <- self$param_increase[[cmt]][grepl(" IF",  self$param_increase[[cmt]])]

  for(a in if_to_handle){

    param <- gsub(" .+", "", a)

    addabove <- paste0("|(line$",param, " >= ref$", param, gsub(paste0(param, " *IF")," &", a) ,")")
    addbelow <- paste0("|(line$",param, " <= ref$", param, gsub(paste0(param, " *IF")," &", a) ,")")

    line_above <- gsub(paste0("line\\$",param), paste0("(line$",param),line_above  ) %>%
      gsub(pattern = paste0("ref\\$",param),replacement =  paste0("ref$",param,addabove, ")") )

    line_below <- gsub(paste0("line\\$",param), paste0("(line$",param),line_below  ) %>%
      gsub(pattern = paste0("ref\\$",param),replacement =  paste0("ref$",param,addbelow, ")") )


  }


  # handle the IF structure for param decrease
  if_to_handle <- self$param_reduce[[cmt]][grepl(" IF",  self$param_reduce[[cmt]])]

  for(a in if_to_handle){

    param <- gsub(" .+", "", a)

    addabove <- paste0("|(line$",param, " <= ref$", param, gsub(paste0(param, " *IF")," &", a) ,")")
    addbelow <- paste0("|(line$",param, " >= ref$", param, gsub(paste0(param, " *IF")," &", a) ,")")

    line_above <- gsub(paste0("line\\$",param), paste0("(line$",param),line_above  ) %>%
      gsub(pattern = paste0("ref\\$",param),replacement =  paste0("ref$",param,addabove, ")") )

    line_below <- gsub(paste0("line\\$",param), paste0("(line$",param),line_below  ) %>%
      gsub(pattern = paste0("ref\\$",param),replacement =  paste0("ref$",param,addbelow, ")") )


  }





  return(c(above = line_above, below = line_below))

})
