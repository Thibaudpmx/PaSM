
library(PaSM)

self <- VP_proj_creator$new()


self$set_targets(filter = cmt == "tumVol",timeforce =  c(10,20,30,40) )

self$targets

