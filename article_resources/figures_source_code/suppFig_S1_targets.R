
library(PaSM)

self <- VP_proj_creator$new()

tiff(width = 2500, height = 1500,filename = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/figS1.tiff", res = 300)

self$set_targets(filter = cmt == "tumVol",timeforce =  c(10,20,30,40) )

dev.off()

shell.exec( "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/figS1.tiff")

self$targets
