
library(QSPVP)

# first producing two patients

self <- VP_proj_creator$new()


VP_df <- tibble(k1 = c(0.5),
                  k2 = c(0.5,1),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 = c(0.1,0.05),
                  lambda1 = c(12),
                  w0 = 50,
                  Vd =  40)



self$set_targets(manual = tibble(protocol = "dose50", cmt = "tumVol", time = 12, min = 0, max = 1E5))


self$add_VP(VP_df, fillatend = F, reducefilteratend = F)

self$poolVP %>%
  unnest() %>%
  ggplot()+
  geom_line(aes(time, tumVol, group = id))+
  geom_segment(aes(x = 40, xend = 40, y = 50, yend = 200), col = "red", size = 2)+
  geom_segment(aes(x = 10, xend = 10, y = 100, yend = 200), col = "red", size = 2)



# With one bloc -------------------------------------------------------------



one <- VP_proj_creator$new()


VP_df <- tibble(k1 = c(0.5),
                k2 = c(0.5,1),
                ke = 1 ,#*  seq(0.6,1.4,0.2),
                lambda0 = c(0.1,0.05),
                lambda1 = c(12),
                w0 = 50,
                Vd =  40)



one$set_targets(manual = tibble(protocol = "dose50", cmt = "tumVol", time = c(10,40), min = c(100,50), max = c(200,200)))
one$add_VP(VP_df, fillatend = F, reducefilteratend = F, npersalve = 1)

set.seed(1320) # change seed until two$plot_2D(...) has one squares
one$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T)



# With two bloc -------------------------------------------------------------



two <- VP_proj_creator$new()


VP_df <- tibble(k1 = c(0.5),
                k2 = c(0.5,1),
                ke = 1 ,#*  seq(0.6,1.4,0.2),
                lambda0 = c(0.1,0.05),
                lambda1 = c(12),
                w0 = 50,
                Vd =  40)



two$set_targets(manual = tibble(protocol = "dose50", cmt = "tumVol", time = c(10,40), min = c(100,50), max = c(200,200)))
two$add_VP(VP_df, fillatend = F, reducefilteratend = F, npersalve = 1)

set.seed(1322) # change seed until two$plot_2D(...) has two squares
two$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T)


# Final plot --------------------------------------------------------------

plotdemo <- self$poolVP %>%
  unnest() %>%
  mutate(id = paste0("\nk2 = ", k2, "\nlbd0 = ", lambda0, "\n" )) %>%
  ggplot()+
  geom_line(aes(time, tumVol, group = id, col = factor(id)), size = 2)+
  geom_segment(aes(x = 40, xend = 40, y = 50, yend = 200, lty= "Targets"), col = "red", size = 2)+
  geom_segment(aes(x = 10, xend = 10, y = 100, yend = 200), col = "red", size = 2)+
  theme_bw()+
  labs( col= "VPs", lty = "", x = "Time (days)", y = "Tumor volume (mm3)")+
  scale_color_manual(values = c("deepskyblue1", "black"))+
  scale_y_log10(); plotdemo



plot1 <- one$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T) +
  geom_hline(data =  self$poolVP, aes(yintercept = lambda0), lty=2)+
  geom_vline(data =  self$poolVP, aes(xintercept = k2), lty=2)+
  geom_text(data = tibble( x = c(0.25,0.75,0.75), y = c(0.15,0.025,0.075), label = c(1, 3, 2)), aes(x=x, y = y, label = label), size = 8)+
  geom_point(data =  self$poolVP %>% filter(id == 1) , aes(k2, lambda0), col = "deepskyblue1", size = 2)+
  geom_point(data =  self$poolVP %>% filter(id == 2) , aes(k2, lambda0), col = "black", size = 2)+

  coord_cartesian(ylim = c(0,0.2), xlim = c(0,1.5)); plot1


plot2 <- two$plot_2D(x = k2, y = lambda0, plotoreturn = 1,add_point = T) +
  geom_text(data = tibble( x = c(0.25,0.75,0.75), y = c(0.15,0.025,0.075), label = c(1, 3, 2)), aes(x=x, y = y, label = label), size = 8)+
  geom_hline(data =  self$poolVP, aes(yintercept = lambda0), lty=2)+
  geom_vline(data =  self$poolVP, aes(xintercept = k2), lty=2)+
  geom_point(data =  self$poolVP %>% filter(id == 1) , aes(k2, lambda0), col = "deepskyblue1", size = 2)+
  geom_point(data =  self$poolVP %>% filter(id == 2) , aes(k2, lambda0), col = "black", size = 2)+
  coord_cartesian(ylim = c(0,0.2), xlim = c(0,1.5)); plot2

cowplot::plot_grid(

  plotdemo

  ,

  plot1+ ggtitle("If black profile sampled"),
  plot2+ ggtitle("If blue profile sampled"),
nrow =1, labels = LETTERS
)
