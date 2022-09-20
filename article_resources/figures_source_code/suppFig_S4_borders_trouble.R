

x <- seq(0,10,2)
y <- seq(0,10,2)

segmv <- tibble(x = x, ymin = min(y), ymax = max(y))
segmh <- tibble(y = y, xmin = min(x), xmax = max(x))

VP4 <- tribble(~x, ~y,
               2,2,
               2,4,
               4,2,
               4,4,
               8,6,
               8,8)

VP2 <-
  bind_rows(
    tibble( x = seq(0,5.5,0.5), y = 2),

    tibble( x = seq(0,5.5,0.5), y = 4),
    tibble( x = seq(6.5,10,0.5), y = 6),
    tibble( x = seq(6.5,10,0.5), y = 8),
    tibble( x = 2, y = seq(0,5.5,0.5)),
    tibble( x = 4, y = seq(0,5.5,0.5)),
    tibble( x = 6, y = seq(4.5,5.5,0.5)),
    tibble( x = 8, y = seq(4.5,10,0.5))

            )%>%
    left_join(VP4 %>% mutate(test = T)) %>%
    filter(is.na(test))

VPuseless <- bind_rows(

  tibble(x = seq(0,6,0.5), y = 6),
  tibble(x = seq(6,10,0.5), y = 4),
  tibble(x =6, y = seq(6,10,0.5)),
  tibble(x =6, y = seq(0,4,0.5))
)



tiff(width = 1800, height = 1500,filename = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/figS4.tiff", res = 300)



ggplot()+
  geom_segment(data = segmv, aes(x = x, xend = x, y = ymin, yend = ymax), lty = 1)+
  geom_segment(data = segmh, aes(x = xmin, xend = xmax, y = y, yend = y), lty = 1)+
  theme_bw()+
  geom_rect(aes(xmin = 0, xmax = 6, ymin = 6, ymax = Inf), fill = "red", alpha = 0.5)+
  geom_rect(aes(xmin = 6, xmax = Inf, ymin = 0, ymax = 4), fill = "chocolate", alpha = 0.5)+
  geom_point( data = VPuseless , aes(x = x, y=y, shape = "\nAlready\nrejected\n(n = 38\n or 47*)"), col = "blue", size = 3)+
  geom_point( data =  VP4, aes(x = x, y=y, fill = "\nfour blocs\n(n = 6)\n"), shape = 21, size = 3)+

  geom_point( data = VP2, aes(x = x, y=y, fill =  "\ntwo blocs\n(n = 67)"), shape = 21, size = 3)+
  scale_shape_manual(values = 8)+
  scale_fill_manual(values = c("red", "darkgreen"))+
  labs(x  = "param X", y = "param Y", shape = "Useless", fill = "Same VPs in")

dev.off()

shell.exec( "D:/these/Second_project/QSP/modeling_work/VT_simeoni/article_QSPVP/figures_300_dpi/figS4.tiff")


crossing(x = seq(0,10,0.5), y = seq(0,10,0.5)) %>%  # 441 en tout
  filter(!(x<=6 & y >=6)) %>%
  filter(!(x>=6 & y <=4)) # 233

# 281 VPs in total
# But what we will have is:

233 + 47 + 3 *6 + 67 # 365 !

# but if we reject the bords

crossing(x = seq(0,10,0.5), y = seq(0,10,0.5)) %>%  # 441 en tout
  filter(!(x<=6 & y >=6)) %>%
  filter(!(x>=6 & y <=4)) %>%
  filter(! x %in% seq(0,10,2), ! y %in% seq(0,10,2) ) # 135


