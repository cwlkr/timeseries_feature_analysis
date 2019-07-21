# animate

pacman::p_load(tidyverse, gganimate, gifski, png)

plot = ggplot(dt.data %>% filter(Image_Metadata_Site == 0 & track_id == 3), aes(x = Image_Metadata_T, y = erk.ratio)) +
  geom_line(aes(group = track_id)) + 
  geom_line(stat = "summary", fun.y = mean, color = "red") +
  transition_reveal(along = Image_Metadata_T) +
  theme_classic() + 
  labs(x = "time (min)" , y = "ERK-KTR C/N")

## doesnt work with the newest ggplot version.
gganimate::animate(plot, fps = 14, nframes = 80 )


anim_save("/home/cedric/Desktop/track_revealpos0.gif")


str(dt.data)
