library(tidyverse)
source("0_prep_data.R")

library(ggmap)
lichen_map <- get_stamenmap(bbox = c(left = -121, bottom = 36,
                                     right=-104, top=47),
                            zoom=5, maptype="terrain-background",
                            color = "bw", force = FALSE)

head(dat_now)
dim(dat_now)

wdth = 4
fntsz = 1.5

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Ca), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_Ca.pdf", height=0.8*wdth, width=wdth)

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=N), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_N.pdf", height=0.8*wdth, width=wdth)

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Fe), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_Fe.pdf", height=0.8*wdth, width=wdth)

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Na), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_Na.pdf", height=0.8*wdth, width=wdth)

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=K), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_K.pdf", height=0.8*wdth, width=wdth)

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Cu), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_Cu.pdf", height=0.8*wdth, width=wdth)

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Si), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_Si.pdf", height=0.8*wdth, width=wdth)

ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Pb), size=fntsz) + 
  scale_color_gradient(low = "blue", high = "red") + xlab("longitude") + ylab("latitude")
ggsave("plots/eda_map_Pb.pdf", height=0.8*wdth, width=wdth)



ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Al), size=3) + 
  scale_color_gradient(low = "blue", high = "red")
ggmap(lichen_map) + geom_point(data=dat_now, aes(x=long, y=lat, color=Ti), size=3) + 
  scale_color_gradient(low = "blue", high = "red")


