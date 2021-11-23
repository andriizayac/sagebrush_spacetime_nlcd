library(magick)
library(raster)

file <- "img_2021_10_20_2220_35"
pathvol <- "~/../../Volumes/My Passport for Mac/"
path <- paste0(pathvol, "DCEW_phenocam/dcew_img_raw/", file, ".jpg")

r <- stack(path)

img <- image_read(path) %>% 
  image_data('rgb') %>% as.integer()

ndvi <- (r[[1]] - 2*r[[3]])/(r[[1]] + 2*r[[3]])
plot(ndvi)


for(i in 1:10){
  if(i > 5) {
    print(i)
  } 
  if(i > 8) {
    print(i)
  }
}



