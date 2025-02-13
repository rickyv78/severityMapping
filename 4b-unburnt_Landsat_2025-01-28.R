library(raster)
library(tidyverse)
library(sf)
library(lubridate)
library(fasterize)
library(here)
library(rgdal)
library(doParallel)

v <- paste0("v", Sys.Date())
#v <- "v2023-02-20"

lshp <- list.files(here::here("inputs", "shpByBurn"), pattern = "shp$")

shpx <- st_read("M:\\Zdrive\\DEC\\Prescribed_Bushfire_Outcomes_2018-134\\DATA\\Working\\Historical\\xIndex\\DBCA_FireHistory_1987to2017_Id.shp", quiet = TRUE) %>%
  st_transform("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

tlist <- as.data.frame(list.files(here::here("dNBR"), pattern = ".tiff$"))
colnames(tlist) <- "tif"
tlist <- mutate(tlist, Burn.Id = str_split_fixed(tif, "_", 2)[,1])

plist <- str_split_fixed(list.dirs(here("all_rgbs"), 
                                   full.names = FALSE, recursive = FALSE), "_", 2)[,2] 

burns <- tlist$Burn.Id[tlist$Burn.Id %in% plist]
#burns <- "BWD-1990-91349817"

w.lst <- list.files("M:\\Zdrive\\DEC\\Prescribed_Bushfire_Outcomes_2018-134\\DATA\\Working\\sw_woodyVeg", pattern = "tif$")

i <- 1
#Define how many cores (memory is limiting factor here)
UseCores <- 10
#Register CoreCluster
cl <- makeCluster(UseCores)
registerDoParallel(cl)













foreach(i = 1:length(burns)) %dopar% {
  library(raster)
  library(tidyverse)
  library(sf)
  library(lubridate)
  library(fasterize)
  library(lwgeom)
  library(here)
  
  t.burn <- filter(tlist, Burn.Id == burns[i])
    
  burn.shp <- st_read(here::here("inputs", "shpByBurn", paste0(burns[i], "_boundry.shp")), 
                        quiet = TRUE)%>%
      st_make_valid()
    plot(burn.shp[,1])
  yr <- year(burn.shp$date)

  burn.shpx <- filter(shpx, FIH_YEAR1 == yr | FIH_YEAR1 == yr-1 | FIH_YEAR1 == yr+1)
  burn.shpx <- filter(burn.shpx, BURNID != burns[i]) %>%
    dplyr::select(BURNID)
  burn.shpx <- rbind(burn.shpx, st_transform(dplyr::select(burn.shp, BURNID), crs(burn.shpx) )) %>%
     st_cast("MULTIPOLYGON")
  ##
  if (dir.exists(here("notVeg"))){
    noVeg.shps.list <- list.files(here("notVeg"), pattern = "shp$", full.names = TRUE)
    noVeg.shps <- st_read(noVeg.shps.list[1])
    plot(noVeg.shps[,1])
    if(("BURNID" %in% colnames(noVeg.shps))==FALSE){
      noVeg.shps$BURNID <- burns[i]
    }else{
      noVeg.shpsF <- filter(noVeg.shps, BURNID == burns[i])
      if(nrow(noVeg.shpsF) != 0){
        noVeg.shps <- filter(noVeg.shps, BURNID == burns[i])
      }else(
        noVeg.shps <- filter(noVeg.shps, is.na(BURNID))
      )
    }
    
    #plot(noVeg.shps[,1])
    burn.shpx <- rbind(burn.shpx, st_transform(dplyr::select(noVeg.shps, BURNID), crs(burn.shpx) )) %>%
      st_cast("MULTIPOLYGON")
  }
  
  burn.shpx$n <- 2
  shp.buf <- st_buffer(burn.shp, dist = 240)
  
  dnbr <- raster(here("dNBR", t.burn$tif[1])) %>%
    mask(mask = shp.buf) %>%
    crop(shp.buf)
  plot(dnbr)
  
  unburnt.rst <- fasterize(burn.shpx, dnbr, field = "n")
  unburnt.rst[is.na(unburnt.rst)] <- 1
  unburnt.rst[unburnt.rst == 2] <- NA
  plot(unburnt.rst) #was burn.rst
  
  dnbr.ub <- dnbr * unburnt.rst # added
  
  #######################################
  # get peremial veg layer
  w.i <- w.lst[str_detect(w.lst, as.character(yr))]
  y <- 1
  while (length(w.i) == 0){
    w.i <- w.lst[str_detect(w.lst, as.character(yr+y))]
    y <- y+1
  }
  
  rst.per <- raster(paste0("M:\\Zdrive\\DEC\\Prescribed_Bushfire_Outcomes_2018-134\\DATA\\Working\\sw_woodyVeg\\", w.i[1]))  
 #########################################
  
  per.i <- crop(rst.per, st_transform(st_buffer(burn.shp, 240), crs(rst.per)))
  per.i <- projectRaster(per.i, dnbr, method = "ngb")
  per.i[per.i >= 1] <- 1

  dnbr.ub <- dnbr.ub * per.i
  plot(dnbr.ub)
  
  burn.rst.buf <- dnbr.ub
  burn.rst.buf[is.na(burn.rst.buf)==FALSE] <- 1
  
  buff.pix.predict <- (sum(as.numeric(st_area(st_buffer(burn.shp, 240)))) - sum(as.numeric(st_area(burn.shp))))/900
  
  if(nrow(freq(burn.rst.buf))==1){
    buff.pix.actual <- 1
  }else{
    buff.pix.actual <- as.numeric(freq(burn.rst.buf)[1,2])
  }
  
  # check if there is more than 5% of the burn perimeter that was used
  if((buff.pix.actual/buff.pix.predict)>0.05){
    q <- quantile(dnbr.ub, probs = seq(0, 1, 0.10))
    threshold <- q[10]
    buf.dif <- q[10]
    if (threshold > 0.05){
      threshold <- 0.05
    }
    unburnt <- dnbr 
    unburnt[unburnt>threshold] <- NA
    unburnt[unburnt<=threshold] <- 1
    cat("buffer difference for", burns[i], " = ", buf.dif, "\n")
  }else{
    unburnt <- dnbr 
    threshold <- 0.05
    buf.dif <- 0.05
    unburnt[unburnt>threshold] <- NA
    unburnt[unburnt<=threshold] <- 1
    cat("unburnt threshold for", burns[i], " = 0.05\n")
  }
  
  ub.stats <- data.frame(BURNID = as.character(), threshold = as.numeric())[1, ]
  ub.stats$BURNID[1] <- burns[i]
  ub.stats$threshold[1] <- buf.dif
  saveRDS(ub.stats, here("tmp", paste0("ub_", burns[i])))
  
  dir.create(here("bufferStats"), showWarnings = FALSE)
  dir.create(here("bufferStats", "figures"), showWarnings = FALSE)
  png(here::here("bufferStats", "figures", paste0(burns[i], ".png")), 
      width = 550, height = 550)
  
  plot(dnbr.ub, main = paste0(burns[i], ": ", round(buf.dif, digits = 2)))
  
  dev.off()
  
  
  unburnt <- mask(unburnt, mask = burn.shp)
  plot(unburnt)
  
  #v <- paste0("v", Sys.Date())
  dir.create(here::here(v), showWarnings = FALSE)
  dir.create(here::here(v, "actual_burnt"), showWarnings = FALSE)
  
 
      writeRaster(unburnt, here(v, "actual_burnt",paste0(burns[i], "_unburnt.tif") ), overwrite=TRUE)
  
  if(T){
  burnt <- dnbr 
  burnt[burnt>=threshold] <- 1
  burnt[burnt<threshold] <- NA
  burnt <- mask(burnt, mask = burn.shp)
  plot(burnt)
  
  #check if there is some burnt area
  df.burnt <- na.omit(as.data.frame(freq(burnt)))
  if (nrow(df.burnt) != 0){
    burnt.ply <- rasterToPolygons(burnt, dissolve = TRUE)
    #plot(burnt.ply)
    st_write(st_as_sf(burnt.ply), here(v, "actual_burnt",paste0(burns[i], "_burnt_area.shp") ), 
           append=FALSE, quiet = TRUE)
  }
  }
#}
}
stopCluster(cl)


ub <- list.files(here("tmp"), pattern = "ub", full.names = TRUE)
ub.stats <- lapply(ub, readRDS) %>% bind_rows()
ub.stats <- filter(ub.stats, BURNID %in% plist)
ggplot(ub.stats, aes(BURNID, threshold))+
  geom_col()+
  geom_hline(yintercept=0.15, linetype="dashed", color = "red")+
  geom_hline(yintercept=0.2, linetype="solid", color = "red")+
  geom_text(aes(label = round(threshold, 2)), hjust = -0.1)+
  #coord_cartesian()+
  coord_flip(ylim=c(0, 0.3))+
  theme_bw()
d <- list.files(here("bufferStats"), pattern = "unburnt_dif")
ggsave(here("bufferStats", paste0("unburnt_dif_",length(d)+1 , ".jpg")), height = nrow(ub.stats)/5, width = 4)
