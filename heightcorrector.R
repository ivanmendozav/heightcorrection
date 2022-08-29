
get_distances <- function(dataset){
  da <- dataset
  da$Distance <- 0
  da$Acc.Distance <- 0
  for(i in 2:nrow(da)){
    pi <- da[i-1,]
    pj <- da[i,]
    d <- sqrt((pi$x-pj$x)^2+(pi$y-pj$y)^2)
    da[i,]$Distance <- round(d,4)
    da[i,]$Acc.Distance <- da[i-1,]$Acc.Distance + d
  }
  return (da)
}

#Via SF library (faster)
transform_to <- function(dataset, target_projection = 6371, source_projection = 4326 ){
  library(sf)
  gps <- dataset
  gps <- gps[!is.na(gps$Latitude),]

  d <- st_multipoint(matrix(c(gps$Longitude,gps$Latitude),nrow=nrow(gps)))
  pt1 <- st_sfc(d,crs = source_projection)
  s.sf.gcs <- st_transform(pt1, target_projection)
  m <- as.matrix(s.sf.gcs[[1]])
  gps$x <- m[,1]
  gps$y <- m[,2]
  return(gps)
}


#Via PostgreSQL
get_xy <- function(g, target_projection = "6371"){
  #6371: Mexico DF
  #32717: Ecuador
  gps <- g
  gps <- gps[!is.na(gps$Latitude),]
  gps$x <- 0; gps$y <- 0
  for(i in 1:nrow(gps)){
    p <- gps[i,]
    t <- transfor2xy(p$Latitude, p$Longitude,target_projection = target_projection)
    gps[i,]$x <- t$x
    gps[i,]$y <- t$y
  }
  return(gps)
}

#sequential density-based discounted clustering
#SDD clustering
sddclust <- function(data,tolerance=20,discount=1,minpoints=2,startclust=1){
  inittol <- tolerance
  sdd_cluster <- startclust
  data$sdd_cluster <- 0
  data[1,]$sdd_cluster <- 1
  for(i in 2:nrow(data)){
    pi <- data[i-1,]
    pj <- data[i,]
    d <- abs(pi$Height-pj$Height)
    if(d>tolerance){
      sdd_cluster <- sdd_cluster+1
      tolerance <-inittol
    }else{
      tolerance <- discount*tolerance
    }
    data[i,]$sdd_cluster <- sdd_cluster  
  }
  t <- as.data.frame(table(data$sdd_cluster))
  noise <- t[t$Freq<minpoints,1]
  return(data)
}


transfor2xy <- function(latitude, longitude, target_projection="31370"){
  library(DBI)
  #31370 Belgica Lambert
  con_t <- connect_pgsql("transmob")
  query <- paste("select round(st_x(st_transform(st_setsrid(st_makepoint(",longitude,",",latitude,"),4326),",target_projection,"))::numeric,3) as x,  round(st_y(st_transform(st_setsrid(st_makepoint(",longitude,",",latitude,"),4326),",target_projection,"))::numeric,3) as y ")              
  df <- dbGetQuery(con_t, query)
  dbDisconnect(con_t)
  return(df[1,])
}


connect_pgsql <- function(database="postgres"){
  require("RPostgreSQL")
  drv <- dbDriver("PostgreSQL")
  pw <- {  "postgres" }
  con_m <- dbConnect(drv, dbname = database, host = "localhost", port = 5432, user = "postgres", password = pw)
  rm(pw) # removes the password
  return(con_m)
}

#fill in gaps in clustered data
height_correction <- function(gps2, max_points=5, tol=1.7, replace_points=10){
  clusters <- length(unique(gps2$sdd_cluster))
  if(clusters>1){
    gps_smooth <- data.frame()
    for(i in 1:(clusters-1)){
      grupo1 <- gps2[gps2$sdd_cluster==i,]
      grupo2 <- gps2[gps2$sdd_cluster==(i+1),]
      grupo1$grupo <- 1; grupo2$grupo <- 2; 
      m <- min(c(min(nrow(grupo1),nrow(grupo2)),max_points))
      #hacer regresión (suavizado de huecos)
      r<-rbind(grupo1[(nrow(grupo1)-(m-1)):nrow(grupo1),],grupo2[1:m,])
      modelo_z <- lm(Height~Acc.Distance,data=r)
      modelo_x <- lm(x~Acc.Distance,data=r)
      modelo_y <- lm(y~Acc.Distance,data=r)
      
      #reemplazar con suavizado
      x_min <- grupo1[(nrow(grupo1)-(m-1)),]$Acc.Distance
      x_max <- grupo2[m,]$Acc.Distance
      gps4 <- data.frame(Acc.Distance=seq(from=x_min,to=x_max,by=((x_max-x_min)/(replace_points-1))))
      gps4$Height <- predict(modelo_z,gps4)
      gps4$x <- predict(modelo_x,gps4)
      gps4$y <- predict(modelo_y,gps4)
      gps4$Latitude<-0
      gps4$Longitude<-0
      gps4$Distance <- 0
      gps4$sdd_cluster <- 0
      gps4$grupo <- 0
      gps4$Acc.Distance <- 0
      gps4$dHeight <- 0
      gps4$group <- "Device M"
      start <- ifelse(i==1,1,m+1)
      if(i>1){
        gps_smooth <- gps_smooth[!gps_smooth$sdd_cluster==i,]
      }
      if(nrow(grupo2)>m){
        gps_smooth <- rbind(gps_smooth,grupo1[start:(nrow(grupo1)-m),],gps4,grupo2[(m+1):nrow(grupo2),])
      }else{
        gps_smooth <- rbind(gps_smooth,grupo1[start:(nrow(grupo1)-m),],gps4)
      }
    }
    gps_smooth <- get_distances(gps_smooth)
    gps_smooth$grupo <- NULL
    gps_smooth$secuencia <- NULL
    gps5 <- sddclust(gps_smooth,tolerance=tol)
    return(gps5)
  }else
    {print("ERROR: No gaps found (clusters number must be longer than 1")}
}


# Google API
API_KEY = "AIzaSyDW0FSCFTdDA_W66u2neDS67K2cQefjCfk"

elevations <- function(locations){
  path <- paste("https://maps.googleapis.com/maps/api/elevation/xml?locations=",locations,"&key=",API_KEY , sep = "")
  fileName <- "elevations.xml"
  download.file(url=path,fileName)
  
  require(XML)
  data <- xmlParse(fileName)
  xml_data <- xmlToDataFrame(data)
  #return (as.numeric(xml_data[["result"]][["elevation"]]))
  return (xml_data)
}

get_elevations <- function(dataset){
  library(stringi)
  coords <- paste(dataset$Latitude,",",dataset$Longitude,sep="")
  locations <- stri_paste(coords,collapse = "|")
  data <- elevations(locations)
  data = data[-1,]
  return (as.numeric(data$elevation))
}

#max locations per request 512
# https://developers.google.com/maps/documentation/elevation/usage-and-billing?hl=es_419#other-usage-limits
chunk_elevations <- function(dataset, chunk_size = 300){
  dataset$galt <- 0
  d <- dataset
  i = 1
  
  while(i < (nrow(d)) ){
    if (i < (nrow(d)-chunk_size)){
      j = i+chunk_size
    }
    else{
      j = nrow(d)
    }
    
    d[i:j,"galt"] <- get_elevations(d[i:j,])
    i=j+1
    
  }
  return(d)
}
