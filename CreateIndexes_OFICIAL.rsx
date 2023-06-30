#Generar INDICES NDVI, NDBI, BSI, LST y UHI
##AnalisisIndices=group
##Folder_correccion_DOS1=Folder
##Banda_Termica=raster
##Layer=vector polygon
##Landsat=selection Landsat5;Landsat7;Landsat8
##Parametros_para_calcular_LST=selection Atmosféricos;Vapor_de_Agua
##Transmisividad_Atmosferica=number
##Radiancia_Atmosferica_Ascendente=number
##Radiancia_Atmosferica_Descendente=number
##Vapor_de_Agua=number
##NDVI=output raster
##NDBI=output raster
##BSI=output raster
##Temperatura_Superficial_Terrestre=output raster
##Isla_de_Calor_Urbano=output raster
##Temperatura_de_brillo_en_sensor_del_satelite=output raster
##Combinacion_Bandas=output raster

#Cargar librerías
library(sp)
library(rgdal)
library(raster)
library(sf)
library(rgeos)
library(gtools)

#Conversión de entradasde selección a vectores para los condicionales
typeLandsat<-c("Landsat5","Landsat7","Landsat8")
typeLandsat<-typeLandsat[Landsat+1]

typeparametro<-c("Atmosféricos","Vapor_de_Agua")
typeparametro<-typeparametro[Parametros_para_calcular_LST+1]

#Funcion Corte 
Corte <- function(B,D){
  RS_Bcrop <-crop(B,D)
  RS_Bmask <- mask(RS_Bcrop,D)
  return(RS_Bmask) }

#Funciones para Indices NDVI, NDBI y BSI
NDVI_f <- function(NIR,RED,LimACRSB){
    NDVI_crsRst<-(Corte(NIR,LimACRSB)-Corte(RED,LimACRSB))/(Corte(NIR,LimACRSB)+Corte(RED,LimACRSB))
    return(NDVI_crsRst)}

NDBI_f <- function(SWIR,NIR,LimACRSB,LayerTemp){
    NDBI_crsRst<-(Corte(SWIR,LimACRSB)-Corte(NIR,LimACRSB))/(Corte(SWIR,LimACRSB)+Corte(NIR,LimACRSB))
    NDBI<-projectRaster(NDBI_crsRst,crs=crs(LayerTemp))
    return(NDBI) }

BSI_f <- function(RED,SWIR,NIR,BLUE,LimACRSB,LayerTemp){
    BSI_crsRst<-((Corte(RED,LimACRSB)+Corte(SWIR,LimACRSB))-(Corte(NIR,LimACRSB)+Corte(BLUE,LimACRSB)))/
    ((Corte(RED,LimACRSB)+Corte(SWIR,LimACRSB))+(Corte(NIR,LimACRSB)+Corte(BLUE,LimACRSB)))
    BSI<-projectRaster(BSI_crsRst,crs=crs(LayerTemp))
    return(BSI) }

#Radiancia del sensor TOA Lsen
Lsen_f<-function(RAD_MULT_BAND,RAD_ADD_BAND,Banda_Termica,Lim){
    LsenB <-RAD_MULT_BAND*Banda_Termica+RAD_ADD_BAND
    LsenLim <- Corte(LsenB,Lim)
    return (LsenLim)}

#Temperatura Brillo sensor
TBS_f<- function(K1_CONSTANT_BAND,K2_CONSTANT_BAND,Lsen,LimACRSB){
    TsenB <- (K2_CONSTANT_BAND)/(log((K1_CONSTANT_BAND/Lsen)+ 1))
    TsenK <- Corte(TsenB,LimACRSB) 
    return (TsenK)}

#Temperatura Superficial Terrestre con Parámetros Atmosféricos Ta,Lu,Ld 
LST_L <- function(Lsen,Tsen,E,Ta,Lu,Ld,b){
    FA1 <- 1/Ta                 #Función atmosférica 1
    FA2 <- -Ld-Lu/Ta            #Función atmosférica 2
    FA3 <- Ld                   #Función atmosférica 3
    gama <- (Tsen^2)/(b*Lsen)
    v2 <- Tsen-(Tsen^2)/b
    Ts <- gama*((1/E)*(FA1*Lsen+FA2)+FA3)+v2 #LST en °C
  return(Ts) }

#Temperatura Superficial Terrestre con Parámetro Vapor de agua 
LST_w <-function(Lsen,Tsen,E,LongOnda,W,b){
    nk1 <- 0.0009*LongOnda^3-0.01638*LongOnda^2+0.04745*LongOnda+0.27436
    ek1 <- 0.00032*LongOnda^3-0.06148*LongOnda^2+1.2021*LongOnda-6.2051
    xk1 <- 0.00986*LongOnda^3-0.23672*LongOnda^2+1.7133*LongOnda-3.2199
    ok1 <- -0.15431*LongOnda^3+5.2757*LongOnda^2-60.117*LongOnda+229.3139
    FA1W <-nk1*W^3+ek1*W^2+xk1*W+ok1        #Función atmosférica 1
    nk2 <- -0.02883*LongOnda^3+0.87181*LongOnda^2-8.82712*LongOnda+29.9092
    ek2 <- 0.13515*LongOnda^3-4.1171*LongOnda^2+41.8295*LongOnda-142.2782
    xk2 <- -0.22765*LongOnda^3+6.8606*LongOnda^2-69.2577*LongOnda+233.0722
    ok2 <- 0.41868*LongOnda^3-14.3299*LongOnda^2+163.6681*LongOnda-623.53
    FA2W <-nk2*W^3+ek2*W^2+xk2*W+ok2        #Función atmosférica 2
    nk3 <- 0.00182*LongOnda^3-0.04519*LongOnda^2+0.32652*LongOnda-0.6003
    ek3 <- -0.00744*LongOnda^3+0.11431*LongOnda^2+0.1756*LongOnda-5.4588
    xk3 <- -0.00269*LongOnda^3+0.31395*LongOnda^2-5.5916*LongOnda+27.9913
    ok3 <- -0.07972*LongOnda^3+2.8396*LongOnda^2-33.6843*LongOnda+132.9798
    FA3W <-nk3*W^3+ek3*W^2+xk3*W+ok3        #Función atmosférica 3
    gama <- (Tsen^2)/(b*Lsen)
    v2 <- Tsen-(Tsen^2)/b
    TsW <- gama*((1/E)*(FA1W*Lsen+FA2W)+FA3W)+v2 #LST en °C
    return(TsW)
}

#CONDICIÓN DE ACUERDO AL SATÉLITE SELECCIONADO Y AL TIPO DE ECUACIÓN LST  
  
#SATÉLITE LANDSAT 5

if (typeLandsat=="Landsat5"){

    setwd(Folder_correccion_DOS1) #Cargamos la tabla
    listabandas <- list.files(path = getwd(),pattern = "[1-7].TIF$",full.names = TRUE) #Cargamos bandas en lista
    #Importación bandas
    RTB1 <- raster(listabandas[1])
    RTB2 <- raster(listabandas[2])
    RTB3 <- raster(listabandas[3])
    RTB4 <- raster(listabandas[4])
    RTB5 <- raster(listabandas[5])
    RTB7 <- raster(listabandas[7]) 
    
    rasterCOMBITOTAL<-stack(RTB1,RTB2,RTB3,RTB4,RTB5,RTB7) #Combinacion de bandas total L5
    LayerTemp <- as(Layer, "Spatial") #Leemos el shapefile limite
    LimACRSB <- spTransform(LayerTemp,CRSobj =crs(rasterCOMBITOTAL)) #Conversion crs
    Combinacion_Bandas <- projectRaster(Corte(rasterCOMBITOTAL,LimACRSB),crs=crs(LayerTemp))

    #Generación índices NDVI, NDBI y BSI
    
    NDVI_crsRst<-NDVI_f(RTB4,RTB3,LimACRSB)
    NDVI<-projectRaster(NDVI_crsRst,crs=crs(LayerTemp))
    NDBI<-NDBI_f(RTB5,RTB4,LimACRSB,LayerTemp)
    BSI<-BSI_f(RTB5,RTB3,RTB4,RTB1,LimACRSB,LayerTemp)
    
    #Proporcion de Vegetacion
    
    PV <- ((NDVI_crsRst-min(getValues(NDVI_crsRst),na.rm = T))/(max(getValues(NDVI_crsRst),na.rm = T)-
                                                                min(getValues(NDVI_crsRst),na.rm = T)))^2
    EL57 <- 0.004*PV+0.986
    
    #Metadatos L5
    
    RADIANCE_MULT_BAND_6 = 5.5375E-02
    RADIANCE_ADD_BAND_6 = 1.18243
    K1_CONSTANT_BAND_6 = 607.76
    K2_CONSTANT_BAND_6 = 1260.56
    
    #Radiancia del sensor TOA cortado en limite
    
    Lsen<-Lsen_f(RADIANCE_MULT_BAND_6,RADIANCE_ADD_BAND_6,Banda_Termica,LimACRSB)
    
    #Temperatura Brillo sensor cortado en límite
    
    TBS <-TBS_f(K1_CONSTANT_BAND_6,K2_CONSTANT_BAND_6,Lsen,LimACRSB)
    Temperatura_de_brillo_en_sensor_del_satelite <- projectRaster(TBS-273.15,crs=crs(LayerTemp))
   
    #LST con Parámetros atmosféricos Ta,Lu,Ld
    
    if (typeparametro=="Atmosféricos") {
        LST_crsRst <-LST_L(Lsen,TBS,EL57,Transmisividad_Atmosferica,Radiancia_Atmosferica_Ascendente,
                           Radiancia_Atmosferica_Descendente,1256)
        LST <- projectRaster(LST_crsRst ,crs=crs(LayerTemp))
        
    #LST con vapor de agua W    
        
    } else if (typeparametro=="Vapor_de_Agua") { 
        LongOndaL5<-11.457  
        LST_crsRst <-LST_w(Lsen,TBS,EL57,LongOndaL5,Vapor_de_Agua,1256)
        LST<-projectRaster(LST_crsRst,crs=crs(LayerTemp))
    }
    Temperatura_Superficial_Terrestre <- LST-273.15
    Isla_de_Calor_Urbano <- (LST-mean(getValues(LST),na.rm=T))/mean(getValues(LST),na.rm=T)

#SATÉLITE LANDSAT 7
    
} else if (typeLandsat=="Landsat7") {

    setwd(Folder_correccion_DOS1) #Cargamos la tabla
    listabandas <- list.files(path = getwd(),pattern = "[1-7].TIF$",full.names = TRUE) #Cargamos bandas en lista
    #Importamos bandas
    RTB1 <- raster(listabandas[1])
    RTB2 <- raster(listabandas[2])
    RTB3 <- raster(listabandas[3])
    RTB4 <- raster(listabandas[4])
    RTB5 <- raster(listabandas[5])
    RTB7 <- raster(listabandas[8])
    
    rasterCOMBITOTAL<-stack(RTB1,RTB2,RTB3,RTB4,RTB5,RTB7) #Combinacion de bandas total L7
    LayerTemp <- as(Layer, "Spatial") #Leemos el shapefile limite
    LimACRSB <- spTransform(LayerTemp,CRSobj =crs(rasterCOMBITOTAL)) #Conversion crs
    Combinacion_Bandas <- projectRaster(Corte(rasterCOMBITOTAL,LimACRSB),crs=crs(LayerTemp))
    
    #Generación índices NDVI, NDBI y BSI

    NDVI_crsRst<-NDVI_f(RTB4,RTB3,LimACRSB)
    NDVI<-projectRaster(NDVI_crsRst,crs=crs(LayerTemp))
    NDBI<-NDBI_f(RTB5,RTB4,LimACRSB,LayerTemp)
    BSI<-BSI_f(RTB5,RTB3,RTB4,RTB1,LimACRSB,LayerTemp)
    
    #Proporcion de Vegetacion
    
    PV <- ((NDVI_crsRst-min(getValues(NDVI_crsRst),na.rm = T))/(max(getValues(NDVI_crsRst),na.rm = T)-
                                                                min(getValues(NDVI_crsRst),na.rm = T)))^2
    EL57 <- 0.004*PV+0.986

    #Metadatos L7
    
    RADIANCE_MULT_BAND_6_VCID_2 = 3.7205E-02
    RADIANCE_ADD_BAND_6_VCID_2 = 3.16280
    K1_CONSTANT_BAND_6_VCID_2 = 666.09
    K2_CONSTANT_BAND_6_VCID_2 = 1282.71

    #Radiancia del sensor TOA
    
    Lsen<-Lsen_f(RADIANCE_MULT_BAND_6_VCID_2,RADIANCE_ADD_BAND_6_VCID_2,Banda_Termica,LimACRSB)

    #Temperatura Brillo sensor
    
    TBS <-TBS_f(K1_CONSTANT_BAND_6_VCID_2,K2_CONSTANT_BAND_6_VCID_2,Lsen,LimACRSB)
    Temperatura_de_brillo_en_sensor_del_satelite <- projectRaster(TBS-273.15,crs=crs(LayerTemp))
       
    #LST con Parámetros atmosféricos Ta,Lu,Ld
    
    if (typeparametro=="Atmosféricos") {
        LST_crsRst<-LST_L(Lsen,TBS,EL57,Transmisividad_Atmosferica,Radiancia_Atmosferica_Ascendente,
                          Radiancia_Atmosferica_Descendente,1277)
        LST<-projectRaster(LST_crsRst,crs=crs(LayerTemp))
        
    #LST con vapor de agua W 
        
    } else if (typeparametro=="Vapor_de_Agua") {
        LongOndaL7<-11.267    
        LST_crsRst <-LST_w(Lsen,TBS,EL57,LongOndaL7,Vapor_de_Agua,1277)
        LST<-projectRaster(LST_crsRst,crs=crs(LayerTemp))
    }
    Temperatura_Superficial_Terrestre <- LST-273.15
    Isla_de_Calor_Urbano <- (LST-mean(getValues(LST),na.rm=T))/mean(getValues(LST),na.rm=T)

#SATÉLITE LANDSAT 8
    
} else if (typeLandsat=="Landsat8"){

    setwd(Folder_correccion_DOS1) #Cargamos la tabla
    listabandas <- list.files(path = getwd(),pattern = "[1-7].TIF$",full.names = TRUE) #Cargamos bandas en lista
    #Importamos bandas
    RTB1 <- raster(listabandas[1])
    RTB2 <- raster(listabandas[3])
    RTB3 <- raster(listabandas[4])
    RTB4 <- raster(listabandas[5])
    RTB5 <- raster(listabandas[6])
    RTB6 <- raster(listabandas[7])
    RTB7 <- raster(listabandas[8])
    rasterCOMBITOTAL<-stack(RTB1,RTB2,RTB3,RTB4,RTB5,RTB6,RTB7) #Combinacion de bandas total L8
    LayerTemp <- as(Layer, "Spatial") #Leemos el shapefile limite
    LimACRSB <- spTransform(LayerTemp,CRSobj =crs(rasterCOMBITOTAL)) #Conversion crs
    Combinacion_Bandas <- projectRaster(Corte(rasterCOMBITOTAL,LimACRSB),crs=crs(LayerTemp))

    #Metadatos L8
    
    RADIANCE_MULT_BAND_10 = 3.3420E-04 
    RADIANCE_ADD_BAND_10 = 0.10000 
    K1_CONSTANT_BAND_10 = 774.8853
    K2_CONSTANT_BAND_10 = 1321.0789
    
    #Generación índices NDVI, NDBI y BSI

    NDVI_crsRst<-NDVI_f(RTB5,RTB4,LimACRSB)
    NDVI<-projectRaster(NDVI_crsRst,crs=crs(LayerTemp))
    NDBI<-NDBI_f(RTB6,RTB5,LimACRSB,LayerTemp)
    BSI<-BSI_f(RTB6,RTB4,RTB5,RTB2,LimACRSB,LayerTemp)
    
    #Proporcion de Vegetacion
    
    PV <- ((NDVI_crsRst-min(getValues(NDVI_crsRst),na.rm = T))/(max(getValues(NDVI_crsRst),na.rm = T)-
                                                                min(getValues(NDVI_crsRst),na.rm = T)))^2
    EL8 <- 0.0015*PV+0.9885

    #Radiancia del sensor TOA
    
    Lsen<-Lsen_f(RADIANCE_MULT_BAND_10,RADIANCE_ADD_BAND_10,Banda_Termica,LimACRSB)

    #Temperatura Brillo Sensor
    
    TBS <-TBS_f( K1_CONSTANT_BAND_10,K2_CONSTANT_BAND_10,Lsen,LimACRSB)
    Temperatura_de_brillo_en_sensor_del_satelite <- projectRaster(TBS-273.15,crs=crs(LayerTemp))
    
    #LST con Parámetros atmosféricos Ta,Lu,Ld
    
    if (typeparametro=="Atmosféricos") {
        LST_crsRst <-LST_L(Lsen,TBS,EL8,Transmisividad_Atmosferica,Radiancia_Atmosferica_Ascendente,
                           Radiancia_Atmosferica_Descendente,1324)
        LST<-projectRaster(LST_crsRst,crs=crs(LayerTemp))
        
    #LST con vapor de agua W 
        
    } else if (typeparametro=="Vapor_de_Agua") {
        LongOndaL8<-10.9    
        LST_crsRst<-LST_w(Lsen,TBS,EL8,LongOndaL8,Vapor_de_Agua,1324)
        LST<-projectRaster(LST_crsRst,crs=crs(LayerTemp))
    }
    Temperatura_Superficial_Terrestre <- LST-273.15
    Isla_de_Calor_Urbano <- (LST-mean(getValues(LST),na.rm=T))/mean(getValues(LST),na.rm=T)
}
