
#####Lasso climate model analysis Canada for BC heatwave project (c) 2022 D. G. Clark

#This script is used to: 1) pull all 19 statistically downscaled GCMs into R (netCDF format); 2) cut the GCM into a region (select region option);
# 3) calculate region SU30 and GSL indicies; 4) plot the SU30 and GSL for each period; 5) create a convex hull to select the outermost points on each plot


#Statistically downscaled GCMs have been downloaded from Environment and Climate Change Canada http://climate-scenarios.canada.ca/?page=bccaqv2-data

##Load libraries
library(RNetCDF)
library(raster)
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(rgdal)
library(dplyr)
library(rgeos)
library(purrr)
library(ggplot2)
library(stringr)
library(sf)
library(ncdf4.helpers)
library(PCICt)
library(reshape2)


######## Pull all 19 netCDFs into R




# read netCDFs [needs to be a loop]
gcmPath<-"D:/GIS/Climate Models/CMIP6/BC/"

# read shapefile polygon for reagional boundaries
boundaryFolder<-"D:/GIS/BC/"
boundaryName<-"BCPopCentresDigitalVan"


####set netcdf variable
varName<-"gslETCCDI"
folder<-"index_gsl/"
#OR
varName<-"su30ETCCDI"
folder<-"index_su30/"
##


#######Function. Install this function and the one below.

## Remember to change the varName file between prcp and temp.

ClimatePlots<-function(folder,varName,boundaryFolder,boundaryName,gcmPath){
  
  ###create table of all statistically downscaled GCMs of the same variable type and note the classes
  oldST<-Sys.time()
  pt<-c("Prince George")
  counter<-0
  cutShp<-readOGR(dsn=boundaryFolder, layer = boundaryName)
  
  gcmTable<-GCMTableFun(gcmPath=gcmPath, var=folder)
  print("wrote GCMTable!")
  

    print(paste0("working on ",pt))
    ###Read shapefile that will be used to cut GCMs for regional analysis
    cutShp<-subset(cutShp, PCNAME==pt)
    cutShp<-gSimplify(cutShp,tol=0.001)
    cutShp<-st_as_sf(spTransform(cutShp,CRS("+init=EPSG:4326")))
    
    ###Iterate through GCMs
    n<-1:nrow(gcmTable)
    for (i in n){
      counter<- counter + 1
      y<-gcmTable[i,]
      gcm<-y[1]
      ssp<-y[2]
      GCMpath<-as.character(paste0(gcmPath,folder,(y[4])))
      
        baseNc<-nc_open(as.character(GCMpath),verbose=F,write=F)
        cal <- as.character(baseNc$dim$time$calendar)
        orig<-as.character(baseNc$dim$time$units)
        tas_time <- nc.get.time.series(baseNc, v="dim",time.dim.name="time")
        tas_time <- as.PCICt(x=tas_time, origin=orig,cal=cal)
        z<-ncvar_get(baseNc,varName)
        layers<-dim(z)[3]
        lon<-ncvar_get(baseNc,"lon")
        lat<-ncvar_get(baseNc,"lat")
        
        for (y in 1:layers){
              print(paste0("working on new layer of ",gcm," - time elapsed ",Sys.time()-oldST))
              
              z1<-z[,,y]
              colnames(z1)<-lat
              rownames(z1)<-lon
              z2<-melt(z1)
              colnames(z2)[3]<-"Value_mean"
              zDate<-as.character(tas_time[y],format="%Y-%m-%d")
              
              ## cut to regions
              newshp<-sf::st_as_sf(z2,coords=c("Var1","Var2"),crs=CRS("+init=EPSG:4326"))
              z2Cut<-as.data.frame(sf::st_join(newshp,cutShp,left=F))
              
              z3<-mean(z2Cut[,1],na.rm=T)
              print(z3)
              z3<-cbind(gcm,ssp,zDate,varName,z3)
              
              if(exists("z3bind")){z3bind<-rbind(z3bind,z3)}else{z3bind<-z3}
              
              remove(z3)        
       
              oldST<-Sys.time()
        }
        
        
        if(exists("zbefore")){zbefore<-rbind(zbefore,z3bind)}else{zbefore<-z3bind}
        saveName<-paste0(gcmPath,"-",varName,".csv")
        write.csv(zbefore,saveName)
        
        remove(z1,z2,zDate,z,baseNc,z3bind,layers)



    }

}


####################GCM Table function
GCMTableFun<-function(gcmPath,var){
  p<-paste0(gcmPath,var)
  gcmname<-list.files(p)
  for(i in gcmname){
    var<-as.character(str_match(i,pattern="(?<=bccaqv2r_)(.*?)(?=_ssp)")[,1])
    ssp<-as.character(str_match(i,pattern="(?<=_ssp)(.*?)(?=_)")[,1])
    gcm<-as.character(str_match(i,pattern="(?<=ssp585_)(.*?)(?=_141)")[,1])
    
    
    
    tbl<-cbind(gcm,ssp,var,i)
      
    if(exists("GCMtab")){GCMtab<-rbind(tbl,GCMtab)}else{GCMtab<-tbl}
      
  }
  return(GCMtab)
}


ClimatePlots(folder=folder,varName = varName,boundaryFolder = boundaryFolder,boundaryName = boundaryName,gcmPath = gcmPath)


#######################Function for file lookup




#####Analysis of the data

Data<-read.csv(file="D:/GIS/Climate Models/CMIP6/BC/VancouverConvex.csv",header=T)

ERA<-levels(as.factor(Data$era))

 for (i in ERA){
        
        pull<-Data %>% 
          dplyr::select(gcm, era, Temp.su30, Temp.gsl) %>%
          filter(era==i)
        
        pull<-pull %>%
          dplyr::group_by(gcm) %>%
          summarise(across(where(is.numeric), ~mean(.x,na.rm=T)))
        
          pull$x<-pull$Temp.su30
          pull$y<-pull$Temp.gsl
          index<-chull(pull$x,pull$y)
          boundaryGCM<-(paste(as.character(gcm[index]), sep=", ", collapse=", "))
          Title<-paste0("ssp585 for Vancouver", i)
          ###Save output
          save<-cbind(i,boundaryGCM)
          if(exists("chullOut")){chullOut<-rbind(save,chullOut)}else{chullOut<-save}
          
          
          ###Draw diagram
          rownames(pull)<-pull$gcm
          xAxis<-"days over 30C"
          yAxis<-"growing degree days"
          hpts <- chull(pull$x, pull$y)
          hpts <- c(hpts, hpts[1])
          tempPlot<-ggplot(d=pull, aes(x,y))+ geom_point(shape=1,color="blue") +
            geom_text(aes(x,y,label=rownames(pull)),
                      position = position_dodge(width=0.5), size=2.5)+
            geom_path(data=pull[hpts,], color="blue")+
            labs(title = Title, x=xAxis,y=yAxis)
        
          ggsave(tempPlot, file=paste0("D:/GIS/Climate Models/CMIP6/BC/plots/", Title,".png"),
                 width = 14, height = 10, units = "cm")
          
}
  


write.csv(chullOut,file="D:/GIS/Climate Models/CMIP6/BC/plots/BC.csv")

remove(chullOut)







