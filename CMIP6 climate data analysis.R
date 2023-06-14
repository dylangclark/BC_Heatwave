### BC Project Climate Data Pull
## 2022, December
# (C) Dylan G. Clark

##### Steps
# 1) Pull CMIP6 models (5) SSP85 for Canada (daily minT and maxT)
# 2) Cut CMIP6 daily projections for each population centre
# 3) Create daily timeseries for each GCM and population centre
# 4) Develop monthly correction constance for each population centre
    # a) pull PRISM monthly 1981 to 2010 minT and maxT data
    # b) pull monthly climate projections from CCCS 1981 to 2010
    # c) cut the monthly Prism and CMIP6 data for each population centre and for each month





####
############################################################ Download data

#HadGEM GCM: gregorian "https://data.pacificclimate.org/data/downscaled_cmip6/tasmin_day_BCCAQv2+ANUSPLIN300_CanESM5_historical+ssp126_r1i1p2f1_gn_19500101-21001231.nc.nc?tasmin[21900:22264][84:352][0:328]"

#EC-Earth GCM: fixed 360 "https://data.pacificclimate.org/data/downscaled_cmip6/tasmax_day_BCCAQv2+ANUSPLIN300_EC-Earth3_historical+ssp585_r4i1p1f1_gr_19500101-21001231.nc.nc?tasmax[18262:55151][0:510][0:1068]"

#FGoals-G3: fixed 365 "https://data.pacificclimate.org/data/downscaled_cmip6/tasmax_day_BCCAQv2+ANUSPLIN300_FGOALS-g3_historical+ssp585_r1i1p1f1_gn_19500101-21001231.nc.nc?tasmax[18250:55113][0:510][0:1068]"

#UKESM1-0-L: fixed 360 "https://data.pacificclimate.org/data/downscaled_cmip6/tasmax_day_BCCAQv2+ANUSPLIN300_UKESM1-0-LL_historical+ssp585_r1i1p1f2_gn_19500101-21001230.nc.nc?tasmax[18000:54359][0:510][0:1068]"

#ACCESS ESM1: gregorian "https://data.pacificclimate.org/data/downscaled_cmip6/tasmax_day_BCCAQv2+ANUSPLIN300_ACCESS-ESM1-5_historical+ssp585_r1i1p1f1_gn_19500101-21001231.nc.nc?tasmax[18262:55151][0:510][0:1068]"


######## Download function Function
websearch<- function (URLbase){
  

  library(rvest)
  library(plyr)
  library(stringr)
  
  options(timeout=1000)
  cnt1<-0

    t1<-seq(18262,55151,by=360)
    cnt2<-0
    for(i in t1){
      t2<-i+359
      cnt2<-cnt2+1
      u<-paste0("https://data.pacificclimate.org/data/downscaled_cmip6/tasmax_day_BCCAQv2+ANUSPLIN300_",URLbase,"_historical+ssp585_r1i1p1f1_gn_19500101-21001231.nc.nc?tasmax[",i,":",t2,"][84:352][0:328]")
      if(cnt2>1.1){U<-c(U,u)}else{U<-u}
    }

  
  for (i in U){
      cnt1<-cnt1+1
      #retry::wait_until(Sys.time() - timeStamp > 120)
      
      name<-as.character(stringr::str_match(i,pattern="(?<=downscaled_cmip6/)(.*?)(?=.nc.nc)")[,1])
      fileDes<-paste0("D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/tasmax/",name,"_",cnt1,".nc")
      tryCatch({
        download.file(url=i,destfile =fileDes,mode="wb")
        print("waiting for download")
      },error =function(e){print("couldn't find file")
        Sys.sleep(30)
      },warning=function(e){print("couldn't find file")}
      )
    
  }
}

######## Script

websearch(URLbase = "ACCESS-ESM1-5")
    




############################################################ NetCDF Analysis 
############################################################ 



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
library(readr)
library(RStoolbox)
library(exactextractr)
library(gdata)
library(epwshiftr)
library(reshape2)
library(svMisc)
library(hdf5r)
library(parallel)
library(profvis)
library(doParallel)


          # read shapefile polygon for reagional boundaries
          boundaryFolder<-"D:/GIS/BC/"
          boundaryName<-"BCPopCentresDigitalVan2"
          
          
          Buffers<-readOGR(dsn=paste0(boundaryFolder,boundaryName,".shp"))
          
          #List all MinT NCDF files in master folder:
          minTPaths<-list.files(path="D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/tasmin/",full.names = T, pattern=".nc",recursive = T)
          write.csv(minTPaths,"D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/minTPaths.csv")
          minTPaths<-read.csv("D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/minTPaths.csv")
          minTPaths<-as.list(minTPaths[,2])
          
          #List all MaxT NCDF files in master folder:
          maxTPaths<-list.files(path="D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/tasmax/",full.names = T, pattern=".nc",recursive = T)
          write.csv(maxTPaths,"D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/maxTPaths.csv")
          maxTPaths<-read.csv("D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/maxTPaths.csv")
          maxTPaths<-as.list(maxTPaths[,2])

          #dirPath<-maxTPaths[7]
          
          

GCMExtract(dirPath = maxTPaths,cutShp = Buffers,variName = "tasmax")


GCMExtract<-function(dirPath,cutShp,variName){
  
  
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
  library(readr)
  library(RStoolbox)
  library(exactextractr)
  library(gdata)
  library(epwshiftr)
  library(reshape2)
  library(svMisc)
  library(hdf5r)
  library(parallel)
  library(profvis)
  library(doParallel)
  
  
  pre1.brick<-brick("D:/GIS/Climate Models/CMIP6/BC/SelectedGCMs/CRS Template.nc")
  
  
  Buffers<-st_as_sf(sp::spTransform(cutShp, crs(pre1.brick[[1]])))
  Region<-as.list(Buffers$PCNAME)
  
  #Create a dataframe with columns: lat, long, buffer, day, rcp, gcm
  cnt<-0
  cnt2<-1
  
  for (i in dirPath){
    
      GCM<-as.character(str_match(i,pattern="(?<=ANUSPLIN300)(.*?)(?=_historical)")[,1])
      RCP<-as.character("ssp585")
    
    fpath<-paste0("D:/R/BC/MasterCMIP6",GCM,"_",variName,".csv")
    
    cnt<-cnt+1
    
    if(cnt2==0){
      
      Master<-try(silent=T,
      expr = {
        Master<-read.csv(fpath)
        Master<-Master[,2:6]
      })
      print("read csv")
    }else{print("didn't read csv")}
    
    cnt2<-cnt2+1
    pre1.brick<-brick(as.character(i))
    layers<-nlayers(pre1.brick)
    pre1.brick<-crop(pre1.brick,extent(Buffers))
    
    
    cores=detectCores()
    my.cluster <- parallel::makeCluster(cores[1]-1, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
      
     M<-foreach (y=1:layers, .combine='rbind', .packages = c("exactextractr","sf","dplyr","rgdal","raster")) %dopar%{
      
      pre2.brick<-pre1.brick[[y]]
      clip2<-mask(pre2.brick,Buffers)
      clip3 <- exactextractr::exact_extract(pre2.brick, Buffers,'mean',progress=F,append_cols=T)
      Date<-as.character(getZ(pre2.brick))
      
      val<-as.numeric(clip3$mean)
      place<-as.character(clip3$PCNAME)
      Observe<-cbind(place,val)
      Observe<-as.data.frame(cbind(RCP,GCM,Date,Observe))
      
      
     }
      #close iterate through day
    
    if(exists("Master")){Master<-rbind(Master,M)
    }else{ print("My first time")
      Master<-M
    }
    print(paste0("Done with ", cnt2, " and ", RCP))
    fileP<-paste0("D:/R/BC/MasterCMIP6",GCM,"_",variName,".csv")
    
    if(cnt %in% c(2,20,39,98,99,100,101,102,108,109,110,111)){
        write.csv(Master,fileP)
        remove(Master)
        cnt2<-0
    }else{print("still rolling")}
    stopCluster(my.cluster)      
    gc()
    
    
    
  }
}







###################################### Merge files from GCMDownExtract and add constants (PRISM offsets)
##

GCMList<-list.files(path="D:/R/BC/",full.names = T, recursive = T,include.dirs = F)

tasminPRISM<-read.csv("D:/R/BC/tasmin difference.csv")
tasmaxPRISM<-read.csv("D:/R/BC/tasmax difference.csv")

  ###Make sure this is maxT
A<-read.csv(GCMList[502])
A<-A[,4:6]
A$Date<-as.Date(strptime(A$Date, format="%Y-%m-%d"))
A$val<-as.numeric(A$val)
gcm<-"ACCESS"
A$place<-as.factor(A$place)

A<-A%>%
  dplyr::group_by(place, Date) %>%
  dplyr::summarise(val=mean(val))

  ###Make sure this is minT
B<-read.csv(GCMList[503])
B<-B[,4:6]
B$Date<-as.Date(strptime(B$Date, format="%Y-%m-%d"))
B$val<-as.numeric(B$val)
B$place<-as.factor(B$place)

B<-B%>%
  dplyr::group_by(place, Date) %>%
  dplyr::summarise(val=mean(val))



C<-merge(A,B,by=c("Date","place"))
    #remove(A,B)
C$GCM<-gcm
colnames(C)<-c("Date","Place","maxT","minT","GCM")
write.csv(C,"D:/R/BC/CMIP6/ACCESS.csv")

C$Date<-as.Date(strptime(C$Date, format="%Y-%m-%d"))
C$Month<-factor(months.Date(C$Date))
C$Year<-as.numeric(format(C$Date,'%Y'))

C$Month<-match(C$Month, month.name)

Places<-as.character(unique(C$Place))

cntrP<-0
for(p in Places){
  temp1<-C %>%
    filter(Place == p)
  
  tminCorrection<-tasminPRISM %>%
      filter(PC==p)
    
  tmaxCorrection<-tasmaxPRISM %>%
      filter(PC==p)
  cntr<-0
        
      for(m in seq(1:12)){
          cntr<-cntr+1
          temp2<-temp1 %>%
            filter(Month==m)
          
          temp2$tasmin<-with(temp2,minT+(tminCorrection[,(m+1)]))
          temp2$tasmax<-with(temp2,maxT+(tmaxCorrection[,(m+1)]))
        
          if(cntr>1.1){new<-rbind(new,temp2)}else{new<-temp2}
      }
  cntrP<-cntrP+1
  
  if(cntrP>1.1){allDF<-rbind(allDF,new)}else{allDF<-new}


}




allDF2<-allDF[,-c(3,4)]



PopCentres<-read.csv("D:/R/BC/Health/HealthModel/HealthModel/data/PopCentres.csv",header=T)
allDF2 <- (merge(x=PopCentres, y=allDF2, by.x = "PCNAME",by.y="Place"))


fn<-paste0("D:/R/BC/CMIP6/ACCESS.csv")
write.csv(allDF2,fn)


############################################ PRISM correction analysis
#################################################

## Steps:
# 1) load buffers for PC (and vancouver splits)
# 2) list of all files (prism and climate baseline)
# 3) cut monthly timeseries (1981 to 2010)
# 4) calculate difference between PRISM and climate baseline results for each month, year, and PC
# 5) aggregate difference by months for each PC



# read shapefile polygon for reagional boundaries
boundaryFolder<-"D:/GIS/BC/"
boundaryName<-"BCPopCentresDigitalVan2"


Buffers<-readOGR(dsn=paste0(boundaryFolder,boundaryName,".shp"))

#List all files in master folder:
    #PrismPaths<-list.files(path="D:/GIS/Climate Models/PRISM/",full.names = T, pattern=".nc",recursive = T)
    #write.csv(PrismPaths,"D:/GIS/Climate Models/CMIP6/BC/PRISM.csv")
PrismPaths<-read.csv("D:/GIS/Climate Models/CMIP6/BC/PRISM.csv")
PrismPaths<-PrismPaths[,-1]


PRISMCalculator(dirPath = PrismPaths,cutShp = Buffers)


PRISMCalculator<-function(dirPath,cutShp){
  
  
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
  library(readr)
  library(RStoolbox)
  library(exactextractr)
  library(gdata)
  library(epwshiftr)
  library(reshape2)
  library(svMisc)
  library(hdf5r)
  library(parallel)
  library(profvis)
  library(doParallel)
  
  
  pre1.brick<-brick(as.character(dirPath[1,1]))
  
  
  Buffers<-st_as_sf(sp::spTransform(cutShp, crs(pre1.brick[[1]])))
  Region<-as.list(Buffers$PCNAME)
  
  #Create a dataframe with columns: lat, long, buffer, day, rcp, gcm
  
  
  
  for (i in 1:nrow(dirPath)){
    
    variName<-dirPath[i,2]
    DataSet<-dirPath[i,3]
    RCP<-as.character("ssp585")
    
    fpath<-paste0("D:/R/BC/MasterPRISM.csv")
    
    Master<-try(silent=T,
                expr = {
                  Master<-read.csv(fpath)
                  Master<-Master[,2:7]
                }
    )
    
    pre1.brick<-brick(as.character(dirPath[i,1]))
    layers<-nlayers(pre1.brick)
    pre1.brick<-crop(pre1.brick,extent(Buffers))
    
    
    cores=detectCores()
    my.cluster <- parallel::makeCluster(cores[1]-1, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
    cntr<-0
    
    M<-foreach (y=1:layers, .combine='rbind', .packages = c("exactextractr","sf","dplyr","rgdal","raster")) %dopar%{
      
      pre2.brick<-pre1.brick[[y]]
      clip2<-mask(pre2.brick,Buffers)
      clip3 <- exactextractr::exact_extract(pre2.brick, Buffers,'mean',progress=F,append_cols=T)
      Date<-as.character(getZ(pre2.brick))
      
      val<-as.numeric(clip3$mean)
      place<-as.character(clip3$PCNAME)
      Observe<-cbind(place,val)
      Observe<-as.data.frame(cbind(RCP,DataSet,variName,Date,Observe))
      
      print(paste0("working on new layer of ",DataSet," - time elapsed ",cntr))
      cntr<-cntr+1
      if(cntr>1.1){M<-rbind(M,Observe)}else{M<-Observe}
      
    }
    #close iterate through day
    
    if(exists("Master")){Master<-rbind(Master,M)
    }else{ print("My first time")
      Master<-M
    }
    print(paste0("Done with ", cntr," for ", DataSet, " and ", RCP))
    fileP<-paste0("D:/R/BC/MasterPRISM.csv")
    write.csv(Master,fileP)
    
    remove(Master)
    stopCluster(my.cluster)      
    gc()
    
    
    
  }
}


############### Data prep for shiny



setwd("D:/R/BC/Health/HealthModel/HealthModel/")

GCM<-"ACCESS"

        if(GCM=="EC-Earth3"){
          DB<-read.csv("D:/R/BC/Health/HealthModel/EC-Earth.csv",header=T)
          DB<-DB[,c(-1,-3,-4,-5,-6,-7,-8,-15)]
          DB$Date<-as.Date(strptime(DB$Date, format="%Y-%m-%d"))
        }else if (GCM=="HadGEM3"){
          DB<-read.csv("D:/R/BC/Health/HealthModel/HadGEM3.csv",header=T)
          DB<-DB[,c(-1,-3,-4,-5,-6,-7,-8,-15)]
          DB$Date<-as.Date(strptime(DB$Date, format="%Y-%m-%d"))
        }else if (GCM=="FGOALS"){
          DB<-read.csv("D:/R/BC/Health/HealthModel/FGOALS.csv",header=T)
          DB<-DB[,c(-1,-3,-4,-5,-6,-7,-8,-15)]
          DB$Date<-as.Date(strptime(DB$Date, format="%Y-%m-%d"))
        }else if (GCM=="UKESM1"){
          DB<-read.csv("D:/R/BC/Health/HealthModel/UKESM1.csv",header=T)
          DB<-DB[,c(-1,-3,-4,-5,-6,-7,-8,-15)]
          DB$Date<-as.Date(strptime(DB$Date, format="%Y-%m-%d"))
        }else if (GCM=="ACCESS"){
          DB<-read.csv("D:/R/BC/Health/HealthModel/ACCESS.csv",header=T)
          DB<-DB[,c(-1,-3,-4,-5,-6,-7,-8,-15)]
          DB$Date<-as.Date(strptime(DB$Date, format="%Y-%m-%d"))
        }

        ## Convert to Health Regions 
        DBScale<-aggregate(DB[,c(10,11)],by=list(DB$Date,DB$Year,DB$Month,DB$HA_Name,DB$Warning_Zone),mean,na.rm=T)
        colnames(DBScale)<-c("Date","year","Month","HA_Name","warningCat","tasmin","tasmax")
        DBScale[,4]<-as.factor(DBScale[,4])
        n<-paste0("data/GCMData/",GCM,"_HR.csv")
        write.csv(DBScale,n)
        
        ## Convert to PC
        DBScale<-DB[,c(7,9,8,1,6,10,11)]
        colnames(DBScale)<-c("Date","year","Month","Place","warningCat","tasmin","tasmax")
        n<-paste0("data/GCMData/",GCM,"_PC.csv")
        write.csv(DBScale,n)
        
        ## Convert to province        
        DB$Province<-"British Columbia"
        DBScale<-aggregate(DB[,c(10,11)],by=list(DB$Date,DB$Year,DB$Month,DB$Province,DB$Warning_Zone),mean,na.rm=T)
        colnames(DBScale)<-c("Date","year","Month","Province","warningCat","tasmin","tasmax")
        n<-paste0("data/GCMData/",GCM,"_BC.csv")
        write.csv(DBScale,n)
        
        # Convert for Excel reading
        DBScale<-DB[,c(7,9,8,1,3,6,10,11)]
        colnames(DBScale)<-c("Date","year","Month","Place","Health_Region","warningCat","tasmin","tasmax")
        
        for(i in unique(DBScale$Health_Region)){
          db<-DBScale %>%
            filter(Health_Region==i)
            
          n<-paste0("data/GCMData/Excel_",GCM,"_",i,".csv")
          write.csv(db,n)
        }
        
########################### Ensemble creation
#####
        ECEarth<-read.csv("D:/R/BC/Health/HealthModel/EC-Earth.csv",header=T)
            GCM<-unique(as.factor(ECEarth$GCM))    
            ECEarth<-ECEarth[,c(14,17,16,2,10,13,18,19)]
            tasmin<-paste0(GCM,"_tasmin")
            tasmax<-paste0(GCM,"_tasmax")
            colnames(ECEarth)<-c("Date","year","Month","Place","Health_Region","warningCat",tasmin,tasmax)
        HadGEM<-read.csv("D:/R/BC/Health/HealthModel/HadGEM3.csv",header=T)
            GCM<-unique(as.factor(HadGEM$GCM))    
            HadGEM<-HadGEM[,c(14,2,18,19)]
            tasmin<-paste0(GCM,"_tasmin")
            tasmax<-paste0(GCM,"_tasmax")
            colnames(HadGEM)<-c("Date","Place",tasmin,tasmax)
        FGOALS<-read.csv("D:/R/BC/Health/HealthModel/FGOALS.csv",header=T)
            GCM<-unique(as.factor(FGOALS$GCM))    
            FGOALS<-FGOALS[,c(14,2,18,19)]
            tasmin<-paste0(GCM,"_tasmin")
            tasmax<-paste0(GCM,"_tasmax")
            colnames(FGOALS)<-c("Date","Place",tasmin,tasmax)
        UKESM<-read.csv("D:/R/BC/Health/HealthModel/UKESM1.csv",header=T)
            GCM<-unique(as.factor(UKESM$GCM))    
            UKESM<-UKESM[,c(14,2,18,19)]
            tasmin<-paste0(GCM,"_tasmin")
            tasmax<-paste0(GCM,"_tasmax")
            colnames(UKESM)<-c("Date","Place",tasmin,tasmax)
        ACCESS<-read.csv("D:/R/BC/Health/HealthModel/ACCESS.csv",header=T)
            GCM<-unique(as.factor(ACCESS$GCM))    
            ACCESS<-ACCESS[,c(14,2,18,19)]
            tasmin<-paste0(GCM,"_tasmin")
            tasmax<-paste0(GCM,"_tasmax")
            colnames(ACCESS)<-c("Date","Place",tasmin,tasmax)
      
            
    Ensemble<-merge(ECEarth,HadGEM,by=c("Date","Place"))
    Ensemble<-merge(Ensemble,FGOALS,by=c("Date","Place"))
    Ensemble<-merge(Ensemble,UKESM,by=c("Date","Place"))
    Ensemble<-merge(Ensemble,ACCESS,by=c("Date","Place"))
    
    write.csv(Ensemble,"D:/R/BC/Health/HealthModel/Ensemble.csv")
        DB<-read.csv("D:/R/BC/Health/HealthModel/Ensemble.csv",header=T)
        
        setwd("D:/R/BC/Health/HealthModel/HealthModel/")
        GCM<-"Ensemble"
        
        DB<-DB[,-1]
          cn<-colnames(DB)
          cn<-cn[7:16]
          DBScale<-aggregate(DB[,c(7:16)],by=list(DB$Date,DB$year,DB$Month,DB$Health_Region,DB$warningCat),mean,na.rm=T)
          colnames(DBScale)<-c("Date","year","Month","HA_Name","warningCat",cn)
          DBScale[,4]<-as.factor(DBScale[,4])
          n<-paste0("data/GCMData/",GCM,"_HR.csv")
          write.csv(DBScale,n)
          
          ## Convert to PC
          DBScale<-DB[,c(1,3,4,2,6,7:16)]
          colnames(DBScale)<-c("Date","year","Month","Place","warningCat",cn)
          n<-paste0("data/GCMData/",GCM,"_PC.csv")
          write.csv(DBScale,n)
          
          ## Convert to province        
          DB$Province<-"British Columbia"
          DBScale<-aggregate(DB[,c(7:16)],by=list(DB$Date,DB$year,DB$Month,DB$Province,DB$warningCat),mean,na.rm=T)
          colnames(DBScale)<-c("Date","year","Month","Province","warningCat",cn)
          n<-paste0("data/GCMData/",GCM,"_BC.csv")
          write.csv(DBScale,n)
          
          # Convert for Excel reading
          DBScale<-DB[,c(1,3,4,2,5,6,7:16)]
          colnames(DBScale)<-c("Date","year","Month","Place","Health_Region","warningCat",cn)
          
          for(i in unique(DBScale$Health_Region)){
            db<-DBScale %>%
              filter(Health_Region==i)
            
            n<-paste0("data/GCMData/Excel_",GCM,"_",i,".csv")
            write.csv(db,n)
          }
        
#####        
############################################ Heat warnings
#################################################


PC<-read.csv("data/PopCentres.csv", header=T)
Warnings<-read.csv("data/WarningCriteria.csv", header=T)


HeatWarning<-function(TS,warnCriteria,PopLookup){
  library(dplyr)
  
  P<-unique(TS$Place)
  TS$Date<-as.Date(strptime(TS$Date, format="%Y-%m-%d"))
  
  for(i in P){
    PC<-PopLookup %>%
      dplyr::filter(PCNAME==P)
    
    wz<-unique(PC$Warning_Zone)
    
    TS<-TS %>%
      dplyr::filter(Place==P)
    
    w<-warnCriteria %>%
      dplyr::filter(Warning_Zone==wz)
    
    minT<-w$Min.temp
    maxT<-w$Max.temp
    
    startD<-min(TS$Date)
    endD<-max(TS$Date)
    rcp<-unique(TS$RCP)
    gcm<-unique(TS$GCM)
    
    for (r in rcp){
      TS %>%
        dplyr::filter(RCP==r)
      print(r)
      for (g in gcm){
        TS %>%
          dplyr::filter(GCM==g)
        
        TS<-TS[order(as.Date(TS$Date, format="%Y-%m-%d")),]
        Lmax<-which(TS$maxT > maxT)
        Lmin<-which(TS$minT > minT)
        print(length(Lmax))
        
        for(d in Lmax){
          d2<-TS[(d+1),7]
          d2night<-TS[d+1,8]
          
          if(d2>maxT & d2night>minT){
            TS[d,9]<-1
            print("found a heatwave")
            
          }else{
            TS[d,9]<-0
          }
          
        }
        
      }
      
    }
  }
  
  
}

HeatWarning(TS=CCSM4, warnCriteria = warnCriteria,PopLookup = PopLookup)



############################################ Check which months can exclude from analysis (Jan 18)
#################################################

DB<-read.csv("D:/R/BC/Health/HealthModel/Ensemble.csv",header=T)
DB<-DB[,-1]
EcoRegions<-read.csv("D:/R/BC/Health/HealthModel/HealthModel/data/PC_EcoZoneMatch.csv",header=T)
EcoRegions<-EcoRegions[,c(1,6)]
colnames(EcoRegions)<-c("Place","Eco_region")

DB2<-merge(x=DB,y=EcoRegions,by="Place")

DB2<-DB2 %>% 
  mutate(flagLow=ifelse(EC.Earth_tasmin<10 & HadGEM3_tasmin<10 & FGOALS_tasmin<10 & UKESM1_tasmin<10 & ACCESS_tasmin<10,1,0))

DB3<-aggregate(DB2[,18],by=list(DB2$Date,DB2$year,DB2$Month),mean,na.rm=T)

write.csv(DB3,"D:/R/BC/Health/HealthModel/MonthExclusion.csv")
