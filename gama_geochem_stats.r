#################################################################### 
#
# process and analyze geochemistry data from California's
# Groundwater Ambient and Monitoring Assessment (GAMA) program
#
# by Walt McNab, Ph.D.
#
# functionality includes quantifying relationships between parameters
# including principal components analysis and analysis of spatial
# correlation
#
####################################################################


### environment ###

setwd("D:/GAMA")
require(compiler)
enableJIT(3)
library(reshape)
library(corrplot)
library(ncf)


### supporting functions ###


Definitions <- function(){
  # assign general definitions in a function ###
	analytes <- c("CA","CL","K","MG","MN","NO3N","SO4","SODIUM")
	attribs <- c("WELL.NAME","APPROXIMATE.LATITUDE","APPROXIMATE.LONGITUDE", "DATE", "COUNTY", "CHEMICAL","CONC", "DATASET_CAT")
	well_features <- attribs[(attribs!="DATE") & (attribs!="CHEMICAL") & (attribs!="CONC")]
	sources <- c("DDW", "ENVIRONMENTAL MONITORING (WELLS)", "WATER SUPPLY (WELLS)")
	loc_dlat <- 0.004           # latitude wiggle (water supply well locations are imprecise and sometimes overlap)
	loc_dlong <- 0.006           # longitude wiggle
	col_names <- c(well_features, c("DATE"), analytes)
	return(list(analytes, attribs, well_features, sources, loc_dlat, loc_dlong, col_names))
	}


WellProps <- function(gama_data, well_features, sources, loc_dlat, loc_dlong){
  # generate a dataframe summarizing well characteristics
  gama_subset <- subset(gama_data, select=well_features)
  well_data <- unique(gama_subset, well_features=FALSE)
  # slightly randomize locations for water supply wells (to make them visible on a map) 
  dlat <- runif(nrow(well_data),-0.5*loc_dlat,loc_dlat)
  dlong <- runif(nrow(well_data),-0.5*loc_dlong,loc_dlong)
  ws_match <- (well_data$DATASET_CAT != "ENVIRONMENTAL MONITORING (WELLS)")
  well_data$APPROXIMATE.LATITUDE <- ws_match*dlat + well_data$APPROXIMATE.LATITUDE
  well_data$APPROXIMATE.LONGITUDE <- ws_match*dlong + well_data$APPROXIMATE.LONGITUDE  
  return (well_data)
  }
	

FillMissing <- function(data_set, full_set_names){
  # fill in missing parameters as NAs so that each county data set is compatible
  missing_names <- setdiff(full_set_names, colnames(data_set))
  data_set[,missing_names] <- NA
  return(data_set)
  }


Correlogram <- function(data_set, z_set, w_set=NULL){
  # compute and plot selected (cross-)correlograms (nothing is returned; plots are simply generated)
  cat("Computing correlogram ...\n")
  correl <- correlog(x=data_set$APPROXIMATE.LONGITUDE,
                     y=data_set$APPROXIMATE.LATITUDE,
                     z=z_set,
                     w=w_set,
                     increment=5,   # a 5-km distance class seems to work with this set
                     resamp=50,
                     latlon=TRUE,
                     quiet=TRUE)  
  plot(correl)
  }   


### input and output functions ###


GetGama <- function(file_name, analytes, attribs, well_features, sources, loc_dlat, loc_dlong, col_names){
  # read in GAMA data set for an indivudal county
  data_set <- read.csv(file=file_name, sep='\t', header=TRUE, stringsAsFactors=FALSE)
  data_set$CHEMICAL[is.na(data_set$CHEMICAL)] <- "SODIUM"               # rename "NA" as "SODIUM"
  data_set$CHEMICAL <- as.character(data_set$CHEMICAL)
  data_set <- subset(data_set, (data_set$CHEMICAL %in% analytes)==TRUE)
  data_set$DATE <- as.Date(data_set$DATE, "%m/%d/%Y")   								# convert date column to date data type
  u_match_1 <- (data_set$UNITS == "UG/L")  												# convert micrograms/L to milligrams/L, if present
  data_set$RESULT <- u_match_1*0.001*data_set$RESULT + (1-u_match_1)*data_set$RESULT
  nd_match_1 <- (data_set$QUALIFIER == "<") 											# process non-detects
  data_set$CONC <- nd_match_1*0.5*data_set$RESULT + (1-nd_match_1)*data_set$RESULT
  data_set <- subset(data_set, CONC > 0) 												# remove rows with reported conc = 0 (or CONC < 0)
  data_set$CONC <-log10(data_set$CONC)                # assume all parameter distributions are lognormal
  data_set <- data_set[attribs]
  well_data <- WellProps(data_set, well_features, sources, loc_dlat, loc_dlong)   # create separate table of unique well properties

  # for water supply wells (assumed to have unique names; need consistent locations)
  ws_data_set <- data_set[(data_set$DATASET_CAT != "ENVIRONMENTAL MONITORING (WELLS)"),]
  ws_data_set <- cast(ws_data_set, WELL.NAME + DATE ~ CHEMICAL, value = "CONC", fun.aggregate=mean)        # pivot the data
  ws_data_set <- merge(ws_data_set, well_data, by = c("WELL.NAME"))   # merge well properties and chemical concentration data
  ws_data_set <- FillMissing(ws_data_set, col_names)

  # for environmental monitoring wells (names may be non-unique, but coordinates are ...)  
  em_data_set <- data_set[(data_set$DATASET_CAT == "ENVIRONMENTAL MONITORING (WELLS)"),]  
  em_data_set <- cast(em_data_set, WELL.NAME + APPROXIMATE.LATITUDE + APPROXIMATE.LONGITUDE + COUNTY + DATASET_CAT + DATE ~ CHEMICAL, value = "CONC", fun.aggregate=mean)        # pivot the data  
  em_data_set <- FillMissing(em_data_set, col_names)

  # combine both sets and return
  full_data_set <- rbind(ws_data_set, em_data_set)
  return(full_data_set) 
  }


GetCounties <- function(){
  # read in list of county GAMA data files and return
  county_data <- read.table("counties.txt", colClasses=c("character"))
  county_files <- county_data$V1
  return(county_files)
  }  


### main script ###


GamaGeochem <- function(d_mode){

  def_output <- Definitions() 			# basic problem definitions (analayte list, data sources, etc.)
  analytes <- def_output[[1]]
  attribs <- def_output[[2]]
  well_features <- def_output[[3]]
  sources <- def_output[[4]]
  loc_dlat <- def_output[[5]]
  loc_dlong <- def_output[[6]]
  col_names <- def_output[[7]]

  if (d_mode==1){
    county_files <- GetCounties()           # read county GAMA data files
    for (i in 1:length(county_files)){ 					# read and process county GAMA data files, turn all of it into one composite pivot table
      # read in data by county
      cat("Processing ", county_files[i], "\n")
      county_data <- GetGama(county_files[i], analytes, attribs, well_features, sources, loc_dlat, loc_dlong, col_names)
      # append to composite data set
      if (i==1) {gama_data <- county_data}
      else {gama_data <- rbind(gama_data, county_data)}
      }
    write.csv(gama_data, file = "all_data.csv", row.names = FALSE)
    }
  else {gama_data <- read.csv(file="all_data.csv", sep=',', header=TRUE)}

  # filter out all NAs and note ratios
  gama_data <- na.omit(gama_data)
     
  # plot histograms for Cl-normalized concentrations
  hist(gama_data$CA-gama_data$CL, main="Calcium/Chloride", xlab="log(Ca/Cl)", col="red", breaks=40)   
  hist(gama_data$K-gama_data$CL, main="Potassium/Chloride", xlab="log(K/Cl)", col="orange", breaks=40)   
  hist(gama_data$MG-gama_data$CL, main="Magnesium/Chloride", xlab="log(Mg/Cl)", col="yellow", breaks=40)   
  hist(gama_data$NO3N-gama_data$CL, main="Nitrate/Chloride", xlab="log(NO3N/Cl)", col="green", breaks=40)      
  hist(gama_data$SO4-gama_data$CL, main="Sulfate/Chloride", xlab="log(SO4/Cl)", col="blue", breaks=40)   
  hist(gama_data$SODIUM-gama_data$CL, main="Sodium/Chloride", xlab="log(Na/Cl)", col="violet", breaks=40)
  hist(gama_data$MN-gama_data$CL, main="Manganese/Chloride", xlab="log(Mn/Cl)", col="red", breaks=40)  
  
  # conduct principal components analysis
  M1 <- cor(gama_data[analytes])
  corrplot.mixed(M1)
  pca_set <- gama_data[analytes]
  pca_object <- prcomp(pca_set)
  print(summary(pca_object))
  print(pca_object$rotation)  
  gama_data <- cbind(gama_data, pca_object$x)
  M2 <- cor(cbind(gama_data[analytes], pca_object$x))
  corrplot.mixed(M2) 

  # write all (time-dependent) results to summary output files
  gama_data$PC2 <--gama_data$PC2
  gama_data$PC3 <--gama_data$PC3
  gama_data$PC5 <--gama_data$PC5  
  write.csv(gama_data, file = "gama_data.csv", na = "", row.names = FALSE)
  
  # note standard deviationsin temporal data in each well 
  sd_data <- aggregate(gama_data[colnames(pca_object$x)], 
                        by=list(gama_data$WELL.NAME,
                                gama_data$APPROXIMATE.LATITUDE,
                                gama_data$APPROXIMATE.LONGITUDE,
                                gama_data$COUNTY,    
                                gama_data$DATASET_CAT),
                        FUN=sd)
  colnames(sd_data)[1:5] = well_features
  sd_data <- na.omit(sd_data)
  write.csv(sd_data, file = "gama_sd.csv", na = "", row.names = FALSE)  
  
  # create time-averaged PCA data frame for a subset of counties
  short_list <- c("SAN JOAQUIN", "STANISLAUS", "MERCED", "MADERA", "FRESNO", "KINGS", "TULARE",
                  "CALAVERAS", "TUOLUMNE", "MARIPOSA")  
  agg_data <- aggregate(gama_data[colnames(pca_object$x)], 
                        by=list(gama_data$WELL.NAME,
                                gama_data$APPROXIMATE.LATITUDE,
                                gama_data$APPROXIMATE.LONGITUDE,
                                gama_data$COUNTY,    
                                gama_data$DATASET_CAT),
                        FUN=sd)
  colnames(agg_data)[1:5] = well_features
  agg_data <- subset(agg_data, (agg_data$COUNTY %in% short_list)==TRUE)
  write.csv(agg_data, file = "agg_data.csv", na = "", row.names = FALSE)

  # plot selected correlograms
  for (i in 1:3){Correlogram(agg_data, agg_data[,colnames(pca_object$x)[i]])}
  Correlogram(agg_data, agg_data$PC1, agg_data$PC2)
  Correlogram(agg_data, agg_data$PC1, agg_data$PC3)
  Correlogram(agg_data, agg_data$PC2, agg_data$PC5)  

  # plot local indicators of spatial association (LISA) for PCs  
  for (i in 1:5){
    cat("LISA calcs ...", i, "\n")
    lisa_plot <- lisa(agg_data$APPROXIMATE.LONGITUDE,
                agg_data$APPROXIMATE.LATITUDE,
                agg_data[,colnames(pca_object$x)[i]],
                neigh=50, resamp=50,  latlon = TRUE, quiet = TRUE)
    plot(lisa_plot)
    }

  print ("Done.")

}


### run script ###

data_mode <- 2 						# data mode=1 --> read GAMA data from individual state files; data_mode=2 --> read from single .csv file
GamaGeochem(data_mode)



