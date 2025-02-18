#### SETUP ####

## LOAD NEEDED PACKAGES

install.packages("ncdf4")
install.packages("openxlsx") 
install.packages("lubridate")
install.packages("timetk")
install.packages("corrr")
install.packages("ggcorrplot")
install.packages("FactoMineR")

library(FactoMineR) # PCA PACKAGE 
library(ggcorrplot) # PCA PACKAGE FOR VISUALIZING
library(corrr) # PCA PACKAGE FOR DATAFRAMES
library(openxlsx) # CREATE XLSX FILES
library(ncdf4) # NETCDF HANDLING
library(dplyr) # DATAFRAME ORGANIZING


## READ IN DATA
bham_current=nc_open("../Data/bellingham_bay_2023.06.01_2024.02.26.nc")
bham_current2=nc_open("../Data/bellingham_bay_2024.02.26_2024.06.08.nc")

#### CURRENT DATA EXTRACTION ####

### RUN THE FOLLOWING EXTRACTION CODE FOR BOTH bham_current AND bham_current2 THEN COMBINE EXCELS BEFORE RUNNING PCA


## OPEN THE NET CDF FILE
ncin <- nc_open(bham_current)
print(bham_current)

## LIST THE NAMES OF THE VARIABLES
names(bham_current$var)

## SET THE VARIABLES WE WANT TO AN OBJECT
ubar <- ncvar_get(bham_current, "ubar")
vbar <- ncvar_get(bham_current, "vbar")
ocean_time <- ncvar_get(bham_current, "ocean_time")

## CREATE A DATA FRAME
current<- data.frame(ubar = ubar, vbar = vbar, ocean_time = ocean_time)

## FIND UNIT OF ocean_time, EPOCH TIME
dlname = ncatt_get(bham_current, "ocean_time", "units")

## TURN EPOCH TIME INTO PST AND PDT TIME
timedate_PDT_PST <- format(timedate_UTC, tz = "America/Los_Angeles", usetz = TRUE)

## ADD timedate  TO CURRENT DATAFRAME
current2.0 <- cbind(current, timedate_PDT_PST, timedate_UTC)

## FIND THE MAGNITUDE OF THE CURRENT (bar) VELOCITY
Vmag <- sqrt((current2.0$ubar)^2 + (current2.0$vbar)^2)

## ADD VELOCITY MAGNITUDE TO current2.0 DATAFRAME
current3.0 <- cbind(current2.0, Vmag)

## TURN current3.0 DATAFRAME INTO A XLSX
write.xlsx(current3.0, file = "second_current_data.xlsx", sheetName = "Sheet1", rowNames = TRUE)

## GET LONGITUDE AND LATITUDE
lonu <- ncvar_get(bham_current,"lon_u")
longv<- ncvar_get(bham_current, "lon_v")
nlon <- dim(lon)


#### PCA ####

## READ IN DATA FOR PCA AND ORGANIZE IT
PCA_data = read.csv("../Data/full_current_data.csv", header=T)
str(PCA_data)
PCA_data= PCA_data[,2:3]

## PERFORM PCA
pca_result <- prcomp(PCA_data, scale. = F)

## GET THE LOADINGS (COEFFICIENTS) OF EACH VARIABLE ON THE PRINCIPAL COMPONENTS, NOT NEEDED FOR FINAL RESULT 
loadings <- pca_result$rotation

## GET THE SCORES (COORDINATES) OF EACH OBSERVATION ON THE PRINCIPAL COMPONENTS, THESE ARE THE VALUES WE NEED. USING THE PC1 VALUES
scores <- pca_result$x

PCA_values <- data.frame(scores)

## EXPORT FINAL PCA RESULTS
write.xlsx(PCA_values, file = "PCA_result_full.xlsx", sheetName = "Sheet1", rowNames = TRUE)

## PRINT THE RESULTS
print(loadings)
print(scores)



## TEST IF THE PCA WORKED ACCURATELY (IF vmag BEFORE AND AFTER IS ABOUT THE SAME), APPEARS TO BE SIMILAR ENOUGH TO MOVE FORWARD WITH IT. 
vmag_PCA = sqrt(scores[,1]^2+scores[,2]^2)
vmag=sqrt(PCA_data$ubar^2+PCA_data$vbar^2)

plot(vmag, vmag_PCA)

## PLOT HODOGRAPH TO SEE HOW MUCH OF A CHANGE THE PCA CREATED IN THE DATASET. NOTHING HUGE
plot(PCA_data$ubar, PCA_data$vbar)
points(scores[,1], scores[,2], col="red")

#### EXTRACT NECESSARY DATA POINTS ####

### WE NEED TO EXTRACT JUST THE SPECIC PC1 VALUES NEEDED FOR OUR OBSERVATIONS

## READ IN DATASET

full_current=read.csv("../Data/full_current_data.csv", header=T)
full_current= full_current[,1:5]

datetime=read.csv("../Data/datetime.csv", header=T)
datetime$date=as.Date(datetime$date)
datetime$time=as.numeric(datetime$time)

## ORGANIZE DATASETS

# CONVERT THE TIME DATE COLUMNS TO POSIXct FORMAT
full_current$timedate <- as.POSIXct(full_current$timedate_PDT_PST, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles")

# COMBINE DATE AND TIME INTO A SINGLE DATETIME STRING
datetime$timedate <- as.POSIXct(paste(datetime$date, sprintf("%02d:00:00", datetime$time)),format="%Y-%m-%d %H:%M:%S", tz="America/Los_Angeles")

## EXTRACTION OF PC1 VALUES FOR OBSERVATIONS

# JOIN THE TWO DATASETS TO FIND MATCHING PC1 VALUES
result <- datetime %>%
  left_join(full_current, by = "timedate") %>%
  select(timedate, PC1)

# PRINT THE RESULT, REMOVE THE RANDOM DUPLICATE OF 2/26/2024
print(result)

# CREATE EXCEL WITH THE PC1 VALUES, EXPORT TO WORKING DIRECTORY
write.xlsx(result, file = "extracted_PC1values.xlsx", sheetName = "Sheet1", rowNames = TRUE)
