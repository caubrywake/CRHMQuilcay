## This script runs the CRHM simulations for the Quilcayhuanca simulations
# created by caroline Aubr-wake, edited 2023-12-12

# obs filG: G:/11_CRHM_cuchi/data/processed/CuchiObs_2014_2020.obs')


# Run main simulations

library(CRHMr)
prjname <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'

# make a copy of the file
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
file.copy(originalFile, copyFile)

## Runs - Basinflow_icemelt
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/RunFlow.txt')

# run2
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Forcings.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/RunForcings.txt')

# run3
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_ET.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/RunET.txt')

# run4
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Routing.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/RunRouting.txt')

# run5
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_glacvol.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/RunGlacVol.txt')

############# Run noice
hru_glacier <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# make a copy of the file
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_noice_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- copyFile
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_glacier, quiet = FALSE)

prjname <- outputPrjFile
## Runs - Basinflow_icemelt
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/NoIce_RunFlow.txt')

#################################################
############# Run no leakage
hru_noleak <- c(-8,-8,-8, -8, -8, -8, -8,0,-8,-8, -8, -8, -8, -8,-8,-8,-8,-8,-8)

# make a copy of the file
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_noleak_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- copyFile
setPrjParameters(outputPrjFile, 'Netroute_D gwwhereto', hru_noleak, quiet = FALSE)

prjname <- outputPrjFile 
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/NoLeakage_RunFlow.txt')

############# Run no gw
hru_nogwmax <- c(9, 9, 9, 9, 9, 9 ,9 ,9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 ,9 )
hru_nogwinit <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2) 
hru_nogwstor <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1)
hru_nogwdelay <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# make a copy of the file
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_nogw_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- copyFile
setPrjParameters(outputPrjFile, 'Netroute_D gwLag', hru_nogwdelay, quiet = FALSE)
setPrjParameters(outputPrjFile, 'Netroute_D gwKstorage', hru_nogwstor, quiet = FALSE)
setPrjParameters(outputPrjFile, 'Soil gw_max', hru_nogwmax, quiet = FALSE)
setPrjParameters(outputPrjFile, 'Soil gw_init', hru_nogwinit, quiet = FALSE)

# set output variables
prjname <- outputPrjFile 
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='G:/11_CRHM_cuchi/CRHM/output/v9/NoGW_RunFlow.txt')

####### Scenarios run
###############################################################################################################################################
library(CRHMr)

# SET SCENARIOS
# 
HRUno = 19
# glacier configuration
hru_noglacier <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
hru_25glacier<- c(0, 5E6, 0, 0, 5E6,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
hru_2glacier<- c(0, 5E6, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
hru_allglacier<- c(0, 5E6,5E6, 0, 5e6,5E6, 0, 0, 5E6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
###############################################################################
##############################################################################

## NO ICE

# make a copy of the file
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

# Temp 0
temp_range <- c(0)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 1 ########################################################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(1)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 2 ##############################################################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(2)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 3 ########################################################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(3)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 4 #####################################################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(4)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 5 ########################################################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(5)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}


###############################################################################
##############################################################################

## Glacier at 2 and 5 ############################################

# Temp 0 ###########################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(0)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 1 ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(1)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp2  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(2)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 3  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(3)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 4  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(4)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 5  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(5)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}



###############################################################################
##############################################################################

## Glacier at 2 ############################################

# Temp 0  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(0)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 1  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(1)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 2  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(2)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 3  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(3)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 4  ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(4)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 5 ############################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}
outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(5)
precip_ranges <- c(0.8, 0.9, 1, 1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}


##########################################################################
### simulation with snow/ice partitining for "likely scenario" for supplementary figure

####### Scenarios run
###############################################################################################################################################
library(CRHMr)

# SET SCENARIOS
# 
HRUno = 19
# glacier configuration
hru_noglacier <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
###############################################################################
##############################################################################

## NO ICE

# make a copy of the file
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_ScenariosSnowRain.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

# Temp + 4 C
temp_range <- c(4)
precip_ranges <- c(1.1, 1.2)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/SnowRain/quilcay_noice_RainSnow_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 1 ########################################################################
originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
copyFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
# Check if the copy file exists
if (file.exists(copyFile)) {
  # Option 1: Overwrite the existing copy file
  file.copy(originalFile, copyFile, overwrite = TRUE)
  
  # Option 2: Delete the existing copy file and create a new one
  file.remove(copyFile)
  file.copy(originalFile, copyFile)
} else {
  # If the copy file doesn't exist, simply create a new copy
  file.copy(originalFile, copyFile)
}

outputPrjFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'G:/11_CRHM_cuchi/CRHM/prjfile/setvariables_ScenariosSnowRain.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(outputPrjFile, variables)# 

temp_range <- c(5)
for (precip_range in precip_ranges) {
  # Set the values of temp_range and precip_range in the corresponding lines
  setPrjParameters(outputPrjFile, 'obs ClimChng_t', rep(temp_range, HRUno), quiet = FALSE)
  setPrjParameters(outputPrjFile, 'obs ClimChng_precip', rep(precip_range, HRUno), quiet = FALSE)
  # automate project
  result <- automatePrj(outputPrjFile)
  # define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/output/scenario/v9/SnowRain/quilcay_noice_RainSnow__t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}
