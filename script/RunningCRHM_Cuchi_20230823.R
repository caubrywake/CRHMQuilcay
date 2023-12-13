## This script runs the CRHM simulations for the Quilcayhuanca simulations
# created by caroline Aubr-wake, edited 2023-12-12

# obs filF: F:/11_CRHM_cuchi/data/processed/CuchiObs_2014_2020.obs')


# Run main simulations

library(CRHMr)
prjname <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'

# make a copy of the file
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
file.copy(originalFile, copyFile)

## Runs - Basinflow_icemelt
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='F:/11_CRHM_cuchi/CRHM/output/v8/RunFlow.txt')

# run2
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Forcings.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='F:/11_CRHM_cuchi/CRHM/output/v8/RunForcings.txt')

# run3
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_ET.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='F:/11_CRHM_cuchi/CRHM/output/v8/RunET.txt')

# run4
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Routing.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='F:/11_CRHM_cuchi/CRHM/output/v8/RunRouting.txt')

############# Run noice
prjname <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823_noice.prj'

## Runs - Basinflow_icemelt
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/Quilcay_withsprings_20230503_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='F:/11_CRHM_cuchi/CRHM/output/v8/NoIce_RunFlow.txt')

############# Run no leakage
prjname <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823_noleakage.prj'

## Runs - Basinflow_icemelt
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/Quilcay_withsprings_20230503_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='F:/11_CRHM_cuchi/CRHM/output/v8/NoLeakage_RunFlow.txt')

############# Run no gw
prjname <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823_nogw.prj'

## Runs - Basinflow_icemelt
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/Quilcay_withsprings_20230503_Basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',prjname, outFile='F:/11_CRHM_cuchi/CRHM/output/v8/NoLeakage_RunFlow.txt')

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
##
# make a copy of the file
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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

outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenarioquilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 1 ########################################################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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

outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 2 ##############################################################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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

outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 3 ########################################################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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

outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 4 #####################################################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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

outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 5 ########################################################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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

outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_noglacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_noice_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}


###############################################################################
##############################################################################

## Glacier at 2 and 5 ############################################

# Temp 0 ###########################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 1 ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp2  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 3  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 4  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 5  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_25glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_25glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}



###############################################################################
##############################################################################

## Glacier at 2 ############################################

# Temp 0  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 1  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 2  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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

outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 3  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 4  ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}

# Temp 5 ############################################
originalFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_20230823.prj'
copyFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
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
outputPrjFile <- 'F:/11_CRHM_cuchi/CRHM/prjfile/quilcay_scenario_20230823.prj'
setPrjParameters(outputPrjFile, 'glacier ice_init', hru_2glacier, quiet = FALSE)
# set output variables
filename <- 'F:/11_CRHM_cuchi/CRHM/prjfile/setvariables_Scenarios.prj'
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
  outputFileName <- sprintf("F:/11_CRHM_cuchi/CRHM/output/scenario/quilcay_2glac_t_%.0f_precip_%.0f.txt", temp_range, precip_range * 100)
  # run
  result<-runCRHM(CRHMfile = 'F:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe',outputPrjFile, outFile=outputFileName)
}
