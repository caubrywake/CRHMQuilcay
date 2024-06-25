###############################################################################################################################################
library(CRHMr)

# SET SCENARIOS

#######################################################
####### for most variable #################
## Variable names
modify_parameters <- function(params, factor) {
  return(params * factor)
}

# Factors to multiply the parameters by
factors <- c(0.1, 0.5, 1, 1.5, 2, 10)

# Full parameter names
varnames <- c( "Netroute_D ssrKstorage", "Netroute_D gwLag", "Netroute_D gwKstorage", 
               "Shared Sdmax", "Soil gw_max",
               "evap_Resist soil_Depth","evap_Resist rcs", "evap_Resist F_Qg",
               "glacier ice_init",  'K_Estimate Ks_gw', "K_Estimate Ks_lower", "K_Estimate Ks_upper")

# Loop through each set of parameters and run simulations
for (i in 1:length(factors)) {
  for (j in 1:length(varnames)) {
    param_full_name <- varnames[j]
    
    # Extract the base parameter name for file naming
    base_param_name <- gsub(".* ", "", param_full_name)
    
    # Copy the original project file for each iteration
    originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
    copyFile <- paste0('G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528_', base_param_name, '_', i, '.prj')
    if (file.exists(copyFile)) {
      file.remove(copyFile)
    }
    file.copy(originalFile, copyFile)
    
    # Output project file
    outputPrjFile <- copyFile
    
    # Set initial parameters and output variables
    filename <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/setvariables_routingsensitivity.prj'
    variables <- readPrjOutputVariables(filename, asDataframe = FALSE)
    setPrjOutputVariables(outputPrjFile, variables)
    
    # Get the current parameter values
    current_params <- readPrjParameters(outputPrjFile, param_full_name)
    
    # Modify the parameter values based on iteration
    new_params <- modify_parameters(current_params, factors[i])
    
    # Set the parameter values
    setPrjParameters(outputPrjFile, param_full_name, new_params, quiet = FALSE)
    
    # Automate project and run the simulation
    result <- automatePrj(outputPrjFile)
    
    # Define output file name
    outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/20240527/output/scenarioquilcay_routing_%s_%d.txt", base_param_name, i)
    
    # Run the model
    result <- runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe', outputPrjFile, outFile = outputFileName)
  }
}

########## export parameter
##########################################
# Define a function to modify parameters
modify_parameters <- function(params, factor) {
  return(params * factor)
}

# Factors to multiply the parameters by
factors <- c(0.1, 0.5, 1, 1.5, 2, 10)

# Full parameter names
varnames <- c("Netroute_D ssrKstorage", "Netroute_D gwLag", "Netroute_D gwKstorage", 
              "Shared Sdmax", "Soil gw_max",
              "evap_Resist soil_Depth", "evap_Resist rcs", "evap_Resist F_Qg",
              "glacier ice_init", 'K_Estimate Ks_gw', "K_Estimate Ks_lower", "K_Estimate Ks_upper")

# Array to store resulting parameter values
resulting_params <- list()

# Loop through each set of factors and parameters
for (i in 1:length(factors)) {
  factor <- factors[i]
  
  for (j in 1:length(varnames)) {
    param_full_name <- varnames[j]
    
    # Extract the base parameter name for file naming
    base_param_name <- gsub(".* ", "", param_full_name)
    
    # Copy the original project file for each iteration
    originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
    copyFile <- paste0('G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528_', base_param_name, '_', i, '.prj')
    if (file.exists(copyFile)) {
      file.remove(copyFile)
    }
    file.copy(originalFile, copyFile)
    
    # Output project file
    outputPrjFile <- copyFile
    
    # Set initial parameters and output variables
    filename <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/setvariables_routingsensitivity.prj'
    variables <- readPrjOutputVariables(filename, asDataframe = FALSE)
    setPrjOutputVariables(outputPrjFile, variables)
    
    # Get the current parameter values
    current_params <- readPrjParameters(outputPrjFile, param_full_name)
    
    # Modify the parameter values based on iteration
    new_params <- modify_parameters(current_params, factor)
    
    # Store the parameter values in the list
    resulting_params[[length(resulting_params) + 1]] <- c(param_full_name, factor, new_params)
  }
}

# Convert the list to a data frame
resulting_params_df <- do.call(rbind, resulting_params)
colnames(resulting_params_df) <- c("Variable", "Factor", paste0("Value", 1:(ncol(resulting_params_df) - 2)))

# Export the table to a CSV file
write.csv(resulting_params_df, "G:/11_CRHM_cuchi/CRHM/20240527/output/resulting_parameters_routing.csv", row.names = FALSE)

####################################################

### GW leakage

#######################################################
####### for most variable #################
## Variable names# Predefined parameter sets
param_sets <- list(
   noleaknosprings =    c(-8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, -8,   -8,  -8, -8, -8, -8,  -8, -8),
   noleak_withsprings = c(-8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, 100,  -8, 100, -8, -8, 100, -8, -8),
   accglacleak =        c(-8, -1, -8, -8, -1, -8, -8, 0, -8, -8, -8, 100,  -8, 100, -8, -8, 100, -8, -8),
   ablglacleak =        c(-8, -8, -1, -8, -8, -1, -8, 0, -9, -8, -8, 100,  -8, 100, -8, -8, 100, -8, -8),
   allglacleak =        c(-8, -1, -1, -8, -1, -1, -8, 0, -9, -8, -8, 100,  -8, 100, -8, -8, 100, -8, -8), 
   onlylakeleak =       c(-8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, 100,  -8, 100, -8, -8, 100, -1, -1),
  onlyupvalleak =      c(-8, -8, -8, -8, -8, -8, -1, 0, -8, -8, -8, 100,  -8, 100, -8, -8, 100, -1, -1),
  rocklakevalleak =    c(-1, -8, -8, -1, -8, -8, -1, 0, -8, -1, -8, 100,  -8, 100, -1, -8, 100, -1, -1),
  rocklakevalglacleak =c(-1, -1, -1, -1, -1, -1, -1, 0, -8, -1, -8, 100,  -8, 100, -1, -8, 100, -1, -1),
  allleak =            c(-1, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, 100,  -1, 100, -1, -1, 100, -1, -1),
   base =               c(-1, -8, -8, -8, -8, -8, -8, 0, -1, -8, -8, 100,  -8, 100, -1, -8, 100, -1, -1))

# Full parameter names
varnames <- c("Netroute_D gwwhereto")

# Loop through each set of parameters and run simulations
for (param_set_name in names(param_sets)) {
  for (j in 1:length(varnames)) {
    param_full_name <- varnames[j]
    
    # Extract the base parameter name for file naming
    base_param_name <- gsub(".* ", "", param_full_name)
    
    # Copy the original project file for each iteration
    originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
    copyFile <- paste0('G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528_', base_param_name, '_', param_set_name, '.prj')
    if (file.exists(copyFile)) {
      file.remove(copyFile)
    }
    file.copy(originalFile, copyFile)
    
    # Output project file
    outputPrjFile <- copyFile
    
    # Set initial parameters and output variables
    filename <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/setvariables_routingsensitivity.prj'
    variables <- readPrjOutputVariables(filename, asDataframe = FALSE)
    setPrjOutputVariables(outputPrjFile, variables)
    
    # Get the current parameter values (this step might be redundant if you are not using current_params)
    current_params <- readPrjParameters(outputPrjFile, param_full_name)
    
    # Use the predefined parameter set
    new_params <- param_sets[[param_set_name]]
    
    # Set the parameter values
    setPrjParameters(outputPrjFile, param_full_name, new_params, quiet = FALSE)
    
    # Automate project and run the simulation
    result <- automatePrj(outputPrjFile)
    
    # Define output file name
    outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/20240527/output/scenarioquilcay_leakage_%s_%s.txt", base_param_name, param_set_name)
    
    # Run the model
    result <- runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe', outputPrjFile, outFile = outputFileName)
  }
}

###########################################################################
################ for soil moisture max - need to vary together
# Define a function to modify parameters
modify_parameters <- function(params, factor) {
  return(params * factor)
}

# Factors to multiply the parameters by
factors <- c(0.1, 0.5, 1, 1.5, 2, 10)

# Full parameter names that need to be modified together
varnames <- c("Shared soil_moist_max", "Shared soil_rechr_max")

# Array to store resulting parameter values
resulting_params <- list()

# Loop through each set of factors and run simulations
for (i in 1:length(factors)) {
  factor <- factors[i]
  
  # Copy the original project file for each iteration
  originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
  copyFile <- paste0('G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528_soilparams_', i, '.prj')
  if (file.exists(copyFile)) {
    file.remove(copyFile)
  }
  file.copy(originalFile, copyFile)
  
  # Output project file
  outputPrjFile <- copyFile
  
  # Set initial parameters and output variables
  filename <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/setvariables_routingsensitivity.prj'
  variables <- readPrjOutputVariables(filename, asDataframe = FALSE)
  setPrjOutputVariables(outputPrjFile, variables)
  
  # Modify both parameters together
  for (j in 1:length(varnames)) {
    param_full_name <- varnames[j]
    
    # Get the current parameter values
    current_params <- readPrjParameters(outputPrjFile, param_full_name)
    
    # Modify the parameter values based on iteration
    new_params <- modify_parameters(current_params, factor)
    
    # Set the parameter values
    setPrjParameters(outputPrjFile, param_full_name, new_params, quiet = FALSE)
  }
  
  # Automate project and run the simulation
  result <- automatePrj(outputPrjFile)
  
  # Define output file name
  outputFileName <- sprintf("G:/11_CRHM_cuchi/CRHM/20240527/output/scenarioquilcay_soilparams_%d.txt", i)
  
  # Run the model
  result <- runCRHM(CRHMfile = 'G:/11_CRHM_cuchi/CRHM/CRHM_042222/CRHM.exe', outputPrjFile, outFile = outputFileName)
}

##########################################################
## Extract param values ##########################################################
##########################################################
################ for soil moisture max - need to vary together
# Define a function to modify parameters
modify_parameters <- function(params, factor) {
  return(params * factor)
}

# Factors to multiply the parameters by
factors <- c(0.1, 0.5, 1, 1.5, 2, 10)

# Full parameter names that need to be modified together
varnames <- c("Shared soil_moist_max", "Shared soil_rechr_max")

# Initialize a list to store resulting parameter values
resulting_params <- list()

# Loop through each set of factors and run simulations
for (i in 1:length(factors)) {
  factor <- factors[i]
  
  # Copy the original project file for each iteration
  originalFile <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_20240528.prj'
  copyFile <- paste0('G:/11_CRHM_cuchi/CRHM/20240527/prjfile/quilcay_scenario_20240528_soilparams_', i, '.prj')
  if (file.exists(copyFile)) {
    file.remove(copyFile)
  }
  file.copy(originalFile, copyFile)
  
  # Output project file
  outputPrjFile <- copyFile
  
  # Set initial parameters and output variables
  filename <- 'G:/11_CRHM_cuchi/CRHM/20240527/prjfile/setvariables_routingsensitivity.prj'
  variables <- readPrjOutputVariables(filename, asDataframe = FALSE)
  setPrjOutputVariables(outputPrjFile, variables)
  
  # Modify both parameters together
  for (j in 1:length(varnames)) {
    param_full_name <- varnames[j]
    
    # Get the current parameter values
    current_params <- readPrjParameters(outputPrjFile, param_full_name)
    
    # Modify the parameter values based on iteration
    new_params <- modify_parameters(current_params, factor)
    
    # Store the parameter values in the list
    resulting_params[[length(resulting_params) + 1]] <- c(param_full_name, factor, new_params)
  }
}

# Convert the list to a data frame
resulting_params_df <- do.call(rbind, resulting_params)
colnames(resulting_params_df) <- c("Variable", "Factor", paste0("HRU", 1:(ncol(resulting_params_df) - 2)))

# Export the table to a CSV file
write.csv(resulting_params_df, "G:/11_CRHM_cuchi/CRHM/20240527/output/resulting_parameters_soilmoist.csv", row.names = FALSE)


