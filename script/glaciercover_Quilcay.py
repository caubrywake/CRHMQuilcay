# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 13:25:49 2023

@author: carol

calculate the glacier % left by Rounce
"""
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np



figdir = 'C:\\CRHMquilcay\\CRHMQuilcay\\fig\\glacierevolution_Rounce\\'
# Open netcdf file from Rounce to look at annual mass balance


# List of RGIId values you want to keep
rgiid_list = [
    'RGI60-16.02197',
    'RGI60-16.02265',
    'RGI60-16.02177',
    'RGI60-16.02183',
    'RGI60-16.02183',
    'RGI60-16.02219',
    'RGI60-16.02225'
]

file_path = 'C:\\CRHMquilcay\\CRHMQuilcay\\data\\glacier_rounce\\R16_glac_area_annual_50sets_2000_2100-ssp126.nc'
dataset = xr.open_dataset(file_path)
glac_mass_annual_data = dataset['glac_area_annual']
mask = glac_mass_annual_data['RGIId'].isin(rgiid_list)
selected_data = glac_mass_annual_data.where(mask, drop=True)
data_values = selected_data.values
gm_ssp126 = np.sum(data_values, axis=1, keepdims=True)


file_path = 'C:\\CRHMquilcay\\CRHMQuilcay\\data\\glacier_rounce\\R16_glac_area_annual_50sets_2000_2100-ssp245.nc'
dataset = xr.open_dataset(file_path)
glac_mass_annual_data = dataset['glac_area_annual']
mask = glac_mass_annual_data['RGIId'].isin(rgiid_list)
selected_data = glac_mass_annual_data.where(mask, drop=True)
data_values = selected_data.values
gm_ssp245 = np.sum(data_values, axis=1, keepdims=True)

file_path = 'C:\\CRHMquilcay\\CRHMQuilcay\\data\\glacier_rounce\\R16_glac_area_annual_50sets_2000_2100-ssp370.nc'
dataset = xr.open_dataset(file_path)
glac_mass_annual_data = dataset['glac_area_annual']
mask = glac_mass_annual_data['RGIId'].isin(rgiid_list)
selected_data = glac_mass_annual_data.where(mask, drop=True)
data_values = selected_data.values
gm_ssp370 = np.sum(data_values, axis=1, keepdims=True)

file_path = 'C:\\CRHMquilcay\\CRHMQuilcay\\data\\glacier_rounce\\R16_glac_area_annual_50sets_2000_2100-ssp585.nc'
dataset = xr.open_dataset(file_path)
glac_mass_annual_data = dataset['glac_area_annual']
mask = glac_mass_annual_data['RGIId'].isin(rgiid_list)
selected_data = glac_mass_annual_data.where(mask, drop=True)
data_values = selected_data.values
gm_ssp585 = np.sum(data_values, axis=1, keepdims=True)
years = selected_data.year.values  # Assuming the 'year' dimension name is 'year'

#%% Plotting mass balance values for 2000-2100

index_2090 = np.where(years == 2100)[0]
mean_change_2090 = np.zeros((4, 2))

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns

# Plot 1: rcp26
data_values = gm_ssp126
average_values = np.mean(data_values[:, 0, :], axis=0)
std_deviation = np.std(data_values[:, 0, :], axis=0)
axs[0].fill_between(years, average_values - std_deviation, average_values + std_deviation, color='blue', alpha=0.2)
axs[0].plot(years, average_values, color='blue', linewidth=2, label='ssp126')


data_values = gm_ssp245
average_values = np.mean(data_values[:, 0, :], axis=0)
std_deviation = np.std(data_values[:, 0, :], axis=0)
axs[0].fill_between(years, average_values - std_deviation, average_values + std_deviation, color='lightgrey', alpha=0.7)
axs[0].plot(years, average_values, color='black', linewidth=2, label='ssp245')

data_values = gm_ssp370
average_values = np.mean(data_values[:, 0, :], axis=0)
std_deviation = np.std(data_values[:, 0, :], axis=0)
axs[0].fill_between(years, average_values - std_deviation, average_values + std_deviation, color='magenta', alpha=0.2)
axs[0].plot(years, average_values, color='magenta', linewidth=2, label='ssp370')
axs[0].set_xlabel('Year')
axs[0].set_ylabel('Mass Balance (kg)')
axs[0].legend()
axs[0].grid(True)  # Add grid

data_values = gm_ssp585
average_values = np.mean(data_values[:, 0, :], axis=0)
std_deviation = np.std(data_values[:, 0, :], axis=0)
axs[0].fill_between(years, average_values - std_deviation, average_values + std_deviation, color='red', alpha=0.2)
axs[0].plot(years, average_values, color='red', linewidth=2, label='ssp585')
axs[0].set_xlabel('Year')
axs[0].set_ylabel('Mass Balance (kg)')
axs[0].legend()
axs[0].grid(True)  # Add grid
axs[0].set_xlim(2015, 2100)

# Plot 2: Mean percentage change with Standard Deviation
data_values = np.squeeze(gm_ssp126, axis =1)
column_15 = data_values[:, 14]  # Note that indexing is zero-based
percentage_change = ((data_values - column_15[:, np.newaxis]) / column_15[:, np.newaxis]) * 100
#percentage_change = ((data_values - data_values[:, :1]) / data_values[:, :1]) * 100
mean_percentage_change = np.mean(percentage_change, axis=0)
std_deviation = np.std(percentage_change, axis=0)
mean_change_2090[0, 0] = mean_percentage_change[index_2090]
mean_change_2090[0, 1] = std_deviation[index_2090]
#axs[1].set_title('Mean Percentage Change (2000-2020)')
axs[1].fill_between(years, mean_percentage_change - std_deviation,
                    mean_percentage_change + std_deviation, color='blue', alpha = 0.2)
axs[1].plot(years, mean_percentage_change, color='blue', linewidth=2, label='ssp126')

data_values = np.squeeze(gm_ssp245, axis =1)
column_15 = data_values[:, 14]  # Note that indexing is zero-based
percentage_change = ((data_values - column_15[:, np.newaxis]) / column_15[:, np.newaxis]) * 100
#percentage_change = ((data_values - data_values[:, :1]) / data_values[:, :1]) * 100
mean_percentage_change = np.mean(percentage_change, axis=0)
std_deviation = np.std(percentage_change, axis=0)
mean_change_2090[1, 0] = mean_percentage_change[index_2090]
mean_change_2090[1, 1] = std_deviation[index_2090]
axs[1].fill_between(years, mean_percentage_change - std_deviation,
                    mean_percentage_change + std_deviation, color='lightgrey', alpha = 0.3)
axs[1].plot(years, mean_percentage_change, color='black', linewidth=2, label='ssp245')


data_values = np.squeeze(gm_ssp370, axis =1)
column_15 = data_values[:, 14]  # Note that indexing is zero-based
percentage_change = ((data_values - column_15[:, np.newaxis]) / column_15[:, np.newaxis]) * 100
#percentage_change = ((data_values - data_values[:, :1]) / data_values[:, :1]) * 100
mean_percentage_change = np.mean(percentage_change, axis=0)
std_deviation = np.std(percentage_change, axis=0)
mean_change_2090[2, 0] = mean_percentage_change[index_2090]
mean_change_2090[2, 1] = std_deviation[index_2090]
axs[1].fill_between(years, mean_percentage_change - std_deviation,
                    mean_percentage_change + std_deviation, color='magenta', alpha = 0.2)
axs[1].plot(years, mean_percentage_change, color='magenta', linewidth=2, label='ssp370')

data_values = np.squeeze(gm_ssp585, axis =1)
column_15 = data_values[:, 14]  # Note that indexing is zero-based
percentage_change = ((data_values - column_15[:, np.newaxis]) / column_15[:, np.newaxis]) * 100
#percentage_change = ((data_values - data_values[:, :1]) / data_values[:, :1]) * 100
mean_percentage_change = np.mean(percentage_change, axis=0)
std_deviation = np.std(percentage_change, axis=0)
mean_change_2090[3, 0] = mean_percentage_change[index_2090]
mean_change_2090[3, 1] = std_deviation[index_2090]
axs[1].fill_between(years, mean_percentage_change - std_deviation,
                    mean_percentage_change + std_deviation, color='red', alpha = 0.2)
axs[1].plot(years, mean_percentage_change, color='red', linewidth=2, label='ssp585')


axs[1].set_xlabel('Year')
axs[1].set_ylabel('Percentage Change in Area (2015-2100) (m2)')
axs[1].legend()
axs[1].grid(True)  # Add grid
axs[1].set_xlim(2015, 2100)
axs[1].set_ylim(-100, 0)

# Adjust layout to prevent overlapping
plt.tight_layout()
figname = 'QuilcayMassbalance_20002100_Rounce_perscenario.png'
plt.savefig(figdir  + figname, dpi=300)  # Save the figure
figname = 'QuilcayMassbalance_20002100_Rounce_perscenario.pdf'
plt.savefig(figdir  + figname, dpi=300)  # Save the figure

# Show the plot
plt.show()


import pandas as pd

# Assuming your array is named 'my_array'

# Create row headers and column headers
row_headers = ['SSP126','SSP245','SSP370', 'SSP585']
column_headers = ['MeanChange', 'StdDev']

# Create a DataFrame with the array and headers
file_path = 'C:/CRHMquilcay/CRHMQuilcay/data/glacier_rounce/ChangeinGlacierVolume_perscenario_2095.csv'
df = pd.DataFrame(mean_change_2090, index=row_headers, columns=column_headers)
df.to_csv(file_path)


#%% For RCP instead
figdir = 'C:\\CRHMquilcay\\CRHMQuilcay\\fig\\glacierevolution_Rounce\\'
# Open netcdf file from Rounce to look at annual mass balance


# List of RGIId values you want to keep
rgiid_list = [
    'RGI60-16.02197',
    'RGI60-16.02265',
    'RGI60-16.02177',
    'RGI60-16.02183',
    'RGI60-16.02183',
    'RGI60-16.02219',
    'RGI60-16.02225'
]

file_path = 'C:\\CRHMquilcay\\CRHMQuilcay\\data\\glacier_rounce\\R16_glac_area_annual_50sets_2000_2100-rcp26.nc'
dataset = xr.open_dataset(file_path)
glac_mass_annual_data = dataset['glac_area_annual']
mask = glac_mass_annual_data['RGIId'].isin(rgiid_list)
selected_data = glac_mass_annual_data.where(mask, drop=True)
data_values = selected_data.values
gm_rcp26= np.sum(data_values, axis=1, keepdims=True)


file_path = 'C:\\CRHMquilcay\\CRHMQuilcay\\data\\glacier_rounce\\R16_glac_area_annual_50sets_2000_2100-rcp45.nc'
dataset = xr.open_dataset(file_path)
glac_mass_annual_data = dataset['glac_area_annual']
mask = glac_mass_annual_data['RGIId'].isin(rgiid_list)
selected_data = glac_mass_annual_data.where(mask, drop=True)
data_values = selected_data.values
gm_rcp45= np.sum(data_values, axis=1, keepdims=True)

file_path = 'C:\\CRHMquilcay\\CRHMQuilcay\\data\\glacier_rounce\\R16_glac_area_annual_50sets_2000_2100-rcp85.nc'
dataset = xr.open_dataset(file_path)
glac_mass_annual_data = dataset['glac_area_annual']
mask = glac_mass_annual_data['RGIId'].isin(rgiid_list)
selected_data = glac_mass_annual_data.where(mask, drop=True)
data_values = selected_data.values
gm_rcp85= np.sum(data_values, axis=1, keepdims=True)

years = selected_data.year.values  # Assuming the 'year' dimension name is 'year'

#%% Plotting mass balance values for 2000-2100

index_2090 = np.where(years == 2100)[0]
mean_change_2090 = np.zeros((3, 2))

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns

# Plot 1: rcp26

data_values = gm_rcp26
average_values = np.mean(data_values[:, 0, :], axis=0)
std_deviation = np.std(data_values[:, 0, :], axis=0)
axs[0].fill_between(years, average_values - std_deviation, average_values + std_deviation, color='lightgrey', alpha=0.7)
axs[0].plot(years, average_values, color='black', linewidth=2, label='rcp26')

data_values = gm_rcp45
average_values = np.mean(data_values[:, 0, :], axis=0)
std_deviation = np.std(data_values[:, 0, :], axis=0)
axs[0].fill_between(years, average_values - std_deviation, average_values + std_deviation, color='magenta', alpha=0.2)
axs[0].plot(years, average_values, color='magenta', linewidth=2, label='rcp45')
axs[0].set_xlabel('Year')
axs[0].set_ylabel('Mass Balance (kg)')
axs[0].legend()
axs[0].grid(True)  # Add grid

data_values = gm_rcp85
average_values = np.mean(data_values[:, 0, :], axis=0)
std_deviation = np.std(data_values[:, 0, :], axis=0)
axs[0].fill_between(years, average_values - std_deviation, average_values + std_deviation, color='red', alpha=0.2)
axs[0].plot(years, average_values, color='red', linewidth=2, label='rcp85')
axs[0].set_xlabel('Year')
axs[0].set_ylabel('Area (km2)')
axs[0].legend()
axs[0].grid(True)  # Add grid
axs[0].set_xlim(2015, 2100)

# Plot 2: Mean percentage change with Standard Deviation

data_values = np.squeeze(gm_rcp26, axis =1)
column_15 = data_values[:, 14]  # Note that indexing is zero-based
percentage_change = ((data_values - column_15[:, np.newaxis]) / column_15[:, np.newaxis]) * 100
#percentage_change = ((data_values - data_values[:, :1]) / data_values[:, :1]) * 100
mean_percentage_change = np.mean(percentage_change, axis=0)
std_deviation = np.std(percentage_change, axis=0)
mean_change_2090[0, 0] = mean_percentage_change[index_2090]
mean_change_2090[0, 1] = std_deviation[index_2090]
axs[1].fill_between(years, mean_percentage_change - std_deviation,
                    mean_percentage_change + std_deviation, color='blue', alpha = 0.3)
axs[1].plot(years, mean_percentage_change, color='blue', linewidth=2, label='rcp26')


data_values = np.squeeze(gm_rcp45, axis =1)
percentage_change = ((data_values - column_15[:, np.newaxis]) / column_15[:, np.newaxis]) * 100
#percentage_change = ((data_values - data_values[:, :1]) / data_values[:, :1]) * 100
mean_percentage_change = np.mean(percentage_change, axis=0)
std_deviation = np.std(percentage_change, axis=0)
mean_change_2090[1, 0] = mean_percentage_change[index_2090]
mean_change_2090[1, 1] = std_deviation[index_2090]
axs[1].fill_between(years, mean_percentage_change - std_deviation,
                    mean_percentage_change + std_deviation, color='magenta', alpha = 0.2)
axs[1].plot(years, mean_percentage_change, color='magenta', linewidth=2, label='rcp45')

data_values = np.squeeze(gm_rcp85, axis =1)
percentage_change = ((data_values - column_15[:, np.newaxis]) / column_15[:, np.newaxis]) * 100
#percentage_change = ((data_values - data_values[:, :1]) / data_values[:, :1]) * 100
mean_percentage_change = np.mean(percentage_change, axis=0)
std_deviation = np.std(percentage_change, axis=0)
mean_change_2090[2, 0] = mean_percentage_change[index_2090]
mean_change_2090[2, 1] = std_deviation[index_2090]
axs[1].fill_between(years, mean_percentage_change - std_deviation,
                    mean_percentage_change + std_deviation, color='red', alpha = 0.2)
axs[1].plot(years, mean_percentage_change, color='red', linewidth=2, label='rcp85')


axs[1].set_xlabel('Year')
axs[1].set_ylabel('Percentage Change in area (2015-2100)')
axs[1].legend()
axs[1].grid(True)  # Add grid
axs[1].set_xlim(2015, 2100)
axs[1].set_ylim(-100, 0)

# Adjust layout to prevent overlapping
plt.tight_layout()
figname = 'QuilcayMassbalance_20002100_Rounce_perscenarioRCP.png'
plt.savefig(figdir  + figname, dpi=300)  # Save the figure
figname = 'QuilcayMassbalance_20002100_Rounce_perscenarioRCP.pdf'
plt.savefig(figdir  + figname, dpi=300)  # Save the figure

# Show the plot
plt.show()

# Create row headers and column headers
row_headers = ['RCP26','RCP45', 'RCP85']
column_headers = ['MeanChange', 'StdDev']

# Create a DataFrame with the array and headers
file_path = 'C:/CRHMquilcay/CRHMQuilcay/data/glacier_rounce/ChangeinGlacierVolume_perscenario_2095RCP.csv'
df = pd.DataFrame(mean_change_2090, index=row_headers, columns=column_headers)
df.to_csv(file_path)

#%% Calculate what would be the basin glacier cover in comparison
# if glacier mass balance goes from 100 to 40%, and initial glacoer cover was 26%, then if we lose 60% of 26%?

init_glaccover = 26;
df['RemainingGlacierCover'] = init_glaccover  * (100 + df['MeanChange']) / 100
# Calculate the upper and lower bounds considering the uncertainty (1 std dev)
df['RemainingGlacierCover_Upper'] = init_glaccover  * (100 +(df['MeanChange'] - df['StdDev']) )/ 100
df['RemainingGlacierCover_Lower'] = init_glaccover  * (100 +(df['MeanChange'] + df['StdDev']) )/ 100

# or
# Calculate the effect of the standard deviation on the remaining glacier cover
df['UncertaintyRem'] = init_glaccover * (df['StdDev'] / 100)

# Summarize the result as RemainingGlacierCover ± Uncertainty
df['RemainingGlacierCover_with_Uncertainty'] = df['RemainingGlacierCover'].astype(str) + ' ± ' + df['UncertaintyRem'].astype(str)
