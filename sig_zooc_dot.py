#! /usr/bin/python
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
#%matplotlib inline
import matplotlib as mpl
from matplotlib import colorbar

import cmocean
import cartopy
import cartopy.util as util
import cartopy.crs as ccrs
import cartopy.feature as cfeature  # features such as the ocean, coastlines rivers, etc

mpl.rcParams['figure.dpi'] = 150
mpl.rcParams['savefig.dpi'] = 200
mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'


path = '/work/bb0820/ISIMIP/ISIMIP3a/SecondaryInputData/climate/ocean/ctrlclim/global/monthly/historical/GFDL-MOM6-COBALT2/gfdl-mom6-cobalt2_ctrlclim_zooc_onedeg_global_monthly_1961_2010.nc'
temp_cntrl = xr.open_mfdataset(path, decode_times=False)

temp_cntrl_data = temp_cntrl['zooc'][:30*12,0,:,:]# get first 30 yeats

stdev_cntrl = np.std(temp_cntrl_data,axis=0) #std of first 30 years
temp_cntrl_mean_First30y = np.mean(temp_cntrl_data, axis = 0)
temp_cntrl_mean_Last10y = np.mean(( temp_cntrl['zooc'][10*12:,0,:,:]) , axis = 0)

stdev_cntrl = np.array(stdev_cntrl)

path = '/work/bb0820/ISIMIP/ISIMIP3a/InputData/climate/ocean/obsclim/global/monthly/historical/GFDL-MOM6-COBALT2/gfdl-mom6-cobalt2_obsclim_zooc_onedeg_global_monthly_1961_2010.nc'
temp_obs = xr.open_mfdataset(path, decode_times=False)

temp_obs_data = temp_obs['zooc'][10*12:,0,:,:]# get last 10 yeats
temp_obs_mean_Last10y = np.mean(temp_obs_data, axis = 0)


last_10_anomaly = temp_obs_mean_Last10y - temp_cntrl_mean_Last10y

L10_anom = temp_obs_mean_Last10y - temp_cntrl_mean_Last10y
L10_anom = np.array(L10_anom)
plot_only_stdev=L10_anom
plot_only_stdev[np.abs(plot_only_stdev) < (stdev_cntrl *1)] = np.nan #the version made without dots #not used here
plot_this = temp_obs_mean_Last10y 
dots_lons=[]
dots_lats=[]
plot_only_stdev = np.array(plot_only_stdev)
s_sig=np.zeros([len(temp_obs_data['lat']),len(temp_obs_data['lon'])]) #inorder to calculate stdev the same way we calculate the dots


dots_lons2=[]
dots_lats2=[]

for lon in range(0,len(temp_obs['lon'])):
    if lon%50 == 0:
        print(lon)
    for lat in range(0,len(temp_obs['lat'])):

        if ((L10_anom[lat,lon])!=np.NaN): #if there is data / (if its not nan)
            if np.abs(L10_anom[lat,lon]) > stdev_cntrl[lat,lon]:
                dots_lons.append(temp_obs_data['lon'][lon].values)
                dots_lats.append(temp_obs_data['lat'][lat].values)
                s_sig[lat,lon] = L10_anom[lat,lon]
#        if plot_only_stdev[lat,lon] != np.nan:
 #               dots_lons2.append(temp_obs_data['lon'][lon].values)
  #              dots_lats2.append(temp_obs_data['lat'][lat].values)

        if s_sig[lat,lon] == 0:
                s_sig[lat,lon] = np.nan

fig,ax= plt.subplots(figsize =(18,13),subplot_kw=dict(projection=ccrs.Robinson()))

ax.add_feature(cfeature.LAND, color = 'lightgray')
ax.add_feature(cfeature.COASTLINE)

#levs=np.arange(0,2.71,0.1)

p = ax.contourf(temp_cntrl['lon'],temp_cntrl['lat'],plot_this,levels = 40, transform=ccrs.PlateCarree(),cmap = cmocean.cm.balance)

plt.scatter(dots_lons,dots_lats,color='k',zorder=1,transform = ccrs.PlateCarree())

cbar = plt.colorbar(p, orientation='horizontal', pad=0.05, fraction=0.05,ax=ax)
cbar.ax.tick_params(labelsize=23)
cbar.set_label('ZOOC' + ' (mol/m$^3$)  ', size = 30)
title = "ZOOC-2000-2010, dots are stat_sig"
plt.title(title,size = 35)
fig.savefig('/scratch/b/b380670/spinUpPeriod/'+title + '.jpg' , bbox_inches='tight') # update with save path
#plt.show()
plt.close()

fig,ax= plt.subplots(figsize =(18,13),subplot_kw=dict(projection=ccrs.Robinson()))

ax.add_feature(cfeature.LAND, color = 'lightgray')
ax.add_feature(cfeature.COASTLINE)

#levs=np.arange(0,2.71,0.1)

p = ax.contourf(temp_cntrl['lon'],temp_cntrl['lat'],s_sig,levels = 40, transform=ccrs.PlateCarree(),cmap = cmocean.cm.balance)

#plt.scatter(dots_lons,dots_lats,color='k',zorder=1,transform = ccrs.PlateCarree())

cbar = plt.colorbar(p, orientation='horizontal', pad=0.05, fraction=0.05,ax=ax)
cbar.ax.tick_params(labelsize=23)
cbar.set_label('ZOOC' + ' (mol/m$^3$) ', size = 30)

title = "ZOOC-2000-2010, only_sig"
plt.title(title,size = 35)
fig.savefig('/scratch/b/b380670/spinUpPeriod/'+title + '.jpg' , bbox_inches='tight') # update with save path

plt.close()
