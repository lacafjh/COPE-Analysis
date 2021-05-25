import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# nc_obj, netCDF object
nc_obj = nc.Dataset("20130703.c1.nc")

# Flight track drawing start
# Create latitude and longitude data https://bit.ly/31LbYDS
# df_ll, DataFrame latlon
# pd.DataFrame(np.array([nc_obj['time'][:], nc_obj['LAT'][:],
#                        nc_obj['LON'][:]]).T,
#              columns=['time', 'LAT', 'LON'])
df_ll = pd.DataFrame({
    'time': np.array(nc_obj['time'][:]),
    'LAT': np.array(nc_obj['LAT'][:]),
    'LON': np.array(nc_obj['LON'][:])
})
vars_ll = ['LAT', 'LON']
for i in vars_ll:
    df_ll = df_ll[df_ll[i] != -32767.0]

# define map extent
extent = [-10, 5, 45, 60]

# define state borders
states_borders = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='50m',
    facecolor='none')

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

# create figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# Add features
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(states_provinces, edgecolor='gray')
ax.add_feature(states_borders, edgecolor='black')
# plot data
ax.plot(df_ll['LON'],
        df_ll['LAT'],
        'o',
        transform=ccrs.PlateCarree(),
        linewidth=0.1)
ax.set_extent(extent)
gl = ax.gridlines(draw_labels=True, linestyle=":", linewidth=0.3, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 5))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, 5))
gl.xlabel_style = {'size': 12}
gl.ylabel_style = {'size': 12}
ax.spines['geo'].set_linewidth(0.5)
ax.set_title('Flight track 2013.7.3', fontsize=14)

plt.show()

# time drawing start
df_T = pd.DataFrame({
    'TIME': np.array(nc_obj['TIME'][:]),
    'lwc100': np.array(nc_obj['lwc100'][:]),
    'trf': np.array(nc_obj['trf'][:]),
    'GALT': np.array(nc_obj['GALT'][:]),
    'hivs': np.array(nc_obj['hivs'][:])
})
df_T = df_T[df_T['TIME'] > 91500][df_T['TIME'] < 92500]

fig, ax = plt.subplots(3, 1, figsize=(8, 12))
plt.subplots_adjust(wspace=0, hspace=0.3)

# TIME - LWC
ax[0].plot(df_T['TIME'], df_T['lwc100'])
# ax[0].set_xlim(1E-7, 1E1)
ax[0].set_xlabel('TIME')
ax[0].set_ylabel('Liquid water content (DMT100) /$(g·m^{-3})$')
ax[0].set_title('(a) TIME - LWC')

# TIME - Temperature & Altitude
ax_alt = ax[1].twinx()
ax[1].plot(df_T['TIME'], df_T['trf'], color='#DC143C')
ax_alt.plot(df_T['TIME'], df_T['GALT'])
# ax[1].set_xlim(1E-7, 1E1)
ax[1].set_xlabel('TIME')
ax[1].set_ylabel('Static temperature /°C')
ax_alt.set_ylabel('Altitude  (Ashtech GPS) /$m$')
ax[1].set_title('(b) TIME - Temperature & Altitude')

# TIME - Vertical speed
ax[2].plot(df_T['TIME'], df_T['hivs'])
# ax[0].set_xlim(1E-7, 1E1)
ax[2].set_xlabel('TIME')
ax[2].set_ylabel('Inertial vertical speed (Honeywell) /$(m·s^{-1})$')
ax[2].set_title('(c) TIME - Vertical speed')

plt.show()

# LWC drawing start
# df_lwc, liquid water content data in DataFrame
df_lwc = pd.DataFrame({
    'lwc100': np.array(nc_obj['lwc100'][:]),
    'GALT': np.array(nc_obj['GALT'][:]),
    'trf': np.array(nc_obj['trf'][:])
})
vars_lwc = ['lwc100', 'GALT', 'trf']
for i in vars_lwc:
    df_lwc = df_lwc[df_lwc[i] != -32767.0]
df_lwc = df_lwc[df_lwc['lwc100'] > 0.01]

# caculate means, asisted by @jifenghi
# df_lwc.sort_values(), ranking; df_lwc.query(), grouping;
# df_lwc.mean(), mean
# ctrf, centre of trf levels; cgalt, centre of galt levels
means_lwc_trf = []
bins = [-7.5, -2.5, 2.5, 6, 8, 10, 12, 14, 17.5]
ctrf = [-5, 0, 5, 7, 9, 11, 13, 15]
for i in range(len(bins) - 1):
    mean = df_lwc.sort_values(by=['trf']).query("trf >= {} & trf <= {}".format(
        bins[i], bins[i + 1])).mean(axis=0)
    means_lwc_trf.append(mean['lwc100'])
means_lwc_galt = []
bins = [0, 700, 1200, 2000, 3350, 4050, 4750, 5250, 6000]
cgalt = [500, 900, 1300, 3000, 3700, 4400, 5100, 5800]
for i in range(len(bins) - 1):
    mean = df_lwc.sort_values(by=['GALT']).query(
        "GALT >= {} & GALT <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_lwc_galt.append(mean['lwc100'])

fig, ax = plt.subplots(1, 2, figsize=(12, 8))

# LWC - Altitude plot
ax[0].scatter(df_lwc['lwc100'], df_lwc['GALT'], s=12, color='#66CCFF')
ax[0].scatter(means_lwc_galt, cgalt, s=64, color='#0066FF')
ax[0].set_xscale('log')
ax[0].set_xlim(1E-2, 1E1)
ax[0].set_xlabel('Liquid water content (DMT100) /$(g·m^{-3})$')
ax[0].set_ylabel('Altitude  (Ashtech GPS) /$m$')
ax[0].set_title('(a) LWC - Altitude')

# LWC - Temperature plot
ax[1].scatter(df_lwc['lwc100'], df_lwc['trf'], s=12, color='#66CCFF')
ax[1].scatter(means_lwc_trf, ctrf, s=64, color='#0066FF')
ax[1].set_xscale('log')
ax[1].set_xlim(1E-2, 1E1)
ax[1].invert_yaxis()
ax[1].set_xlabel('Liquid water content (DMT100) /$(g·m^{-3})$')
ax[1].set_ylabel('Static temperature /°C')
ax[1].set_title('(b) LWC - Temperature')

plt.show()

# pd.pivot_table(df_lwc, values='lwc100',index=[u'trf'],
#                aggfunc={'lwc100': np.mean})
# pd.cut(df_lwc['trf'], 6, labels=[-7.5, -2.5, 2.5, 7.5, 12.5, 17.5])
# Temperature bins
# tbins = [-7.5, -2.5, 2.5, 7.5, 12.5, 17.5]

# FSSP droplet concentration
# df_DSSP_conc, FSSP droplet concentration data in DataFrame
df_DSSP_conc = pd.DataFrame({
    'lwc100': np.array(nc_obj['lwc100'][:]),
    'jlb_conc2_IBL': np.array(nc_obj['jlb_conc2_IBL'][:]),
    'GALT': np.array(nc_obj['GALT'][:]),
    'trf': np.array(nc_obj['trf'][:])
})
vars_DSSP_conc = ['lwc100', 'jlb_conc2_IBL', 'GALT', 'trf']
for i in vars_DSSP_conc:
    df_DSSP_conc = df_DSSP_conc[df_DSSP_conc[i] != -32767.0]
df_DSSP_conc = df_DSSP_conc[df_DSSP_conc['jlb_conc2_IBL'] != 0.0]
df_DSSP_conc = df_DSSP_conc[df_DSSP_conc['lwc100'] > 0.01]
df_DSSP_conc = df_DSSP_conc[df_DSSP_conc['trf'] < 12]

# caculate means, asisted by @jifenghi
# df_lwc.sort_values(), ranking; df_lwc.query(), grouping;
# df_lwc.mean(), mean
# ctrf, centre of trf levels; cgalt, centre of galt levels
means_DSSP_conc_trf = []
bins = [5, 6.5, 7.5, 8.5, 9.5, 10.5]
ctrf = [6, 7, 8, 9, 10]
for i in range(len(bins) - 1):
    mean = df_DSSP_conc.sort_values(by=['trf']).query(
        "trf >= {} & trf <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_DSSP_conc_trf.append(mean['jlb_conc2_IBL'])
means_DSSP_conc_galt = []
bins = [600, 775, 925, 1075, 1225, 1400]
cgalt = [700, 850, 1000, 1150, 1300]
for i in range(len(bins) - 1):
    mean = df_DSSP_conc.sort_values(by=['GALT']).query(
        "GALT >= {} & GALT <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_DSSP_conc_galt.append(mean['jlb_conc2_IBL'])

fig, ax = plt.subplots(1, 2, figsize=(12, 8))

# FSSP droplet concentration (JLB) method 2 - Altitude plot
ax[0].scatter(df_DSSP_conc['jlb_conc2_IBL'],
              df_DSSP_conc['GALT'],
              s=12,
              color='#66CCFF')
ax[0].scatter(means_DSSP_conc_galt, cgalt, s=64, color='#0066FF')
ax[0].set_xscale('log')
ax[0].set_xlim(1E-1, 1E4)
ax[0].set_xlabel('FSSP droplet concentration (JLB) method 2 /$(cm^{-3})$')
ax[0].set_ylabel('Altitude  (Ashtech GPS) /$m$')
ax[0].set_title('(a) Droplet concentration - Altitude')

# FSSP droplet concentration (JLB) method 2 - Temperature plot
ax[1].scatter(df_DSSP_conc['jlb_conc2_IBL'],
              df_DSSP_conc['trf'],
              s=12,
              color='#66CCFF')
ax[1].scatter(means_DSSP_conc_trf, ctrf, s=64, color='#0066FF')
ax[1].set_xscale('log')
ax[1].set_xlim(1E-1, 1E4)
ax[1].invert_yaxis()
ax[1].set_xlabel('FSSP droplet concentration (JLB) method 2 /$(cm^{-3})$')
ax[1].set_ylabel('Static temperature /°C')
ax[1].set_title('(b) Droplet concentration - Temperature')

plt.show()

# Fast2DC concentration
# df_Fast2DC_conc, Fast2DC concentration data in DataFrame
df_Fast2DC_conc = pd.DataFrame({
    'lwc100': np.array(nc_obj['lwc100'][:]),
    'jlb_conc2_IBL': np.array(nc_obj['jlb_conc2_IBL'][:]),
    'CONC0_cip_IBR': np.array(nc_obj['CONC0_cip_IBR'][:]),
    'GALT': np.array(nc_obj['GALT'][:]),
    'trf': np.array(nc_obj['trf'][:])
})
vars_Fast2DC_conc = ['lwc100', 'CONC0_cip_IBR', 'GALT', 'trf']
for i in vars_Fast2DC_conc:
    df_Fast2DC_conc = df_Fast2DC_conc[df_Fast2DC_conc[i] != -32767.0]
df_Fast2DC_conc = df_Fast2DC_conc[df_Fast2DC_conc['jlb_conc2_IBL'] > 2.0]
df_Fast2DC_conc = df_Fast2DC_conc[df_Fast2DC_conc['lwc100'] > 0.01]
df_Fast2DC_conc = df_Fast2DC_conc[df_Fast2DC_conc['GALT'] < 1500]
df_Fast2DC_conc = df_Fast2DC_conc[df_Fast2DC_conc['trf'] < 12][
    df_Fast2DC_conc['trf'] > 5]

# caculate means, asisted by @jifenghi
# df_lwc.sort_values(), ranking; df_lwc.query(), grouping;
# df_lwc.mean(), mean
# ctrf, centre of trf levels; cgalt, centre of galt levels
means_Fast2DC_conc_galt = []
bins = [600, 850, 950, 1050, 1150, 1250, 1400]
cgalt = [700, 900, 1000, 1100, 1200, 1300]
for i in range(len(bins) - 1):
    mean = df_Fast2DC_conc.sort_values(by=['GALT']).query(
        "GALT >= {} & GALT <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_Fast2DC_conc_galt.append(mean['CONC0_cip_IBR'])
means_Fast2DC_conc_trf = []
bins = [5, 6.25, 6.75, 7.25, 7.75, 8.25, 12]
ctrf = [6, 6.5, 7, 7.5, 8, 9.5]
for i in range(len(bins) - 1):
    mean = df_Fast2DC_conc.sort_values(by=['trf']).query(
        "trf >= {} & trf <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_Fast2DC_conc_trf.append(mean['CONC0_cip_IBR'])

fig, ax = plt.subplots(1, 2, figsize=(12, 8))

# Fast2DC (aka CIP) concentration - Altitude plot
ax[0].scatter(df_Fast2DC_conc['CONC0_cip_IBR'],
              df_Fast2DC_conc['GALT'],
              s=12,
              color='#66CCFF')
ax[0].scatter(means_Fast2DC_conc_galt, cgalt, s=64, color='#0066FF')
ax[0].set_xscale('log')
ax[0].set_xlim(1E-2, 1E4)
ax[0].set_xlabel('Fast2DC (aka CIP) concentration /$(L^{-1})$')
ax[0].set_ylabel('Altitude  (Ashtech GPS) /$m$')
ax[0].set_title('(a) Droplet concentration - Altitude')

# Fast2DC (aka CIP) concentration - Temperature plot
ax[1].scatter(df_Fast2DC_conc['CONC0_cip_IBR'],
              df_Fast2DC_conc['trf'],
              s=12,
              color='#66CCFF')
ax[1].scatter(means_Fast2DC_conc_trf, ctrf, s=64, color='#0066FF')
ax[1].set_xscale('log')
ax[1].set_xlim(1E-2, 1E4)
ax[1].invert_yaxis()
ax[1].set_xlabel('Fast2DC (aka CIP) concentration /$(L^{-1})$')
ax[1].set_ylabel('Static temperature /°C')
ax[1].set_title('(b) Droplet concentration - Temperature')

plt.show()

# LWC - Wind vertical component
# df_LWC_hw, LWC - Vertical speed data in DataFrame
df_LWC_hw = pd.DataFrame({
    'lwc100': np.array(nc_obj['lwc100'][:]),
    'jlb_conc2_IBL': np.array(nc_obj['jlb_conc2_IBL'][:]),
    'CONC0_cip_IBR': np.array(nc_obj['CONC0_cip_IBR'][:]),
    'hw': np.array(nc_obj['hw'][:])
})
vars_LWC_hw = ['lwc100', 'jlb_conc2_IBL', 'CONC0_cip_IBR', 'hw']
for i in vars_LWC_hw:
    df_LWC_hw = df_LWC_hw[df_LWC_hw[i] != -32767.0]
df_LWC_hw = df_LWC_hw[df_LWC_hw['lwc100'] > 0.01]
df_LWC_hw = df_LWC_hw[df_LWC_hw['jlb_conc2_IBL'] > 2.0]

# caculate means, asisted by @jifenghi
# df_lwc.sort_values(), ranking; df_lwc.query(), grouping;
# df_lwc.mean(), mean
# ctrf, centre of trf levels; cgalt, centre of galt levels
means_LWC_hw = []
bins = [-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75]
chw = [-1, -0.5, 0, 0.5, 1, 1.5]
for i in range(len(bins) - 1):
    mean = df_LWC_hw.sort_values(by=['hw']).query(
        "hw >= {} & hw <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_LWC_hw.append(mean['lwc100'])

fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# LWC - Vertical speed plot
ax.scatter(df_LWC_hw['lwc100'], df_LWC_hw['hw'], s=12, color='#66CCFF')
ax.scatter(means_LWC_hw, chw, s=64, color='#0066FF')
ax.set_xscale('log')
ax.set_xlim(1E-2, 1E1)
ax.set_xlabel('Liquid water content (DMT100) /$(g·m^{-3})$')
ax.set_ylabel('Wind vertical component (Honeywell) /$(m·s^{-1})$')
ax.set_title('LWC - Wind vertical component')

plt.show()

# concentration - Vertical speed
# caculate means, asisted by @jifenghi
# df_lwc.sort_values(), ranking; df_lwc.query(), grouping;
# df_lwc.mean(), mean
# ctrf, centre of trf levels; cgalt, centre of galt levels
means_LWC_hw_FSSP = []
bins = [-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75]
cFSSP = [-1, -0.5, 0, 0.5, 1, 1.5]
for i in range(len(bins) - 1):
    mean = df_LWC_hw.sort_values(by=['hw']).query(
        "hw >= {} & hw <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_LWC_hw_FSSP.append(mean['jlb_conc2_IBL'])
means_LWC_hw_Fast2DC = []
bins = [-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75]
cFast2DC = [-1, -0.5, 0, 0.5, 1, 1.5]
for i in range(len(bins) - 1):
    mean = df_LWC_hw.sort_values(by=['hw']).query(
        "hw >= {} & hw <= {}".format(bins[i], bins[i + 1])).mean(axis=0)
    means_LWC_hw_Fast2DC.append(mean['CONC0_cip_IBR'])

fig, ax = plt.subplots(1, 2, figsize=(12, 8))

# FSSP droplet concentration (JLB) method 2 - Vertical speed plot
ax[0].scatter(df_LWC_hw['jlb_conc2_IBL'],
              df_LWC_hw['hw'],
              s=12,
              color='#66CCFF')
ax[0].scatter(means_LWC_hw_FSSP, cFSSP, s=64, color='#0066FF')
ax[0].set_xscale('log')
ax[0].set_xlim(1E0, 1E3)
ax[0].set_xlabel('FSSP droplet concentration (JLB) method 2 /$(cm^{-3})$')
ax[0].set_ylabel('Wind vertical component (Honeywell) /$(m·s^{-1})$')
ax[0].set_title('(a) FSSP droplet concentration - Vertical speed')

# Fast2DC (aka CIP) concentration - Vertical speed plot
ax[1].scatter(df_LWC_hw['CONC0_cip_IBR'],
              df_LWC_hw['hw'],
              s=12,
              color='#66CCFF')
ax[1].scatter(means_LWC_hw_Fast2DC, cFast2DC, s=64, color='#0066FF')
ax[1].set_xscale('log')
ax[1].set_xlim(1E-2, 1E4)
ax[1].set_xlabel('Fast2DC (aka CIP) concentration /$(L^{-1})$')
ax[1].set_ylabel('Wind vertical component (Honeywell) /$(m·s^{-1})$')
ax[1].set_title('(b) Fast2DC concentration - Vertical speed')

plt.show()
