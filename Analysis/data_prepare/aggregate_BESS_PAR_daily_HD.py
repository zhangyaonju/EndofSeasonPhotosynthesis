import multiprocessing
#import h5py
import os
import numpy as np
import fnmatch
import sys
from netCDF4 import Dataset
from skimage.measure import block_reduce

lat = np.arange(-90 + 0.5 / 2., (90 + 1e-8) - 0.5 / 2., 0.5)
lon = np.arange(-180 + 0.5 / 2., (180 + 1e-8) - 0.5 / 2., 0.5)
y = np.arange(2001,2017)
def exportpar(predY,outfile):
    nday = predY.shape[0]
    print('creating: '+outfile)
    f = Dataset(outfile, 'w', format='NETCDF4')
    lati = f.createDimension('lat', len(lat))
    loni = f.createDimension('lon', len(lon))
    yearday = f.createDimension('yday', nday)
    latitudes = f.createVariable('lat', 'f4', 'lat')
    longitudes = f.createVariable('lon', 'f4', 'lon')
    yearmonth = f.createVariable('yday', 'f4', 'yday')

    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    latitudes.standard_name = 'latitude'
    longitudes.standard_name = 'longitude'
    latitudes.axis = 'Y'
    longitudes.axis = 'X'
    latitudes.long_name = 'latitude'
    longitudes.long_name = 'longitude'
    latitudes[:] = lat
    longitudes[:] = lon
    yearmonth[:] = np.arange(nday)+1
    bess_par = f.createVariable('PAR', 'f4', ('yday','lat','lon'),zlib=True)
    bess_par.missing_value=-999.9
    bess_par [:] = predY.astype('f4').reshape(nday,360,720)
    f.close()

def get_bess_par(y):
    bess_year = fnmatch.filter(bessfiles,'*A'+str(y)+'*.nc')
    if bess_year is None:
        return None
    #print('bess file: ')
    #print(bess_year)
    bess_hd_year = []
    for ibess in bess_year:
        with Dataset(ibess,mode='r') as fh:
            try:
                parday = fh.variables['surface_downwelling_photosynthetic_radiative_flux_in_air'][:,:]*1.0
                print('bess radiation shape: '+str(parday.shape))
                parday[parday<=-999]=np.nan
            except:
                print('cannot read bess file: '+ibess)
        bessrad = np.flipud(np.array(parday.T))
        bess_hd_year += [block_reduce(bessrad,block_size=(10,10),func = np.nanmean)]
        #print('radiaiton shape: '+str(bess_hd_mon.shape))
    bess_hd_year = np.array(bess_hd_year)
    exportpar(bess_hd_year,"/rigel/glab/users/zy2309/DATA/BESS_daily_HD/"+str(y)+".bess_hd_daily_PAR.nc")
    #return bess_hd_mon

def getallfile (directory):
    fileList=[]
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if ".nc4" == name[-4:] or '.nc' == name[-3:]: fileList.append(os.path.join(path,name))
    #print fileList
    fileList.sort()
    return fileList

bess_dir = '/rigel/glab/users/Datasets/BESS/BESSRadiation/realtimedata.snu.ac.kr/BESSRadiation/BESS_PAR_Daily/'

bessfiles = getallfile(bess_dir)

#bess_HD = np.zeros((360,720,17*12))
#for i in range(17*12):
#    bess_HD[:,:,i] = get_bess_par(ym[i])

#exportpar(bess_HD,"/rigel/glab/users/zy2309/DATA/bess_hd_monthly_PAR.nc")


pool = multiprocessing.Pool(8)
pool.map(get_bess_par,y)

