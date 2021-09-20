# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 07:06:46 2021

@author: Rathliob
"""
# ----------------------------------------------------------------------------
# Program timer:
# https://stackoverflow.com/questions/1557571/how-do-i-get-time-of-a-python-programs-execution
import timeit

start = timeit.default_timer()
# ----------------------------------------------------------------------------

#### Skript to prepare all the input data for using in the analysis of managment units 
    
    # List of all Input data:
        # dhm 2m - tif
        # wetlands
        # conservation
        # rockfall zome
        # verdrisk
        # forsteignung
        
        # transport border
        # gemeindegrenzen
        # ecosystem services
        # laubnadel - tif
        # forest
        # perimeter
        # roads - line shape

   
# ----------------------------------------------------------------------------

# Import the data folders and workign folders to the system path
import sys
sys.path.insert(0, r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/mhh')
sys.path.insert(0, r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath')
sys.path.insert(0, r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/Daten/temp')
sys.path.insert(0, r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/Python/python_code_Leo')
sys.path.insert(0, r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_preprocessing/_mhh')
sys.path.insert(0, r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_preprocessing/_ZG_Stw')
sys.path.insert(0, r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Output_preprocessing')


## Load packages and Perimeter

import numpy as np
import geopandas as gpd
from osgeo import gdal, gdalconst
import shutil
import subprocess

from functions_geoprocessing import CreateAndWriteRaster, RasterizePolygon
from functions_geoprocessing import TransferHeaderGeoTiffToHeaderASCII
from functions_geoprocessing import Vector2Raster, WriteASCII
from functions_geometrical_analysis_BWEs import Vector2Raster_with_attribute
from Forest_Inventory_Func import ReadHeaderAndLoadRaster

# ----------------------------------------------------------------------------
# Variables that can be changed to be able to change some simple parameters

delta = 2 # new cellsize

test_area = '_ZG_Stw' # _mhh or '_ZG_Stw' depending on which testing site is used

# in this folder all output data will be written
basic_path_out = r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Output_preprocessing' + '/' + test_area + '_8'

basic_path_in = r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_preprocessing' + '/' + test_area
# ----------------------------------------------------------------------------

## Input Files and paths:

name_conservation     = 'conservation' + test_area + '.shp'
name_wetlands         = 'wetlands' + test_area + '.shp'
name_dhm              = 'dhm_2m' + test_area + '.tif'
name_dhm_copy         = 'dhm_2m' + test_area + '_test.tif'
name_verdrisk         = 'verdrisk' + test_area + '.shp'
name_rockfall         = 'rockfall_zone' + test_area + '.shp'

name_adminborders     = 'gemeindegrenzen_4m' + test_area + '.shp'
name_transportborders = 'transport_borders_4m' + test_area + '.shp'
name_es               = 'ecosystemservices_4m' + test_area + '.shp'
name_forsteignung     = 'forsteignung' + test_area + '.shp'
name_wuchs            = 'wuchs' + test_area + '.shp'
name_laubnadel        = 'laubnadel' + test_area + '.tif'

name_forest           = 'forest_minus' + test_area + '.shp'
name_area             = test_area + '.shp' 

name_roads            = 'roads' + test_area + '_split.shp'

# ----------------------------------------------------------------------------

path_conservation     = basic_path_in + '/' + name_conservation
path_wetlands         = basic_path_in + '/' + name_wetlands
path_dhm              = basic_path_in + '/' + name_dhm
path_dhm_copy         = basic_path_in + '/' + name_dhm_copy
path_verdrisk         = basic_path_in + '/' + name_verdrisk
path_rockfall         = basic_path_in + '/' + name_rockfall

path_adminborders     = basic_path_in + '/' + name_adminborders
path_transportborders = basic_path_in + '/' + name_transportborders
path_es               = basic_path_in + '/' + name_es
path_forsteignung     = basic_path_in + '/' + name_forsteignung
path_wuchs            = basic_path_in + '/' + name_wuchs
path_laubnadel        = basic_path_in + '/' + name_laubnadel

path_forest           = basic_path_in + '/' + name_forest
path_area             = basic_path_in + '/' + name_area

path_roads            = basic_path_in + '/' + name_roads

print('Name and pathgiving for inputfiles done')

# ----------------------------------------------------------------------------
# Einlesen des Strassennetzwerks:


print('Reading of my Road Network')
road_network = gpd.read_file(path_roads)
print('Reading of my Road Network Done')

# ----------------------------------------------------------------------------
# defining the names of the output data all in raster

name_NULLRASTER              = 'NULL_Raster.tif' # raster all zeros to get Nodata value for others

name_RASTER_conservation     = 'conservation' + test_area + '.tif'
name_RASTER_wetlands         = 'wetlands' + test_area + '.tif'
name_RASTER_dhm              = 'dhm2m' + test_area + '.tif'
name_RASTER_slopework        = 'slopework' + test_area + '.tif'
name_RASTER_slope            = 'slope' + test_area + '.tif'
name_ASCII_verdrisk          = 'verdrisk' + test_area + '.asc'
name_RASTER_rockfall         = 'rockfall' + test_area + '.tif'

name_RASTER_adminborders     = 'adminborders' + test_area + '.tif'
name_RASTER_transportborders = 'transportborders' + test_area + '.tif'
name_RASTER_es               = 'es' + test_area + '.tif'
name_ASCII_forsteignung      = 'forsteignung' + test_area + '.asc'
name_ASCII_wuchs             = 'wuchs' + test_area + '.asc'
name_RASTER_laubnadel        = 'laubnadel' + test_area + '.tif'

name_RASTER_area             = 'perimeter' + test_area + '.tif'
name_RASTER_forest           = 'forest' + test_area + '.tif'

name_RASTER_roads            = 'forestroads' + test_area + '.tif'
# ----------------------------------------------------------------------------
path_NULLRASTER              =  basic_path_out + '/' + name_NULLRASTER

path_RASTER_conservation     =  basic_path_out + '/' + name_RASTER_conservation
path_RASTER_wetlands         =  basic_path_out + '/' + name_RASTER_wetlands
path_RASTER_dhm              =  basic_path_out + '/' + name_RASTER_dhm
path_RASTER_slopework        =  basic_path_out + '/' + name_RASTER_slopework
path_RASTER_slope            =  basic_path_out + '/' + name_RASTER_slope
path_ASCII_verdrisk          =  basic_path_out + '/' + name_ASCII_verdrisk
path_RASTER_rockfall         =  basic_path_out + '/' + name_RASTER_rockfall

path_RASTER_adminborders     =  basic_path_out + '/' + name_RASTER_adminborders
path_RASTER_transportborders =  basic_path_out + '/' + name_RASTER_transportborders
path_RASTER_es               =  basic_path_out + '/' + name_RASTER_es
path_ASCII_forsteignung      =  basic_path_out + '/' + name_ASCII_forsteignung
path_ASCII_wuchs             =  basic_path_out + '/' + name_ASCII_wuchs
path_RASTER_laubnadel        =  basic_path_out + '/' + name_RASTER_laubnadel
 
path_RASTER_area             =  basic_path_out + '/' + name_RASTER_area
path_RASTER_forest           =  basic_path_out + '/' + name_RASTER_forest

path_RASTER_roads            =  basic_path_out + '/' + name_RASTER_roads
print('Name and path giving for outputfiles done')
# ----------------------------------------------------------------------------

#### write the important output files as txt files

ins = '.tif'
outs = '.txt'
path_ASCII_conservation     = path_RASTER_conservation.replace(ins, outs)
path_ASCII_wetlands         = path_RASTER_wetlands.replace(ins, outs)
path_ASCII_dhm              = path_RASTER_dhm.replace(ins, outs)
path_ASCII_slope            = path_RASTER_slope.replace(ins, outs)
path_ASCII_rockfall         = path_RASTER_rockfall.replace(ins, outs)

path_ASCII_adminborders     = path_RASTER_adminborders.replace(ins, outs)
path_ASCII_transportborders = path_RASTER_transportborders.replace(ins, outs)
path_ASCII_es               = path_RASTER_es.replace(ins, outs)
path_ASCII_laubnadel        = path_RASTER_laubnadel.replace(ins, outs)

path_ASCII_area             = path_RASTER_area.replace(ins, outs)
path_ASCII_forest           = path_RASTER_forest.replace(ins, outs)

path_ASCII_roads            = path_RASTER_roads.replace(ins, outs)

print('path giving to output asciis done')
##-----------------------------------------------------------------------------
# Open the raster files

raster_dhm_original = gdal.Open(path_dhm)
print('opening dhm raster done')

raster_laubnadel_original = gdal.Open(path_laubnadel)
print('opening laubnadel raster done')
##-----------------------------------------------------------------------------

## Create and Write Empty Reference Raster    

# create an empty raster with the perimeter extent
area_df = gpd.read_file(path_area)

# defining the upper and lower bounds as well as the cell size
bbox = np.array(area_df.total_bounds)
bbox_f = np.floor(bbox / delta)*delta-delta/2.
bbox_c = np.ceil(bbox / delta)*delta+delta/2.

X = np.array((bbox_f[0],bbox_c[2]))
Y = np.array((bbox_f[1],bbox_c[3]))

xmin, ymin, xmax, ymax = [min(X)-delta/2, min(Y)-delta/2, max(X)+delta/2, max(Y)+delta/2]         

anz_x = (xmax - xmin) / delta  
anz_y = (ymax - ymin) / delta           
image_size = (int(anz_y),int(anz_x))
        
v_pixels = np.zeros(image_size).astype(int)                   
x_index = np.arange(anz_x).astype(int)
y_index = np.arange(anz_y).astype(int)
x_values = np.arange(xmin,xmax+delta,delta).astype(int)
y_values = np.arange(ymin,ymax+delta,delta).astype(int)

xres = delta
yres = delta

geotransform = (xmin, xres, 0, ymin, 0, yres)
referenceProj = raster_dhm_original.GetProjection() # defining coordinate system


# create the 1-band raster file
number_of_bands = 1
reference = gdal.GetDriverByName('GTiff').Create(path_NULLRASTER, int(anz_x), int(anz_y), number_of_bands, gdal.GDT_Float32)
reference.SetGeoTransform(geotransform)    # specify coords
reference.GetRasterBand(1).WriteArray(v_pixels)   # write r-band to the raster
reference.GetRasterBand(1).SetNoDataValue(-9999)
reference.SetProjection(referenceProj)    # define coordinate system
reference.FlushCache()                     # write to disk

# =============================================================================
## Resampling: 
# eg turning all the rasters into same extent, same cell size, overlapping cells
# writing -9999 as noData value
## https://gis.stackexchange.com/questions/234022/resampling-a-raster-from-python-without-using-gdalwarp

# define a function that turns the input data raster into an output data raster
# with the right cell size and arrangement, extent and coordinate system

def resampling(raster_original,reference):
    

    inputf = raster_original # gdal.Open(inputfile, gdalconst.GA_ReadOnly)
    inputProj = inputf.GetProjection()
    inputTrans = inputf.GetGeoTransform()
    
    referenceProj = reference.GetProjection()
    referenceTrans = reference.GetGeoTransform()
    bandreference = reference.GetRasterBand(1)  
      
    x = reference.RasterXSize 
    y = reference.RasterYSize
  
    driver= gdal.GetDriverByName('GTiff') # defining format
    output = driver.Create(path_dhm_copy,x,y,1,bandreference.DataType) # creating file path, x,y,z, extent, and band number
    output.SetGeoTransform(referenceTrans) # set transformation
    output.SetProjection(referenceProj) # set coordinate system
    
      # write file: inputfile, output file, input coord, output coord, data type
    gdal.ReprojectImage(inputf,output,inputProj,referenceProj,gdalconst.GRA_Bilinear)
    
    return(output)
# =============================================================================

#### Actually resample the tif data--------------------------------------------
                                 
# Resampling dhm
output = resampling(raster_dhm_original, reference)
# Turning into array
arraydhm = np.array(output.GetRasterBand(1).ReadAsArray())
arraydhm[arraydhm==0] = -9999
# create the variable that represents the geotiff dhm
dhm = CreateAndWriteRaster(path_RASTER_dhm, reference, arraydhm)

print('Resampling dhm finished')

# Resampling laubnadel
output = resampling(raster_laubnadel_original, reference)
# Turning into array
arraylaubnadel = np.array(output.GetRasterBand(1).ReadAsArray())
arraylaubnadel[arraylaubnadel==0] = -9999
# create the variable that represents the geotiff laubnadel
laubnadel = CreateAndWriteRaster(path_RASTER_laubnadel, reference, arraylaubnadel)

print('Resampling laubandel finished')
                                
##-----------------------------------------------------------------------------

## https://gis.stackexchange.com/questions/286055/densify-shapely-polygon-using-ogr-segmentize
## https://gis.stackexchange.com/questions/257109/gdaldem-alpha-flag-in-python-bindings-is-missing/257118

## Compute Slope, turn DEM into slope
                              
slope_as_percent = gdal.DEMProcessing(path_RASTER_slopework, dhm, 'slope', slopeFormat = 'percent')
slope_as_percent.GetRasterBand(1).ReadAsArray()

slope_raster = CreateAndWriteRaster(path_RASTER_slope, reference, slope_as_percent.GetRasterBand(1).ReadAsArray())
array_slope = np.array(slope_raster.GetRasterBand(1).ReadAsArray())
print ('Computing slope from dhm done')
                     

##-----------------------------------------------------------------------------

# Create all the rasters, verdrisk, forsteignung and wuchs will be created in the end

conservation_raster     = RasterizePolygon(path_conservation, reference, path_RASTER_conservation)
wetlands_raster         = RasterizePolygon(path_wetlands, reference, path_RASTER_wetlands)
rockfall_raster         = RasterizePolygon(path_rockfall, reference, path_RASTER_rockfall)

adminborders_raster     = RasterizePolygon(path_adminborders, reference, path_RASTER_adminborders)
transportborders_raster = RasterizePolygon(path_transportborders, reference, path_RASTER_transportborders)
es_raster               = RasterizePolygon(path_es, reference, path_RASTER_es)

# Create rasters for test area and forest and turn into numpy arrays
area_raster             = RasterizePolygon(path_area, reference, path_RASTER_area)
array_area              = np.array(area_raster.GetRasterBand(1).ReadAsArray(),  dtype = 'uint8')

forest_raster           = RasterizePolygon(path_forest, reference, path_RASTER_forest)
array_forest            = np.array(forest_raster.GetRasterBand(1).ReadAsArray(),  dtype = 'uint8')

forestroads_raster      = RasterizePolygon(path_roads, reference, path_RASTER_roads)

print('creation of simple rasters done')

## ----------------------------------------------------------------------------
## Writing of ASCII Files:

print ('Write ASCII Files')

header_info, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value = TransferHeaderGeoTiffToHeaderASCII(reference)

                                
# dhm
pathout = path_ASCII_dhm
arrfl = np.flipud(arraydhm)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# slope
pathout = path_ASCII_slope
arrfl = np.flipud(array_slope)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# laubnadel
pathout = path_ASCII_laubnadel
arrfl = np.flipud(arraylaubnadel)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
                          
# area
pathout = path_ASCII_area
arrfl = np.flipud(array_area)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# forest
pathout = path_ASCII_forest
array_forest[array_area==0] = 0 # reclass to zero outside perimeter
arrfl = np.flipud(array_forest)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# conservation
pathout = path_ASCII_conservation
conservation_raster_0_1 = conservation_raster.GetRasterBand(1).ReadAsArray()
conservation_raster_0_1[conservation_raster_0_1>0] = 1 # reclass if other values than 1 exist

conservation_raster_0_1[array_forest==0] = 0 # reclass to zero outside forest
arrfl = np.flipud(conservation_raster_0_1)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# wetlands
pathout = path_ASCII_wetlands
wetlands_raster_0_1 = wetlands_raster.GetRasterBand(1).ReadAsArray()
wetlands_raster_0_1[wetlands_raster_0_1>0] = 1 # reclass if other values than 1 exist

wetlands_raster_0_1[array_forest==0] = 0  # reclass to zero outside forest
arrfl = np.flipud(wetlands_raster_0_1)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# rockfall
pathout = path_ASCII_rockfall
rockfall_raster_0_1 = rockfall_raster.GetRasterBand(1).ReadAsArray()
rockfall_raster_0_1[rockfall_raster_0_1>0] = 1 # reclass if other values than 1 exist

rockfall_raster_0_1[array_forest==0] = 0 # reclass to zero outside forest
arrfl = np.flipud(rockfall_raster_0_1)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# adminborders
pathout = path_ASCII_adminborders
adminborders_raster_wo = adminborders_raster.GetRasterBand(1).ReadAsArray()

adminborders_raster_wo[array_forest==0] = 0 # reclass to zero outside forest
arrfl = np.flipud(adminborders_raster_wo)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# transport borders
pathout = path_ASCII_transportborders
transportborders_raster_wo = transportborders_raster.GetRasterBand(1).ReadAsArray()

transportborders_raster_wo[array_forest==0] = 0 # reclass to zero outside forest
arrfl = np.flipud(transportborders_raster_wo)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# ecosystem services
pathout = path_ASCII_es
es_raster_wo = es_raster.GetRasterBand(1).ReadAsArray()

es_raster_wo[array_forest==0] = 0 # reclass to zero outside forest
arrfl = np.flipud(es_raster_wo)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# forest roads
pathout = path_ASCII_roads
roads_raster_wo = forestroads_raster.GetRasterBand(1).ReadAsArray()

roads_raster_wo[array_forest==0] = 0 # reclass to zero outside forest
arrfl = np.flipud(roads_raster_wo)
WriteASCII(pathout, arrfl, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

print ('writing of ASCII files done, start QGIS shape to raster')

# =============================================================================
# subprocess.run([r"C:\OSGeo4W64\bin\qgis_process-qgis.bat", 
#                 "run", "gdal:rasterize",
#                 "INPUT="+path_forsteignung,
#                 "FIELD=forstEig",
#                 "BURN=0", "UNITS:1", "WIDTH:2", "HEIGHT:2",
#                 "EXTENT:'2644116.000000000,2649070.000000000,1239858.000000000,1244240.000000000 [EPSG:2056]'",
#                 "NODATA:-9999", "OPTIONS:", "DATA_TYPE:5", "INIT:None",
#                 "INVERT:False", "EXTRA:",
#                 "OUTPUT:"+path_ASCII_forsteignung])
# =============================================================================

# =============================================================================
# subprocess.run([r"C:\OSGeo4W64\bin\qgis_process-qgis.bat", 
#                 "run", "grass7:v.to.rast",
#                 "input:"+path_forsteignung,
#                 "type:[0,1,3]", "where:", "use:0", 
#                 "attribute_column:forstEig",
#                 "rgb_column:", "label_column:", "value:1", "memory:300",
#                 "output:"+ path_ASCII_forsteignung,
#                 "GRASS_REGION_PARAMETER:2644116.000000000,2649070.000000000,1239858.000000000,1244240.000000000 [EPSG:2056]",
#                 "GRASS_REGION_CELLSIZE_PARAMETER:2", "GRASS_RASTER_FORMAT_OPT:",
#                 "GRASS_RASTER_FORMAT_META:", "GRASS_SNAP_TOLERANCE_PARAMETER:-1",
#                 "GRASS_MIN_AREA_PARAMETER:0.0001"])
# =============================================================================

print("preprocessing data done")

# ----------------------------------------------------------------------------
# copy the dhm copy file back from the original, cause the copy is broken during the run

source = path_dhm
destination = path_dhm_copy
shutil.copy(source, destination)
print("DHM_test file restored for next run")

# copy all the output txt files into the BWE input folder
folder_BWEin = r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_BWEs' + '/' + test_area

source = path_ASCII_conservation
destination = folder_BWEin + '/' + name_RASTER_conservation.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_wetlands
destination = folder_BWEin + '/' + name_RASTER_wetlands.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_dhm
destination = folder_BWEin + '/' + name_RASTER_dhm.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_slope
destination = folder_BWEin + '/' + name_RASTER_slope.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_rockfall
destination = folder_BWEin + '/' + name_RASTER_rockfall.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_adminborders
destination = folder_BWEin + '/' + name_RASTER_adminborders.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_transportborders
destination = folder_BWEin + '/' + name_RASTER_transportborders.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_es
destination = folder_BWEin + '/' + name_RASTER_es.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_laubnadel
destination = folder_BWEin + '/' + name_RASTER_laubnadel.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_area
destination = folder_BWEin + '/' + name_RASTER_area.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_forest
destination = folder_BWEin + '/' + name_RASTER_forest.replace(ins, outs)
shutil.copy(source, destination)

source = path_ASCII_roads
destination = folder_BWEin + '/' + name_RASTER_roads.replace(ins, outs)
shutil.copy(source, destination)

source = path_NULLRASTER
destination = folder_BWEin + '/' + name_NULLRASTER
shutil.copy(source, destination)


# =============================================================================
# source = path_ASCII_verdrisk
# destination = folder_BWEin + '/' + name_ASCII_verdrisk
# shutil.copy(source, destination)
# 
# source = path_ASCII_wuchs
# destination = folder_BWEin + '/' + name_ASCII_wuchs
# shutil.copy(source, destination)
# 
# source = path_ASCII_forsteignung
# destination = folder_BWEin + '/' + name_ASCII_forsteignung
# shutil.copy(source, destination)
# =============================================================================
print("Output ASCIIs copied into BWE Input folder")

# ----------------------------------------------------------------------------
# Program timer end
stop = timeit.default_timer()
execution_time = stop - start

print("Execution time: " + str(round(execution_time,2)) + "s") # It returns time in seconds
# ----------------------------------------------------------------------------