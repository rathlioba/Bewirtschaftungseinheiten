# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 14:02:57 2021

## Geometrical Analysis of Forest Road Networks:

@author: Rathliob


"""
# ----------------------------------------------------------------------------
# Program timer:
# https://stackoverflow.com/questions/1557571/how-do-i-get-time-of-a-python-programs-execution
import timeit

start = timeit.default_timer()
# ----------------------------------------------------------------------------

# Input files used:
    # for defining nodes of the network:
        # raster of forest area
    # for defining trafficability raster case:
        # slope
        # conservation
        # wetlands
        # rockfall 
        # verdrisk
    # for defining transport borders case
        # transport_borders
    # for defining admin borders case
        # admin_borders
    # for defining ecosystem services case
        # es
        
    # for aggregating the units based on forest mixture type
        # laubnadel
    # for aggregating the units based on growth
        # wuchs / forsteignung
    # for aggregating the units based on transport to forest edge
        # landing points forest edge

## Remark: All Input Data must have the same Coordinate System! (EPSG:2056)

####--------------------------------------------------------------------------

import sys
sys.path
sys.path.append('N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/Python/Skripte_MA/GitLab_clone')
sys.path.append('N:/forema/FPS/Projekte_der_Gruppe/Erschliessung/Aargau/Data/Vector_Input')
sys.path.append('N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/Python/python_code_Leo')
sys.path.append('N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_preprocessing/_mhh')
sys.path.append('N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_preprocessing/_ZG_Stw')
sys.path.append('N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Output_preprocessing')
sys.path.append('C:/OSGeo4W64/apps/Python37/lib')
sys.path.append('C:/OSGeo4W64/apps/qgis/bin')
sys.path.append('C:/OSGeo4W64/bin')
# sys.path.append ('path to environment')
# print(sys.path)

import geopandas as gpd
import numpy as np,sys
from Forest_Inventory_Func import ReadHeaderAndLoadRaster
from functions_geoprocessing import WriteASCII# , CreateAndWriteRaster
import igraph
import subprocess
# from osgeo import gdal
import functions_geometrical_analysis_BWEs as fga

####--------------------------------------------------------------------------
# Methods evaluated:
## igraph
## network x

## https://graph-tool.skewed.de/performance
## https://www.timlrx.com/2019/05/05/benchmark-of-popular-graph-network-packages/

# --> igraph is used! 

#### XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# CHANGE THESE VARIABLES IN THE BEGINNING
    
# path to run qgis subprocesses: Liobas Laptop
qgis_path = r"C:\Program Files\QGIS 3.16\bin\qgis_process-qgis.bat" 

maxDistanzBefahrbar = 300 # <- distance in m, as threshold how far the skidding is allowed 
maxSlope = 30
maxZuzug = '35' # maximum distance for Bodenzug. defines buffer around groundbased and forestroads
maxSlopeZuzug = 300

test_area = '_ZG_Stw'  # or 'ZG_Stw'

out_bem = test_area + '_admin+transport'  # tag in the name of the output files to know the version


#### XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Decision variables for network analysis criteria
# !!!!Attention!!!!
# Do not set any variable to 'True' that would produce an empty raster. 
# These cases for now are: 
    # es for _mhh
    # 

# Criterion 1: Trafficability
all_terrain       = False
ground_based      = False
winch_assist      = False

# Criterion 2: Transport borders
transport_borders = True

# Criterion 3: Administrative borders
admin_borders     = True

# Criterion 4: Borders between ecosystem services
es_borders        = False


# Decision variables for aggregation

# Criterion 5: based on timber transport
forest_edge       = True

# Criterion 6: based on forest growth or production suitability
production_suit   = False # this is the standard
forest_growth     = False

# Criterion 7a: based on separating forest mixture types
separate_mixture  = False
# Criterion 7b: based on integrating forest mixture types 
integrate_mixture = False

# Criterion 8: based on number of BWEs
max_amount        = False
max_size          = False

#### XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# =============================================================================

## Define In & Output Directories:
folder_out = r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Output_BWEs' + '/' + test_area
folder_in = r'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_BWEs' + '/' + test_area

## Define Names of Input files:
filename_forest           = 'forest' + test_area + '.txt'
name_NULLRASTER           = 'NULL_raster.tif'
filename_perimeter        = 'perimeter' + test_area + '.txt'

filename_slope            = 'slope' + test_area + '.txt'
filename_conservation     = 'conservation' + test_area + '.txt'
filename_wetlands         = 'wetlands' + test_area + '.txt'
filename_rockfall         = 'rockfall' + test_area + '.txt'
filename_verdrisk         = 'verdrisk' + test_area + '.asc'

filename_adminborders     = 'adminborders' + test_area + '.txt'
filename_transportborders = 'transportborders' + test_area + '.txt'
filename_es               = 'es' + test_area + '.txt'

filename_edgepoints       = 'edgepoints' + test_area + '.shp'
filename_growth           = 'forsteignung' + test_area + '.asc'
if forest_growth:
    filename_growth       = 'wuchs' + test_area + '.asc'
filename_mixture          = 'laubnadel' + test_area + '.txt'

filename_forestroads      = 'roads' + test_area + '_split.shp'
# filename_RASTER_forestroads='forestroads' + test_area + '.txt'
# --------------------------------------------------------------------------
# Define file paths to all input files

path_forest           = folder_in + '/' + filename_forest
path_NULLRASTER       = folder_in + '/' + name_NULLRASTER
path_perimeter        = folder_in + '/' + filename_perimeter

path_slope            = folder_in + '/' + filename_slope
path_conservation     = folder_in + '/' + filename_conservation
path_wetlands         = folder_in + '/' + filename_wetlands
path_rockfall         = folder_in + '/' + filename_rockfall
path_verdrisk         = folder_in + '/' + filename_verdrisk

path_adminborders     = folder_in + '/' + filename_adminborders
path_transportborders = folder_in + '/' + filename_transportborders
path_es               = folder_in + '/' + filename_es

path_edgepoints       = folder_in + '/' + filename_edgepoints
path_growth           = folder_in + '/' + filename_growth
path_mixture          = folder_in + '/' + filename_mixture

path_forestroads      = folder_in + '/' + filename_forestroads
# path_RASTER_forestroads=folder_in + '/' + filename_RASTER_forestroads
# --------------------------------------------------------------------------
# Define paths for output files

path_ASCII_groundbased     = folder_out + '/' + 'groundbased' + out_bem + '.txt'
path_ASCII_BWEarea         = folder_out + '/' + 'BWE_area'+out_bem+'.txt'
path_ASCII_networkroads    = folder_out + '/' + 'network_roads'+out_bem+'.txt'
path_ASCII_trafficability  = folder_out + '/' + 'trafficability'+out_bem+'.txt'
path_ASCII_connected       = folder_out + '/' + 'map_connected'+out_bem+'.txt'
path_ASCII_skiddist        = folder_out + '/' + 'skid_dist'+out_bem+'.txt'
path_ASCII_mapunits        = folder_out + '/' + 'map_units'+out_bem+'.txt'

path_ASCII_roaddist        = folder_out + '/' + 'roaddist' + out_bem + '.txt'
path_ASCII_bweunits        = folder_out + '/' + 'BWEunits' + out_bem + '.txt'

# intermediate results
path_ASCII_forestroads_buffer80 = folder_out + '/' + 'forestroads_buffer80' + out_bem + '.asc'
path_ASCII_groundbased_buffer80 = folder_out + '/' + 'groundbased_buffer80' + out_bem + '.asc'

# intermediate rasters written for control
path_ASCII_forestroads_withID = folder_out + '/raster_forestroads_withID' + out_bem + '.txt'
path_ASCII_networkinput       = folder_out + '/network_input' + out_bem + '.txt'
path_ASCII_edgepoints         = folder_out + '/edgepoints' + out_bem + '.txt'

# --------------------------------------------------------------------------

# Read in all the files

print('Load Raster Data')    
header_forest, raster_forest                     = ReadHeaderAndLoadRaster(path_forest,True, True)
# header_forestroads, raster_forestroads           = ReadHeaderAndLoadRaster(path_RASTER_forestroads, True, True)
header_perimeter, raster_perimeter               = ReadHeaderAndLoadRaster(path_perimeter, True, True)

header_slope, raster_slope                       = ReadHeaderAndLoadRaster(path_slope,True, True)
header_conservation, raster_conservation         = ReadHeaderAndLoadRaster(path_conservation,True, True)
header_wetlands, raster_wetlands                 = ReadHeaderAndLoadRaster(path_wetlands,True, True)
header_rockfall, raster_rockfall                 = ReadHeaderAndLoadRaster(path_rockfall,True, True)
header_verdrisk, raster_verdrisk                 = ReadHeaderAndLoadRaster(path_verdrisk,True, True)

header_adminborders, raster_adminborders         = ReadHeaderAndLoadRaster(path_adminborders,True, True)
header_transportborders, raster_transportborders = ReadHeaderAndLoadRaster(path_transportborders,True, True)
header_es, raster_es                             = ReadHeaderAndLoadRaster(path_es,True, True)

header_growth, raster_growth                     = ReadHeaderAndLoadRaster(path_growth,True, True)
header_mixture, raster_mixture                   = ReadHeaderAndLoadRaster(path_mixture,True, True)
print('Load Raster Data finished')    

## Get Raster Shape:
raster_shape = raster_forest.shape

## Get Raster Cellsize:
cellsize = header_forest['cellsize']

##----------------------------------------------------------------------------
## Read in Road Network
print('Load Road Network')  
gdf_forestroads = gpd.read_file(path_forestroads)
print('Loading Road Network finished')  

## Org Strassen Index
    # again turns the road network into a raster, this time with the raster cells
    # keeping the OBJECTID attribute value in the cell
    # is used in the map of landings, as the nodes that fall on the road are the 
    # end-nodes for the network analysis

attribute = 'OBJECTID'
raster_forestroads_withid = fga.Vector2Raster_with_attribute(gdf_forestroads, header_forest[0],attribute)

# write forestroads with ID
pathout = path_ASCII_forestroads_withID
ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value = header_forest[0]
WriteASCII(pathout, raster_forestroads_withid, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

##----------------------------------------------------------------------------
## Read in network of Forest Edgepoints
print('Load Forest Edgepoints')  
gdf_edgepoints = gpd.read_file(path_edgepoints)
print('Loading Forest Edgepoints finished')  

# turn the edgepoints into a raster that has the same format aseverythign else
# keep the ID value in the raster
# use as end nodes in the second network aggregation analysis

attribute = 'OBJECTID'
raster_edgepoints_withid = fga.Vector2Raster_with_attribute_polygon(gdf_edgepoints, header_forest[0],attribute)

# write edgepoint raster with IDs
pathout = path_ASCII_edgepoints
ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value = header_forest[0]
WriteASCII(pathout, raster_edgepoints_withid, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

####===========================================================================

#### change input for the geometrical analysis with cases.
# geometrical analysis input that is targeted -> 'network_input'

# while going through the cases the network_input is updated constantly. 
# if one variable is true, the network_input used in the next case already includes that previous variable
# if no variable is true then network_input stays the baseline of all forest area

# !! from original script
# !! changed 'befahrbar' to 'network_input'
# !! changed 'befahrbare_ind' to 'input_ind'

print("Create input raster for network")
####-----------------------------------
## baseline case 0: everything within the forest
network_input = raster_forest

####-----------------------------------
## case 1: BWEs according to method of harvesting

    # create rasters for trafficability analysis.
    # three options: 
        # all terrain without conservation and wetlands area is included
        # only terrain harvestable with ground based machines
        # terrain reachable with a winch or ground based methods
        
raster_allterrain = np.logical_and(raster_forest == 255, 
                                   np.logical_and(raster_wetlands != 1, raster_conservation != 1)) 

raster_groundbased_1 = np.logical_and(raster_allterrain == 1, 
                                      np.logical_and(raster_slope <= maxSlope,
                                                     raster_slope >= 0))  
raster_groundbased_2 = np.logical_and(raster_verdrisk != 0,
                                      np.logical_and(raster_verdrisk != 4, raster_verdrisk != 5))
raster_groundbased_3 = np.logical_and(raster_groundbased_1 == 1, raster_groundbased_2 == 1)
raster_groundbased   = np.logical_or(raster_groundbased_3 == 1, raster_forestroads_withid != 0)

####----
print("buffer roads and groundbased")
# Write groundbased as ASCII file for buffer input
# have to turn bool array into 0_1 int array first!! 
raster_groundbased_0_1 = raster_groundbased.astype(int)

# Write area where groundbased transport is possible
pathout = path_ASCII_groundbased
ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value = header_forest[0]
WriteASCII(pathout, raster_groundbased_0_1, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)  

subprocess.run([qgis_path, 
                "run", "grass7:r.buffer", "input="+path_ASCII_forestroads_withID, "distances="+maxZuzug, "units=0", 
                "-z=True", "output="+path_ASCII_forestroads_buffer80, "GRASS_REGION_PARAMETER=", 
                "GRASS_REGION_CELLSIZE_PARAMETER=0", "GRASS_RASTER_FORMAT_OPT=", 
                "GRASS_RASTER_FORMAT_META="], shell=True)
subprocess.run([qgis_path, 
                "run", "grass7:r.buffer", "input="+path_ASCII_groundbased, "distances="+maxZuzug, "units=0", 
                "-z=True", "output="+path_ASCII_groundbased_buffer80, "GRASS_REGION_PARAMETER=", 
                "GRASS_REGION_CELLSIZE_PARAMETER=0", "GRASS_RASTER_FORMAT_OPT=", 
                "GRASS_RASTER_FORMAT_META="], shell=True)

# read .asc files back as arrays
header_forestroads_buffer80, raster_forestroads_buffer80 = ReadHeaderAndLoadRaster(path_ASCII_forestroads_buffer80,True, True)
header_groundbased_buffer80, raster_groundbased_buffer80 = ReadHeaderAndLoadRaster(path_ASCII_groundbased_buffer80,True, True)
####----

raster_winchbuffer_1 = np.logical_or(raster_forestroads_buffer80 == 2, raster_groundbased_buffer80 == 2)
raster_winchbuffer_2 = np.logical_or(raster_winchbuffer_1 == 1, raster_groundbased_0_1 == 1)
raster_winchbuffer_3 = np.logical_and(raster_allterrain == 1, 
                                   np.logical_and(raster_rockfall != 1, 
                                                  np.logical_and(raster_slope >= maxSlope, 
                                                                 raster_slope <= maxSlopeZuzug)))
raster_winchassist = np.logical_and(raster_winchbuffer_2 == 1, raster_winchbuffer_3 == 0)

if all_terrain:
   network_input = np.logical_and(network_input != 0, raster_allterrain == 1)
   pathout = path_ASCII_BWEarea
   WriteASCII(pathout, network_input +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

                      
if ground_based:
   network_input = np.logical_and(network_input != 0, raster_groundbased == 1)
   pathout = path_ASCII_BWEarea
   WriteASCII(pathout, network_input +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

    
if winch_assist:
   network_input = np.logical_and(network_input != 0, raster_winchassist == 1)
   pathout = path_ASCII_BWEarea
   WriteASCII(pathout, network_input +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

####-----------------------------------
# case 2: transport borders shouldnt be crossed by BWEs

if transport_borders:
    raster_transportborders = raster_transportborders.astype(int)
    network_input = np.logical_and(network_input != 0, raster_transportborders != 255)
    pathout = path_ASCII_BWEarea
    WriteASCII(pathout, network_input +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
    
####-----------------------------------    
# case 3: administrative borders shouldnt be crossed by BWEs

if admin_borders:
    raster_adminborders = raster_adminborders.astype(int)
    network_input = np.logical_and(network_input != 0, raster_adminborders != 255)
    pathout = path_ASCII_BWEarea
    WriteASCII(pathout, network_input +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

####-----------------------------------    
# case 4: borders between ecosystem services shouldnt be crossed by BWEs

if es_borders:
    raster_es = raster_es.astype(int)
    network_input = np.logical_and(network_input != 0, raster_es != 255)  
    pathout = path_ASCII_BWEarea
    WriteASCII(pathout, network_input +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
    

# Turn network input raster into a boolean raster :D 
network_input = network_input!= 0

print('preparation end****************************')

####--------------------------------------------------------------------------
print("Create Road split based on forest edgepoints")
## Aggregation of the map_of_landings units based on further cases

# case 0: timber transport: second network analysis
if forest_edge:
    
    # find shortest path from every point of the road to a set on forst edge points
    # network input (up until now is the forest) will be the road network as raster/array
    # end points (up until now the road network) will be set of forest edge points

    # replace network_input with the forest road network
    header_forestroads, raster_forestroads = ReadHeaderAndLoadRaster(path_ASCII_forestroads_withID, True, True)
    # turn forestroad raster into boolean raster
    raster_forestroads = raster_forestroads != 0

    network_input_roads = raster_forestroads
    
    # if transport border, es borders or admin borders are chosen, these also have to cut the roads! 
    # so that there is no transport movement across a border by the road. 
    # open case still for allterrain, groundbased and winchassist! 
    
# =============================================================================
#     if transport_borders:
#         network_input_roads = np.logical_and(network_input_roads != 0, raster_transportborders != 255)
# =============================================================================
    
    if admin_borders:
        network_input_roads = np.logical_and(network_input_roads != 0, raster_adminborders != 255)
        
    if es_borders:
        network_input_roads = np.logical_and(network_input_roads != 0, raster_es != 255)

    # replace the road network with the edgepoints

    network_points = raster_edgepoints_withid != 0

# run the network analysis
    print('Start edgepoint run of network analysis')
    dist_raster_roads, bweunits_raster = fga.RoadToForestEdgeNetwork(network_input_roads, network_points, raster_forest, raster_edgepoints_withid, raster_shape, cellsize)
    
    print('Finish edgepoint run of network analysis')
    
    # write txt of the distances from each node to the end node
    pathout = path_ASCII_roaddist
    WriteASCII(pathout, dist_raster_roads, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

    # write txt of the end node OBJECTID (from road network) for each forest node
    pathout = path_ASCII_bweunits
    WriteASCII(pathout, bweunits_raster, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)


####===========================================================================
print("prepare network for forest to road analysis")

# Prepare base dataset for the nodes in the network

## Outputs the linear Index of the map that is to be used as base for network analysis
# find(argument) returns a list of linear indices where argument is true
input_ind = fga.find(network_input)

# define input for network, only forest roads! 
# should be the same as raster_forestroads in boolean
network_roads = np.logical_and(raster_forest != 0 , raster_forestroads_withid != 0)
# turn to linear coordinates
networkroads_ind = fga.find(network_roads)

unique_bwe_ind = np.unique(input_ind)
red_ind_bwe = np.arange(len(unique_bwe_ind))
tf = fga.ismember(input_ind, networkroads_ind)
red_ind_roads = red_ind_bwe[tf]

# prepare the position and weights of the network edges

L = len(input_ind)
lll = len(red_ind_roads)
lplus = np.ones([lll]).astype('int')*(L) 


V1_N,V2_N,L_N = fga.defineNeighbourhood(network_input,input_ind,N8 = True)

V1 = np.concatenate((lplus,red_ind_roads,V1_N))   
V2 = np.concatenate((red_ind_roads,lplus,V2_N));


weights = np.ones(V1.shape)
weights[:lll*2]=0  #0.001
weights[lll*2:]=L_N

Edges =np.concatenate(([[V1],[V2]])).T
Edges_List = Edges.tolist()
##-----------------------------------------------------------------------------
print ("build network for forest to road analysis")
# build the nodes and edges of the network

g = igraph.Graph()
g.add_vertices(len(Edges_List))
g.add_edges(Edges_List)
g.es["weight"] = weights
#m = g.get_adjacency()
#exit()
#g.write_pickle(r"N:\forema\FPS\Projekte_der_Gruppe\Masterarbeit_Lioba_Rath\testdaten\Output_BWEs\_ZG_Stw\g.pkl")
dist_list = g.shortest_paths_dijkstra(source=L, target=None, weights="weight", mode = 'OUT')
path_list = g.get_shortest_paths(L)

##-----------------------------------------------------------------------------    
## shortest path

"""
shortest_paths_dijkstra(source=None, target=None, weights=None, mode=OUT)
source code 
Calculates shortest path lengths for given vertices in a graph.

The algorithm used for the calculations is selected automatically: a simple BFS is used for unweighted graphs, Dijkstra's algorithm is used when all the weights are positive. Otherwise, the Bellman-Ford algorithm is used if the number of requested source vertices is larger than 100 and Johnson's algorithm is used otherwise.

Parameters:
source - a list containing the source vertex IDs which should be included in the result. If None, all vertices will be considered.
target - a list containing the target vertex IDs which should be included in the result. If None, all vertices will be considered.
weights - a list containing the edge weights. It can also be an attribute name (edge weights are retrieved from the given attribute) or None (all edges have equal weight).
mode - the type of shortest paths to be used for the calculation in directed graphs. OUT means only outgoing, IN means only incoming paths. ALL means to consider the directed graph as an undirected one.
Returns:
the shortest path lengths for given vertices in a matrix

"""

##-----------------------------------------------------------------------------
print("Compute results")
## Distance Computations:
    # results in an array (txt) of the length of the shortest path from every node to the end node

dist =  np.array(dist_list)   

dist_finite = np.isfinite(dist)
dist_finite_values = dist[dist_finite]

dist_raster = np.zeros(raster_shape)-9999

ind_sol_1 = fga.find(dist_finite)
Lin_Index_Distance_Raster = input_ind[ind_sol_1[:-1]]
dist_raster.ravel()[Lin_Index_Distance_Raster] = dist_finite_values[:-1]*cellsize

##-----------------------------------------------------------------------------
## Map Trafficable:
    # returns an array (txt) of all nodes in the trafficable area
    # that have a connection to an end node (excludes the isolated ones)

ind_sol = fga.find(dist < maxDistanzBefahrbar/cellsize)
linCoord_sol = input_ind[ind_sol[0:-1]]

map_connected  = np.zeros(raster_shape)
map_connected.ravel()[linCoord_sol]=1

##-----------------------------------------------------------------------------
## Map of Landings:
    # returns an array (txt) based on the nodes of the network, where every node
    # when turned into a cell retains the value of the end node (on the road) 
    # it reaches with the shortest path

Org_Number_Landing = np.zeros(len(ind_sol)-1).astype('int')
 
for q in range(len(ind_sol)-1):
    Org_Number_Landing[q] =  path_list[ind_sol[q]][1]
    
LinIndexLanding = input_ind[Org_Number_Landing]
    
landing_raster = np.zeros(raster_shape)-9999

landing_raster.ravel()[linCoord_sol] = Org_Number_Landing  ## Oder LinIndexLanding

Org_Number_Landing_ID = raster_forestroads_withid.ravel()[input_ind[Org_Number_Landing]]
if forest_edge:
    Org_Number_Landing_ID = bweunits_raster.ravel()[input_ind[Org_Number_Landing]]
    
landing_raster.ravel()[linCoord_sol] = Org_Number_Landing_ID  ## Oder LinIndexLanding

##-----------------------------------------------------------------------------
# create raster that has allt he classes of trafficable are in the forest
# 0 = outside perimeter, outside forest
# 1 = roads
# 2 = terrain reachable with groundbased methods
# 3 = terrain reachable with a winch from groundbased
# 4 = forest terrain outside winch reach or groundbased area 

trafficability = network_roads.astype(int)
trafficability[raster_allterrain == 1] = 4
trafficability[raster_winchassist == 1] = 3
trafficability[raster_groundbased_0_1 == 1] = 2
trafficability[raster_forestroads_withid != 0] = 1
trafficability[raster_forest == 0] = 0
trafficability[raster_perimeter == 0] = 0

##-----------------------------------------------------------------------------
print ("Write ASCII files")
## write all output for control:

# Write area where groundbased transport is possible
pathout = path_ASCII_groundbased
WriteASCII(pathout, raster_groundbased_0_1, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)    

# write BWE area area as txt
pathout = path_ASCII_BWEarea
WriteASCII(pathout, network_input +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# # write road raster that is used as network input
network_roads = network_roads.astype(int)
pathout = path_ASCII_networkroads
WriteASCII(pathout, network_roads, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# write txt of trafficable area that has a connection to an end node
pathout = path_ASCII_connected
WriteASCII(pathout, map_connected +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# write raster that has classes of trafficbility
pathout = path_ASCII_trafficability
WriteASCII(pathout, trafficability +0, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# write txt of the distances from each node to the end node
pathout = path_ASCII_skiddist
WriteASCII(pathout, dist_raster, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

# write txt of the end node OBJECTID (from road network) for each forest node
pathout = path_ASCII_mapunits
WriteASCII(pathout, landing_raster, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)


# ----------------------------------------------------------------------------
# Program timer end
stop = timeit.default_timer()
execution_time = stop - start

print("Execution time: " + str(round(execution_time,2)) + "s") # It returns time in seconds
# ----------------------------------------------------------------------------


####--------------------------------------------------------------------------                
    # case 0: timber transport: second network analysis
        # connect all road-nodes to a shapefile set of points on the border of the forest
            # connect those end points to a virtual end point again, run shortest path
            # same as before, create a map which forest road node ends up at which
            # border of the forest point
            # aggregate the previous smaller units (based on same road segments) 
            # into bigger units (based on same border of the forest point)  
            
            # replace network and roads input
            # netork_input = np.logical_and(raster_forestroads != 0, raster_forest != 0)
            # network_roads = RasterizePolygon(path_to_point_layer, header_forest[0])
            
            # only compute the skid dist raster and then the raster of landings
            # attribute = 'ID'
            # raster_points_on_forest_edge = fga.Vector2Raster_with_attribute(gdf_points_on_forest_edge, header_forest[0],attribute)
            # 
            # make the base for the trafficbale
            # ind_sol = fga.find(dist < maxDistanzBefahrbar/cellsize)
            # linCoord_sol = input_ind[ind_sol[0:-1]]
            # landings map
            # Org_Number_Landing = np.zeros(len(ind_sol)-1).astype('int')
 
            # for q in range(len(ind_sol)-1):
            #     Org_Number_Landing[q] =  path_list[ind_sol[q]][1]
    
            # LinIndexLanding = input_ind[Org_Number_Landing]
    
            # landing_raster = np.zeros(raster_shape)-9999

            # landing_raster.ravel()[linCoord_sol] = Org_Number_Landing  ## Oder LinIndexLanding

            # Org_Number_Landing_ID = raster_forestroads_withid.ravel()[input_ind[Org_Number_Landing]]
            # landing_raster.ravel()[linCoord_sol] = Org_Number_Landing_ID  ## Oder LinIndexLanding

            
            
            # reshape input map to represent the points including ID on forest roads
            # reshape the output map to contain only certain points at the forest edge
            # run analysis
            # map of skid dist
            # map of landings as in: all road points going to the same end point get
                # end point ID, which is then written in the map of landings from before
# further option to aggregate            
    # case 1: tree growth: 
        # input raster of forest area classified by growth strength (wuchs/eignung) (4 classes)
        # aggregate adjacent units to bigger units based on an even mixture 
        # between growth classes
        
        # overlay growth map with the current landings
        # calculate area percent of growth categories per landing
        # combine neighbors based on equal bad growth percentages
        
    # case 2a: deciduous vs coniferous balanced
        # input raster of forest area classified into forest mixture type
        # aggregate adjacent units to bigger units based on an even mixture 
        # between forest types
        
        # overlay forest mix map with landings
        # calculate percentage of deciduous per landing
        # combine neighbors bason on getting as close to 50 percent as possible
        
    # case 2b: deciduous vs coniferous separated
        # input raster of forest area classified into forest mixture type
        # aggregate adjacent units to bigger units based on a separation 
        # between forest types     
        
        # combine neighbors either maximizing or minimizing the percentage
    
    # case 3: number of units
        # a) aggregate adjacent units up until a certain size
        # b) or aggregate them to reach a certain number in the area, 
            # all with +/- similar sizes
        # c) number of units is based on the rotation period of interventions
        
        # calculate size of all landings
        # a) combine neighbors based up to max size
        # b) calculate average size from max number, aggregate neighbors up to that size
        # c) calculate max number from rotation period, calculate max size, then aggregate

####--------------------------------------------------------------------------














