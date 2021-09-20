# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 08:40:53 2021

@author: Rathliob
"""

# This script contains the functions that are used in the geometrical network
# analysis for the BWEs


# =============================================================================
import geopandas as gpd
import numpy as np,sys
from osgeo import gdal


# =============================================================================

## Required functions are listed below:

    
# https://gis.stackexchange.com/questions/250555/buffering-around-raster-using-gdal-and-numpy

def raster_buffer(raster_filepath, dist=80):
     """This function creates a distance buffer around the given raster file with non-zero values.
     The value in output raster will have value of the cell to which it is close to."""
     d=gdal.Open(raster_filepath)
     if d is None:
         print("Error: Could not open image ")
         sys.exit(1)
     global proj,geotrans,row,col
     proj=d.GetProjection()
     geotrans=d.GetGeoTransform()
     row=d.RasterYSize
     col=d.RasterXSize
     inband=d.GetRasterBand(1)
     in_array = inband.ReadAsArray(0,0,col,row).astype(int)
     Xcell_size=int(abs(geotrans[1]))
     Ycell_size=int(abs(geotrans[5]))
     cell_size = (Xcell_size+Ycell_size)/2
     cell_dist=dist/cell_size
     in_array[in_array == (inband.GetNoDataValue() or 0 or -999)]=0
     out_array=np.zeros_like(in_array)
     temp_array=np.zeros_like(in_array)
     i,j,h,k=0,0,0,0
     print("Running distance buffer...")
     while(h<col):
         k=0
         while(k<row): 
             if(in_array[k][h]>=1):
                 i=int(h-cell_dist)
                 while((i<cell_dist+h) and i<col):
                     j=int(k-cell_dist)
                     while(j<(cell_dist+k) and j<row):
                         if(((i-h)**2+(j-k)**2)<=cell_dist**2):
                             if(temp_array[j][i]==0 or temp_array[j][i]>((i-h)**2+(j-k)**2)):
                                 out_array[j][i]= in_array[k][h]
                                 temp_array[j][i]=(i-h)**2+(j-k)**2
                         j+=1
                     i+=1
             k+=1
         h+=1
     d,temp_array,in_array=None,None, None
     return out_array

def export_array(in_array,output_path):
    """This function is used to produce output of array as a tif map."""
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(output_path,col,row,1)
    outband=outdata.GetRasterBand(1)  
    outband.SetNoDataValue(np.nan)
    outband.WriteArray(in_array)
    # Georeference the image
    outdata.SetGeoTransform(geotrans)
    # Write projection information
    outdata.SetProjection(proj)       
    outdata.FlushCache()
    outdata = None

## change to linear coordinates:
def sub2ind(array_shape, rows, cols):
    return rows*array_shape[1] + cols

## change back to raster
def ind2sub(array_shape, ind):
    ## rows
    #rows = (ind.astype('int') / array_shape[1])
    rows = np.floor(ind.astype('int') / array_shape[1]).astype('int')
    ## columns
    cols = (ind.astype('int') % array_shape[1]) # or numpy.mod(ind.astype('int'), array_shape[1])
    return (rows, cols)


def RealCoordinatesToGridIndex( coordX,coordY, IG_info):

    #nCols = IG.info[1);
    nRows = IG_info[1];
    xllcorner = IG_info[2];
    yllcorner = IG_info[3];
    cellsize = IG_info[4];
    
    # Python
    XGridCoord = (coordX - xllcorner)/cellsize -0.5 ;
    YGridCoord = (-1*coordY + yllcorner + cellsize * (nRows-1))/ cellsize +0.5;

    return(XGridCoord,YGridCoord)


def GiveLinearCoords( cellsize_new,IG_info,X,Y ):
    # GIVELINEARCOORDS Summary of this function goes here
    #   Detailed explanation goes here
    
    # Calculates the linear coordinatees based on the IG.Info
    # cellsize can be adjusted at will

    # INPUT:
    # cellsize_new:     resolution Raster [M]
    # IG                Info about terrain characteristics (DEM?) 
    # X                 X Koordinate , Real
    # Y                 Y Koordinate , Real
    
    # OUTPUT:
    # LinCoord          linear coordinate, relative to IG with cellsize new
    
    cellsize_org = IG_info[5];
    fact = cellsize_org / cellsize_new;
    IG_1m_info = np.zeros(5)
    IG_1m_info[0] = IG_info[0]*fact;   # nCols
    IG_1m_info[1] = IG_info[1]*fact;   # nRows
    IG_1m_info[2] = IG_info[2];      # xllcorner
    IG_1m_info[3] = IG_info[4];      # yllcorner
    IG_1m_info[4] = cellsize_new;    # cellsize
    
    nRows_1m = IG_1m_info[2];
    nCols_1m = IG_1m_info[1];
    size_1m = [nRows_1m,nCols_1m];
    
    [XGridCoord,YGridCoord] = RealCoordinatesToGridIndex(X,Y,IG_1m_info);
    
    LinCoord = sub2ind(size_1m, round(YGridCoord), round(XGridCoord));
    
    return (LinCoord)


def GridIndexToRealCoordinates( XGridCoord,YGridCoord,IG_info ):

    nRows = IG_info[1];
    xllcorner = IG_info[2];
    yllcorner = IG_info[3];
    cellsize = IG_info[4];

    # Matlab
    coordX = xllcorner + cellsize * (XGridCoord + 0.5);
    coordY = yllcorner + cellsize * (nRows-1) - cellsize * (YGridCoord - 0.5);

    return(coordX,coordY)


def find(raster_tf):
    res = np.where(raster_tf)
    lin_ind = sub2ind(raster_tf.shape, res[0], res[1])
    return(lin_ind)


def segmentize(geom):
    
    from osgeo import ogr
    from shapely.wkt import loads

    
    wkt = geom.wkt  # shapely Polygon to wkt
    geom = ogr.CreateGeometryFromWkt(wkt)  # create ogr geometry
    geom.Segmentize(2)  # densify geometry
    wkt2 = geom.ExportToWkt()  # ogr geometry to wkt
    new = loads(wkt2)  # wkt to shapely Polygon
    return new

def Vector2Raster_II(gdf, header):

        
    gdf['geometry'] = gdf['geometry'].map(segmentize)
    
    for j in range(len(gdf)):
        x_coords = np.array(gdf['geometry'][j].coords.xy[0])
        y_coords = np.array(gdf['geometry'][j].coords.xy[1])
        if j == 0:
            x_collector = x_coords
            y_collector = y_coords
        else:
            x_collector = np.concatenate((x_collector,x_coords))
            y_collector = np.concatenate((y_collector,y_coords))
                
    points = np.concatenate(([x_collector],[y_collector]),axis = 0).T
    
    cellsize = header[4]
    xmin = header[2]
    ymin = header[3]
    ncols = header[0]
    nrows = header[1]
    xmax = xmin+ncols*cellsize
    ymax = ymin+nrows*cellsize
    anz_y = nrows
    anz_x = ncols
    
    values_within_boundaries = np.logical_and(np.logical_and(points[:,0] > xmin, points[:,0] < xmax),np.logical_and(points[:,1] > ymin, points[:,1] < ymax))
        
    points = points[values_within_boundaries,:]
        
    points_ind_x=np.floor((points[:,0]-xmin)/cellsize).astype(int)
    points_ind_y=np.floor((points[:,1]-ymin)/cellsize).astype(int)
    
    points_ind_lin_ind = sub2ind([anz_x,anz_y], points_ind_x, points_ind_y)
    ## Unique
    points_ind_lin_ind_unique = np.unique(points_ind_lin_ind)
    ## Write to Raster
    rows, cols = ind2sub([anz_x,anz_y], points_ind_lin_ind_unique) 
    image_size = (int(anz_y),int(anz_x))     
    v_pixels = np.zeros(image_size).astype(int)  

    v_pixels[cols,rows] = 1
    
    v_pixels = np.flipud(v_pixels)
    
    return(v_pixels)
                

def Vector2Raster_with_attribute(gdf, header,attribute):

    # attribute = 'OBJECTID'
    # gdf = gdf_forest_road
    # header = header_forest[0]
        
    gdf['geometry'] = gdf['geometry'].map(segmentize)
    
    for j in range(len(gdf)):
        # for use in Polygons instead of lines:
        # x_coords = np.array(gdf['geometry'][j].exterior.coords.xy[0])
        # y_coords = np.array(gdf['geometry'][j].exterior.coords.xy[1])
        x_coords = np.array(gdf['geometry'][j].coords.xy[0])
        y_coords = np.array(gdf['geometry'][j].coords.xy[1])
        #attr = np.array([gdf[attribute][j]])
        attr =np.ones(x_coords.shape,dtype='int8')*gdf[attribute][j]
        if j == 0:
            x_collector = x_coords
            y_collector = y_coords
            attr_col = attr
        else:
            x_collector = np.concatenate((x_collector,x_coords))
            y_collector = np.concatenate((y_collector,y_coords))
            attr_col =    np.concatenate((attr_col,attr))
            
                
    points = np.concatenate(([x_collector],[y_collector]),axis = 0).T
    
    cellsize = header[4]
    xmin = header[2]
    ymin = header[3]
    ncols = header[0]
    nrows = header[1]
    xmax = xmin+ncols*cellsize
    ymax = ymin+nrows*cellsize
    anz_y = nrows
    anz_x = ncols
    
    values_within_boundaries = np.logical_and(np.logical_and(points[:,0] > xmin, points[:,0] < xmax),np.logical_and(points[:,1] > ymin, points[:,1] < ymax))
        
    points = points[values_within_boundaries,:]
    attr_col = attr_col[values_within_boundaries]
        
    points_ind_x=np.floor((points[:,0]-xmin)/cellsize).astype(int)
    points_ind_y=np.floor((points[:,1]-ymin)/cellsize).astype(int)
    
    points_ind_lin_ind = sub2ind([anz_x,anz_y], points_ind_x, points_ind_y)
    ## Unique
    #points_ind_lin_ind_unique = np.unique(points_ind_lin_ind)
    
    points_ind_lin_ind_unique, indices = np.unique(points_ind_lin_ind, return_index=True)
    attr_col_unique = attr_col[indices]
    ## Write to Raster
    rows, cols = ind2sub([anz_x,anz_y], points_ind_lin_ind_unique) 
    image_size = (int(anz_y),int(anz_x))     
    v_pixels = np.zeros(image_size).astype(int)  

    #v_pixels[cols,rows] = 1
    v_pixels[cols,rows] = attr_col_unique
    
    v_pixels = np.flipud(v_pixels)
    
    return(v_pixels)
        
def Vector2Raster_with_attribute_polygon(gdf, header,attribute):

    # attribute = 'OBJECTID'
    # gdf = gdf_forest_road
    # header = header_forest[0]
        
    gdf['geometry'] = gdf['geometry'].map(segmentize)
    
    for j in range(len(gdf)):
        # for use in Polygons instead of lines:
        # x_coords = np.array(gdf['geometry'][j].exterior.coords.xy[0])
        # y_coords = np.array(gdf['geometry'][j].exterior.coords.xy[1])
        x_coords = np.array(gdf['geometry'][j].exterior.coords.xy[0])
        y_coords = np.array(gdf['geometry'][j].exterior.coords.xy[1])
        #attr = np.array([gdf[attribute][j]])
        attr =np.ones(x_coords.shape,dtype='int8')*gdf[attribute][j]
        if j == 0:
            x_collector = x_coords
            y_collector = y_coords
            attr_col = attr
        else:
            x_collector = np.concatenate((x_collector,x_coords))
            y_collector = np.concatenate((y_collector,y_coords))
            attr_col =    np.concatenate((attr_col,attr))
            
                
    points = np.concatenate(([x_collector],[y_collector]),axis = 0).T
    
    cellsize = header[4]
    xmin = header[2]
    ymin = header[3]
    ncols = header[0]
    nrows = header[1]
    xmax = xmin+ncols*cellsize
    ymax = ymin+nrows*cellsize
    anz_y = nrows
    anz_x = ncols
    
    values_within_boundaries = np.logical_and(np.logical_and(points[:,0] > xmin, points[:,0] < xmax),np.logical_and(points[:,1] > ymin, points[:,1] < ymax))
        
    points = points[values_within_boundaries,:]
    attr_col = attr_col[values_within_boundaries]
        
    points_ind_x=np.floor((points[:,0]-xmin)/cellsize).astype(int)
    points_ind_y=np.floor((points[:,1]-ymin)/cellsize).astype(int)
    
    points_ind_lin_ind = sub2ind([anz_x,anz_y], points_ind_x, points_ind_y)
    ## Unique
    #points_ind_lin_ind_unique = np.unique(points_ind_lin_ind)
    
    points_ind_lin_ind_unique, indices = np.unique(points_ind_lin_ind, return_index=True)
    attr_col_unique = attr_col[indices]
    ## Write to Raster
    rows, cols = ind2sub([anz_x,anz_y], points_ind_lin_ind_unique) 
    image_size = (int(anz_y),int(anz_x))     
    v_pixels = np.zeros(image_size).astype(int)  

    #v_pixels[cols,rows] = 1
    v_pixels[cols,rows] = attr_col_unique
    
    v_pixels = np.flipud(v_pixels)
    
    return(v_pixels)
        

def ismember(A, B):
    tf = np.in1d(A, B)
    return (tf)


def defineNeighbourhood(network_input,input_ind,N8 = True):

    # if N8 --> 8 neighbourhood, else 4 neighbourhood
    
    s = network_input.shape
    
    M = np.zeros(s).astype('int')
    M.ravel()[input_ind] = np.arange(len(input_ind))
    
    ##-----------------------------------------------------------------------------
    # neighbourhood: 4 neighbourhood!
    
    B1 = np.zeros(s).astype('int')
    B1[:-1,:] = network_input[1:,:]
    
    B2 = np.zeros(s).astype('int')
    B2[1:,:] = network_input[:-1,:]
    
    B3 = np.zeros(s).astype('int')
    B3[:,:-1] = network_input[:,1:]
    
    B4 = np.zeros(s).astype('int')
    B4[:,1:] = network_input[:,:-1]
    

    M1 = np.zeros(s).astype('int')
    M1[:-1,:] = M[1:,:]
    
    M2 = np.zeros(s).astype('int')
    M2[1:,:] = M[:-1,:]
    
    M3 = np.zeros(s).astype('int')
    M3[:,:-1] = M[:,1:]
    
    M4 = np.zeros(s).astype('int')
    M4[:,1:] = M[:,:-1]
    
    N1 = find(np.logical_and(network_input,B1))
    N2 = find(np.logical_and(network_input,B2))
    N3 = find(np.logical_and(network_input,B3))
    N4 = find(np.logical_and(network_input,B4))
    
    V1_4 = np.concatenate((M.ravel()[N1],M1.ravel()[N1],M.ravel()[N2],M2.ravel()[N2],M.ravel()[N3],M3.ravel()[N3],M.ravel()[N4],M4.ravel()[N4]))
    V2_4 = np.concatenate((M1.ravel()[N1],M.ravel()[N1],M2.ravel()[N2],M.ravel()[N2],M3.ravel()[N3],M.ravel()[N3],M4.ravel()[N4],M.ravel()[N4]));
    L_4 = np.zeros(V1_4.shape)+1
        ##-------------------------------------------------------------------
        # neighbourhood: 8 neighbourhood!
    if N8:
        
        B5 = np.zeros(s).astype('int')
        B5[:-1,:-1] = network_input[1:,1:]
        
        B6 = np.zeros(s).astype('int')
        B6[1:,:-1] = network_input[:-1,1:]
        
        B7 = np.zeros(s).astype('int')
        B7[1:,1:] = network_input[:-1,:-1]
        
        B8 = np.zeros(s).astype('int')
        B8[:-1,1:] = network_input[1:,:-1]
        
        
        
        M5 = np.zeros(s).astype('int')
        M5[:-1,:-1] = M[1:,1:]
        
        M6 = np.zeros(s).astype('int')
        M6[1:,:-1] = M[:-1,1:]
        
        M7 = np.zeros(s).astype('int')
        M7[1:,1:] = M[:-1,:-1]
        
        M8 = np.zeros(s).astype('int')
        M8[:-1,1:] = M[1:,:-1]
        
        N5 = find(np.logical_and(network_input,B5))
        N6 = find(np.logical_and(network_input,B6))
        N7 = find(np.logical_and(network_input,B7))
        N8 = find(np.logical_and(network_input,B8))
        
        V1_8 = np.concatenate((M.ravel()[N5],M5.ravel()[N5],M.ravel()[N6],M6.ravel()[N6],M.ravel()[N7],M7.ravel()[N7],M.ravel()[N8],M8.ravel()[N8]))
        V2_8 = np.concatenate((M5.ravel()[N5],M.ravel()[N5],M6.ravel()[N6],M.ravel()[N6],M7.ravel()[N7],M.ravel()[N7],M8.ravel()[N8],M.ravel()[N8]));
        L_8 = np.zeros(V1_8.shape)+np.sqrt(2)
    
        V1_N = np.concatenate((V1_4,V1_8))
        V2_N = np.concatenate((V2_4,V2_8))
        L_N = np.concatenate((L_4,L_8))
        
    else:
        V1_N = V1_4
        V2_N = V2_4
        L_N = L_4


    return(V1_N,V2_N,L_N)

# =============================================================================




####===========================================================================
# function to run the network analysis a second time
# basically calculating the distance form a road network inside the forest to
    # a fix set of storage site at the forest edge. 
    # calculate the shortest path to the forest edge, along the forest network

def RoadToForestEdgeNetwork(network_input, network_roads, raster_forest, raster_edgepoints_withid, raster_shape, cellsize):
    
    import igraph    

    print("prepare second network run for aggregation")
    
    # network_inout = roads in the forest
    # network_roads = network of collection points at the forest edge
    
    # Prepare base dataset for the nodes in the network
    
    ## Outputs the linear Index of the road entwork that is to be used as base for network analysis
    # find(argument) returns a list of linear indices where argument is true
    input_ind = find(network_input)
    
    # define input for network, only forest edgepoints! 
    # turn to linear coordinates
    networkroads_ind = find(network_roads)
    
    unique_bwe_ind = np.unique(input_ind)
    red_ind_bwe = np.arange(len(unique_bwe_ind))
    tf = ismember(input_ind, networkroads_ind)
    red_ind_roads = red_ind_bwe[tf]
    
    # prepare the position and weights of the network edges
    
    L = len(input_ind)
    lll = len(red_ind_roads)
    lplus = np.ones([lll]).astype('int')*(L) 
    
    
    V1_N,V2_N,L_N = defineNeighbourhood(network_input,input_ind,N8 = True)
    
    V1 = np.concatenate((lplus,red_ind_roads,V1_N))   
    V2 = np.concatenate((red_ind_roads,lplus,V2_N));
    
    
    weights = np.ones(V1.shape)
    weights[:lll*2]=0  #0.001
    weights[lll*2:]=L_N
    
    Edges =np.concatenate(([[V1],[V2]])).T
    Edges_List = Edges.tolist()
    ##-----------------------------------------------------------------------------
    print ("build network")
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
    
    dist_raster_roads = np.zeros(raster_shape)-9999
    
    ind_sol_1 = find(dist_finite)
    Lin_Index_Distance_Raster = input_ind[ind_sol_1[:-1]]
    dist_raster_roads.ravel()[Lin_Index_Distance_Raster] = dist_finite_values[:-1]*cellsize
    
    ##-----------------------------------------------------------------------------
    ## Map Trafficable:
        # returns an array (txt) of all nodes in the trafficable area
        # that have a connection to an end node (excludes the isolated ones)
    
    ind_sol = find(dist_finite)
    linCoord_sol = input_ind[ind_sol[:-1]]
    
    ##-----------------------------------------------------------------------------
    ## Map of Landings:
        # returns an array (txt) based on the nodes of the network, where every node
        # when turned into a cell retains the value of the end node (on the road) 
        # it reaches with the shortest path
    
    Org_Number_Landing = np.zeros(len(ind_sol)-1).astype('int')
     
    for q in range(len(ind_sol)-1):
        Org_Number_Landing[q] =  path_list[ind_sol[q]][1]
        
    # LinIndexLanding = input_ind[Org_Number_Landing]
        
    bweunit_raster = np.zeros(raster_shape)-9999
    
    bweunit_raster.ravel()[linCoord_sol] = Org_Number_Landing  ## Oder LinIndexLanding

    Org_Number_Landing_ID = raster_edgepoints_withid.ravel()[input_ind[Org_Number_Landing]]
    bweunit_raster.ravel()[linCoord_sol] = Org_Number_Landing_ID  ## Oder LinIndexLanding

    ##-----------------------------------------------------------------------------
    
    return(dist_raster_roads, bweunit_raster)

####===========================================================================
























