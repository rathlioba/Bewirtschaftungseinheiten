# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 11:02:52 2021

@author: Rathliob
"""

# This script is preparing the data to be used in the assessment of the drivability of forest soil
# List of raw input data for criteria:
    # Slope: DHM25
    # Gewässer, Strasse: TLMregio
    # Schutzwald: kantonaler Perimeter
    # Verdichtungsrisiko: Karte Aargau
    # Waldstandorte: Kartierung Staatswald Menzingen
    # Bodeneigenschaften: nationale Bodeneignungskarte
    # Feuchtgebiete: Bundesinventare, TLMregio
    # Schutzgebiete: Bundesinventare, TLMregio, kantonale Daten
    # Forest roads: National and cantonal inventory
    # forest area: cantonal perimeter
    
# Desired Format: Raster with a 2m grid, coded according to the classification of the original values into drivability classes$
    # Slope
    # Obstacles (Roads, Water, Trainlines, dry stone walls)
    # Verdichtungsrisikokarte
    # Karte Waldgesellschaften
    # Karte Bodeneigenschaften (sceleton, permeability, retention, waterlogging)
    # Wetlands (Hochmoor, Flachmoor, Auen)
    # conservation areas (conservation national and cantonal)
    # rockfall protection forest
    # forest roads
    # forest area
    
    # perimeter of observation

# -----------------------------------------------------------------------------
#### Slope

    # Input: swissALTI3D dhm 2m
    # load raster tiles: ZG: 86-95&22-24 AG: 44-49& 39-43
    # GDAL merge the tiles for mhh and zg
            # leave liek this!
            
    # calculate slope in percent and degree
    # clip to forest are in ZG and AG
    
    # Output: slope_degree_mhh/ZG_Stw.tif, slope_percent_mhh/ZG_Stw.tif
        # Output: slope_degree_mhh/ZG_Stw_forest.tif, slope_percent_mhh/ZG_Stw_forest.tif
        # Output: dhm_2m_mhh/ZG_Stw_test.tif
        
    # turn into ascii
    
    
# -----------------------------------------------------------------------------
#### Obstacles

    # Input: tlm: TLM_BB_FELS, TLM_einzelobjekt, TLM_Verbauung, TLM_Versorgungsbaute, 
    # clip to mhh and ZG
        # Versorgungsbaute MHH empty, eisenbahn zg empty, übrige bahn zg and AG empty
    # select fels from tlm bb
    # clip to forest in both areas
    # merge Versorgungsbaute und Verbauung
    # create buffer around all obstacle data, apart from fels with 5m
    # merge fels, verbausorgung and einzelobjekt
    
    # Input: TLM:_Gebäude, 
    # clip to mhh and ZG, then clip by forest in both areas
    # create buffer of 5m around houses
    
    # Input: TLM_Bahn, TLM_übrigebahn
    # clip to mhh and ZG, then clip by forest in both areas
    # merge datasets and create buffer of 5m
    
    # merge all obstacles
    
    # create transport borders by merging the tlm obstacles with bahn and water
    
    # Output: obstacles_tlm_mhh/ZG_Stw.shp, obstacles_buildings_mhh/ZG_Stw.shp
        # obstacles_rail_mhh/ZG_Stw.shp, obstacles_mhh/ZG_Stw.shp
        # transport_borders_mhh/ZG_Stw.shp

# -----------------------------------------------------------------------------
#### soil properties

    # Input: AG only: Verdichtungskarte, ZG only: BEK / forest types Stw
    # clip BEK and forest types ZG to Stw ZG 
    # clip AG Verd to mmh
    # duplicate one of the bek layers
    # classify ZG BEK by attributes: Vernässung = 4 und Skelett = 3, 4 -> klasse 3, hoch
         # Vernässung < 4 und Skelett = 3, 4 -> Klasse 2, mittel
         # delete all other attributes
    # reclassify the forest types ZG layer, add new attribute verdrisk, 
        # reclassify attribute by table based on AG and ZG type description
        # delete all other attributes, but verdrisk

    # duplicate the classified verdrisk shapes twice for both mhh and ZG
    # first set: create new attribute wuchs
        # reclassify based on forest type based on list in ZG classification
        # delete all attributes but wuchs
    # second set: create new attribute eignung    
        # reclassify based on forest type based on list in ZG classification
        # delete all attributes but eignung
        
    
    # Output: verdrisk_mhh, bek_Wasserspeicher/Wasserdurchlass/Vernaessung/Skelett_ZG_Stw.shp
        # verdrisk_bek_ZG_Stw.shp, verdrisk_foresttype_ZG_Stw.shp
        # wuchs_mhh/ZG_Stw.shp, forsteignung_mhh/ZG_Stw.shp
        
        
        
    # Input: verdrisk_mhh/ZG_Stw.shp, forsteignung_mhh/ZG_Stw.shp, wuchs_mhh/ZG_Stw.shp
        # turn into raster with respective burn in field
        # using extent of previously created output
        # processing.run("gdal:rasterize", {'INPUT':'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/testdaten/Input_preprocessing/_mhh/verdrisk_mhh.shp','FIELD':'VERD_RISK','BURN':0,'UNITS':1,'WIDTH':2,'HEIGHT':2,'EXTENT':'2644116.000000000,2649070.000000000,1239858.000000000,1244240.000000000 [EPSG:2056]','NODATA':-9999,'OPTIONS':'','DATA_TYPE':5,'INIT':None,'INVERT':False,'EXTRA':'','OUTPUT':'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/Daten/temp/verdrisk_mhh.tif'})
        # turn into ascii using extent of previously reated output
        # processing.run("gdal:translate", {'INPUT':'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/Daten/temp/verdrisk_mhh.tif','TARGET_CRS':QgsCoordinateReferenceSystem('EPSG:2056'),'NODATA':-9999,'COPY_SUBDATASETS':False,'OPTIONS':'','EXTRA':'','DATA_TYPE':0,'OUTPUT':'N:/forema/FPS/Projekte_der_Gruppe/Masterarbeit_Lioba_Rath/Daten/temp/verdrisk_mhh.asc'})
        
    # Output: verdrisk_mhh/ZG_Stw.asc, forsteignung_mhh/ZG_Stw.asc, wuchs_mhh/ZG_Stw.asc
        
# -----------------------------------------------------------------------------
#### wetlands

   # Input: bundesinventare: Hochmoor, Flachmoor, Moorlandschaft, Auen (alle, alpin, ausserhalb)
   # combine all, then clip by mhh and Stw_perimeter
   # MHH is empty!
   # tlm_Bodenbedeckung clip by ZG Stw, then filter by object type wetlands (11)
   # merge with wetlands zg and dissolve, multi to single
   # clip to ZG Stw and AG
   
   # Output: wetlands_mhh/ZG_Stw.shp


# -----------------------------------------------------------------------------
#### conservation areas

    # Input: waldreservate swiss geodienste, Aargau Vertragsflächen, Waldnaturschutzgebiete
    # select keine gezielten eingriffe from waldreservate and waldnaturschutzgebiete
    # select from langfrsitige Verträge LV TYP: Nutzungsverzicht
    # merge all and then clip to mhh and ZG Stw
    # conservation ZG is empty! 
    
    # Output: conservation_mhh.shp, conservation_ZG_StW.shp

# -----------------------------------------------------------------------------
#### rockfall protection forest

    # Input: Schutzwaldtypen Zug and felssturz ausbruch und transitfläche, NONE for AG
    # aus Schutzwald Zieltyp select Prio 1, and type F, G, H, I, J, K , L
    # merge with felssturz and clip to Stw perimeter
    # dissolve, create spatial index andmultipart to singlepart
    # clip to ZG Stw
    
    # Output: rockfall_protection_ZG_Stw.shp

# -----------------------------------------------------------------------------
#### forest roads

    # Input: TLM_Strasse, Feinerschliessung ZG, Forsterschliessung ZG
    # clip TLM to mhh and ZG
    # select from TLM_Strasse: Objektart: 0, 4, 8, 9, 10 ,11, 15, 16, 20
    # clip AG selection by AG forest
    # tlm ZG selection, select Befahrbarkeit Falsch and delete
    # select wanderweg (and befahrbarkeit) k.W. and delete
    # Select Kunstbaute Treppe and delete
    # select 1m Weg and delete
    # clip Feinerschliessung and FOrsterschliessung by perimeter ZG Stw
    # merge tlm selection with ZG feinerschliessung, forsterschliessung
    # make separate layer with ZG roads clipped to actual Stw perimeter
    
    # Output: roads_mhh.shp, roads_ZG_Stw.shp, roads_ZG_Stw_forest.shp
    

# -----------------------------------------------------------------------------
#### forest area
    
    # Input: for ZG: ZG Staatswaldperimeter, for AG: AG wald
    # clip AG Wald to mhh 
    
    # exclude obstacels_buildings from the forest area
    # Output: forest_nobuild_mhh.shp, forest_minus_ZG_Stw.shp

# -----------------------------------------------------------------------------
#### water

    # Input: TLM_Fliessgewässer, TLM_stehendeGewässer, ZG_gewässernetz
    # clip TLMs to both perimeter MHH and perimeter ZG Stw, clip ZG gewässer to ZGStw
    # combine tlm_fliess+ tlm_stehen for MHH
    # combine all three for ZG
    # delete unterirdisch kuenstlich features from zg dataset
    # dissolve both datasets
    # multipart to singlepart
    # clip AG and ZG water to forest
    
    # Output: water_mhh.shp and water_ZG_Stw.shp


# -----------------------------------------------------------------------------
#### borders

   # Input: swissBOUNDARIES_kantonsgebeit
   # select and extract 'Aargau', then 'Zug'
   # Input: swissBOUNDARIES_hoheitsgebiet
   # select and extract 'Muhen, Hirschthal, Holziken'
   # ZG: Reproject Staatswaldperimeter_Berg to EPSG:2056
   
   # Output: Zug.shp, Aargau.shp, mhh.shp, Staatswaldperimeter_Berg_LV95
   # created: ZG_Stw_perimeter.shp

# -----------------------------------------------------------------------------
#### conifers / broadleaves
    # turn into ascii of specific format
    
# -----------------------------------------------------------------------------
#### Step II
# Database for Criteria

#### trafficability
    # load slope_percent, wetlands, rockfall_zone, verdrisk, conservation
        # between 0% and 30% slope 
        # AND where verdrisk is not 0, 4 or 5. (IS  1, 2 or 3)
        # AND where IS NOT wetlands -----------------------> trafficable 1
        # > 30% slope
        # OR verdrisk 0, 4 or 5 is non-trafficable
        # OR if conservation = False in the conservation areas
        # OR wetlands -------------------------------------> non-trafficable 3
        # between 30% and 100% slope
        # AND within 80m of a forest road downhill and 35m uphill
        # OR within 80m of trafficable downhill and 35m uphill
        # IS NOT rockfall protection forest ---------------> winch_assist 2
        
        
        
#### transport borders
    # load water, obstacles_tlm, obstacles_rail
    # buffer of 1m around water
    # merge all (create spatial index, multipart to singlepart)
    
    # Output: transport_borders_mhh/ZG_Stw.shp
    
#### borders
    # load swissBoundaries, Hoheitsgebiet
    # extract: Muhen, Hirschthal, Holziken, and then Oberägeri, Unterägeri, Menzingen
    # polygons to lines, convert the boundaries to lines
    # buffer of 0.5 m on the outside of the boundaries
    # clip the buffers by the perimeter extent, dissolve
      
    
    # Ouput: gemeinden_mhh.shp, gemeinden_ZG_Stw.shp
    # gemeindegrenzen_buffer_ZG_Stw/mhh.shp
    
#### ecosystem services
    # load Waldnaturschutzgebiete, Besondere Lebensräume, Schutzwald Zieltypen, Erholungswald, forest_ZG_Stw
    # none for AG
    # clip BL, Waldnaturschutz, Schutzwald und Erholungswald by forest ZG Stw
    # clip Stw by BL and recreation, merge all three back
    # polygons to lines, convert the boundaries to lines
    # buffer of 0.5 m on the outside of the boundaries
    # clip the buffers by the perimeter extent, multipart to singelpart
    # manually delete all the buffers that are in the border odf the forest
    
    
    # output: ecosystemservices_buffer_ZG_Stw.shp
    

#### roads in segments
    # Input: roads, as created above
    # use GRASS v.split with max length 5000m (to make it split at intersections)
    
    # Output: roads_mhh_split.shp
        
    
    
    


