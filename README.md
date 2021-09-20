# Bewirtschaftungseinheiten

This code is a result of my master thesis "Waldplanung im stufigen Wald - Individuell anpassbare Bewirtschaftungseinheiten basierend auf einem Set an objektiven Kriterien". It is used to group an area of forest into management units (BWEs) based on several criteria. These criteria can be chosen at the beginning of the code to be included or not. 

The grouping of forest into management units is a method originally developed by Leo Bont (WSL), only applied and slightly changed by me. It uses the method of finding the shortest path from the forest to a forest road or the forest edge. 
As this is my first coding, the programming can be optimized and will be further developed. 

This repository contains several files:
- Grouping_Forest_In_BWEs:
containing the actual selection of criteria and process of grouping area into managment units
- datapreprocessing_BWEs:
containing preprocessing stepts to bring all input data for the Grouping_Forest_In_BWEs script into the correct format
- functions_GroupingForestInBWEs:
Authored entirely by Leo Bont! Contains functions used in the Grouping_Forest_In_BWEs script
- description_DatasourcesAndProcessing:
NOT containing code, but the detailed description of all required datasources and their filtering and processing before         using them in the datapreprocessing_BWEs and Grouping_Forest_In_BWEs scripts

Lioba Rath, Forschungsgruppe Nachhaltige Forstwirtschaft, WSL, 20.09.2021

DOI: https://zenodo.org/badge/latestdoi/408371437
