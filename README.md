# ClineDescription

* This document provides information about the datasets and scripts of the manuscript entitled "Selection maintains repeated covariation between color and morphology in a spatially recurrent, temporally stable cline" 

## Outline of contents
````
├──── *README.md
├──── Data
│ ├── Data_IndividualLevel.csv
│ ├── Data_PopulationLevel.csv
│ ├── Data_ReflectanceSpectrometry.csv
│ ├── Data_Simulations_ConsecutiveGenerations.tsv
│ ├── Data_Simulations_ExpectedDefault.tsv
│ ├── Data_Simulations_ReplicatesDefault.tsv
│ ├── Data_Simulations_ReplicatesDrift.tsv
│ ├── Data_Simulations_ReplicatesMigration.tsv
│ └── Data_TemporalSeries.tsv
└──── Scripts
  ├── Scripts_01_EcosystemEquivalence.R
  ├── Scripts_02_ShellColor.R
  ├── Scripts_03_Covariation.R
  └── Scripts_04_SelectivePressures.R
````

## `Data` folder

* This folder contains the raw datasets used in the different scripts. For each of the datasets, a description of the variables and values they can take is provided. 

* Here are hyperlinks to each of them:

1. [Data_IndividualLevel] (#Data_IndividualLevel)
2. [Data_PopulationLevel] (#Data_PopulationLevel)
3. [Data_ReflectanceSpectrometry] (#Data_ReflectanceSpectrometry)
4. [Data_Simulations\_ConsecutiveGenerations](#Data_Simulations\_ConsecutiveGenerations) 
5. [Data_Simulations\_ExpectedDefault] (#Data_Simulations\_ExpectedDefault) 
6. [Data_Simulations\_ReplicatesDefault] (#Data_Simulations\_ReplicatesDefault)
7. [Data_Simulations\_ReplicatesDrift] (#Data_Simulations\_ReplicatesDrift)
8. [Data_Simulations\_ReplicatesMigration] (#Data_Simulations\_ReplicatesMigration)
9. [Data_TemporalSeries] (#Data\_TemporalSeries)

 
## <a id='Data_IndividualLevel'></a> 1.- Data_IndividualLevel

* This dataset contains the dataset for analyses at the individual snail level. 

* The locality-level variables were taken from the `PopulationLevel` dataset, assigning individual snails the locality-level variables of the locality from which they were sampled.  

* The following variables were not used in the analyses of the present paper: `SaxatilisDensity_Mean`, `SaxatilisDensity_SD`, or `RW2-3`. Nevertheless, we keep them if they can be of use. 

### Variables

#### Locality-level variables (Taken from the `PopulationLevel` dataset)

* These variables refer to the locality in which the individuals live, so they are shared between all individual snails inhabiting a particular locality.

**Ria**: Ría to which the sample individual belongs. It can be `Vigo`, `Pontevedra`, `Arousa`, and `MurosNoia`. 

**Coast**: Coast of the ría to which the sample individual belongs. It can be `North` or `South`. 

**RiaCoast**: Ria and coast to which the sample belongs. It constitutes the variables `Ria` and `Coast` fused. 

**LocNumRiaCoast**: Within a given ria and coast, the number of the locality to which the sample individual belongs; it can go from 1 to 6, 7, 8, or 9 (depending on ria and coast), with 1 being the closest to the river mouth (innermost locality) and 6, 7, 8, or 9 the farthest (outermost locality). 

	1=Innermost locality

	...
 
	9=Outermost locality

**LocNumCumulative**: Locality number starting from Vigo, northern coast, innermost to outermost locality, then southern coast. Then Pontevedra, northern coast, innermost to outermost locality, then southern coast. And continuing with the same rule for Arousa, and then Muros-Noia. 

	1=First northern Vigo locality (PuntaCabalo)
	
	...
	
	55=Last southern Muros-Noia locality (Corrubedo)

**LocName**: Name of the locality where the sample individual was collected.

**Latitude**: Latitude of the locality where the sample individual was collected (decimal degrees). Nearest to the sampling points of *L. saxatilis*.

**Longitude**: Longitude of the locality where the sample individual was collected (decimal degrees). Nearest to the sampling points of *L. saxatilis*.

**DistanceRiver_Absolute**: Distance of each locality from the river mouth, bordering the coast (in meters). The measurement  was taken using Google Earth’s path option and following the following criteria: 

	-Eye altitude: 2.85 km
	-Starting point for measurement: the bridge that crosses each river mouth. 
	-Shipyards, ports, and/or any other human constructions on the coast were included in the measurements (because they are part of the coast on which the Littorinids can settle). 
	-Localities at Illa de Arousa (IllaDeArousaE and IllaDeArousaW) were excluded as they are not directly connected to the coast, and no meaningful distance from the river mouth can be assessed.
	-Salt marshes were also included in the measurements, measuring the path from the low tide point. 
	-Salt mine in Vigo (northern coast; 42.351096, -8.631387) was avoided as it has a wall bordering it (measurement was done through that wall).
	-The end of the path was established on the completely exposed coast out of each ria.

**DistanceRiver_Normalized**: Normalized distance of the sampled locality to the outermost point of each ría, as measured in meters (Vigo: North = 54586 m, South = 92171; Pontevedra: North = 53814, South = 36472; Arousa: North = 88131, South = 104121; Muros-Noia: North = 46368, South = 57628). Localities on Illa de Arousa were designated NAs because they are not directly connected to the coast.

**Rock_Type**: Type of rock that makes up the locality, assessed in situ and corroborated with the Mapa Litolóxico de Galicia ([link](http://descargas.xunta.es/844f77ea-5168-4112-8e08-4114a0d0a2bf1561704726980)). 

	AlkalineGranitoids = Alkaline granitoids.
	CalcAlkalineGranitoids = Calc-alkaline granitoids.
	Gneisses = Gneisses
	Schists = Schists

**Wave_Mean**: Mean wave height estimated for each locality. 

**Ascophyllum_Coverage**: Relative coverage of rocks by *Ascophyllum nudosum* in a 2 m (width) x 10 m (long) (20 m2) sample of the lower intertidal zone (usually below Fucus sp.) of the sampled locality, closest to the sampling points of *L. saxatilis*.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering of <20% of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the rocks.
	5 = SuperAbundant = Covering >81% of the rocks.
	
**Barnacles_Coverage**: Relative coverage of the surface of the rocks by barnacles (mostly *Chthamalus sp.* and *Semibalanus balanoides*) in a 2 m (width) x 10 m (length) (20 m^2) sample of the low and mid-upper intertidal zones of the sampled locality, closest to the sampling points of *L. saxatilis*. Barnacles tend to show a uniform distribution, increasing up to the point of covering the whole rock.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering <20% of the surface of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the surface of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the surface of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the surface of the rocks.
	5 = SuperAbundant = Covering >81% of the surface of the rocks.

**Corallina_Presence**: Whether *Corallina officinalis* is found in the lower intertidal zone after 15 minutes of sampling. 

	0 = Absent = Not found.
	1 = Present = Found.
	
**FabalisObtusata_Presence**: Whether *Littorina fabalis* and/or *L. obtusata* is found in the lower intertidal zone of the sampled locality after 15 minutes of sampling. (These species are grouped together as they are very difficult to distinguish in situ.)

	0 = Absent = Not found.
	1 = Present = Found.
	
**Fucus_Coverage**: Relative coverage of rocks by *Fucus sp.* (mostly *F. vesiculosus*) in a 2 m (width) x 10 m (long) (20 m^2) sample of the lower intertidal zone (usually above *Ascophyllum nudosum* and/or below *Pelvetia canaliculata*) of the sampled locality, closest to the sampling points of *L. saxatilis*.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering <20% of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the rocks.
	5 = SuperAbundant = Covering >81% of the rocks.

**Gibbula_Presence**: Whether *Gibbula umbilicalis* is found in the mid intertidal zone of the sampled locality after 15 minutes of sampling. 

	0 = Absent = Not found.
	1 = Present = Found.
	
**Lichina_Coverage**: Relative coverage of rocks by patches of the black lichen *Lichina pygmaea* in a 2 m (width) x 10 m (long) (20 m^2) sample of the upper intertidal zone of the sampled locality, closest to the sampling points of *L. saxatilis*. Patches must be >20x20cm (approx.) to count.

	0 = Absent = No patches.
	1 = Occasional = 1-5 patches.
	2 = Frequent = 6-10 patches.
	3 = Common = 11-16 patches.
	4 = Abundant = 17-21 patches.
	5 = SuperAbundant = Most rocks in the upper zone of the intertidal covered by patches.
	
**Littorea_Presence**: Whether *Littorina littorea* is found in the mid intertidal zone of the sampled locality after 15 minutes of sampling. 

	0 = Absent = Not found.
	1 = Present = Found.
	
**Melarhaphe_Presence**: Whether *Melarhaphe neritoides* is found in the upper intertidal zone of the sampled locality after 15 minutes of sampling.

	0 = Absent = Not found.
	1 = Present = Found.
	
**Mytilus_Coverage**: presence of patches of *Mytilus galloprovincialis* in a 2 m (width) x 10 m (long) (20 m^2) sample of the mid-low intertidal zone of the sampled locality, closest to the sampling points of *L. saxatilis*. Patches should contain groups of at least 10 (most of them mature) individuals.

	0 = Absent = No patches.
	1 = Occasional = 1-5 patches; usually of small mussels.
	2 = Frequent = 6-10 patches.
	3 = Common = 11-16 patches.
	4 = Abundant = 17-21 patches.
	5 = SuperAbundant = Most rocks in the mid-low zone covered by patches.
	
**Nucella_Presence**: Whether *Nucella lapillus* is found in the mid-low intertidal zone of the sampled locality after 15 minutes of sampling.

	0 = Absent = Not found.
	1 = Present = Found.
	
**Pelvetia_Coverage**: Relative coverage of rocks by *Pelvetia canaliculata* in a 2 m (width) x 10 m (long) (20 m^2) sample of the lower intertidal zone of the sampled locality (above *Fucus sp.*), closest to the sampling points of *L. saxatilis*.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering <20% of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the rocks.
	5 = SuperAbundant = Covering >81% of the rocks.
	
**Phorcus_Presence**: Whether *Phorcus lineatus* is found in the mid-upper intertidal zone of the sampled locality after 15 minutes of sampling.

	0 = Absent = Not found.
	1 = Present = Found.

**EcolPC1**: Value of Ecological PC1 for the locality where the snail was caught. Imported from the `PopulationLevel`dataset.

**EcolPC2**: Value of Ecological PC2 for the locality where the snail was caught. Imported from the `PopulationLevel`dataset.

**Cluster**: Cluster to which the locality was assigned when performing hierarchical clustering on a Bray-Curtis distance matrix of species composition (with presence/absence data). The variable was imported from the `PopulationLevel` dataset. It has the following categories, confirmed after inspecting their relationship to PC1:

	WaveSheltered = Cluster grouping all variables that have low EcolPC1 values (wave-sheltered profile).
	WaveIntermediate = Cluster grouping all variables that have intermediate EcolPC1 values (wave-intermediate profile).
	WaveExposed = Cluster grouping all variables that have high EcolPC1 values (wave-exposed profile).
 
**SaxatilisDensity_Mean**: Average counts of *L. saxatilis* in three samples of a quadrat of Ø = 60 cm placed randomly in the sampled locality's mid-upper intertidal zone.

**SaxatilisDensity_SD**: Standard deviation of counts of *L. saxatilis* in three samples of a quadrat of Ø = 60 cm placed randomly in the sampled locality's mid-upper intertidal zone.

#### Individual-level variables

-**IndividualNumCumulative**: Cumulative number of the individual snail, ignoring the locality to which it belongs.

	1 = First snail of the sample.
	
	...
	
	3085 = Last snail of the sample.

-**IndividualNumLoc**: Number of the individual snail within a given locality.

	1 = First snail of the locality
	
	...
	
	n = last snail of the locality

-**IDCodeName**: Identification code for each snail, constructed as the concatenation of the following variables: `Ria`, `Coast`, `LocName`, `IndividualNumLoc`.

-**IDCodeNameColorMorph**: Identification code for each snail, constructed as the concatenation of the following variables: `Ria`, `Coast`, `LocName`, `IndividualNumLoc`, `ColorMorph`. 

**ShellSculpture**: Shell sculpture as described by Fischer-Piette and collaborators (cited in Reid 1996: 318). It ranks the approximate distance between the tip of the ribs and the bottom of the grooves in an ordinal manner. The hypothesis is that this difference is caused by variation in the grooves’ height. This parameter is determined through the inspection of the shell along the last whorl and outer lip.

	Rudis = Smoothed shells, with ribs and grooves very close in height.
	Rudissima = Grooves with intermediate height, so they create a visible sculptured pattern but not too extreme.
	Jugosa = Strongly ribbed shells, with a clear distance between the grooves and the ribs, as seen in the outer lip of the shell.
	
**ColorMorph**: Color morph using a color code system inspired by Sacchi (1979).

	Albida = White. Described by Sacchi (1979:12) as "candida, or cream-white (…)".
	Fulva = Fawn. Described by Sacchi (1979:12) as: "a wide class of colorations, ranging from greyish to hazel, but always with yellow-straw tones".
	Lutea = Yellow. Described by Sacchi (1979:12) as: "yellow (…), from citrine to golden-yellow".
	Aurantia = Orange. Described by Sacchi (1979:12) as: "the orange (…)".
	Fusca = Brown. Described by Sacchi (1979:12) as: "brown-blackish".
	Nigra = Black. Described by Sacchi (1979:12) as: "brilliant-black phenotype with sometimes reddish shades".
	Lineata = Light ribs color with lines in grooves. Described by Sacchi (1979:12) as individuals in which there are "thin dark lines coinciding with the grooves of the sculpture, but with lots of modifications".
	Bande = Two bands. Described by Sacchi (1979:12) as: "shells with two large dark main bands on a lighter background, orange on citrine, red on yellow".
	Tessellata = Tessellations along the shell. Described by Sacchi (1979:12) as individuals “with a reticulated or ‘checkered’ tint, with many intermediate nuances".

**Scars**: Marks on the shell reflecting past attempts of predation by crabs. These are small breaks in the line of growth of the shell followed by new lines of growth. Breaks in the outer lip are only computed as scars if they are large and distinguishable, to avoid computing small breaks caused by the handling of the samples. The snail can present scars in each of the whorls of the shell.

	Absent = No scars.
	Present = One or more scars.

**ShellLength**: Shell length (mm) of the snail from apex to posterior end of body whorl, measured with a digital caliper to the nearest 0.1 mm.

**ShellThickness**: Shell thickness (mm) at the central point of the first whorl, obtained with a gauge to the nearest 0.001 from the aperture of the shell.

**CS**: Centroid size of the individual snail; the square root of the sum of squared distances of the different landmarks from their centroid. 

**RW1**: Relative warp 1 extracted from the snail, explaining the most variance (39.06%).

**RW2**: Relative warp 2 extracted from the snail, explaining 15.47% of the variance.

**RW3**: Relative warp 3 extracted from the snail, explaining 12.05% of the variance.

**SculpturedShell**: Dichotomic variable created during the analysis, where `if_else(ShellSculpture == "Jugosa", "1", "0")`, so:

	0 = Rudis or rudissima
	1 = Jugosa

**Lineata**: Dichotomic variable created during the analysis, where `if_else(ColorMorph == "Lineata", "1", "0")`.

**Fulva**: Dichotomic variable created during the analysis, where `if_else(ColorMorph == "Fulva", "1", "0")`. 

**RW1_Z, RW2_Z, ShellLength_Z**: Variables created during the analysis that consist of Z-scaled versions of the original variables.

**ShellThickness_Rel**: Variable created during the analysis consisting of size-relativized shell thickness values obtained from the residuals of lm(log(ShellThickness) ~ log(ShellLength)). 

## <a id='Data_PopulationLevel'></a> 2.- Data_PopulationLevel

* This dataset contains the population summaries of the geophysical, ecological, and snail data for each of the 55 studied localities. 

* The data concerning snail population means and proportions were extracted from the `IndividualLevel` files, summarizing the variables with R (Script available upon request).

* Several variables of the present dataset were not used or studied only in exploratory analyses. They include: `SaxatilisDensity_1-3`, `Sheltered_Frequency`, or `RW2-3_Mean`. 

### Variables

**Ria**: Ría to which the sample belongs. It can be `Vigo`, `Pontevedra`, `Arousa`, and `MurosNoia`. 

**Coast**: Coast of the ría to which the sample belongs. It can be `North` or `South`. 

**RiaCoast**: Ria and coast to which the sample belongs. It constitutes the variables `Ria` and `Coast` fused. 

**LocNumRiaCoast**: Within a given ria and coast, the number of the locality to which the sample individual belongs; it can go from 1 to 6, 7, 8, or 9 (depending on ria and coast), with 1 being the closest to the river mouth (innermost locality) and 6, 7, 8, or 9 the farthest (outermost locality).  

	1=Innermost locality

	...
 
	9=Outermost locality

**LocNumCumulative**: Locality number starting from Vigo, northern coast, innermost to outermost locality, then southern coast. Then Pontevedra, northern coast, innermost to outermost locality, then southern coast. And continuing with the same rule for Arousa, and then Muros-Noia. 

	1=First northern Vigo locality (PuntaCabalo)
	
	...
	
	55=Last southern Muros-Noia locality (Corrubedo)

**LocName**: Name of the locality.

**Latitude**: Latitude of the locality (decimal degrees). Nearest to the sampling points of *L. saxatilis*.

**Longitude**: Longitude of the locality (decimal degrees). Nearest to the sampling points of *L. saxatilis*.

**DistanceRiver_Absolute**: Distance of each locality from the river mouth, bordering the coast (in meters). The measurement  was taken using Google Earth’s path option and following the following criteria: 

	-Eye altitude: 2.85 km
	-Starting point for measurement: the bridge that crosses each river mouth. 
	-Shipyards, ports, and/or any other human constructions on the coast were included in the measurements (because they are part of the coast on which the Littorinids can settle). 
	-Localities at Illa de Arousa (IllaDeArousaE and IllaDeArousaW) were excluded as they are not directly connected to the coast, and no meaningful distance from the river mouth can be assessed.
	-Salt marshes were also included in the measurements, measuring the path from the low tide point. 
	-Salt mine in Vigo (northern coast; 42.351096, -8.631387) was avoided as it has a wall bordering it (measurement was done through that wall).
	-The end of the path was established on the completely exposed coast out of each ria.

**DistanceRiver_Normalized**: Normalized distance of the sampled locality to the outermost point of each ría, as measured in meters (Vigo: North = 54586 m, South = 92171; Pontevedra: North = 53814, South = 36472; Arousa: North = 88131, South = 104121; Muros-Noia: North = 46368, South = 57628). Localities on Illa de Arousa were designated NAs because they are not directly connected to the coast.

**Rock_Type**: Type of rock that makes up the locality, assessed in situ and corroborated with the Mapa Litolóxico de Galicia ([link](http://descargas.xunta.es/844f77ea-5168-4112-8e08-4114a0d0a2bf1561704726980)). 

	AlkalineGranitoids = Alkaline granitoids.
	CalcAlkalineGranitoids = Calc-alkaline granitoids.
	Gneisses = Gneisses
	Schists = Schists

**Wave_Mean**: Mean wave height estimated for each locality. 

**Ascophyllum_Coverage**: Relative coverage of rocks by *Ascophyllum nudosum* in a 2 m (width) x 10 m (long) (20 m2) sample of the lower intertidal zone (usually below Fucus sp.) of the sampled locality, closest to the sampling points of *L. saxatilis*.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering of <20% of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the rocks.
	5 = SuperAbundant = Covering >81% of the rocks.

**Ascophyllum_Presence**: whether Ascophyllum nudosum is found in the lower intertidal zone (usually below Fucus sp.) after 15 minutes of sampling. (NOTE: Created in Excel with function '=IF(Ascophyllum_Coverage>0, "1", "0").

	0 = Absent = Not found.
	1 = Present = Found.
	
**Barnacles_Coverage**: Relative coverage of the surface of the rocks by barnacles (mostly *Chthamalus sp.* and *Semibalanus balanoides*) in a 2 m (width) x 10 m (length) (20 m^2) sample of the low and mid-upper intertidal zones of the sampled locality, closest to the sampling points of *L. saxatilis*. Barnacles tend to show a uniform distribution, increasing up to the point of covering the whole rock.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering <20% of the surface of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the surface of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the surface of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the surface of the rocks.
	5 = SuperAbundant = Covering >81% of the surface of the rocks.

**Barnacles_Presence**:  Whether barnacles (mostly *Chthamalus sp.* and *Semibalanus balanoides*) are found at a relative abundance in the lower and/or mid-upper intertidal zone after 15 minutes of sampling. (NOTE: Created in Excel with function '=IF(Barnacles_Coverage>2, "1", "0". The reason for this is that Barnacles_Coverage = 0 is extremely rare, and most of the sheltered localities have = 1 or = 2, so to cover variation in this parameter, it is better to increase the threshold for presence.) 
 
	0 = Absent = Not found
	1 = Present = Found
	
**Corallina_Presence**: Whether *Corallina officinalis* is found in the lower intertidal zone after 15 minutes of sampling. 

	0 = Absent = Not found.
	1 = Present = Found.
	
**FabalisObtusata_Presence**: Whether *Littorina fabalis* and/or *L. obtusata* is found in the lower intertidal zone of the sampled locality after 15 minutes of sampling. (These species are grouped together as they are very difficult to distinguish in situ.)

	0 = Absent = Not found.
	1 = Present = Found.
	
**Fucus_Coverage**: Relative coverage of rocks by *Fucus sp.* (mostly *F. vesiculosus*) in a 2 m (width) x 10 m (long) (20 m^2) sample of the lower intertidal zone (usually above *Ascophyllum nudosum* and/or below *Pelvetia canaliculata*) of the sampled locality, closest to the sampling points of *L. saxatilis*.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering <20% of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the rocks.
	5 = SuperAbundant = Covering >81% of the rocks.

**Fucus_Presence**: whether *Fucus sp.* (mostly *F. vesiculosus*) is found in the lower intertidal zone (usually above *Ascophyllum nudosum* and/or below *Pelvetia canaliculata*) of the sampled locality after 15 minutes of sampling. (NOTE: Created in Excel with function '=IF(Fucus_Coverage>0, "1", "0").

	0 = Absent = Not found.
	1 = Present = Found.
	
**Gibbula_Presence**: Whether *Gibbula umbilicalis* is found in the mid intertidal zone of the sampled locality after 15 minutes of sampling. 

	0 = Absent = Not found.
	1 = Present = Found.
	
**Lichina_Coverage**: Relative coverage of rocks by patches of the black lichen *Lichina pygmaea* in a 2 m (width) x 10 m (long) (20 m^2) sample of the upper intertidal zone of the sampled locality, closest to the sampling points of *L. saxatilis*. Patches must be >20x20cm (approx.) to count.

	0 = Absent = No patches.
	1 = Occasional = 1-5 patches.
	2 = Frequent = 6-10 patches.
	3 = Common = 11-16 patches.
	4 = Abundant = 17-21 patches.
	5 = SuperAbundant = Most rocks in the upper zone of the intertidal covered by patches.
	
**Lichina_Presence**: Whether any patch of *Lichina pygmaea* is found in the upper intertidal zone. (NOTE: Created in Excel with function '=IF(Lichina_Coverage>0, "1", "0").
 
	0 = Absent = Not found
	1 = Present = Found
	
**Littorea_Presence**: Whether *Littorina littorea* is found in the mid intertidal zone of the sampled locality after 15 minutes of sampling. 

	0 = Absent = Not found.
	1 = Present = Found.
	
**Melarhaphe_Presence**: Whether *Melarhaphe neritoides* is found in the upper intertidal zone of the sampled locality after 15 minutes of sampling.

	0 = Absent = Not found.
	1 = Present = Found.
	
**Mytilus_Coverage**: presence of patches of *Mytilus galloprovincialis* in a 2 m (width) x 10 m (long) (20 m^2) sample of the mid-low intertidal zone of the sampled locality, closest to the sampling points of *L. saxatilis*. Patches should contain groups of at least 10 (most of them mature) individuals.

	0 = Absent = No patches.
	1 = Occasional = 1-5 patches; usually of small mussels.
	2 = Frequent = 6-10 patches.
	3 = Common = 11-16 patches.
	4 = Abundant = 17-21 patches.
	5 = SuperAbundant = Most rocks in the mid-low zone covered by patches.

**Mytilus_Presence**: whether patches of *Mytilus galloprovincialis* are found in the mid-low intertidal zone of the sampled locality.(NOTE: Created in Excel with function '=IF(Mytilus_Coverage>0, "1", "0").


	0 = Absent = Not found.
	1 = Present = Found.

**Nucella_Presence**: Whether *Nucella lapillus* is found in the mid-low intertidal zone of the sampled locality after 15 minutes of sampling.

	0 = Absent = Not found.
	1 = Present = Found.
	
**Pelvetia_Coverage**: Relative coverage of rocks by *Pelvetia canaliculata* in a 2 m (width) x 10 m (long) (20 m^2) sample of the lower intertidal zone of the sampled locality (above *Fucus sp.*), closest to the sampling points of *L. saxatilis*.

	0 = Absent = No exemplars covering the rocks.
	1 = Occasional = Covering <20% of the rocks.
	2 = Frequent = Covering approximately between 21% and 40% of the rocks.
	3 = Common = Covering approximately between 41% and 60% of the rocks.
	4 = Abundant = Covering approximately between 61% and 80% of the rocks.
	5 = SuperAbundant = Covering >81% of the rocks.

**Pelvetia_Presence**: whether *Pelvetia canaliculata* is found in the lower intertidal zone of the sampled locality (above *Fucus sp.*).(NOTE: Created in Excel with the function '=IF(Pelvetia_Coverage>0, "1", "0").

	0 = Absence = Not found.
	1 = Presence = Found.

**Phorcus_Presence**: Whether *Phorcus lineatus* is found in the mid-upper intertidal zone of the sampled locality after 15 minutes of sampling.

	0 = Absent = Not found.
	1 = Present = Found.
	
**SaxatilisDensity1**: Absolute number of *L. saxatilis* in the first sample of a quadrat of Ø = 60 cm placed randomly in the sampled locality's mid-upper intertidal zone.

**SaxatilisDensity2**: Absolute number of *L. saxatilis* in the second sample of a quadrat of Ø = 60 cm placed randomly in the sampled locality's mid-upper intertidal zone.

**SaxatilisDensity3**: Absolute number of *L. saxatilis* in the third sample of a quadrat of Ø = 60 cm placed randomly in the sampled locality's mid-upper intertidal zone.

**SaxatilisDensity_Mean**: Average counts of *L. saxatilis* in three samples of a quadrat of Ø = 60 cm placed randomly in the sampled locality's mid-upper intertidal zone.

**SaxatilisDensity_SD**: standard deviation of counts of *L. saxatilis* in three samples of a quadrat of Ø = 60 cm placed randomly in the sampled locality's mid-upper intertidal zone.

**n**: Number of *L. saxatilis* snails collected. 

**Albida_Counts**: Absolute number of albida (white) *L. saxatilis* snails in the sample of the population.

**Fulva_Counts**: Absolute number of fulva (cream, greyish) *L. saxatilis* snails in the sample of the population. 

**Lutea_Counts**: Absolute number of lutea (yellow) *L. saxatilis* snails in the sample of the population. 

**Aurantia_Counts**: Absolute number of aurantia (orange) *L. saxatilis* snails in the sample of the population. 

**Fusca_Counts**: Absolute number of fusca (brown) *L. saxatilis* snails in the sample of the population. 

**Nigra_Counts**: Absolute number of nigra (black) *L. saxatilis* snails in the sample of the population. 

**Lineata_Counts**: Absolute number of lineata (lineated with paler background) *L. saxatilis* snails in the sample of the population. 

**Bande_Counts**: Absolute number of bande (orange, brown, or black bands on white or cream background) *L. saxatilis* snails in the sample of the population. 

**Tessellata_Counts**: Absolute number of tessellata (with tessellations) *L. saxatilis* snails in the sample of the population. 

**ColorMorphNA_Counts**: Absolute number of unidentified color morph *L. saxatilis* snails in the sample of the population. 

**Sheltered_Counts**: Number of snails of the morphs that tend to live in sheltered ecosystems (fulva, lutea, and fusca). Obtained in Excel as '=SUM(Fulva\_Counts, Lutea\_Counts, Fusca\_Counts)'.

**Rest_Counts**: Number of snails of all classes except for fulva and lineata). Obtained in Excel as '=SUM(Albida\_Counts, Lutea\_Counts, Aurantia\_Counts, Fusca\_Counts, Nigra\_Counts, Bande\_Counts, Tessellata\Counts, ColorMorph\_NAs)'.

**Minority_Counts**: Number of snails of the minority morphs (all but fulva, fusca, lutea, i.e., the sheltered morphs, and lineata). Obtained in Excel as '=SUM(Albida\_Counts, Aurantia\_Counts, Nigra\_Counts, Bande\_Counts, Tessellata\Counts, ColorMorph\_NAs)'.

**Fulva_Frequency**: Relative frequency of fulva snails within the sample. Obtained in Excel as '=Fulva\_Counts/n\_Snails'. 

**Lutea_Frequency**: Relative frequency of lutea snails within the sample. Obtained in Excel as '=Lutea\_Counts/n\_Snails'. 

**Fusca_Frequency**: Relative frequency of fusca snails within the sample. Obtained in Excel as '=Fusca\_Counts/n\_Snails'. 

**Lineata_Frequency**: Relative frequency of lineata snails within the sample. Obtained in Excel as '=Lineata\_Counts/n\_Snails'. 

**Sheltered_Frequency**: Relative frequency of snails of the morphs that tend to live in sheltered ecosystems (fulva, lutea, and fusca). Obtained in Excel as '=Sheltered\_Counts/n\_Snails'. 

**Rest_Frequency**: Relative frequency of snails of all morphs but fulva and lineata within the sample. Obtained in Excel as '=Rest\_Counts/n\_Snails'. 

**Minority_Frequency**: Relative frequency of snails of the minority morphs (all except for fulva, lutea, fusca, and lineata) within the sample. Obtained in Excel as '=Minority\_Counts/n\_Snails'.

**Rudis\_Counts**: Absolute number of rudis (smoothed shells, with ribs and grooves very close in height) *L. saxatilis* snails in the sample of the population.

**Rudissima\_Counts**: Absolute number of rudissima (grooves with intermediate height, so they create a visible sculptured pattern but not too extreme) *L. saxatilis* snails in the sample of the population.

**Jugosa\_Counts**: Absolute number of jugosa (strongly ribbed shells, with a clear distance between the grooves and the ribs, as seen in the outer lip of the shell) *L. saxatilis* snails in the sample of the population.

**Jugosa\_Frequency**: Relative frequency of jugosa snails within the sample. Obtained in Excel as '=Jugosa\_Counts/n' (NOTE: in the population with ShellSculptureNA\_Counts=1, the formula was '=Jugosa\_Counts/(n-1)' to account for the NA case in the relative frequency). 

**ShellSculptureNA_Counts**: Absolute number of *L. saxatilis* snails of unidentified shell sculpture type in the sample of the population.

**ScarsAbsent_Counts**: Absolute number of *L. saxatilis* snails in the sample of the population without scars or marks on the shell reflecting past attempts of predation by crabs.

**ScarsPresent_Counts**: Absolute number of *L. saxatilis* snails in the sample of the population with scars or marks on the shell reflecting past attempts of predation by crabs.

**ScarsPresent_Frequency**: Relative frequency of snails in the sample with one or more visible predation marks (i.e., scars) in their shells. Obtained in Excel as `=Scars\_Counts/n` (NOTE: in the population with ScarsNA\_Counts=2, the formula was `=ScarsPresent\_Counts/(n-2)` to account for the NA case in the relative frequency).

**ScarsNA_Counts**: Absolute number of *L. saxatilis* snails in the sample of the population in which it was not possible to determine the presence of scars.

**ShellLength_Mean**: Mean value of shell length (mm) of the *L. saxatilis* snails from the sample of the population, from apex to posterior end of body whorl.

**ShellLength_SD**: Standard deviation of shell length (mm) of the *L. saxatilis* snails from the sample of the population, from apex to posterior end of body whorl.

**ShellThickness_Mean**: Mean value of shell thickness (mm) of the *L. saxatilis* snails from the sample of the population, at the central point of the first whorl.

**ShellThickness_SD**: Standard deviation of shell thickness (mm) of the *L. saxatilis* snails from the sample of the population, at the central point of the first whorl.

**ShellThickness_Rel_Mean**: Mean value of relative shell thickness (weighted by length) of the *L. saxatilis* snails from the sample of the population, at the central point of the first whorl. It was obtained by extracting the residuals of a linear model using shell length as a predictor of shell thickness.

**ShellThickness_Rel_SD**: Standard deviation of relative shell thickness (weighted by length) of the *L. saxatilis* snails from the sample of the population, at the central point of the first whorl. It was obtained by extracting the residuals of a linear model using shell length as a predictor of shell thickness.

**CS_Mean**: Mean value of the centroid size of the *L. saxatilis* snails from the sample of the population.

**CS_SD**: Standard deviation of the centroid size of the *L. saxatilis* snails from the sample of the population.

**RW1_Mean**: Mean value of the relative warp 1 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (39.06%).

**RW1_SD**: Standard deviation of the relative warp 1 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (39.06%).

**RW2_Mean**: Mean value of the relative warp 2 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (15.47%).

**RW2_SD**: Standard deviation of the relative warp 2 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (15.47%).

**RW3_Mean**: Mean value of the relative warp 3 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (12.05%).

**RW3_SD**: Standard deviation of the relative warp 3 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (12.05%).

**RW4_Mean**: Mean value of the relative warp 4 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (10.45%).

**RW4_SD**: Standard deviation of the relative warp 4 of the *L. saxatilis* snails from the sample of the population, explaining the most variance (10.45%).

**EcolPC1**: Ecological PC1 calculated from the species variables recorded. It explains 49.2% of the variance.

**EcolPC2**: Ecological PC2 calculated from the species variables recorded. It explains 13.2% of the variance.

**Cluster**: Cluster to which the locality was assigned when performing hierarchical clustering on a Bray-Curtis distance matrix of species composition (with presence/absence data). It has the following categories, confirmed after inspecting their relationship to PC1:

	WaveSheltered = Cluster grouping all variables that have low EcolPC1 values (wave-sheltered profile).
	WaveIntermediate = Cluster grouping all variables that have intermediate EcolPC1 values (wave-intermediate profile).
	WaveExposed = Cluster grouping all variables that have high EcolPC1 values (wave-exposed profile).

## <a id='Data_ReflectanceSpectrometry'></a>  3.- Data_ReflectanceSpectrometry

* This dataset corresponds to spectrometry measurements from *lineata* snails and barnacles from Punta Igreixiña (Vigo N-5), as well as from *fulva* snails and rocks from Punta Cabalo (Vigo N-1). 

* The dataset is formatted for analyses in `pavo`, with each measured object in a column.

### Variables

**wl**: wavelength values as recorded by the spectrometer (`pavo` package format).

**PuntaIgreixinha\_Lineata\_01**: First lineata snail, corresponding to IDCodeName `Vigo_North_PuntaIgreixinha_1` in the IndividualLevel dataset.

...

**PuntaIgreixinha\_Lineata\_46**: Last lineata snail, corresponding to IDCodeName `Vigo_North_PuntaIgreixinha_46` in the IndividualLevel dataset.

**PuntaIgreixinha\_Barnacles\_01**: First random measurement of barnacle coverage in a big rock collected from Punta Igreixinha.

...

**PuntaIgreixinha\_Barnacles\_20**: Last random measurement of barnacle coverage in a big rock collected from Punta Igreixinha.

**PuntaCabalo\_Fulva\_01**: First fulva snail, corresponding to IDCodeName `Vigo_North_PuntaCabalo_1` in the IndividualLevel dataset.

...

**PuntaCabalo\_Fulva\_35**: Last fulva snail, corresponding to IDCodeName `Vigo_North_PuntaCabalo_1` in the IndividualLevel dataset.

**PuntaCabalo\_Rock\_01**: First random measurement of rocks collected from Punta Cabalo.

...

**PuntaCabalo\_Rock\_20**: First random measurement of rocks collected from Punta Cabalo.

## <a id='Data_Simulations_ConsecutiveGenerations'></a>  4.- Data_Simulations\_ConsecutiveGenerations

* This dataset contains the distribution of consecutive generations with significant associations from beta-binomial models across 100 simulation replicates. 

* This code corresponds only to that needed to generate the figures for the manuscript and supplementary materials. All raw data extracted directly from SLiM is available at https://gitlab.com/elcortegano/rias_gallegas.

### Variables

**T_type** Type of concurrent significance condition. 

	T1 = Significance in a single transect (transect 1, i.e., one coast).  
	T2 = Simultaneous significance in transects 1 and 2. 
	T3 = Simultaneous significance in transects 1, 2, and 3. 
	T4 = Simultaneous significance in all four transects. 

**run_length**: Number of consecutive generations with significant associations. Ranges from 1 to 11, with 11 representing `>10` (i.e., more than 10 generations). 

**score**: Median weighted run total across the 100 simulation replicates. For each replicate, it is calculated as the sum of `count × run_length` for runs of consecutive generations with significant associations (`P < 0.05`), where `count` is the number of runs of a given length. 

**sl**: Lower bound (2.5th percentile) of the 95% CI for the score value.

**sh**: Upper bound (97.5th percentile) of the 95% CI for the score value.

## <a id='Data_Simulations_ExpectedDefault'></a>  5.- Data_Simulations\_ExpectedDefault

* This dataset contains the slope coefficients (b) of the beta-binomial regression of the morph frequencies ~ position along the transect, calculated on a per-transect and per-generation basis (i.e., using pooled data from all 100 simulation replicates). 

### Variables

**generation**: Generation number for which the regression was evaluated. It ranges from 1 to 1210 (10 burn-in generations, 600 expansion, 600 end of expansion). 

**transect**: Transect to which the values of the regression correspond. It goes from 1 to 4 (i.e., two full rias). 

**b**: Slope coefficient from the beta-binomial regression.  

**p_beta**: Probability value of the beta-binomial regression slope.

**lrt_p_value**: P-value of the likelihood ratio test comparing the full (with position along the coast as predictor) and null (intercept-only) models. 

## <a id='Data_Simulations_ReplicatesDefault'></a>  6.- Data_Simulations\_ReplicatesDefault

* This dataset contains summary data from the beta-binomial regression of the morph frequencies ~ position along the transect, calculated on a per-transect, per-generation, and per-replicate basis. 

### Variables:

**generation**: Generation number for which the regression was evaluated. It ranges from 1 to 1210 (10 burn-in generations, 600 expansion, 600 end of expansion). 

**transect**: Transect to which the values of the regression correspond. It goes from 1 to 4 (i.e., two full rias). 

**bl**: Lower bound (2.5th percentile) of the 95% CI for the beta-binomial regression slope.

**bh**: Upper bound (97.5th percentile) of the 95% CI for the beta-binomial regression slope.

**b**: Median value of the slope coefficient from the beta-binomial regression, calculated from 100 simulation replicates.

**N**: Number of simulations that were run. 

**nsig**: Number of simulation replicates that returned a statistically significant slope.

**p**: Percentage of simulation replicates that returned a statistically significant slope.

## <a id='Data_Simulations_ReplicatesDrift'></a>  7.- Data_Simulations\_ReplicatesDrift

* This dataset contains the summary results of the SLiM simulations, but with enhanced drift settings. It contains the same variables as `ReplicatesDefault`. 
 
## <a id='Data_Simulations_ReplicatesMigration'></a>  8.- Data_Simulations\_ReplicatesMigration

* This dataset contains the summary results of the SLiM simulations, but with enhanced migration settings. It contains the same variables as `ReplicatesDefault`. 

## <a id='Data_TemporalSeries'></a> 9.- Data_TemporalSeries

* This dataset is a simplified version of the one used for Gefaell et al. (2024), Current Zoology ([link] (https://doi.org/10.6084/m9.figshare.21771152)).

**IMPORTANT NOTE**: The distance from the river mouth data (`DistanceRiver` and `DistanceRiver_Normalized`) were selected from the 2022 sampling points and used for the 1979 sampling points. The distances from the river mouth in the `PopulationLevel` dataset and this one were calculated using slightly different methods (see `DistanceRiver_Absolute` in `IndividualLevel` and `PopulationLevel` above), so there is not a perfect correspondence between the localities both datasets share. 

### Variables

**Year**: Year when the sample was taken. 

	1979 = Sacchi's study (Sacchi, 1979, Tables I-V, pp. 14-16). 
	2022 = Our samples.

**LocName**: Name of locality as provided by Sacchi (1979: 10-11) with some corrections based on contemporary toponimy.

**Coast**: Coast of the ría to which the sample individual belongs. It can be `North` or `South`. 

**LocNumCumulative**: Locality number starting from the northern coast, innermost to outermost locality, then the southern coast.

**DistanceRiver_Absolute**: Distance of each location to Verdugo's rivermouth (start of the estuary) bordering the coast. Measured in meters and estimated using Google Earth. 

**n**: Sample size gathered by Sacchi for his paper and by us in 2022. 

**Fulva_Counts**: Absolute number of fulva *L. saxatilis* snails in the sample of the population.

**Lineata_Counts**: Absolute number of lineata *L. saxatilis* snails in the sample of the population. 

**Fulva_Frequency**: Relative frequency of fulva snails within the sample. Obtained in Excel as '=Fulva\_Counts/n\_Snails'. 

**Lineata_Frequency**: Relative frequency of lineata snails within the sample. Obtained in Excel as '=Fulva\_Counts/n\_Snails'. 


## `Scripts` folder

* This folder contains the scripts used for the data analyses. Scripts normally use more than one dataset. They are structured based on the main questions addressed in the manuscript, as well as in the same order as the main narrative. 

* All scripts have the same basic structure:

		-BASIC INFO: Explains the aim, the researcher in charge, and the date of last update.
		-SETUP SECTION: Loads the packages and datasets. It handles the data and sets themes for figures.
		-ANALYSIS SECTION: Features the code for each of the analyses. It is divided between `Main analyses`(Analyses in the main manuscript) and `Supplementary analyses`(for supplementary material). 
	
* Here are hyperlinks to each of them:

1. [Scripts_EcosystemEquivalence] (#Scripts_01_EcosystemEquivalence)
2. [Scripts_ShellColor] (#Scripts_02_ShellColor)
3. [Scripts_Covariation] (#Scripts_03_Covariation)
4. [Scripts_SelectivePressures](#Scripts_04_SelectivePressures) 

## <a id='Scripts_01_EcosystemEquivalence'></a> 1.- Scripts\_01_EcosystemEquivalence

* This script aims to answer the question of whether the ecosystems of the different Rías Baixas are ecologically equivalent. 

* In the `Main analyses` section, it includes:
	*  **1.-** The ecological PCA (**Section 1**).
	*  **2.-** An exploration of the relationship between the Ecological PC1 and wave exposure (**Section 2**).
	*  **3.-** An exploration of how the Ecological PC1 varies across the geographic space (**Section 3**).
	*  **4.-** A study of ecological differences between ecosystem types (wave-sheltered vs. wave-exposed) across the rias (**Section 4**).
	*  **5.-** The code for creating the composite figure number 2 (**Section 5**). 

## <a id='Scripts_02_ShellColor'></a> 2.- Scripts\_02_ShellColor

* This script aims to answer the question of whether the main colors of the cline are associated with specific ecosystem types, and if so, whether non-selective mechanisms can explain it. 

* In the `Main analyses` section, it includes:

	* **1.-** An exploration of how colors vary with the Ecological PC1 across rias (**Section 1**).
	* **2.-** A study of the temporal stability of the pattern in the Ría de Vigo (with data from Gefaell et al., 2024, Curr Zool) (**Section 2**). 
	* **3.-** A test of isolation by distance as explanation of morph composition differences along the ria (**Section 3**). 
	* **4.-** Figures summarizing the individual-based modeling data (SLiM) (**Section 4**).  
	* **5.-** The code for creating the composite figure number 4 (**Section 5**). 

## <a id='Scripts_03_Covariation'></a> 3.- Scripts\_03_Covariation

* This script aims to answer the question of whether there is covariation between color and other shell traits and whether this covariation is equivalent across the Rías Baixas. 

* In the `Main analyses` section, it includes:

	* **1.-** An analysis of whether the main morphs occupy different positions in the trait multivariate space (**Section 1**).
	* **2.-** Whether the covariation is equivalent across evolutionary units (coasts) (**Section 2**). 
	* **3.-** Whether lineage or morph is more important in explaining patterns of trait covariation (**Section 3**).
	* **4.-** A study of the extent to which phenotypic differences between morphs across rias are of the same magnitude and direction (**Section 4**)
	* **5.-** The code for creating the composite figure number 5 (**Section 5**). 

## <a id='Scripts_04_SelectivePressures'></a> 4.- Scripts\_04_SelectivePressures

* This script aims to answer the question of what the selective pressures are that maintain the traits (color, shell shape, and shell length). 

* In the `Main analyses` section, it includes:

	* **1.-** A test of the hypothesis that colors are maintained by selection for crypsis (**Section 1**). 
	* **2.-** A test of the hypothesis that shell shape and length respond to predation by crabs and wave action, which operate in opposite directions (**Section 2**). 
	* **3.-** The code for creating the composite figure number 6 (**Section 3**).
