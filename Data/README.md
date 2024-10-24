# Calixto et al. Oikos Monarda phytochemistry paper
Repo for Monarda fistulosa phytochemistry and herbivory common garden project. The common garden experiment was located in Wisconsin, USA using seeds collected from M. fistulosa populations originating from Montana, USA and Wisconsin, USA. 
Contact Phil Hahn (hahnp@ufl.edu) with any questions. This manuscript has been accepted for publication (10/21/2024): 

Calixto, E.S., J.L. Maron, K. Keefover-Ring, J.H. Cammarano, and P.G. Hahn. Phytochemical diversity increases with resources availability but has mixed effects on herbivory. Oikos (accepted). 10.1111/oik.10914


Data files include:

## Monarda_Garden_Traits_2019.csv
This file contains data on monarda leaf traits and phytochemistry.

| variables | units | description |
| :--- | :---: | :--- |
| Region | categorical | Geographical region, MT = Montana, WI = Wisconsin |
| Site | categorical | Code for site location, see file 'WIMT_Sites_2019.csv' or Table S1 for more details about locations | 
| Plant | number | numerical code for each maternal plant for where seed was collected |
| Rep | categorical | A-E, code for indivdual plants grown from the maternal plant ID (ie. half siblings) |
| Row | number | row in common garden, used to locate plants in the garden |
| OrderGard | number | order planted in garden, used to locate plants in the garden |
| Unique | number | unique number for each plant |
| Biomass_g | g | dry biomass of plants harvested on 13 August 2019 |
| codedup | categorical | Plant and Rep column pasted together |
| Chemo | categorical | chemotype, T = thymol, C = carvacrol |
| a_thujene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| a_pinene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| sabinene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| b_pinene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| octen_3_ol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| myrcene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| a_phellandrene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| carene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| a_terpinene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| p_cymene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| limonene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| beta_ocimene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| gama_terpinene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| cis_sab_hydrate | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| trans_sab_hydrate | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| borneol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| terp_4_ol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| cym_8_ol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| a_terpineol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| nerol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| isoborneol_formate | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| thymoquinone | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| geraniol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| geranial | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| thymol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| carvacrol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| eugenol | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| caryophyllene | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| germ_D | mg / g dry leaf tissue | concentration of chemical compound, see Table S2 for full compound name additional info |
| total | mg / g dry leaf tissue | concentration of total chemical compounds |
| totTC | mg / g dry leaf tissue | concentration of only thymol and carvacrol |


## MonardaGarden_Data2021.csv 
This file contains data on Monarda plant measurements in the field.

| variables | units | description |
| :--- | :---: | :--- |
| Region      |    categorical     | Geographical region, MT = Montana, WI = Wisconsin  |
| Site        |    categorical     | Code for site location, see file 'WIMT_Sites_2019.csv' or Table S1 for more details about locations   |
| Plant       |    number     | Identifier for plant sample          |
| Rep         |   number      | Replication number                   |
| Row         |   number      | Row position in the field            |
| Order       |    number     | order planted in garden, used to locate plants in the garden  |
| hg_Tallest  |    cm   | Height of the tallest plant part     |
| volD1       |   cm    | Diameter of plant at widest point, used to estimate volume     |
| volD2       |   cm    | Diameter of plant orthogonal to volD1, used to estimate volume |
| b           |  presence  | Presence (1) or absence (0) of buds  |
| fl          |  presence  | Presence (1) or absence (0) of flowers                    |
| fr          |  presence  | Presence (1) or absence (0) of fruits                     |
| infl        |  count  | Number of inflorescences             |
| ch          |    %    | Percent of chewing damage estimated on whole plant            |
| prc         |    %    | Percent of piercing/sucking damage estimated on whole plant |
| notes       |         | Additional notes or observations     |


## MonardaGarden_SeedDamage2021 (5).csv 
This file contains data on seed head damage to Monarda plants in the common garden.

| variables | units | description |
| :--- | :---: | :--- |
| Region     |    categorical     | Geographical region, MT = Montana, WI = Wisconsin   |
| Site       |    categorical     | Code for site location, see file 'WIMT_Sites_2019.csv' or Table S1 for more details about locations   |
| Plant      |    number     | Identifier for plant sample           |
| Rep        |    number     | Replication number                    |
| Row        |     number    | Row position in the field             |
| Flowered   |  Presence | Whether the plant flowered            |
| HeadSize   |   cm    | Size of the flower head               |
| HeadMass   |    g    | Mass of the flower head               |
| N_tubes    |  count  | Number of tubes                       |
| HerbDmg    |    %    | Percentage of herbivory damage        |
| Col_grub   |  count  | Number of Coleoptera grubs            |
| Lep_grub   |  count  | Number of Lepidoptera grubs           |
| Lep_Phy    |  count  | Number of Phytophagous Lepidoptera    |
| Fly        |  count  | Number of flies                       |
| Unknwn     |  count  | Number of unknown organisms           |
| Seeds_No   |  count  | Number of seeds                       |
| Seed_wt    |    g    | Weight of the seeds                   |
| Notes      |         | Additional notes or observations      |
| DateP      |  Date   | Date of the recorded measurement      |



## MonardaGarden2020.csv
This file contains data collected in the common garden site in Wisconsin in the summer of 2020.

| variables | units | description |
| :--- | :---: | :--- |
| Region | categorical | MT = Montana, WI = Wisconsin |
| Site | categorical | Code for site location, see file 'WIMT_Sites_2019.csv' or Table S1 for more details about locations  | 
| Plant | number | numerical code for each maternal plant for where seed was collected |
| Rep | categorical | A-E, code for indivdual plants grown from the maternal plant ID (ie. half siblings) |
| Row | number | row in common garden, used to locate plants in the garden |
| Order | number | order planted in garden, used to locate plants in the garden |
| survey | category	| month of survey, July or Aug=August |
| height | cm	| stem height, ground to top of plant |
| leaves_total | count | number of leaves on plant |	
| leaves_damaged | count | number of leaves with damage |	
| herbpercentage | percent | overall estimate of percent damage to the whole plant |	
| leaf_1 | percent | percent damage visually estimated on a single leaf |	
| leaf_2 | percent | percent damage visually estimated on a single leaf |	
| leaf_3 | percent | percent damage visually estimated on a single leaf |	
| leaf_4 | percent | percent damage visually estimated on a single leaf |	
| leaf_5 | percent | percent damage visually estimated on a single leaf |	
| leaf_6 | percent | percent damage visually estimated on a single leaf |	
| leaf_7 | percent | percent damage visually estimated on a single leaf |	
| leaf_8 | percent | percent damage visually estimated on a single leaf |	
| leaf_9 | percent | percent damage visually estimated on a single leaf |	
| leaf_10 | percent | percent damage visually estimated on a single leaf |
| mildew | presence absence | presence of powdery mildew, 0=absent, 1=present |
| mines_blotch | count | number of powdery mine blotches |	
| mines_linear | count | number of leaf mines |	
| galls_leaf | count | number of leaf galls by Ametrodiplosis fistulosae |	
| galls_stem | count | number of stem galls |	
| galls_bud | count | number of bud galls |	
| leps | count | number of lepidoptera caterpillars counted on plant |	
| aphids_parasitized | count |	number of parasitized aphid mummies on plant |
| aphids_unparasitized | count | number of live aphids (Aphis monardae) on plant |	
| ants | count | number of ants on plant |	
| ant_genus | category | genus of ant |	
| flowerbuds | count | number of flower buds |	
| bloom_date | date | date of bloom |	
| damage_chewing | percent | percent chewing damage to whole plant |	
| damage_spot | presence | presence of spotting damage |	
| agromyzid | count | number of agromyzid flies on plant |	


## WIMT_Sites_2019.csv 
This file contains information on site characteristics for locations where seed was originally collected.

| variables | units | description |
| :--- | :---: | :--- |
| Region     |     categorical    | Geographical region, MT = Montana, WI = Wisconsin                  |
| Site       |     categorical    | Code for site location  |
| Lat        |  degrees | Latitude coordinate                  |
| Long       |  degrees | Longitude coordinate                 |
| pH         |         | Soil pH level                        |
| orgmat     |    %    | Organic matter content               |
| nit        |    %    | Nitrogen content                     |
| K          |    %    | Potassium content                    |
| cec        |  meq/100g | Cation exchange capacity            |
| P          |    %    | Phosphorus content                   |
| mat1       |     cm    | bioclim1, mean annual precipitation    |
| iso3       |         | bioclim3, isothermatic               |
| tsd4       |         | bioclim4, standard deviation of monthly temperatures  |
| map12      |  mm/year | bioclim12, Mean annual precipitation |
| pwrq18     |         | bioclim18, precipication in the warmest quarter      |

