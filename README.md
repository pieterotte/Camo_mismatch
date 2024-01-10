# Camo_mismatch
Repository for data and code for the manuscript: Otte, Pieter J., Joris P.G.M. Cromsigt, Christian Smit, and Tim R. Hofmeester. “Snow cover-related camouflage mismatch increases detection by predators.” Journal of Experimental Zoology Part A: Ecological and Integrative Physiology, 2024.

**Abstract**

Camouflage expressed by animals is an adaptation to local environments that certain animals express to maximise survival and fitness. Animals at higher latitudes change their coat colour according to a seasonally changing environment, expressing a white coat in winter and a darker coat in summer. The timing of moulting is tightly linked to the appearance and disappearance of snow and is mainly regulated by photoperiod. However, due to climate change an increasing mismatch is observed between the coat colour of these species and their environment. Here, we conducted an experiment in northern Sweden, with white and brown decoys to study how camouflage (mis)-match influenced 1) predator attraction to decoys, and 2) predation events. Using camera trap data, we showed that mismatching decoys attracted more predators and experienced a higher likelihood of predation events in comparison to matching decoys, suggesting that camouflage mismatched animals experience increased detection by predators. These results provide insight into the function of a seasonal colour coat and the need for this adaptation to maximise fitness in an environment that is exposed to high seasonality. Thus, our results suggest that, with increasing climate change and reduced snow cover, animals expressing a seasonal colour coat will experience a decrease in survival. 

**Data**

The analysis is split into two parts, the GLMM and the Time-to-event analyis. 
The GLMM.R file contains the models as used and described in the manuscript. This file requires the dat.pred.csv, dat.pred.type.csv, and dat.interaction.csv files. These files contained processed data that has been retrieved from anotated camera trap data. 
The time-to-event.R file contains the survival models as used and described in the manuscript. This file requires the dat.surv.int.csv, and dat.surv.dec.csv files. These files contain the processed data needed to create survival histories, retrieved from anotated camera trap data.

**Repository**

A secured version of the data and code as used in the published article can be found on Zenodo: 10.5281/zenodo.10478713
