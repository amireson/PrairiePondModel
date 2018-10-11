# Prairie Pond Model

This model is applied to partition the loses of water from a prairie wetland pond between evaporation and infiltration, using observations of water level and stable isotope composition of the pond water.

Edward Bam(1,2)<br>
Andrew Ireson(1,3)
 
(1) Global Institute for Water Security, University of Saskatchewan www.usask.ca/GIWS <br>
(2) e.bam@usask.ca <br>
(3) andrew.ireson@usask.ca   http://usask.ca/~andrew.ireson

# The model

The prairie pond model is written in Python 2.7. To run the model, it is necessary to configure the inputs in the file Input.xlsx and to place the script in the same folder as the input. The script can run in three modes: 1) simulate partitioning of evaporation and infiltration with fixed parameters; 2) optimize the proportion of evaporation and infiltration to fit water level and isotopic composition of the pond; 3) run a sensitivity analysis. The script can also be adapted to other run cases.

Also provided are bash scripts and input files with parameters and driving data for four ponds and various years.To reproduce the analysis in *Bam et al. (Journal of Hydrology, under review)* clone this repository and run the bash script "batchrun.sh".
