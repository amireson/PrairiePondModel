#!/bin/bash

Pond=$1
Year=$2

if [[ $Pond = "" ]]; then
    echo Do nothing
    
elif [[ $Pond = "Sensi" ]]; then
    echo Sensitivity analysis
    mkdir TempRun
    cp src/PrairiePondModel.py TempRun/.
    cp InputFiles/SensitivityAnalysis/Input.xlsx TempRun/.

    cd TempRun
    python PrairiePondModel.py
    cd ..

    rm -rf OutputFiles/Sensitivity
    mkdir -p OutputFiles/Sensitivity
    cp TempRun/*.png OutputFiles/Sensitivity/.
    cp TempRun/*.txt OutputFiles/Sensitivity/.
    rm -rf TempRun

else
    echo Running pond $Pond, year $Year
    mkdir TempRun
    cp src/*.py TempRun/.
    cp InputFiles/Pond_$Pond/$Year/Input.xlsx TempRun/.

    cd TempRun
    python PrairiePondModel.py
    python PPM_Plot.py
    cd ..

    rm -rf OutputFiles/Pond_$Pond/$Year
    mkdir -p OutputFiles

    echo $Pond
    cp TempRun/Figure.png OutputFiles/Pond${Pond}_$Year.png
    cp TempRun/Output.xlsx OutputFiles/Pond${Pond}_$Year.xlsx
    rm -rf TempRun
fi



