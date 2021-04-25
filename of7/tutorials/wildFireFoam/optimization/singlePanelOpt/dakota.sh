#!/bin/bash

#Adapted from tutorial created by Tobais Holzmann
#http://www.holzmann-cfd.de/index.php/en/

# $1 change from DAKOTA to OF file
# $2 pipe OF output to DAKOTA

# Update OF through DAKOTA
# ------------------------------------------------------------------------------
dprepro $1 system/dakotaParameters.orig system/dakotaParameters

# Run simulation with new parameter set
# ------------------------------------------------------------------------------

    # Get parameters and remove possible scientific notation
    #---------------------------------------------------------------------------
    filmFlowRate=`head -5 system/dakotaParameters | tail -1 | cut -d'=' -f2`
    filmFlowRate=`echo $filmFlowRate | sed -e 's/[eE]+*/\*10\^/'`

    # Rescale parameters to real values
    #---------------------------------------------------------------------------
    filmFlowRateScaled=`echo "$filmFlowRate*0.02 + 0.0" | bc -l` #0 to 0.02 m/s

    # Print current loop number and parameter values
    #---------------------------------------------------------------------------
    funcIter=`cat .funcIter`

    >&2 echo -e ""
    >&2 echo -e "Function Iteration : $funcIter"
    >&2 echo "Film inflow velocity : $filmFlowRateScaled"
    >&2 echo -e ""

    # Set new parameter values
    #---------------------------------------------------------------------------
    sed -i "52s/.*/        value           uniform (0 -$filmFlowRateScaled 0);/" 0/filmRegion/Uf

    # Run OpenFOAM
    #---------------------------------------------------------------------------
    wildFireFoam > log.wildFireFoam 2>&1

    # Get desired quantitites and update cost function
    #---------------------------------------------------------------------------
    solidMassLoss=`grep "Total gas mass produced " log.wildFireFoam | tail -1 | awk '{print $7}'`
    solidMassLoss=`echo $solidMassLoss | sed -e 's/[eE]+*/\*10\^/'`
    func=`echo "0.5*$solidMassLoss/0.00000823 + 0.5*$filmFlowRateScaled/0.01" | bc -l`

    echo $func > .dakotaInput.dak

    >&2 echo -e ""
    >&2 echo "Solid mass loss : $solidMassLoss"
    >&2 echo "Cost Function : $func"
    >&2 echo -e ""

    # Clean Case
    #---------------------------------------------------------------------------
    mkdir runFiles/Optimization$funcIter
    cp -rf 5 constant system log.wildFireFoam runFiles/Optimization$funcIter

    foamListTimes -rm
    rm -rf postProcessing

    # Increase the loop number and store in dummy file
    #--------------------------------------------------------------------------
    echo $((funcIter+1)) > .funcIter

# Pipe output file to DAKOTA
#------------------------------------------------------------------------------
cp .dakotaInput.dak $2

sleep 1

#------------------------------------------------------------------------------
