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
    injectionAngle1=`head -5 system/dakotaParameters | tail -1 | cut -d'=' -f2`
    injectionMass=`head -8 system/dakotaParameters | tail -1 | cut -d'=' -f2`

    injectionAngle1=`echo $injectionAngle1 | sed -e 's/[eE]+*/\*10\^/'`
    injectionMass=`echo $injectionMass | sed -e 's/[eE]+*/\*10\^/'`

    # Rescale parameters to real values
    #---------------------------------------------------------------------------
    injectionAngle1Scaled=`echo "tmp=$injectionAngle1*90; scale=2; tmp/1" | bc -l`
    injectionMassScaled=`echo "tmp=$injectionMass*2 + 0.0; scale=2; tmp/1" | bc -l`

    # Print current loop number and parameter values
    #---------------------------------------------------------------------------
    funcIter=`cat .funcIter`

    >&2 echo -e ""
    >&2 echo -e "Function Iteration : $funcIter"
    >&2 echo "Injection cone angle : $injectionAngle1Scaled"
    >&2 echo "Injection mass : $injectionMassScaled"
    >&2 echo -e ""

    # Set new parameter values
    #---------------------------------------------------------------------------
    sed -i "107s/.*/           thetaOuter       $injectionAngle1Scaled;/" constant/reactingCloud1Properties
    sed -i "88s/.*/            massTotal        $injectionMassScaled;/" constant/reactingCloud1Properties

    # Run OpenFOAM
    #---------------------------------------------------------------------------
    # Source tutorial run functions
    . $WM_PROJECT_DIR/bin/tools/RunFunctions

    application=wildFireFoam
    nProcs=$(getNumberOfProcessors)
    
    if [ "$nProcs" -gt "1" ];
    then
        echo "Running $application in parallel on $nProcs processors"
        ( mpirun -np $nProcs $application -parallel > log.$application 2>&1 )
    else
        echo "Running $application in serial"
        ( $application > log.$application 2>&1 )
    fi

    # Get desired quantitites and update cost function
    #---------------------------------------------------------------------------
    solidMassLoss=`grep "Total gas mass produced " log.wildFireFoam | tail -1 | awk '{print $7}'`
    solidMassLoss=`echo $solidMassLoss | sed -e 's/[eE]+*/\*10\^/'`
    func=`echo "0.5*$solidMassLoss/0.00108471 + 0.5*$injectionMassScaled/1" | bc -l`

    echo $func > .dakotaInput.dak

    >&2 echo -e ""
    #>&2 echo "Mean Integrated HRR : $meanIntHRR"
    >&2 echo "Solid mass loss : $solidMassLoss"
    >&2 echo "Cost Function : $func"
    >&2 echo -e ""

    # Clean Case
    #---------------------------------------------------------------------------
    mkdir runFiles/Optimization$funcIter
    cp -rf 5 constant system log.wildFireFoam runFiles/Optimization$funcIter

    foamListTimes -rm

    # Increase the loop number and store in dummy file
    #--------------------------------------------------------------------------
    echo $((funcIter+1)) > .funcIter

# Pipe output file to DAKOTA
#------------------------------------------------------------------------------
cp .dakotaInput.dak $2

sleep 1

#------------------------------------------------------------------------------
