#!/usr/bin/python3

import numpy as np
import ddpy
import os
import uuid
import pickle
import pandas as pd
import matplotlib.pyplot as plt

import statsmodels.api as sm

from multiprocessing import Process
import signal
#import subprocess

def main():
    simulationCount = 1
    #for ii in range(2):
    #    subprocess.run( [singleModelib, ''] )
    for ii in range( simulationCount):
        print(f"loopIdx: {ii}")
        p = Process( target=singleModelib)
        p.start()
        p.join()

    #signal.signal( signal.SIGABRT, signal_handler)
    #os.kill( os.getpid(), signal.SIGABRT)

    return

def singleModelib():
    #hiRateSteps = 300
    #lowRateSteps = 1000
    hiRateSteps = 300
    zeroRateSteps = 300
    lowRateSteps = 100000
    #lowRateSteps = 3000
    stepsBetweenMeasurements = 1000 # steps
    measurementPeriod = 1
    evlOutputPeriod1 = 100
    #evlOutputPeriod1 = 300
    evlOutputPeriod2 = 100
    #smoothingWindowTimePeriod = 1.0
    smoothingWindowTimePeriod = 5.0e-9
    #smoothingWindowTimePeriod = 1.0e-14
    #smoothingWindowTimePeriod = 700.0
    dddFolderPath = "./testSerialRuns-allSlipSystems"#contains inputFiles/{DD.txt,mesh,load,material}
    lattice = 'bcc' # TODO: read crystalStructure from {material}.txt
    #lattice = 'fcc'
    # TODO: create uniformExternalLoadController.txt via python
    # TODO: create DD.txt via python
    # TODO: create mesh, material, and load files via python

    meshFilePath = "./unitCube.msh" # relative to inputFiles directory
    material = 'Fe_320' # used to find file {material}.txt
    #material = 'Cu'

    # simulation box dimensions 10^{-10}[m]
    xlo = 0
    #xhi = 1000
    xhi = 15000 # angstrom
    ylo = 0
    #yhi = 1000
    yhi = 15000
    zlo = 0
    #zhi = 1000
    zhi = 15000

    #alignToSlipSystem0 = 1 # 1 : true, 0 : false ## implemented for 'true'

    # user must identify slip system 0 slip direction ss and plane normal nn
    if lattice == 'bcc':
        ss=np.array([5.773502691896258e-01,  5.773502691896258e-01, -5.773502691896258e-01])
        nn=np.array([8.660254037844385e-01, 0.000000000000000e+00, 8.660254037844385e-01])
        nn=nn/np.linalg.norm(nn)
    elif lattice == 'fcc':
        ss=np.array([ 0.000000000000000e+00, 7.071067811865475e-01, 7.071067811865475e-01])
        nn=np.array([-7.071067811865476e-01, 7.071067811865476e-01, -7.071067811865476e-01])
        nn=nn/np.linalg.norm(nn)

    # create C2G matrix
    print(f"createC2G( lattice, ss, nn)")
    c2g = createC2G( lattice, ss, nn) #alignToSlipSystem0)

    print(f"creating DDInterface({dddFolderPath})")
    dc = ddpy.DDInterface( dddFolderPath)
    print(f"reading Burger's vector magnitude from {dddFolderPath}/inputFiles/{material}.txt")
    dc.readBurgersMagnitude( dddFolderPath
                                + "/inputFiles/"
                                + material
                                + ".txt") # required for readLammpsBounds

    print(f"setting simulation box bounds: [{xlo},{xhi}],[{ylo},{yhi}],[{zlo},{zhi}]") # prerequisite: burgersMagnitude
    dc.setBoxBounds( xlo, xhi, ylo, yhi, zlo, zhi) # [m]
    print(f"regenerating polycrystal.txt") # prerequisite: burgersMagnitude
    dc.regeneratePolycrystalFile(
            c2g=c2g,
            lattice=lattice,
            material=material,
            meshFilePath=meshFilePath)

    parameterLog = dict()
    ###

    #for loopIdx in range(simulationCount):
    #print(f"loop {loopIdx}")
    tag = str(uuid.uuid1())

    outputFolder = dddFolderPath + "/output/" + tag # to be created within dddFolderPath

    # create output directories
    if not os.path.exists( outputFolder):
        print(f"creating directory: {outputFolder}")
        try:
            os.makedirs( outputFolder)
        except:
            print(f"failed to create: {outputFolder}")
    if not os.path.exists( outputFolder + "/evl"):
        print(f"creating directory: {outputFolder}/evl")
        try:
            os.makedirs( outputFolder + "/evl")
        except:
            print(f"failed to create: {outputFolder}/evl")
    if not os.path.exists( outputFolder + "/F"):
        print(f"creating directory: {outputFolder}/F")
        try:
            os.makedirs( outputFolder + "/F")
        except:
            print(f"failed to create: {outputFolder}/F")

    #periodicLoopTargetDensity = 1e15
    #periodicLoopSegmentCount = 50
    #periodicLoopRadiusDistributionMean = 4.0e-8
    #periodicLoopRadiusDistributionStd = 2.0e-8
    #periodicLoopTargetDensity = 1e13
    periodicLoopSegmentCount = 50
    periodicLoopRadiusDistributionMean = 1.0e-6
    periodicLoopRadiusDistributionStd = 1.0e-8
    #periodicLoopRadiusDistributionMean = 1.0e-7
    #periodicLoopRadiusDistributionStd = 1.0e-8
    loopDensitiesPerSlipSystem \
            = {
                    0:1.0e0,
                    1:1.0e0,
                    2:1.0e0,
                    3:1.0e0,
                    4:1.0e0,
                    5:1.0e0,
                    6:1.0e0,
                    7:1.0e0,
                    8:1.0e0,
                    9:1.0e0,
                    10:1.0e0,
                    11:1.0e0,
                    12:1.0e0,
                    13:1.0e0,
                    14:1.0e0,
                    15:1.0e0,
                    16:1.0e0,
                    17:1.0e0,
                    18:1.0e0,
                    19:1.0e0,
                    20:1.0e0,
                    21:1.0e0,
                    22:1.0e0,
                    23:1.0e0
            }

    parameterLog[ tag] = {
            'loopDensitiesPerSlipSystem':loopDensitiesPerSlipSystem,
            #'periodicLoopTargetDensity':periodicLoopTargetDensity,
            'periodicLoopSegmentCount':periodicLoopSegmentCount,
            'periodicLoopRadiusDistributionMean':periodicLoopRadiusDistributionMean,
            'periodicLoopRadiusDistributionStd':periodicLoopRadiusDistributionStd
            }

    #print(f"loopIdx:{loopIdx}, parameterLog:\n{parameterLog}")
    print("dc.clearMicrostructureSpecifications()")
    dc.clearMicrostructureSpecifications()

    #print("dc.specifyPrismaticLoopDensityPerSlipSystem()")
    #dc.specifyPrismaticLoopDensitiesPerSlipSystem(
    #        tag=tag,
    #        prismaticLoopDensitiesPerSlipSystem=loopDensitiesPerSlipSystem,
    #        #loopSegmentCount=periodicLoopSegmentCount,
    #        prismaticLoopRadiusMean=periodicLoopRadiusDistributionMean,
    #        prismaticLoopRadiusStd=periodicLoopRadiusDistributionStd
    #        )

    print("dc.specifyLoopDensityPerSlipSystem()")
    dc.specifyLoopDensitiesPerSlipSystem(
            tag=tag,
            loopDensitiesPerSlipSystem=loopDensitiesPerSlipSystem,
            loopSegmentCount=periodicLoopSegmentCount,
            loopRadiusMean=periodicLoopRadiusDistributionMean,
            loopRadiusStd=periodicLoopRadiusDistributionStd
            )


    #print("dc.specifyLoopDensity()")
    #dc.specifyLoopDensity(
    #        tag=tag,
    #        loopDensity=periodicLoopTargetDensity,
    #        loopSegmentCount=periodicLoopSegmentCount,
    #        loopRadiusMean=periodicLoopRadiusDistributionMean,
    #        loopRadiusStd=periodicLoopRadiusDistributionStd
    #        )
    print(f"trying to run in {outputFolder}")
    try:
        #runInNewDir( dc, nSteps, outputFolder)
        dc.setOutputPath( outputFolder)
        print(f"dc.generateMicrostructure()")
        dc.generateMicrostructure() # will write configuration to file
        dc.readDefectiveCrystal()

        print(f"dc.setOutputFrequency({evlOutputPeriod1})")
        dc.setOutputFrequency( evlOutputPeriod1)

        # attempt to change the external load applied
        #strain \
        #    = np.array([
        #        [0,0,0],
        #        [0,0.0,0],
        #        [0,0,0]],dtype=np.float64)
        #strainRate \
        #    = np.array([
        #        [0,0,0],
        #        #[0,1.0e-5,0],
        #        [0,1.0e-12,0],
        #        [0,0,0]],dtype=np.float64)
        #stress \
        #    = np.array([
        #        [0,0,0],
        #        [0,0,0],
        #        [0,0,0]],dtype=np.float64)
        #stressRate \
        #    = np.array([
        #        [0.0,0.0,0.0],
        #        [0.0,0.0,0.0],
        #        [0.0,0.0,0.0]],dtype=np.float64)
        #stiffnessRatio \
        #    = np.array([0,1.0e20,0,0,0,0],dtype=np.float64)

        #dc.setExternalLoad(
        #        strain=strain,
        #        strainRate=strainRate,
        #        stress=stress,
        #        stressRate=stressRate,
        #        stiffnessRatio=stiffnessRatio
        #        )

        #parameterLog[ tag]['strainRate1'] = strainRate
        #try:
        #    print(f"trying to write parameterLog[{tag}]")
        #    with open(outputFolder + "/parameterLog.pkl", 'wb') as handle:
        #        pickle.dump(
        #                parameterLog,
        #                handle,
        #                protocol=pickle.HIGHEST_PROTOCOL
        #                )
        #except:
        #    print(f"error: failed to write parameterLog.pkl")


        #dc.disableMechanicalMeasurements()
        #dc.runGlideSteps( hiRateSteps)

        #print(f"dc.writeConfigToTxt()")
        #dc.writeConfigToTxt()

        #####
        #strainRate \
        #    = np.array([
        #        [0,0,0],
        #        [0,0,0],
        #        [0,0,0]],dtype=np.float64)
        #dc.setExternalLoad(
        #        strainRate=strainRate,
        #        )
        #print(f"dc.enableMechanicalMeasurements({measurementPeriod})")
        #dc.enableMechanicalMeasurements( measurementPeriod)

        #print(f"dc.setOutputFrequency({evlOutputPeriod2})")
        #dc.setOutputFrequency( evlOutputPeriod2)

        #print(f"dc.runGlideSteps({zeroRateSteps})")
        #dc.runGlideSteps( zeroRateSteps)

        #print(f"dc.writeConfigToTxt()")
        #dc.writeConfigToTxt()
        #####
        #print(f"dc.writeConfigToTxt()")
        #dc.writeConfigToTxt()

        strainRate \
            = np.array([
                [0,0,0],
                [0,1.0e+3,0], # [%/s]
                [0,0,0]],dtype=np.float64)

        parameterLog[ tag]['strainRate'] = strainRate
        try:
            print(f"trying to write parameterLog[{tag}]")
            with open(outputFolder + "/parameterLog.pkl", 'wb') as handle:
                pickle.dump(
                        parameterLog,
                        handle,
                        protocol=pickle.HIGHEST_PROTOCOL
                        )
        except:
            print(f"error: failed to write parameterLog.pkl")

        dc.setExternalLoad(
                strainRate=strainRate,
                )
        print(f"dc.enableMechanicalMeasurements({measurementPeriod})")
        dc.enableMechanicalMeasurements( measurementPeriod)

        print(f"dc.setOutputFrequency({evlOutputPeriod2})")
        dc.setOutputFrequency( evlOutputPeriod2)

        totalStepsRun = 0
        while totalStepsRun < lowRateSteps:

            print(f"dc.runGlideSteps({stepsBetweenMeasurements})")
            dc.runGlideSteps( stepsBetweenMeasurements )
            totalStepsRun += stepsBetweenMeasurements 
            print(f"totalStepsRun: {totalStepsRun}")

            print(f"dc.writeConfigToTxt()")
            dc.writeConfigToTxt()

            print(f"measurementTimes, measurementIDs, measuredProperties = dc.getMechanicalMeasurements()")
            measurementTimes, measurementIDs, measuredProperties \
                    = dc.getMechanicalMeasurements()
            measurementTimes = measurementTimes.squeeze()
            measurementIDs = measurementIDs.squeeze()

            try:
                print(f"trying to write measurements.pkl")
                with open(outputFolder + "/measurements.pkl", 'wb') as handle:
                    pickle.dump(
                            measuredProperties,
                            handle,
                            protocol=pickle.HIGHEST_PROTOCOL
                            )
                with open(outputFolder + "/times.pkl", 'wb') as handle:
                    pickle.dump(
                            measurementTimes,
                            handle,
                            protocol=pickle.HIGHEST_PROTOCOL
                            )
                with open(outputFolder + "/runIDs.pkl", 'wb') as handle:
                    pickle.dump(
                            measurementIDs,
                            handle,
                            protocol=pickle.HIGHEST_PROTOCOL
                            )
            except:
                print(f"error: failed to write measurements.pkl and times.pkl")

            # calculate total density time series
            totalDensity = np.zeros( ( len( measuredProperties['density'][0])))
            for ssID in measuredProperties['density'].keys():
                for timeStep in range( len(measuredProperties[ 'density'][ ssID])):
                    totalDensity[ timeStep] += measuredProperties[ 'density'][ ssID][ timeStep]

            #print(f"totalDensity.shape: {totalDensity.shape}")
            #print(f"measurementTimes.shape: {measurementTimes.shape}")
            #print(f"totalDensity: {totalDensity}")

            # calculate strainRates
            strainRates = dict() # dict of np.array, one per slip system
            smoothedPlasticDistortion = dict()
            smoothedStrainRates = dict()
            #for ii in measuredProperties['slipSystemPlasticDistortion'].keys(): # keys are ssID
            #    strainRates[ii] = np.ndarray(measurementTimes.shape[0]-1)
            #    if (measurementTimes.shape[0] != measuredProperties['slipSystemPlasticDistortion'][ii].shape[0]):
            #        print(f"error: measurementTimes.shape()[0] != measuredProperties['slipSystemPlasticDistortion'][ii].shape[0]")
            #        return
            #    # iterate over time step and calculate rates
            #    for tt in range( measurementTimes.shape[0] -1):
            #        delta_t = (measurementTimes[tt+1] - measurementTimes[tt])
            #        strainRates[ii][tt] = (measuredProperties['slipSystemPlasticDistortion'][ii][tt+1] - measuredProperties['slipSystemPlasticDistortion'][ii][tt])/delta_t
            #for ssID in measuredProperties['slipSystemPlasticDistortion'].keys():
            #    if (measurementTimes.shape[0] != measuredProperties['slipSystemPlasticDistortion'][ ssID].shape[0]):

            for ssID in measuredProperties['slipSystemPlasticDistortion'].keys():
                strainRates[ ssID] = calculate_rate( measuredProperties['slipSystemPlasticDistortion'][ ssID], measurementTimes)

                smoothedPlasticDistortion[ ssID] = smooth_measurement( measuredProperties['slipSystemPlasticDistortion'][ ssID], measurementTimes, smoothingWindowTimePeriod)
                smoothedStrainRates[ ssID] = calculate_rate( smoothedPlasticDistortion[ ssID][:,1], smoothedPlasticDistortion[ ssID][:,0])

            create_measurement_plots(
                    measuredProperties,
                    measurementTimes,
                    measurementIDs,
                    strainRates,
                    smoothedStrainRates,
                    totalDensity,
                    outputFolder)

        ################################################################
        ### end loop
        #### resume run

        #print(f"dc.runGlideSteps({lowRateSteps})")
        #dc.runGlideSteps( lowRateSteps)

        #print(f"dc.writeConfigToTxt()")
        #dc.writeConfigToTxt()

        #print(f"measurementTimes, measurementIDs, measuredProperties = dc.getMechanicalMeasurements()")
        #measurementTimes, measurementIDs, measuredProperties \
        #        = dc.getMechanicalMeasurements()
        #measurementTimes = measurementTimes.squeeze()
        #measurementIDs = measurementIDs.squeeze()

        #try:
        #    print(f"trying to write measurements2.pkl")
        #    with open(outputFolder + "/measurements.pkl", 'wb') as handle:
        #        pickle.dump(
        #                measuredProperties,
        #                handle,
        #                protocol=pickle.HIGHEST_PROTOCOL
        #                )
        #    with open(outputFolder + "/times.pkl", 'wb') as handle:
        #        pickle.dump(
        #                measurementTimes,
        #                handle,
        #                protocol=pickle.HIGHEST_PROTOCOL
        #                )
        #    with open(outputFolder + "/runIDs.pkl", 'wb') as handle:
        #        pickle.dump(
        #                measurementIDs,
        #                handle,
        #                protocol=pickle.HIGHEST_PROTOCOL
        #                )
        #except:
        #    print(f"error: failed to write measurements2.pkl and times2.pkl")

        ## calculate total density time series
        #totalDensity = np.zeros( ( len( measuredProperties['density'][0])))
        #for ssID in measuredProperties['density'].keys():
        #    for timeStep in range( len(measuredProperties[ 'density'][ ssID])):
        #        totalDensity[ timeStep] += measuredProperties[ 'density'][ ssID][ timeStep]

        ##print(f"totalDensity.shape: {totalDensity.shape}")
        ##print(f"measurementTimes.shape: {measurementTimes.shape}")
        ##print(f"totalDensity: {totalDensity}")

        ## calculate strainRates
        #strainRates = dict() # dict of np.array, one per slip system
        #smoothedPlasticDistortion = dict()
        #smoothedStrainRates = dict()
        ##for ii in measuredProperties['slipSystemPlasticDistortion'].keys(): # keys are ssID
        ##    strainRates[ii] = np.ndarray(measurementTimes.shape[0]-1)
        ##    if (measurementTimes.shape[0] != measuredProperties['slipSystemPlasticDistortion'][ii].shape[0]):
        ##        print(f"error: measurementTimes.shape()[0] != measuredProperties['slipSystemPlasticDistortion'][ii].shape[0]")
        ##        return
        ##    # iterate over time step and calculate rates
        ##    for tt in range( measurementTimes.shape[0] -1):
        ##        delta_t = (measurementTimes[tt+1] - measurementTimes[tt])
        ##        strainRates[ii][tt] = (measuredProperties['slipSystemPlasticDistortion'][ii][tt+1] - measuredProperties['slipSystemPlasticDistortion'][ii][tt])/delta_t
        ##for ssID in measuredProperties['slipSystemPlasticDistortion'].keys():
        ##    if (measurementTimes.shape[0] != measuredProperties['slipSystemPlasticDistortion'][ ssID].shape[0]):

        #for ssID in measuredProperties['slipSystemPlasticDistortion'].keys():
        #    strainRates[ ssID] = calculate_rate( measuredProperties['slipSystemPlasticDistortion'][ ssID], measurementTimes)

        #    smoothedPlasticDistortion[ ssID] = smooth_measurement( measuredProperties['slipSystemPlasticDistortion'][ ssID], measurementTimes, smoothingWindowTimePeriod)
        #    smoothedStrainRates[ ssID] = calculate_rate( smoothedPlasticDistortion[ ssID][:,1], smoothedPlasticDistortion[ ssID][:,0])

        #create_measurement_plots(
        #        measuredProperties,
        #        measurementTimes,
        #        measurementIDs,
        #        strainRates,
        #        smoothedStrainRates,
        #        totalDensity,
        #        outputFolder)

        ##### resume a 2nd time
        #print(f"dc.runGlideSteps({lowRateSteps})")
        #dc.runGlideSteps( lowRateSteps)

        #print(f"dc.writeConfigToTxt()")
        #dc.writeConfigToTxt()

        #print(f"measurementTimes, measurementIDs, measuredProperties = dc.getMechanicalMeasurements()")
        #measurementTimes, measurementIDs, measuredProperties \
        #        = dc.getMechanicalMeasurements()
        #measurementTimes = measurementTimes.squeeze()
        #measurementIDs = measurementIDs.squeeze()

        #try:
        #    print(f"trying to write measurements.pkl")
        #    with open(outputFolder + "/measurements3.pkl", 'wb') as handle:
        #        pickle.dump(
        #                measuredProperties,
        #                handle,
        #                protocol=pickle.HIGHEST_PROTOCOL
        #                )
        #    with open(outputFolder + "/times3.pkl", 'wb') as handle:
        #        pickle.dump(
        #                measurementTimes,
        #                handle,
        #                protocol=pickle.HIGHEST_PROTOCOL
        #                )
        #except:
        #    print(f"error: failed to write measurements3.pkl and times3.pkl")

        ## calculate total density time series
        #totalDensity = np.zeros( ( len( measuredProperties['density'][0])))
        #for ssID in measuredProperties['density'].keys():
        #    for timeStep in range( len(measuredProperties[ 'density'][ ssID])):
        #        totalDensity[ timeStep] += measuredProperties[ 'density'][ ssID][ timeStep]

        ##print(f"totalDensity.shape: {totalDensity.shape}")
        ##print(f"measurementTimes.shape: {measurementTimes.shape}")
        ##print(f"totalDensity: {totalDensity}")

        ## calculate strainRates
        #strainRates = dict() # dict of np.array, one per slip system
        #smoothedPlasticDistortion = dict()
        #smoothedStrainRates = dict()
        ##for ii in measuredProperties['slipSystemPlasticDistortion'].keys(): # keys are ssID
        ##    strainRates[ii] = np.ndarray(measurementTimes.shape[0]-1)
        ##    if (measurementTimes.shape[0] != measuredProperties['slipSystemPlasticDistortion'][ii].shape[0]):
        ##        print(f"error: measurementTimes.shape()[0] != measuredProperties['slipSystemPlasticDistortion'][ii].shape[0]")
        ##        return
        ##    # iterate over time step and calculate rates
        ##    for tt in range( measurementTimes.shape[0] -1):


        ##print(f"df:{df}")
        ##print(f"ax = df.plot()")
        ##ax = df.plot('time','totalDensity')
        ##print(f" fig.savefig(outputFolder + '/time_totalDensity.png')")
        ##ax.figure.savefig( outputFolder + "/time_totalDensity.png")
        ##df.to_pickle( outputFolder + "/time_totalDensity.pkl")
        ##except:
        ##    print(f"exception while saving plot {outputFolder}/time_totalDensity.png")
        ##try:
        ##    with open(outputFolder + "/densities.pkl", 'wb') as handle:
        ##        pickle.dump(
        ##                measuredProperties['density'],
        ##                handle,
        ##                protocol=pickle.HIGHEST_PROTOCOL
        ##                )
        ##    with open(outputFolder + "/times.pkl", 'wb') as handle:
        ##        pickle.dump(
        ##                measurementTimes,
        ##                handle,
        ##                protocol=pickle.HIGHEST_PROTOCOL
        ##                )
        ##    #print(f"trying to create dataframe time_vs_density")
        ##    #myDf = pd.DataFrame(
        ##    #    {
        ##    #        'time':measurementTimes,
        ##    #        'density':measuredProperties['density']
        ##    #        }
        ##    #    ).to_pickle( outputFolder + "/time_vs_density.pkl")
        ##    #print(f"trying to create {outputFolder}/time_vs_density.pkl")
        ##    #myDf.to_pickle( outputFolder + "/time_vs_density.pkl")
        ##except:
        ##    raise Exception("could not create pickle file "
        ##            + outputFolder + "/time_vs_density.pkl")

        ##create_density_plot(
        ##        measuredProperties['times'],
        ##        measuredProperties['density'],
        ##        outputFolder
        ##        )

        ##print(f"dc.getDensityPerSlipSystem():\n{dc.getDensityPerSlipSystem()}")
        ##print(f"measurementTimes:\n{measurementTimes}")
        ##print(f"measurementIDs:\n{measurementIDs}")
        ##print(f"measuredProperties:\n{measuredProperties}")

        ##measuredProperties['density'][1][2] = 3.3
        ##measuredDensities = measuredProperties['density']
        ##measuredDensities[1][1] = 100.0
        ##print(f"measuredDensities: {np.array(measuredDensities)}")
        ##print(f"np.shape(measuredDensities[23]): {np.shape(measuredDensities[23])}")
        ##print(f"len(measurementTimes):\n{len(measurementTimes)}")
        ##print(f"len(measurementIDs):\n{len(measurementIDs)}")
        ##print(f"len(measuredProperties['density'][0]:\n{len(measuredProperties['density'][0])}")
        ##print(f"len(measuredProperties['stressTensorComponent'][11]:\n{len(measuredProperties['stressTensorComponent'][11])}")

        ##dc.disableMechanicalMeasurements()
        ##dc.runGlideSteps( 10)
        ##dc.enableMechanicalMeasurements(2)
        ##dc.runGlideSteps( 10)
        ##measurementTimes, measurementIDs, measuredProperties \
        ##        = dc.getMechanicalMeasurements()
        ##print(f"measuredProperties['density'][0].flags:\n{measuredProperties['density'][0].flags}")
        ##print(f"measurementTimes:\n{measurementTimes}")
        ##print(f"measurementIDs:\n{measurementIDs}")
        ##print(f"measuredProperties:\n{measuredProperties}")
        ##
        ##print(f"len(measurementTimes):\n{len(measurementTimes)}")
        ##print(f"len(measurementIDs):\n{len(measurementIDs)}")
        ##print(f"len(measuredProperties['density'][0]:\n{len(measuredProperties['density'][0])}")
        ##print(f"len(measuredProperties['stressTensorComponent'][11]:\n{len(measuredProperties['stressTensorComponent'][11])}")
        ##dc.disableMechanicalMeasurements()
    except OSError:
        print(f"error: runInNewDir " + tag + " failed using parameters: {parameterLog[tag]}")
    except AssertionError:
        print(f"error: runInNewDir " + tag + " failed using parameters: {parameterLog[tag]}")
    print(f"outside of loop")
    print(f"parameterLog:\n{parameterLog}")

    return

def runInNewDir( dc, Nsteps, outputPath):
    print(f"output path: {outputPath}")
    if not os.path.exists( outputPath):
        print(f"creating directory: {outputPath}")
        os.makedirs( outputPath)
    if not os.path.exists( outputPath + "/evl"):
        print(f"creating directory: {outputPath + '/evl'}")
        os.makedirs( outputPath + "/evl")
    if not os.path.exists( outputPath + "/F"):
        print(f"creating directory: {outputPath + '/F'}")
        os.makedirs( outputPath + "/F")
    #print(f"dc.setCurrentStep(0) { dc.setCurrentStep(0)}")
    print(f"dc.setOutputPath( {outputPath})")
    dc.setOutputPath( outputPath)
    print(f"dc.generateMicrostructure()")
    dc.generateMicrostructure() # will write configuration to file

    # read the new configuration
    print(f"dc.readDefectiveCrystal()")
    dc.readDefectiveCrystal()

    print(f"dc.setOutputFrequency({Nsteps-1})")
    #dc.setOutputFrequency( Nsteps)
    #dc.setOutputFrequency( 1000)

    # attempt to change the external load applied
    strain \
        = np.array([
            [0,0,0],
            [0,0.0,0],
            [0,0,0]],dtype=np.float64)
    strainRate \
        = np.array([
            [0,0,0],
            [0,1.0e-7,0],
            [0,0,0]],dtype=np.float64)
    stress \
        = np.array([
            [0,0,0],
            [0,0,0],
            [0,0,0]],dtype=np.float64)
    stressRate \
        = np.array([
            [0.0,0.0,0.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0]],dtype=np.float64)
    stiffnessRatio \
        = np.array([0,1.0e20,0,0,0,0],dtype=np.float64)

    dc.setExternalLoad(
            strain=strain,
            strainRate=strainRate,
            stress=stress,
            stressRate=stressRate,
            stiffnessRatio=stiffnessRatio
            )

    #dc.enableMechanicalMeasurements( 10)

    print(f"dc.runGlideSteps({Nsteps})")
    dc.runGlideSteps(Nsteps)
    #print(f"dc.writeConfigToTxt()")
    dc.writeConfigToTxt()

    return

def createC2G( lattice, ss, nn ): #alignToSlipSystem0):
    #alignToSlipSystem0 = 1 # 1 : true, 0 : false
    ## create C2G matrix
    #if ( alignToSlipSystem0 ):
    #    if lattice == 'bcc':
    #        ss=np.array([5.773502691896258e-01,  5.773502691896258e-01, -5.773502691896258e-01])
    #        nn=np.array([8.660254037844385e-01, 0.000000000000000e+00, 8.660254037844385e-01])
    #        nn=nn/np.linalg.norm(nn)
    #    elif lattice == 'fcc':
    #        ss=np.array([ 0.000000000000000e+00, 7.071067811865475e-01, 7.071067811865475e-01])
    #        nn=np.array([-7.071067811865476e-01, 7.071067811865476e-01, -7.071067811865476e-01])
    #        nn=nn/np.linalg.norm(nn)
    #    else:
    #        return np.eye(3)
    c2g = np.array( [ss, np.cross(nn,ss), nn])
    #else:
    #    c2g = np.eye(3)
    return c2g

def create_measurement_plots(
    measuredProperties,
    measurementTimes,
    measurementIDs,
    strainRates,
    smoothedStrainRates,
    totalDensity,
    outputFolder):

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('time [s]')
    ax.set_ylabel(r'strainRates [s$^{-1}$]')
    for ssID in strainRates.keys():
        ax.plot( measurementTimes[1:], strainRates[ ssID], alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/time_strainRates.png")
    plt.close()

    for ssID in strainRates.keys():
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.set_xlabel('time [s]')
        ax.set_ylabel(r'strainRates [s$^{-1}$]')
        ax.plot( measurementTimes[1:], smoothedStrainRates[ ssID], alpha=0.5, label='smoothed')
        ax.plot( measurementTimes[1:], strainRates[ ssID], alpha=0.5, label='raw')
        ax.legend( loc="best", framealpha=0.0)
        plt.tight_layout()
        plt.savefig(outputFolder + "/time_strainRates_" + str(ssID) + ".png")
        plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('time [s]')
    ax.set_ylabel(r'smoothedStrainRates [s$^{-1}$]')
    for ssID in strainRates.keys():
        ax.plot( measurementTimes[1:], smoothedStrainRates[ ssID], alpha=0.5, label='smoothed')
    plt.tight_layout()
    plt.savefig(outputFolder + "/time_smoothedStrainRates" + ".png")
    plt.close()

    for ssID in measuredProperties['density'].keys():
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.set_xlabel('resolvedShearStress [Pa]')
        ax.set_ylabel(r'smoothedStrainRate [s$^{-1}$]')
        ax.plot(
                measuredProperties['resolvedShearStress'][ ssID][1:],
                smoothedStrainRates[ ssID],
                alpha=0.5)
        plt.tight_layout()
        plt.savefig(outputFolder + "/rss_plasticDeformationRate_" + str(ssID) + ".png")
        plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('resolvedShearStress [Pa]')
    ax.set_ylabel(r'smoothedStrainRate [s$^{-1}$]')
    for ssID in measuredProperties['density'].keys():
        ax.plot(
                measuredProperties['resolvedShearStress'][ ssID][1:],
                smoothedStrainRates[ ssID],
                alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/rss_plasticDeformationRates" + ".png")
    plt.close()


    for ssID in measuredProperties['density'].keys():
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.set_xlabel('time [s]')
        ax.set_ylabel(r'density [m$^{-2}$]')
        ax.plot( measurementTimes, measuredProperties['density'][ ssID], alpha=0.5)
        plt.tight_layout()
        plt.savefig(outputFolder + "/time_density_" + str(ssID) + ".png")
        plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('time [s]')
    ax.set_ylabel(r'density [m$^{-2}$]')
    for ssID in measuredProperties['density'].keys():
        ax.plot( measurementTimes, measuredProperties['density'][ ssID], alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/time_densities" + ".png")
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel(r'density [m$^{-2}$]')
    ax.set_ylabel(r'strainRate [s$^{-1}$]')
    for ssID in strainRates.keys():
        ax.plot( measuredProperties['density'][ ssID][1:],
                strainRates[ ssID], alpha=0.5, label='raw')
    plt.tight_layout()
    plt.savefig(outputFolder + "/density_strainRates" + ".png")
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel(r'density [m$^{-2}$]')
    ax.set_ylabel(r'smoothedStrainRate [s$^{-1}$]')
    for ssID in strainRates.keys():
        ax.plot( measuredProperties['density'][ ssID][1:],
                smoothedStrainRates[ ssID], alpha=0.5, label='smoothed')
    plt.tight_layout()
    plt.savefig(outputFolder + "/density_smoothStrainRates" + ".png")
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel(r'totalDensity [m$^{-2}$]')
    ax.set_ylabel(r'smoothedStrainRate [s$^{-1}$]')
    for ssID in strainRates.keys():
        ax.plot( totalDensity[1:],
                smoothedStrainRates[ ssID], alpha=0.5, label='smoothed')
    plt.tight_layout()
    plt.savefig(outputFolder + "/totalDensity_smoothStrainRates" + ".png")
    plt.close()

    for ssID in strainRates.keys():
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.set_xlabel(r'density [m$^{-2}$]')
        ax.set_ylabel(r'strainRate [s$^{-1}$]')
        ax.plot( measuredProperties['density'][ ssID][1:],
                smoothedStrainRates[ ssID], alpha=0.5, label='smoothed')
        ax.plot( measuredProperties['density'][ ssID][1:], strainRates[ ssID], alpha=0.5, label='raw')
        ax.legend( loc="best", framealpha=0.0)
        plt.tight_layout()
        plt.savefig(outputFolder + "/density_strainRate_" + str(ssID) + ".png")
        plt.close()


    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('time [s]')
    ax.set_ylabel(r'totalDensity [m$^{-2}$]')
    ax.plot( measurementTimes, totalDensity)
    plt.tight_layout()
    plt.savefig( outputFolder + "/time_totalDensity.png")

    #df = pd.DataFrame(
    #        {
    #            'time':measurementTimes,
    #            'totalDensity' :totalDensity
    #            #'plasticDistortion':measuredProperties['slipSystemPlasticDistortion']
    #            }
    #    )

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('time [s]')
    ax.set_ylabel('plasticDistortion')
    for ii in measuredProperties['slipSystemPlasticDistortion'].keys(): # keys are ssID
        ax.plot( measurementTimes, measuredProperties['slipSystemPlasticDistortion'][ii], alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/time_plasticDistortion.png")
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('time [s]')
    ax.set_ylabel('resolvedShearStress [Pa]')
    for ii in measuredProperties['resolvedShearStress'].keys(): # keys are ssID
        ax.plot( measurementTimes, measuredProperties['resolvedShearStress'][ii], alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/time_resolvedShearStresses.png")
    plt.close()

    for ii in measuredProperties['resolvedShearStress'].keys(): # keys are ssID
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.set_xlabel('time [s]')
        ax.set_ylabel('resolvedShearStress [Pa]')
        ax.plot( measurementTimes, measuredProperties['resolvedShearStress'][ii], alpha=0.5)
        plt.tight_layout()
        plt.savefig(outputFolder + "/time_resolvedShearStress_" + str(ii) + ".png")
        plt.close()


    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel('plasticDistortion')
    ax.set_ylabel('resolvedShearStress [Pa]')
    for ii in measuredProperties['resolvedShearStress'].keys(): # keys are ssID
        ax.plot(  measuredProperties['slipSystemPlasticDistortion'][ii], measuredProperties['resolvedShearStress'][ii], alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/plasticDistortion_resolvedShearStresses.png")
    plt.close()

    for ii in measuredProperties['resolvedShearStress'].keys(): # keys are ssID
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.set_xlabel('|plasticDistortion|')
        ax.set_ylabel('|resolvedShearStress| [Pa]')
        ax.plot(
            np.abs(measuredProperties['slipSystemPlasticDistortion'][ii]),
            np.abs(measuredProperties['resolvedShearStress'][ii]),
            alpha=0.5)
        plt.tight_layout()
        plt.savefig(outputFolder + "/plasticDistortion_resolvedShearStress_" + str(ii) + ".png")
        plt.close()


    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel(r'totalDensity [m$^{-2}$]')
    ax.set_ylabel('plasticDistortion')
    for ii in measuredProperties['slipSystemPlasticDistortion'].keys():
        ax.plot(
                totalDensity,
                measuredProperties['slipSystemPlasticDistortion'][ii],
                alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/plasticDistortion_totalDensity.png")
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.set_xlabel(r'density [m$^{-2}$]')
    ax.set_ylabel('plasticDistortion')
    for ii in measuredProperties['slipSystemPlasticDistortion'].keys():
        ax.plot(
                measuredProperties['density'][ii],
                measuredProperties['slipSystemPlasticDistortion'][ii],
                alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputFolder + "/plasticDistortion_density.png")
    plt.close()

    #rss = dict()
    #for ti in range(len(measuredTimes)):
    #    # calculate rss
    #    stress = np.array([
    #        [measuredProperties['stressTensorComponents'][11][ti],
    #         measuredProperties['stressTensorComponents'][12][ti],
    #         measuredProperties['stressTensorComponents'][13][ti]],
    #        [measuredProperties['stressTensorComponents'][12][ti],
    #         measuredProperties['stressTensorComponents'][22][ti],
    #         measuredProperties['stressTensorComponents'][23][ti]],
    #        [measuredProperties['stressTensorComponents'][13][ti],
    #         measuredProperties['stressTensorComponents'][23][ti],
    #         measuredProperties['stressTensorComponents'][33][ti]
    #        ])
    #    for ssid in measuredProperties['slipSystemPlasticDistortion'].keys():
    #        rss[ssid] = np.ndarray(measuredTimes.shape[0])
    #        print(f"plotting ss{ti}")
    #        ax.plot( measuredProperties['slipSystemPlasticDistortion'][ti], stress)
    #plt.tight_layout()
    #plt.savefig(outputFolder + "/time_plasticDistortion.png")
    #plt.close()

    return

def smooth_measurement( measurementSeries, times, windowTimePeriod):
    #print(f"smoothing data with window {windowTimePeriod}")
    smoothed = sm.nonparametric.lowess( exog=times, endog=measurementSeries, frac=(windowTimePeriod/(times[-1] - times[0])))
    return smoothed

def calculate_rate( measurementSeries, times):
    rates = np.ndarray( times.shape[0]-1)
    if ( times.shape[0] != measurementSeries.shape[0]):
        print(f"error: times.shape()[0] != measurementSeries.shape[0]")
    # iterate over time step and calculate rates
    for tt in range( times.shape[0] -1):
        delta_t = (times[tt+1] - times[tt])
        rates[ tt] = ( measurementSeries[tt+1] - measurementSeries[tt])/delta_t
    return rates


if __name__ == "__main__": main()
