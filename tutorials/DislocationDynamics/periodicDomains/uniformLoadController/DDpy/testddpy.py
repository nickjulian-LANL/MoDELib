#!/usr/bin/python3

import numpy as np
#import pickle
import ddpy
import os

def main():

    dddFolderPath = "./testDDpy"#contains inputFiles/{DD.txt,mesh,load,material}
    lattice = 'bcc' # TODO: read crystalStructure from {material}.txt
    #lattice = 'fcc'
    # TODO: create uniformExternalLoadController.txt via python
    # TODO: create DD.txt via python
    # TODO: create mesh, material, and load files via python

    meshFilePath = "./unitCube.msh" # relative to inputFiles directory
    material = 'Fe_320' # used to find file {material}.txt
    #material = 'Cu'

    # simulation box dimensions [m]
    xlo = 0
    xhi = 1000
    ylo = 0
    yhi = 1000
    zlo = 0
    zhi = 1000

    #alignToSlipSystem0 = 1 # 1 : true, 0 : false

    nSteps = 11

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

    # runInNewDir( dc, Nsteps, c2g, tagIn, outputPath):
    tag = "pdl0"
    outputFolder = dddFolderPath + "/testrun1" # to be created within dddFolderPath
    
    # specify individual dipoles
    ssIDs = np.array([0,1])
    exitFaceIDs = np.array([1,1])
    points = np.array([[30,30,30],[40,40,40]])
    heights = np.array([3,4])
    dipoleNodes = np.array([0,0])
    dipoleGlideSteps = np.array([10,20])


    dc.clearMicrostructureSpecifications()
    dc.specifyDipoles(
            tag=tag,
            slipSystemIDs=ssIDs,
            exitFaceIDs=exitFaceIDs,
            points=points,
            heights=heights,
            dipoleNodes=dipoleNodes,
            dipoleGlideSteps=dipoleGlideSteps 
            )

    print(f"trying to run in {outputFolder}")
    runInNewDir( dc, 2, c2g, tag, outputFolder)

    print(f"dc.setCurrentStep(0)")
    dc.setCurrentStep(0)

    print("calling dc.getPlasticStrains()")
    plasticStrains = dc.getPlasticStrains()
    # plasticStrains now contains one 1x3 vector for each strained slip system
    print(f"calling dc.getDensityPerSlipSystem()")
    densities = dc.getDensityPerSlipSystem()
    print(f"densities:\n{densities}")

    print("calling dc.getResolvedShearStrains()")
    shearStrains = dc.getResolvedShearStrains()
    print("calling dc.getResolvedShearStresses()")
    shearStresses = dc.getResolvedShearStresses()

    print(f"plasticStrains:\n{plasticStrains}")
    print(f"shearStrains:\n{shearStrains}")
    print(f"shearStresses:\n{shearStresses}")

    outputFolder = dddFolderPath + "/testrun2" # to be created within dddFolderPath
    tag = "pll0"
    # specify individual dipoles
    ssIDs = np.array([ 4, 0])
    radii = np.array([ 1.0e-8, 3.0e-8])
    loopSegmentCounts = np.array([ 100, 100]) # periodicLoopSides
    loopCenters = np.array([[ 50, 100, 150],
                            [ 150, 100, 50]])# [m]


    dc.clearMicrostructureSpecifications()
    dc.specifyLoops(
            tag=tag,
            slipSystemIDs=ssIDs,
            loopRadii=radii,
            loopSegmentCounts=loopSegmentCounts,
            loopCenters=loopCenters
            )
    print(f"trying to run in {outputFolder}")
    runInNewDir( dc, 3, c2g, tag, outputFolder)

    print(f"dc.setCurrentStep(0)")
    dc.setCurrentStep(0)

    print("calling dc.getPlasticStrains()")
    plasticStrains = dc.getPlasticStrains()
    print("calling dc.getResolvedShearStrains()")
    shearStrains = dc.getResolvedShearStrains()
    print("calling dc.getResolvedShearStresses()")
    shearStresses = dc.getResolvedShearStresses()
    print(f"plasticStrains:\n{plasticStrains}")
    print(f"shearStrains:\n{shearStrains}")
    print(f"shearStresses:\n{shearStresses}")

    densities = dc.getDensityPerSlipSystem()
    print(f"densities:\n{densities}")

    outputFolder = dddFolderPath + "/testrun3" # to be created within dddFolderPath
    # specify a density of loops
    tag = "pld0"
    periodicLoopTargetDensity = 1e15;
    periodicLoopSegmentCount = 100;
    periodicLoopRadiusDistributionMean = 4e-8;
    periodicLoopRadiusDistributionStd = 1e-8;

    dc.clearMicrostructureSpecifications()
    dc.specifyLoopDensity(
            tag=tag,
            loopDensity=periodicLoopTargetDensity,
            loopSegmentCount=periodicLoopSegmentCount,
            loopRadiusMean=periodicLoopRadiusDistributionMean,
            loopRadiusStd=periodicLoopRadiusDistributionStd
            )
    print(f"trying to run in {outputFolder}")
    runInNewDir( dc, 3, c2g, tag, outputFolder)

    print(f"dc.setCurrentStep(0)")
    dc.setCurrentStep(0)

    print("calling dc.getPlasticStrains()")
    plasticStrains = dc.getPlasticStrains()
    print("calling dc.getResolvedShearStrains()")
    shearStrains = dc.getResolvedShearStrains()
    print("calling dc.getResolvedShearStresses()")
    shearStresses = dc.getResolvedShearStresses()
    print(f"plasticStrains:\n{plasticStrains}")
    print(f"shearStrains:\n{shearStrains}")
    print(f"shearStresses:\n{shearStresses}")

    densities = dc.getDensityPerSlipSystem()
    print(f"densities:\n{densities}")

    outputFolder = dddFolderPath + "/testrun4" # to be created within dddFolderPath
    # specify a density of loops
    tag = "pdd0"
    periodicDipoleTargetDensity = 1e15;

    dc.clearMicrostructureSpecifications()
    dc.specifyDipoleDensity(
            tag=tag,
            dipoleDensity=periodicDipoleTargetDensity
            )
    print(f"trying to run in {outputFolder}")
    runInNewDir( dc, 10, c2g, tag, outputFolder)

    print("calling dc.getPlasticStrains()")
    plasticStrains = dc.getPlasticStrains()
    print("calling dc.getResolvedShearStrains()")
    shearStrains = dc.getResolvedShearStrains()
    print("calling dc.getResolvedShearStresses()")
    shearStresses = dc.getResolvedShearStresses()
    print(f"plasticStrains:\n{plasticStrains}")
    print(f"shearStrains:\n{shearStrains}")
    print(f"shearStresses:\n{shearStresses}")

    densities = dc.getDensityPerSlipSystem()
    print(f"densities:\n{densities}")

    # attempt to change the external load applied

    strain \
        = np.array([
            [0,0,0],
            [0,0.0,0],
            [0,0,0]],dtype=np.float64)
    strainRate \
        = np.array([
            [0,0,0],
            [0,0,0],
            [0,0,0]],dtype=np.float64)
    stress \
        = np.array([
            [0,0,0],
            [0,2e-6,0],
            [0,0,0]],dtype=np.float64)
    stressRate \
        = np.array([
            [0.0,0.0,0.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0]],dtype=np.float64)
    stiffnessRatio \
        = np.array([0,0.0,0,0,0,0],dtype=np.float64)

    dc.setExternalLoad(
            #strain=strain,
            strainRate=strainRate,
            stress=stress,
            #stressRate=stressRate,
            stiffnessRatio=stiffnessRatio
            )

    #print("after setting load: dc.getPlasticStrains()")
    plasticStrains = dc.getPlasticStrains()
    #print("after setting load: dc.getResolvedShearStrains()")
    shearStrains = dc.getResolvedShearStrains()
    #print("after setting load: dc.getResolvedShearStresses()")
    shearStresses = dc.getResolvedShearStresses()
    print(f"after setting load, plasticStrains:\n{plasticStrains}")
    print(f"after setting load, shearStrains:\n{shearStrains}")
    print(f"after setting load, shearStresses:\n{shearStresses}")

    print(f"dc.runGlideSteps({10})")
    dc.runGlideSteps(10)
    dc.writeConfigToTxt()
    #print("after setting load: dc.getPlasticStrains()")
    plasticStrains = dc.getPlasticStrains()
    #print("after setting load: dc.getResolvedShearStrains()")
    shearStrains = dc.getResolvedShearStrains()
    #print("after setting load: dc.getResolvedShearStresses()")
    shearStresses = dc.getResolvedShearStresses()
    print(f"after setting load and running 2 steps, plasticStrains:\n{plasticStrains}")
    print(f"after setting load and running 2 steps, shearStrains:\n{shearStrains}")
    print(f"after setting load and running 2 steps, shearStresses:\n{shearStresses}")
    return

def runInNewDir( dc, Nsteps, c2g, tagIn, outputPath):

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
    print(f"dc.setOutputPath( {outputPath})")
    dc.setOutputPath( outputPath)
    print(f"dc.generateMicrostructure()")
    dc.generateMicrostructure() # will write configuration to file

    # read the new configuration
    print(f"dc.readDefectiveCrystal()")
    dc.readDefectiveCrystal()

    print(f"dc.setOutputFrequency(1)")
    dc.setOutputFrequency(1)
    print(f"dc.runGlideSteps({Nsteps})")
    dc.runGlideSteps(Nsteps)
    #print(f"dc.writeConfigToTxt()")
    #dc.writeConfigToTxt()

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

if __name__ == "__main__": main()
