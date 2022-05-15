/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_cpp_
#define model_MicrostructureGenerator_cpp_

#include <fstream>


#include <MicrostructureGenerator.h>
#include <PeriodicDipoleGenerator.h>
#include <PeriodicLoopGenerator.h>

namespace model
{

    /**********************************************************************/
    MicrostructureGenerator::MicrostructureGenerator(const std::string& folderName) :
    /* init*/ traitsIO(folderName)
    /* init*/,configIO(traitsIO.evlFolder)
    /* init*/,auxIO(traitsIO.auxFolder)
    /* init*/,outputBinary(TextFileParser(traitsIO.ddFile).readScalar<int>("outputBinary",true))
    /* init */,periodicFaceIDs(TextFileParser(traitsIO.polyFile).template readSet<int>("periodicFaceIDs",true))
//    /* init */,meshFilename(folderName+"/inputFiles/"+TextFileParser(folderName+"/inputFiles/polycrystal.txt").readString("meshFile",true))
    /* init */,mesh(traitsIO.meshFile,TextFileParser(traitsIO.polyFile).readMatrix<double>("A",3,3,true),
                    TextFileParser(traitsIO.polyFile).readMatrix<double>("x0",1,3,true).transpose(),periodicFaceIDs)
    /* init*/,minSize(0.1*std::min(mesh.xMax(0)-mesh.xMin(0),std::min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
    /* init*/,maxSize(std::max(mesh.xMax(0)-mesh.xMin(0),std::max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
    /* init*/,poly(traitsIO.polyFile,mesh)
    /* init*/,glidePlaneFactory(poly)
    /* init*/,periodicGlidePlaneFactory(poly, glidePlaneFactory)
    {
        
        std::cout<<greenBoldColor<<"Generating microstructure for "<<folderName<<defaultColor<<std::endl;
        
        
        // Some sanity checks
        if(mesh.volume()<FLT_EPSILON)
        {
            throw std::runtime_error("mesh "+traitsIO.meshFile+" is empty.");
        }
        
        
        std::ifstream polyFile(folderName+"/inputFiles/initialMicrostructure.txt");
        if(polyFile.is_open())
        {
            std::string line;
            while (std::getline(polyFile, line))
            {
                const std::string microstructureFileName(folderName+"/inputFiles/"+line);
                const std::string microstructureType(TextFileParser(microstructureFileName).readString("type",false));
                const std::string tag(TextFileParser(microstructureFileName).readString("tag",false));
                bool success(false);
                if(microstructureType=="PeriodicDipole")
                {
                    success=this->emplace(tag,new PeriodicDipoleGenerator(microstructureFileName)).second;
//                    PeriodicDipoleGenerator generator(microstructureFileName);
//                    generator.generate(*this);
                }
                else if(microstructureType=="PeriodicLoop")
                {
                    success=this->emplace(tag,new PeriodicLoopGenerator(microstructureFileName)).second;

//                    PeriodicLoopGenerator generator(microstructureFileName);
//                    generator.generate(*this);
                }
                else
                {
                    std::cout<<"unkown microstructure type "<<microstructureType<<std::endl;
                }
                if(!success)
                {
                    throw std::runtime_error("Duplicate microstructure tag "+tag+".");
                }
            }
        }
        else
        {
            throw std::runtime_error("Cannot open file "+folderName+"/inputFiles/initialMicrostructure.txt");
        }
        
        for(auto& gen : *this)
        {
            if(gen.second->style=="individual")
            {
                gen.second->generateIndividual(*this);
            }
            else if(gen.second->style=="density")
            {
                gen.second->generateDensity(*this);
            }
            else
            {
                throw std::runtime_error("Uknown style for generator "+gen.second->tag);
            }
        }
        
        // Call individual generators
        //            addStraightDislocations();
        //            addFrankReadSources();
        //            addSingleArmDislocations();
        //            addPrismaticLoops();
        //            addIndividualStraightDislocations();
        //            addFrankLoops();
        //            addNonPlanarLoops();
        //            // addPeriodicLoops();
        //            addStatisticallyHomegeneousPeriodicLoops();
        //    addIndividualPeriodicStraighDislocation();
        //            addStatisticallyHomegeneousPlanarDipolarLoops();
        //            addPeriodicJunctionLoops();
        //            addIrradiationLoops();
        //            addStackingFaultTetrahedra();
        //            addEshelbyInclusions();
        writeConfigFiles(0);
        
    }


    size_t MicrostructureGenerator::insertLoop(const VectorDimD& b,const VectorDimD& unitNormal,const VectorDimD& P0,const size_t& grainID,const DislocationLoopType& loopType)
    {
        const size_t loopID(configIO.loops().size());
        configIO.loops().emplace_back(loopID, b,unitNormal,P0,grainID,loopType);
        return loopID;
    }

    size_t MicrostructureGenerator::insertLoopNode(const size_t& loopID,const VectorDimD& loopNodePos,const size_t& networkNodeID,const VectorDimD& loopNodeShift,const std::pair<short int,short int>& periodicEdgeIDs)
    {
        const size_t loopNodeID(configIO.loopNodes().size());
        configIO.loopNodes().emplace_back(loopNodeID,loopID,loopNodePos,networkNodeID,loopNodeShift,periodicEdgeIDs);
        return loopNodeID;
    }

    std::vector<size_t> MicrostructureGenerator::insertLoopLinks(const size_t& loopID,const std::vector<size_t>& loopNodeIDs)
    {
        std::vector<size_t> temp;
        for(size_t k=0;k<loopNodeIDs.size();++k)
        {
            const size_t k1=(k+1)<loopNodeIDs.size()? k+1 : 0;
            temp.push_back(configIO.loopLinks().size());
            const size_t sourceNodeID(loopNodeIDs[k ]);
            const size_t   sinkNodeID(loopNodeIDs[k1]);
            const auto sourceNode(configIO.loopNodes()[sourceNodeID]);
            const auto   sinkNode(configIO.loopNodes()[sinkNodeID]);
            configIO.loopLinks().emplace_back(loopID,sourceNodeID,sinkNodeID,(sourceNode.P-sinkNode.P).norm()>FLT_EPSILON,0);
        }
        return temp;
    }


    size_t MicrostructureGenerator::insertNetworkNode(const VectorDimD& networkNodePos)
    {
        const size_t networkNodeID(configIO.nodes().size());
        configIO.nodes().emplace_back(networkNodeID,networkNodePos,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
        return networkNodeID;
    }




    /**********************************************************************/
    void MicrostructureGenerator::writeConfigFiles(const size_t& fileID)
    {
        
        //    auxIO.setGlidePlaneBoundaries(glidePlaneFactory); // change this function to take a GlidePlaneFactory during write
        
        if(outputBinary)
        {
            std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getBinFilename(fileID)<<defaultColor<<std::endl;
            configIO.writeBin(fileID);
            auxIO.writeBin(fileID);
        }
        else
        {
            std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getTxtFilename(fileID)<<defaultColor<<std::endl;
            configIO.writeTxt(fileID);
            auxIO.writeTxt(fileID);
        }
    }


}
#endif