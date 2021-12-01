/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetwork_CPP_
#define model_DislocationNetwork_CPP_


#include <DislocationNetwork.h>


namespace model
{
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    DislocationNetwork<dim,corder>::DislocationNetwork(int& argc, char* argv[],
                                                                         const DefectiveCrystalParameters& _simulationParameters,
                                                                         const SimplicialMesh<dim>& _mesh,
                                                                         const Polycrystal<dim>& _poly,
                                                                         const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                                                                         const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                                                                         const std::vector<VectorDim>& _periodicShifts,
                                                                         long int& runID) :
    //    /* init */ LoopNetworkType(_simulationParameters.isPeriodicSimulation()? std::shared_ptr<NetworkComponentType>(new NetworkComponentType()) : nullptr)
    /* init */ simulationParameters(_simulationParameters)
    /* init */,mesh(_mesh)
    /* init */,poly(_poly)
    /* init */,glidePlaneFactory(poly)
    /* init */,periodicGlidePlaneFactory(simulationParameters.isPeriodicSimulation()? new PeriodicGlidePlaneFactory<dim>(poly, glidePlaneFactory) : nullptr)
    //    /* init */,periodicDislocationLoopFactory(simulationParameters.isPeriodicSimulation()? new PeriodicDislocationLoopFactory<DislocationNetworkType>(poly,glidePlaneFactory) : nullptr)
    /* init */,bvpSolver(_bvpSolver)
    /* init */,externalLoadController(_externalLoadController)
    /* init */,periodicShifts(_periodicShifts)
    /* init */,networkRemesher(*this)
    /* init */,junctionsMaker(*this)
    /* init */,crossSlipMaker(*this)
    /* init */,nodeContractor(*this)
    /* init */,timeIntegrator("inputFiles/DD.txt")
    //    /* init */,gbTransmission(*this)
    //        /* init */,timeIntegrationMethod(TextFileParser("inputFiles/DD.txt").readScalar<int>("timeIntegrationMethod",true))
    ///* init */,maxJunctionIterations(TextFileParser("inputFiles/DD.txt").readScalar<int>("maxJunctionIterations",true))
    //        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true)),
    //        /* init */,totalTime(0.0),
    //        /* init */ dt(0.0),
    //        /* init */ vMax(0.0),
    //        /* init */ Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true)),
    //        /* init */,_plasticDistortionFromVelocities(MatrixDim::Zero())
    // /* init */,oldPlasticDistortionFromAreas(std::make_pair(0.0,MatrixDim::Zero()))
    // /* init */,_plasticDistortionRateFromVelocities(MatrixDim::Zero())
    // /* init */,_plasticDistortionRateFromAreas(MatrixDim::Zero())
    /* init */,ddSolverType(TextFileParser("inputFiles/DD.txt").readScalar<int>("ddSolverType",true))
    /* init */,computeDDinteractions(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeDDinteractions",true))
    /* init */,crossSlipModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("crossSlipModel",true))
    /* init */,outputFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFrequency",true))
    /* init */,outputBinary(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputBinary",true))
    /* init */,outputGlidePlanes(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputGlidePlanes",true))
    /* init */,outputElasticEnergy(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputElasticEnergy",true))
    /* init */,outputMeshDisplacement(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputMeshDisplacement",true))
    /* init */,outputFEMsolution(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFEMsolution",true))
    /* init */,outputDislocationLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputDislocationLength",true))
    //        /* init */,outputPlasticDistortion(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortion",true))
    /* init */,outputPlasticDistortionRate(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortionRate",true))
    /* init */,outputQuadraturePoints(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputQuadraturePoints",true))
    /* init */,outputLinkingNumbers(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLinkingNumbers",true))
    /* init */,outputLoopLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLoopLength",true))
    /* init */,outputSegmentPairDistances(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSegmentPairDistances",true))
    /* init */,computeElasticEnergyPerLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeElasticEnergyPerLength",true))
    //    /* init */,outputPeriodicConfiguration(simulationParameters.isPeriodicSimulation()? TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPeriodicConfiguration",true) : false)
    //        /* init */ _userOutputColumn(3)
    /* init */,use_stochasticForce(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_stochasticForce",true))
    /* init */,surfaceAttractionDistance(TextFileParser("inputFiles/DD.txt").readScalar<double>("surfaceAttractionDistance",true))
    //        /* init */,computePlasticDistortionRateFromVelocities(TextFileParser("inputFiles/DD.txt").readScalar<int>("computePlasticDistortionRateFromVelocities",true))
    /* init */,folderSuffix("")
    /* init */,use_velocityFilter(TextFileParser("inputFiles/DD.txt").readScalar<double>("use_velocityFilter",true))
    /* init */,velocityReductionFactor(TextFileParser("inputFiles/DD.txt").readScalar<double>("velocityReductionFactor",true))
    /* init */,verboseDislocationNode(TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationNode",true))
    {

        assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
        
        // Some sanity checks
        //            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
        
        // Initialize static variables
        LoopNetworkType::verboseLevel=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseLoopNetwork",true);
        verboseDislocationNetwork=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationNetwork",true);
        LoopType::initFromFile("inputFiles/DD.txt");
        LoopNodeType::initFromFile("inputFiles/DD.txt");
        LoopLinkType::initFromFile("inputFiles/DD.txt");
        NetworkLinkType::initFromFile("inputFiles/DD.txt");
        // NetworkNodeType::initFromFile("inputFiles/DD.txt");
        //        PeriodicDislocationBase::initFromFile("inputFiles/DD.txt");
        //        DislocationNetworkComponentType::initFromFile("inputFiles/DD.txt");
        DislocationStressBase<dim>::initFromFile("inputFiles/DD.txt");
        
//        DislocationCrossSlip<LoopNetworkType>::initFromFile("inputFiles/DD.txt");
//        DislocationJunctionFormation<LoopNetworkType>::initFromFile("inputFiles/DD.txt");
        /*
        int stochasticForceSeed=TextFileParser("inputFiles/DD.txt").readScalar<int>("stochasticForceSeed",true);
        if(stochasticForceSeed<0)
        {
            StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
        }
        else
        {
            StochasticForceGenerator::init(stochasticForceSeed);
        }
        */
        
        if(argc>1)
        {
            folderSuffix=argv[1];
            //                std::cout<<"folderSuffix="<<folderSuffix<<std::endl;
        }
        
        DDconfigIO<dim> evl(folderSuffix);
        evl.read(runID);
        setConfiguration(evl);
        createEshelbyInclusions();
    }
    
    /**********************************************************************/
    //Giacomo Version
    // template <int dim, short unsigned int corder>
    // void DislocationNetwork<dim,corder>::setConfiguration(const DDconfigIO<dim>& evl)
    // {
    //     this->loopLinks().clear(); // erase base network to clear current config
    
    //     // Create Loops
    //     std::deque<std::shared_ptr<LoopType>> tempLoops; // keep loops alive during setConfiguration
    //     size_t loopNumber=1;
    //     for(const auto& loop : evl.loops())
    //     {
    //         const bool faulted= poly.grain(loop.grainID).rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false;
    //         std::cout<<"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
    //         const size_t loopIDinFile(loop.sID);
    //         LoopType::set_count(loopIDinFile);
    //         GlidePlaneKey<dim> loopPlaneKey(loop.P,poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
    //         tempLoops.push_back(this->loops().create(loop.B,glidePlaneFactory.getFromKey(loopPlaneKey)));
    //         assert(this->loops().get(loopIDinFile)->sID==loopIDinFile);
    //         loopNumber++;
    //     }
    
    //     // Create NetworkNodes
    //     std::deque<std::shared_ptr<NetworkNodeType>> tempNetNodes; // keep loops alive during setConfiguration
    //     size_t netNodeNumber=1;
    //     for(const auto& node : evl.nodes())
    //     {
    //         std::cout<<"Creating DislocationNode "<<node.sID<<" ("<<netNodeNumber<<" of "<<evl.nodes().size()<<")"<<std::endl;
    //         const size_t nodeIDinFile(node.sID);
    //         NetworkNodeType::set_count(nodeIDinFile);
    //         tempNetNodes.push_back(this->networkNodes().create(node.P,node.V,node.velocityReduction));
    //         assert(this->networkNodes().get(nodeIDinFile)->sID==nodeIDinFile);
    //         netNodeNumber++;
    //     }
    
    //     // Create LoopNodes
    //     std::deque<std::shared_ptr<LoopNodeType>> tempLoopNodes; // keep loops alive during setConfiguration
    //     size_t loopNodeNumber=1;
    //     for(const auto& node : evl.loopNodes())
    //     {
    //         std::cout<<"Creating DislocationLoopNode "<<node.sID<<" ("<<loopNodeNumber<<" of "<<evl.loopNodes().size()<<")"<<std::flush;
    //         const size_t nodeIDinFile(node.sID);
    //         LoopNodeType::set_count(nodeIDinFile);
    //         const auto loop(this->loops().get(node.loopID));
    //         const auto netNode(this->networkNodes().get(node.networkNodeID));
    //         assert(loop && "Loop does not exist");
    //         assert(netNode && "NetworkNode does not exist");
    //         const auto periodicPatch(loop->periodicGlidePlane? loop->periodicGlidePlane->patches().getFromKey(node.periodicShift) : nullptr);
    //         const auto periodicPatchEdge((periodicPatch && node.edgeID>=0)? periodicPatch->edges()[node.edgeID] : nullptr);
    //         if(periodicPatch)
    //         {
    //             std::cout<<", on patch "<<periodicPatch->shift.transpose()<<std::flush;
    //         }
    //         if(periodicPatchEdge)
    //         {
    //             std::cout<<" on edge "<<periodicPatchEdge->edgeID<<std::flush;
    //         }
    //         std::cout<<std::endl;
    //         tempLoopNodes.push_back(this->loopNodes().create(loop,netNode,node.P,periodicPatch,periodicPatchEdge));
    //         assert(this->loopNodes().get(nodeIDinFile)->sID==nodeIDinFile);
    //         loopNodeNumber++;
    //     }
    
    //     // Insert Loops
    //     std::map<size_t,std::map<size_t,size_t>> loopMap;
    //     for(const auto& looplink : evl.loopLinks())
    //     {// Collect LoopLinks by loop IDs
    //         loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
    //     }
    //     assert(loopMap.size()==evl.loops().size());
    
    //     for(const auto& loop : evl.loops())
    //     {// for each loop in the DDconfigIO<dim> object
    
    //         const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
    //         assert(loopFound!=loopMap.end());
    //         std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
    //         loopNodes.push_back(this->loopNodes().get(loopFound->second.begin()->first));
    //         for(size_t k=0;k<loopFound->second.size();++k)
    //         {
    //             const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
    //             if(k<loopFound->second.size()-1)
    //             {
    //                 loopNodes.push_back(this->loopNodes().get(nodeFound->second));
    //             }
    //             else
    //             {
    //                 assert(nodeFound->second==loopNodes[0]->sID);
    //             }
    //         }
    //         this->insertLoop(this->loops().get(loop.sID),loopNodes);
    //     }
    //     updateGeometry();
    // }
    
    //New Version
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::setConfiguration(const DDconfigIO<dim>& evl)
    {
        this->loopLinks().clear(); // erase base network to clear current config
        
        // Create Loops
        std::deque<std::shared_ptr<LoopType>> tempLoops; // keep loops alive during setConfiguration
        size_t loopNumber=1;
        for(const auto& loop : evl.loops())
        {
            const bool faulted= poly.grain(loop.grainID).rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false;
            model::cout<<"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
            const size_t loopIDinFile(loop.sID);
            LoopType::set_count(loopIDinFile);

            
            switch (loop.loopType)
            {
                case DislocationLoopIO<dim>::GLISSILELOOP:
                {
                    GlidePlaneKey<dim> loopPlaneKey(loop.P, poly.grain(loop.grainID).reciprocalLatticeDirection(loop.N));
                    tempLoops.push_back(this->loops().create(loop.B, glidePlaneFactory.getFromKey(loopPlaneKey)));
                    assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
                    loopNumber++;
                    break;
                }
                case DislocationLoopIO<dim>::SESSILELOOP:
                {
                    tempLoops.push_back(this->loops().create(loop.B,loop.grainID,loop.loopType ));
                    assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
                    loopNumber++;
                    break;
                }
                default:
                    assert(false && "Unknown DislocationLoop type");
                    break;
            }
        }
        
        // Create NetworkNodes
        std::deque<std::shared_ptr<NetworkNodeType>> tempNetNodes; // keep loops alive during setConfiguration
        size_t netNodeNumber=1;
        for(const auto& node : evl.nodes())
        {
            model::cout<<"Creating DislocationNode "<<node.sID<<" ("<<netNodeNumber<<" of "<<evl.nodes().size()<<")"<<std::endl;
            const size_t nodeIDinFile(node.sID);
            NetworkNodeType::set_count(nodeIDinFile);
            tempNetNodes.push_back(this->networkNodes().create(node.P,node.V,node.velocityReduction));
            assert(this->networkNodes().get(nodeIDinFile)->sID==nodeIDinFile);
            netNodeNumber++;
        }
        
        // Create LoopNodes
        std::deque<std::shared_ptr<LoopNodeType>> tempLoopNodes; // keep loops alive during setConfiguration
        size_t loopNodeNumber=1;
        for(const auto& node : evl.loopNodes())
        {
            model::cout<<"Creating DislocationLoopNode "<<node.sID<<" ("<<loopNodeNumber<<" of "<<evl.loopNodes().size()<<")"<<std::flush;
            // std::cout<<"Printing stuff for loopNodes "<<std::endl;

            // std::cout << node.loopID << std::endl;
            // std::cout << node.sID << std::endl;
            // std::cout << std::setprecision(15) << std::scientific << node.P.transpose() << std::endl;
            // std::cout << node.networkNodeID << std::endl;
            // std::cout << std::setprecision(15) << std::scientific << node.periodicShift.transpose() << std::endl;
            // std::cout << node.edgeID << std::endl;

            const size_t nodeIDinFile(node.sID);
            LoopNodeType::set_count(nodeIDinFile);
            const auto loop(this->loops().get(node.loopID));
            const auto netNode(this->networkNodes().get(node.networkNodeID));
            assert(loop && "Loop does not exist");
            assert(netNode && "NetworkNode does not exist");
            const auto periodicPatch(loop->periodicGlidePlane? loop->periodicGlidePlane->patches().getFromKey(node.periodicShift) : nullptr);
            const auto periodicPatchEdge((periodicPatch && node.edgeIDs.first>=0)? (node.edgeIDs.second>=0 ? std::make_pair(periodicPatch->edges()[node.edgeIDs.first],
            periodicPatch->edges()[node.edgeIDs.second]): std::make_pair(periodicPatch->edges()[node.edgeIDs.first],nullptr)):std::make_pair(nullptr,nullptr)); 
            // std::cout<<"PeriodicPlane edge created "<<std::endl;
            if(periodicPatch)
            {
                model::cout<<", on patch "<<periodicPatch->shift.transpose()<<std::flush;
            }
            if(periodicPatchEdge.first)
            {
                model::cout<<" on edge "<<periodicPatchEdge.first->edgeID<<std::flush;
            }
            if(periodicPatchEdge.second)
            {
                model::cout<<" and "<<periodicPatchEdge.second->edgeID<<std::flush;
            }
            // std::cout<<"Trying to create the loop node with loopID "<<loop->sID<<" if loop has GP "<<(loop->glidePlane!=nullptr)<<" networkID is "
            // <<netNode->sID<<std::endl;
            model::cout<<std::endl;
            tempLoopNodes.push_back(this->loopNodes().create(loop,netNode,node.P,periodicPatch,periodicPatchEdge));
            assert(this->loopNodes().get(nodeIDinFile)->sID==nodeIDinFile);
            loopNodeNumber++;
        }
        
        // Insert Loops
        std::map<size_t,std::map<size_t,size_t>> loopMap;
        for(const auto& looplink : evl.loopLinks())
        {// Collect LoopLinks by loop IDs
            loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
        }
        assert(loopMap.size()==evl.loops().size());
        
        for(const auto& loop : evl.loops())
        {// for each loop in the DDconfigIO<dim> object
            
            const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
            assert(loopFound!=loopMap.end());
            std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
            loopNodes.push_back(this->loopNodes().get(loopFound->second.begin()->first));
            for(size_t k=0;k<loopFound->second.size();++k)
            {
                const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
                if(k<loopFound->second.size()-1)
                {
                    loopNodes.push_back(this->loopNodes().get(nodeFound->second));
                }
                else
                {
                    assert(nodeFound->second==loopNodes[0]->sID);
                }
            }
            std::cout<<" Inserting loop "<<loop.sID<<std::endl;
            this->insertLoop(this->loops().get(loop.sID),loopNodes);
        }
        updateGeometry();
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    const typename DislocationNetwork<dim,corder>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder>::eshelbyInclusions() const
    {
        return *this;
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder>::eshelbyInclusions()
    {
        return *this;
    }
    
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::createEshelbyInclusions()
    {
        for(const auto& grain : poly.grains())
        {
            EshelbyInclusion<dim>::addSlipSystems(grain.second.slipSystems());
        }
        
        
        IDreader<'E',1,14,double> inclusionsReader;
        inclusionsReader.read(0,true);
        
        const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
        for(const auto& pair : inclusionsReader)
        {
            
            const size_t& inclusionID(pair.first);
            Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
            
            const VectorDim C(row.template segment<dim>(0));
            const double a(row(dim+0));
            MatrixDim eT(MatrixDim::Zero());
            const int typeID(row(13));
            int k=dim+1;
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    eT(i,j)=row(k);
                    k++;
                }
            }
            
            
            
            EshelbyInclusion<dim>::set_count(inclusionID);
            eshelbyInclusions().emplace(std::piecewise_construct,
                                        std::make_tuple(inclusionID),
                                        std::make_tuple(C,a,eT,poly.nu,poly.mu,inclusionsMobilityReduction[typeID],typeID) );
        }
    }
    
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::updateGeometry()
    {
        VerboseDislocationNetwork(2,"DislocationNetwork::updateGeometry"<<std::endl;);
        for(auto& loop : this->loops())
        {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
            loop.second.lock()->updateGeometry();
        }
        // updatePlasticDistortionRateFromAreas();
        VerboseDislocationNetwork(3,"DislocationNetwork::updateGeometry DONE"<<std::endl;);
    }
    
    // template <int dim, short unsigned int corder>
    // void DislocationNetwork<dim,corder>::updatePlasticDistortionRateFromAreas()
    // {
    //     VerboseDislocationNetwork(2,"DislocationNetwork::updatePlasticDistortionRateFromAreas"<<std::endl;);
    //     const double dt(simulationParameters.totalTime-oldPlasticDistortionFromAreas.first);
    //     if(dt>=0.0)
    //     {
    //         const MatrixDim pd(plasticDistortion());
    //         if(dt>0.0)
    //         {
    //             _plasticDistortionRateFromAreas=(pd-oldPlasticDistortionFromAreas.second)/dt;
    //         }
    //         oldPlasticDistortionFromAreas=std::make_pair(simulationParameters.totalTime,pd); // update old values
    //     }
    //     VerboseDislocationNetwork(3,"DislocationNetwork::updatePlasticDistortionRateFromAreas DONE"<<std::endl;);
    // }
    
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::DislocationNetworkIOType DislocationNetwork<dim,corder>::io()
    {
        return DislocationNetworkIOType(*this,folderSuffix);
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::DislocationNetworkIOType DislocationNetwork<dim,corder>::io() const
    {
        return DislocationNetworkIOType(*this,folderSuffix);
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::plasticDistortion() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& loop : this->loops())
        {
            temp+= loop.second.lock()->plasticDistortion();
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::plasticDistortionRate() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& loop : this->loops())
        {
            temp+= loop.second.lock()->plasticDistortionRate();
        }
        return temp;
    }
    
    // template <int dim, short unsigned int corder>
    // const typename DislocationNetwork<dim,corder>::MatrixDim& DislocationNetwork<dim,corder>::plasticDistortionRate() const
    // {
    //     return  _plasticDistortionRateFromAreas;
    // }
    
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::plasticStrainRate() const
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
    }
    
    template <int dim, short unsigned int corder>
    std::tuple<double,double,double,double> DislocationNetwork<dim,corder>::networkLength() const
    {/*!\returns the total line length of the DislocationNetwork. The return
      * value is a tuple, where the first value is the length of bulk glissile
      * dislocations, the second value is the length of bulk sessile
      * dislocations, and the third value is the length accumulated on
      * the mesh boundary.
      */
        double bulkGlissileLength(0.0);
        double bulkSessileLength(0.0);
        double boundaryLength(0.0);
        double grainBoundaryLength(0.0);
        
        for(auto& loop : this->loops())
        {
            for(const auto& loopLink : loop.second.lock()->loopLinks())
            {
                if(loopLink->networkLink())
                {
                    if(!loopLink->networkLink()->hasZeroBurgers())
                    {
                        if(loopLink->networkLink()->isBoundarySegment())
                        {
                            boundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else if(loopLink->networkLink()->isGrainBoundarySegment())
                        {
                            grainBoundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else
                        {
                            if(loopLink->networkLink()->isSessile())
                            {
                                bulkSessileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                            else
                            {
                                bulkGlissileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                        }
                    }
                }
            }
        }
        return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
    }
    
    
    template <int dim, short unsigned int corder>
    bool DislocationNetwork<dim,corder>::contract(std::shared_ptr<NetworkNodeType> nA,
                                                                    std::shared_ptr<NetworkNodeType> nB)
    {
        return nodeContractor.contract(nA,nB);
    }
    
    

    
    
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::VectorDim DislocationNetwork<dim,corder>::displacement(const VectorDim& x) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        
        assert(!simulationParameters.isPeriodicSimulation());
        
        VectorDim temp(VectorDim::Zero());
        
        for(const auto& loop : this->loops())
        {// sum solid angle of each loop
            for(const auto& shift : periodicShifts)
            {
                temp-=loop.second.lock()->solidAngle(x+shift)/4.0/M_PI*loop.second.lock()->burgers();
            }
            //                if(!(loop.second->isVirtualBoundaryLoop() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES))
            //                {
            //                    for(const auto& shift : periodicShifts)
            //                    {
            //                        temp-=loop.second->solidAngle(x+shift)/4.0/M_PI*loop.second->burgers();
            //                    }
            //                }
        }
        
        for(const auto& link : this->networkLinks())
        {// sum line-integral part of displacement field per segment
            if(   !link.second.lock()->hasZeroBurgers()
               //                   && (!link.second->isBoundarySegment() || simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM)
               //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
               )
            {
                for(const auto& shift : periodicShifts)
                {
                    temp+=link.second.lock()->straight.displacement(x+shift);
                }
            }
        }
        
        return temp;
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(size_t k=0;k<fieldPoints.size();++k)
        {
            fieldPoints[k]=displacement(fieldPoints[k].P);
        }
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    typename DislocationNetwork<dim,corder>::MatrixDim DislocationNetwork<dim,corder>::stress(const VectorDim& x) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& link : this->networkLinks())
        {// sum stress field per segment
            if(   !link.second.lock()->hasZeroBurgers()
               && !(link.second.lock()->isBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_NO_FEM) // exclude boundary segments even if they are non-zero Burgers
               //                   && !(link.second->isVirtualBoundarySegment() && simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
               )
            {
                for(const auto& shift : periodicShifts)
                {
                    temp+=link.second.lock()->straight.stress(x+shift);
                }
            }
        }
        return temp;
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>>& fieldPoints) const
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(size_t k=0;k<fieldPoints.size();++k)
        {
            fieldPoints[k]=stress(fieldPoints[k].P);
        }
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::assembleAndSolveGlide()
    {/*! Performs the following operatons:
      */
#ifdef _OPENMP
        const size_t nThreads = omp_get_max_threads();
#else
        const size_t nThreads = 1;
#endif
        
        //! -1 Compute the interaction StressField between dislocation particles
        
        
        if(corder==0)
        {// For straight segments use analytical expression of stress field
            const auto t1= std::chrono::system_clock::now();
            std::cout<<"Computing analytical stress field at quadrature points ("<<nThreads<<" threads) "<<std::flush;
            for (const auto& links : this->networkLinks())
            {
                links.second.lock()->updateQuadraturePointsSeg();
            }
#ifdef _OPENMP
#pragma omp parallel for
            for (size_t k = 0; k < this->networkLinks().size(); ++k)
            {
                auto linkIter(this->networkLinks().begin());
                std::advance(linkIter, k);
                linkIter->second.lock()->assembleGlide();
            }
#else
            for (auto& linkIter : this->networkLinks())
            {
                linkIter.second.lock()->assembleGlide();
            }
#endif
            std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        else
        {// For curved segments use quandrature integration of stress field
            assert(0 && "ALL THIS MUST BE RE-IMPLEMENTED FOR CURVED SEGMENTS");
        }
        
        //! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
        DislocationGlideSolver<LoopNetworkType>(*this).solve();
    }
    
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::updateBoundaryNodes()
    {

        if (danglingBoundaryLoopNodes.size())
        {
            std::cout << "Removing bnd Nodes" << std::endl;
            for (const auto &node : danglingBoundaryLoopNodes)
            {
                this->removeLoopNode(node->sID);
            }
        }
            std::cout<<"Inserting new boundary nodes"<<std::endl;
            std::map<std::tuple<const NetworkNodeType*,const NetworkNodeType*,std::set<const PlanarMeshFace<dim>*>>,std::shared_ptr<NetworkNodeType>> newNetworkNodesMap;
            for (const auto &weakLoop : this->loops())
            {
                
                const auto loop(weakLoop.second.lock());
                if (loop->periodicGlidePlane)
                {
                    VerboseDislocationNetwork(1, " DisloationNetwork::DislocationLoop "<<loop->sID<<" Updating boundary nodes " << std::endl;);

                    std::vector<std::pair<VectorLowerDim, const LoopNodeType *const>> loopNodesPos;
                    std::map<std::tuple<const LoopNodeType *, const LoopNodeType *, const std::pair<const PeriodicPlaneEdge<dim> *,const PeriodicPlaneEdge<dim> *>>, const LoopNodeType *> bndNodesMap;
                    for (const auto &loopLink : loop->linkSequence())
                    {
                        if (!loopLink->source->periodicPlaneEdge.first)
                        {
                            loopNodesPos.emplace_back(loop->periodicGlidePlane->referencePlane->localPosition(loopLink->source->get_P()), loopLink->source.get());
                        }
                        else
                        {
                            bndNodesMap.emplace(std::make_tuple(loopLink->source->periodicPrev(), loopLink->source->periodicNext(), std::make_pair(loopLink->source->periodicPlaneEdge.first.get(),loopLink->source->periodicPlaneEdge.second.get())), loopLink->source.get());
                        }
                    }


                    if (loopNodesPos.size())
                    {
                        const auto polyInt(loop->periodicGlidePlane->polygonPatchIntersection(loopNodesPos)); // 2dPos,shift,edgeID,LoopNodeType*
                        std::map<const LoopNodeType *, size_t> polyIntMap;
                        for (size_t p = 0; p < polyInt.size(); ++p)
                        {
                            if (std::get<3>(polyInt[p]))
                            {
                                polyIntMap.emplace(std::get<3>(polyInt[p]), p);
                            }
                        }

                        for (size_t k = 0; k < loopNodesPos.size(); ++k)
                        {
                            const auto periodicPrev(loopNodesPos[k].second);
                            const auto periodicPrevNetwork(periodicPrev->networkNode.get());
                            const auto polyIter(polyIntMap.find(periodicPrev));
                            assert(polyIter != polyIntMap.end());
                            const size_t p(polyIter->second);

                            const size_t k1(k < loopNodesPos.size() - 1 ? k + 1 : 0);
                            const auto periodicNext(loopNodesPos[k1].second);
                            const auto periodicNextNetwork(periodicNext->networkNode.get());
                            const auto polyIter1(polyIntMap.find(periodicNext));
                            assert(polyIter1 != polyIntMap.end());
                            const size_t p1(polyIter1->second);

                            const auto periodicNetworkSource(periodicPrevNetwork->sID < periodicNextNetwork->sID ? periodicPrevNetwork : periodicNextNetwork);
                            const auto periodicNetworkSink(periodicPrevNetwork->sID < periodicNextNetwork->sID ? periodicNextNetwork : periodicPrevNetwork);

                            // std::cout<<"periodicPrev->periodicNext"<<periodicPrev->tag()<<"==>"<<periodicNext->tag()<<std::endl;

                            const LoopNodeType *currentSource(periodicPrev);
                            size_t p2 = (p + 1) % polyInt.size();
                            while (p2 < p1)
                            {
                                const auto periodicPatch(loop->periodicGlidePlane->getPatch(std::get<1>(polyInt[p2])));
                                // const auto periodicPatchEdge(periodicPatch->edges()[std::get<2>(polyInt[p2])]);
                                const auto periodicPatchEdge(std::get<2>(polyInt[p2]).second < 0 ? std::make_pair(periodicPatch->edges()[std::get<2>(polyInt[p2]).first], nullptr) : 
                                std::make_pair(periodicPatch->edges()[std::get<2>(polyInt[p2]).first], periodicPatch->edges()[std::get<2>(polyInt[p2]).second]));
                                const auto bndIter(bndNodesMap.find(std::make_tuple(periodicPrev, periodicNext, std::make_pair(periodicPatchEdge.first.get(),periodicPatchEdge.second.get()))));

                                if (bndIter != bndNodesMap.end())
                                { // exising bnd node found
                                    // std::cout<<" Using boundary node "<<bndIter->second->tag()<<std::endl;
                                    currentSource = bndIter->second;
                                }
                                else
                                {
                                    //Original

                                    // const VectorDim loopNodePos(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                                    // const VectorDim networkNodePos(loopNodePos + std::get<1>(polyInt[p2]));

                                    //There may be some error due to 2D usage of ssd.dMin
                                    const VectorDim loopNodePostemp(loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(polyInt[p2])));
                                    
                                    VectorDim networkNodePos(VectorDim::Zero());
                                    if (periodicPatchEdge.second)
                                    {
                                        SegmentSegmentDistance<dim> ssd (periodicPatchEdge.first->meshIntersection->P0,periodicPatchEdge.first->meshIntersection->P1,
                                        periodicPatchEdge.second->meshIntersection->P0,periodicPatchEdge.first->meshIntersection->P1);
                                        assert(ssd.dMin<FLT_EPSILON && "Two edges must intersect");
                                        networkNodePos=0.5*(ssd.x0+ssd.x1);
                                        assert((networkNodePos-(loopNodePostemp + std::get<1>(polyInt[p2]))).norm()<FLT_EPSILON && "Position mismatch");
                                    }
                                    else
                                    {
                                        networkNodePos=periodicPatchEdge.first->meshIntersection->snap(loopNodePostemp + std::get<1>(polyInt[p2]));
                                    }
                                    const VectorDim loopNodePos(networkNodePos - std::get<1>(polyInt[p2]));

                                    const auto currentLoopLink(currentSource->next.second);

                                    const auto currentNetworkLink(currentLoopLink->networkLink());

                                    // float u(0.0);
                                    // if (periodicPrevNetwork == periodicNetworkSource)
                                    // {
                                    //     const auto chord(periodicNext->get_P() - periodicPrev->get_P());
                                    //     const auto chordNorm2(chord.squaredNorm());
                                    //     u = (loopNodePos - periodicPrev->get_P()).dot(chord) / chordNorm2;
                                    //     // key = std::make_tuple(periodicNetworkSource, periodicNetworkSink, periodicPatchEdge->meshIntersection->faces, u);
                                    // }
                                    // else if (periodicNextNetwork == periodicNetworkSource)
                                    // {
                                    //     const auto chord(periodicPrev->get_P() - periodicNext->get_P());
                                    //     const auto chordNorm2(chord.squaredNorm());
                                    //     u=(loopNodePos - periodicNext->get_P()).dot(chord) / chordNorm2;
                                    //     // key = std::make_tuple(periodicNetworkSource, periodicNetworkSink, periodicPatchEdge->meshIntersection->faces, u);
                                    // }
                                    // else
                                    // {
                                    //     assert(false && "periodicNetworkSource must be periodicPrevNetwork or periodicNextNetwork");
                                    // }
                                    // if (u<FLT_EPSILON)
                                    // {
                                    //     u=0.0;
                                    // }
                                    // if (u>1.0-FLT_EPSILON)
                                    // {
                                    //     u=1.0;
                                    // }
                                    std::set <const PlanarMeshFace<dim>*>tmpMeshFaces;
                                    for (const auto& pmface : periodicPatchEdge.first->meshIntersection->faces)
                                    {
                                        tmpMeshFaces.emplace(pmface);
                                    }
                                    if (periodicPatchEdge.second)
                                    {
                                        for (const auto &pmface : periodicPatchEdge.second->meshIntersection->faces)
                                        {
                                            tmpMeshFaces.emplace(pmface);
                                        }
                                    }
                                    const auto key(std::make_tuple(periodicNetworkSource, periodicNetworkSink, tmpMeshFaces));

                                    const auto networkNodeIter(newNetworkNodesMap.find(key));
                                    if (networkNodeIter != newNetworkNodesMap.end())
                                    {
                                        // std::cout<<" CaseA key is "<<std::get<0>(key)->tag()<<" "<<std::get<1>(key)->tag()<<" and mesh faces "<<std::flush;
                                        // for (const auto& mf : std::get<2>(key))
                                        // {
                                        //     std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush; 
                                        // }
                                        // std::cout<<" and patch shift "<<std::get<1>(polyInt[p2]).transpose()<<std::endl;

                                        if ((networkNodeIter->second->get_P() - networkNodePos).norm() > 10000*FLT_EPSILON) //THis condition is just to check for widely different positions
                                        {
                                            std::cout<<std::scientific<<std::setprecision(15) << " PeriodicNetwork Source "<<periodicNetworkSource->sID <<"P= "<<periodicNetworkSource->get_P().transpose()<<
                                            "\n PeriodicNetwork Sink "<<periodicNetworkSink->sID<<" P = "<<periodicNetworkSink->get_P().transpose()<<std::endl; 
                                            std::cout << " Position of the network node "<<std::scientific<<std::setprecision(15) << networkNodeIter->second->get_P().transpose()<<std::endl;
                                            std::cout << " Actual Position of the network node should be " <<std::scientific<<std::setprecision(15) << networkNodePos.transpose()<<std::endl;
                                            std::cout << " Position difference "<<std::setprecision(15)<< (networkNodeIter->second->get_P()- networkNodePos).transpose()<<std::endl;
                                            std::cout << " Position difference norm "<<std::setprecision(15)<< (networkNodeIter->second->get_P()- networkNodePos).norm()<<std::endl;
                                            assert(false && "Loop node position and network node position mismatch");
                                        }
                                        //Here we want the loop node position to be commensurate with the network node (This is important for minimizing accumulating error)
                                        const VectorDim commensurateLoopNodePos(networkNodeIter->second->get_P()-std::get<1>(polyInt[p2]));
                                        const auto newLoopNode(this->loopNodes().create(loop, networkNodeIter->second, commensurateLoopNodePos, periodicPatch, periodicPatchEdge));
                                        currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                        
                                        //Loop into this this may be important
                                        // if((networkNodeIter->second->get_P()-networkNodePos).norm()> FLT_EPSILON )
                                        // {
                                        //     /*This is a very important definition...This means a new network node must be created...This condition may arise when a loop link exits from the edge of the cube
                                        //     Since for a loop node we can have just one shared pointer to the periodicplaneedge this means that the loop link exiting must
                                        //     from the edge of the cube must create nodes on all equivalen edges even through the actual connection will take place on diagnonally
                                        //     opposite ends. For such a case the network node position picked here may be different than the actual one... In such a case, we must
                                        //     create a new dislocation. If not satisfied enable the assert condition and let the code run... It will come into this conditio
                                        //     */
                                        //     const auto newNetNode(this->networkNodes().create(networkNodePos, VectorDim::Zero(), 1.0)); // TODO compute velocity and velocityReduction by interpolation
                                        //     std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                            
                                        //     newNetworkNodesMap.emplace(key, newNetNode);
                                            
                                        //     const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                        //     currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                        //     // std::cout<<" Position of the network node "<<networkNodeIter->second->get_P().transpose();
                                        //     // std::cout<<" Actual Position of the network node should be "<<networkNodePos.transpose();
                                        //     // assert(false && "Loop node position and network node position mismatch");
                                        // }
                                        // else
                                        // {
                                        
                                            
                                        // }
                                    }
                                    else
                                    {
                                        // std::cout << " CaseB key is " << std::get<0>(key)->tag() << " " << std::get<1>(key)->tag() << " and mesh faces " << std::flush;
                                        // for (const auto& mf : std::get<2>(key))
                                        // { 
                                        //     std::cout<<mf->sID<<" with outnormal"<<mf->outNormal().transpose()<< "\t "<<std::flush; 
                                        // }
                                        // std::cout<<" and patch shift "<<std::get<1>(polyInt[p2]).transpose()<<std::endl;

                                        const auto newNetNode(this->networkNodes().create(networkNodePos, VectorDim::Zero(), 1.0)); // TODO compute velocity and velocityReduction by interpolation
                                        // std::cout << "emplacing " << currentNetworkLink->tag() << "@" << std::setprecision(15) << std::scientific  << ", newNetNode=" << newNetNode->tag() << std::endl;
                                        newNetworkNodesMap.emplace(key, newNetNode);
                                        const auto newLoopNode(this->loopNodes().create(loop, newNetNode, loopNodePos, periodicPatch, periodicPatchEdge));
                                        currentSource = this->expandLoopLink(*currentLoopLink, newLoopNode).get();
                                    }
                                }
                                p2 = (p2 + 1) % polyInt.size();
                            }
                        }
                    }

                }
            }

            danglingBoundaryLoopNodes.clear();
        
    }
    
    
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::moveGlide(const double & dt_in)
    {/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
      */
        
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
        
        if(simulationParameters.isPeriodicSimulation())
        {
            
            
            danglingBoundaryLoopNodes.clear();
            
            // std::cout<<"Printing network connectivity info network Nodes "<<std::endl;
            
            // for (const auto& node : this->networkNodes())
            // {
            //     std::cout<<"Node ID "<<node.second.lock()->tag()<<"=>"<<node.second.lock()->glidePlanes().size()<<" neighbors size "
            //     <<node.second.lock()->neighbors().size()<<" loopNodes size "<<node.second.lock()->loopNodes().size()<<std::endl;
            //     std::cout<<" Printing the loop node info "<<std::endl;
            //     for (const auto& loopN : node.second.lock()->loopNodes())
            //     {
            //         std::cout<<"LoopNode Under consideration is "<<loopN->tag()<<std::endl;
            //         if (loopN->prev.first)
            //         {
            //             std::cout<<" Prev "<<loopN->prev.first->tag()<<std::endl;
            //             if (loopN->prev.second->networkLink())
            //             {
            //                 std::cout<<" Burgers of link "<<loopN->prev.second->networkLink()->burgers().transpose()<<std::endl;
            //             }
            //         }
            //         if (loopN->next.first)
            //         {
            //             std::cout<<" Next "<<loopN->next.first->tag()<<std::endl;
            //             if (loopN->next.second->networkLink())
            //             {
            //                 std::cout << " Burgers of link " << loopN->next.second->networkLink()->burgers().transpose() << std::endl;
            //             }
            //         }
            //     }
            //     std::cout<<" Printing neighbor links "<<std::endl;
            //     for (const auto& neigh : node.second.lock()->neighbors())
            //     {
            //         std::cout<<" Neighbor link tag "<<std::get<1>(neigh.second)->tag()<<" loopLinks size "<<std::get<1>(neigh.second)->loopLinks().size()
            //         <<" glidePlane size "<<std::get<1>(neigh.second)->glidePlanes().size()<<std::endl;
            //     }
            // }
            
            // std::cout<<"Printing GlidePlane Size for the network Nodes "<<std::endl;
            
            // for (const auto& node : this->networkNodes())
            // {
            //     std::cout<<"Node ID"<<node.second.lock()->tag()<<"=>"<<node.second.lock()->glidePlanes().size()<<std::endl;
            // }
            
            // std::cout<<"Printing Slip Systems of the network Links "<<std::endl;
            
            // for (const auto& link : this->networkLinks())
            // {
            //     std::cout<<" Link "<<link.second.lock()->tag()<<" [ "<<link.second.lock()->loopLinks().size()<<" ] "<<" has slipSystem Compare to nullPtr "<<(link.second.lock()->slipSystem()==nullptr)<<std::endl;
            //     std::cout<<" Glide Plane size "<<link.second.lock()->glidePlanes().size()<<" Burgers "<<link.second.lock()->burgers().transpose()<<" GlidePlaneNormal "
            //     <<link.second.lock()->glidePlaneNormal().transpose()<<std::endl;
            //     std::cout<<" Printing the slip ssytem of loopLinks "<<std::endl;
            //     if (link.second.lock()->source->sID==903 && link.second.lock()->sink->sID==79259)
            //     {
            //         for (const auto &looplink : link.second.lock()->loopLinks())
            //         {
            //             if (looplink->loop->slipSystem())
            //             {
            //                 std::cout << " For loopLink " << looplink->tag() << " slip system ID is " << looplink->loop->slipSystem()->sID << " Burgers "
            //                           << looplink->loop->slipSystem()->s.cartesian().transpose() << " and normal " << looplink->loop->slipSystem()->n.cartesian().transpose() << std::endl;
            //             }
            //         }
            //     }
            
            // }
            
            
            
            for(auto& node : this->networkNodes())
            {// Expansion
                // std::cout<<"Trying set p for "<<node.second.lock()->sID<<" LoopNdoe size ==>"<<node.second.lock()->loopNodes().size()<<" GlidePlane size ==>"
                // <<node.second.lock()->glidePlanes().size()<<std::endl;
                // std::cout<<std::setprecision(15)<<" Old Position "<<node.second.lock()->get_P().transpose()<<std::endl;
                // std::cout<<std::setprecision(15)<<" New Position "<<(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in).transpose()<<std::endl;
                // std::cout<<"Checking if the loop contains the position "<<std::endl;
                //                for (const auto& loopN : node.second.lock()->loopNodes())
                //                {
                //                    const VectorDim patchShift(loopN->periodicPlanePatch()? loopN->periodicPlanePatch()->shift : VectorDim::Zero());
                //                    const VectorDim oldPosition(loopN->get_P());
                //                    const VectorDim oldPositionN(node.second.lock()->get_P()-patchShift);
                //                    const VectorDim newPosition(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in-patchShift);
                //
                //                    // std::cout<<" Loop "<<loopN->loop()->sID<<" contains position OLD Position"<<loopN->loop()->glidePlane->contains(oldPosition)<<std::endl;
                //                    // std::cout<<" Loop "<<loopN->loop()->sID<<" contains position New Position"<<loopN->loop()->glidePlane->contains(newPosition)<<std::endl;
                //                }
                node.second.lock()->trySet_P(node.second.lock()->get_P()+node.second.lock()->get_V()*dt_in);
            }
            std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

            updateBoundaryNodes();
            
            
        }
        else
        {
            std::cout<<"        Moving DislocationNodes by glide (dt="<<dt_in<< ")... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            for (auto& nodeIter : this->networkNodes())
            {
                //                nodeIter.second.lock()->moveGlide(dt_in);
                nodeIter.second.lock()->set_P(nodeIter.second.lock()->get_P()+nodeIter.second.lock()->get_V()*dt_in);
            }
            std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        }
        
        
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder>
    void DislocationNetwork<dim,corder>::singleGlideStepDiscreteEvents(const long int& runID)
    {
        
        //#ifdef DislocationNucleationFile
        //            if(use_bvp && !(runID%use_bvp))
        //            {
        //                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
        //                updateQuadraturePoints();
        //            }
        //#endif
        //            assembleAndSolve(runID,straightSegmentsDeq);
        //            computeNodaVelocities(runID);
        
        //! 10- Cross Slip (needs upated PK force)
        //        DislocationCrossSlip<DislocationNetworkType> crossSlip(*this);
        //        crossSlip.execute();
        //
        //
        //        //                        gbTransmission.transmit();
        //        //gbTransmission.directTransmit();
        //        //            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();
        //        //            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations
        //
        //
        //        //! 12- Form Junctions
        // junctionsMaker.formJunctions(DDtimeIntegrator<0>::dxMax);  //Length of the junction is changed
        junctionsMaker.formJunctions(networkRemesher.Lmin);
        // io().output(simulationParameters.runID);
        //
        //        //            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
        //        //            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();
        //
        //
        //
        //        //! 13- Node redistribution
        networkRemesher.remesh(runID);
        //
        //        updateVirtualBoundaryLoops();
        
    }
    
    
    template <int dim, short unsigned int corder>
    int DislocationNetwork<dim,corder>::verboseDislocationNetwork=0;
    
    template class DislocationNetwork<3,0>;
    // template class DislocationNetwork<3,1>;
    
}
#endif
