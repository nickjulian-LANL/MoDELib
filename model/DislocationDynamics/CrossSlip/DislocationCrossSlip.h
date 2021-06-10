/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_H_
#define model_DislocationCrossSlip_H_

#include <utility> // for std::pair
#include <vector>
#include <Eigen/Dense>
#include <SegmentSegmentDistance.h>
//#include <DislocationSegmentIntersection.h>
// #include <DislocationNetworkRemesh.h>
#include <CrossSlipModels.h>

#include <PlanePlaneIntersection.h>

#include <MPIcout.h>
#include <EqualIteratorRange.h>
#include <N2IteratorRange.h>
#include <BCClattice.h>
#include <FCClattice.h>


#ifndef NDEBUG
#define VerboseCrossSlip(N,x) if(verboseCrossSlip>=N){model::cout<<x;}
#else
#define VerboseCrossSlip(N,x)
#endif


namespace model
{
    
    template <typename DislocationNetworkType>
    class DislocationCrossSlip
    {
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        typedef typename DislocationNetworkType::NetworkLinkType NetworkLinkType;
        typedef typename DislocationNetworkType::NetworkNodeType NetworkNodeType;
        typedef typename DislocationNetworkType::LoopNodeType LoopNodeType;
        // typedef typename DislocationNetworkType::IsNetworkEdgeType IsNetworkLinkType;
        // typedef typename DislocationNetworkType::IsNodeType IsNodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        // typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        
//        typedef std::tuple<size_t,size_t,size_t,size_t> CrossSlipTupleType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;

        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        //! A reference to the DislocationNetwork
        DislocationNetworkType& DN;
        CrossSlipContainerType crossSlipDeq;
        
    public:
        
        const int verboseCrossSlip;
        const double crossSlipDeg;
        
        /**********************************************************************/
        DislocationCrossSlip(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,verboseCrossSlip(TextFileParser("inputFiles/DD.txt").readScalar<double>("crossSlipDeg",true))
        /* init */,crossSlipDeg(TextFileParser("inputFiles/DD.txt").readScalar<double>("crossSlipDeg",true))
        {
            assert(crossSlipDeg>=0.0 && DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");

//            if(DN.crossSlipModel)
//            {
//                const auto t0= std::chrono::system_clock::now();
//                model::cout<<"Finding CrossSlip segments: "<<std::flush;
//                crossSlipDeq=findCrossSlipSegments(DN.poly,DN.crossSlipModel);
//                VerboseCrossSlip(1,crossSlipDeq.size()<<" found"<<std::endl;);
//                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
            
        }
        
        /**********************************************************************/
        void findCrossSlipSegments()
        {

            if(DN.crossSlipModel)
            {
                model::cout<<"Finding CrossSlip segments... "<<std::flush;
                const auto t0= std::chrono::system_clock::now();

             crossSlipDeq.clear();
            const double sinCrossSlipRad(std::sin(crossSlipDeg*M_PI/180.0));
//            CrossSlipContainerType crossSlipDeq;

            
            for(const auto& linkIter : DN.networkLinks())
            {
                const auto link(linkIter.second.lock());
                if(   !link->isBoundarySegment()
                   && !link->source->isBoundaryNode()
                   && !link->sink->isBoundaryNode()
                   && !link->isGrainBoundarySegment()
                   && !link->source->isGrainBoundaryNode()
                   && !link->sink->isGrainBoundaryNode()
                   && !link->hasZeroBurgers()
                   && link->isGlissile()
                   && link->chord().normalized().cross(link->burgers().normalized()).norm()<=sinCrossSlipRad
                   && link->chord().norm()>2.0*DN.networkRemesher.Lmin
                   )
                {
                    //                    const auto& grain(**link.second->grains().begin());
                    
                    
                    if(DN.poly.crystalStructure=="BCC")
                    {
                        CrossSlipModels<BCClattice<dim>>::addToCrossSlip(*link,crossSlipDeq,DN.crossSlipModel);
                    }
                    else if(DN.poly.crystalStructure=="FCC")
                    {
                        CrossSlipModels<FCClattice<dim>>::addToCrossSlip(*link,crossSlipDeq,DN.crossSlipModel);
                    }
                    else
                    {
                        std::cout<<"Unknown cross-slip model for crystal structure '"<<DN.poly.crystalStructure<<"'. Exiting."<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                    
                }
            }
                model::cout<<crossSlipDeq.size()<<" found "<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            
//            return crossSlipDeq;
        }
        
//        /******************************************************************/
//        static void initFromFile(const std::string& fileName)
//        {
//            crossSlipDeg=TextFileParser(fileName).readScalar<double>("crossSlipDeg",true);
//            assert(crossSlipDeg>=0.0 && DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
//            verboseCrossSlip=TextFileParser(fileName).readScalar<int>("verboseCrossSlip",true);
//        }
        
        /**********************************************************************/
        void execute()
        {
            if(DN.crossSlipModel)
            {
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"Executing cross slip "<<std::flush;
            size_t executed(0);
            for(const auto& tup : crossSlipDeq)
            {
                const std::shared_ptr<NetworkNodeType>& source(std::get<0>(tup));
                const std::shared_ptr<NetworkNodeType>& sink(std::get<1>(tup));
                const size_t& sourceID(source->sID);
                const size_t& sinkID(sink->sID);
                const size_t& grainID(std::get<2>(tup));
                const size_t& slipID(std::get<3>(tup));
                
                const std::shared_ptr<NetworkNodeType> isSource(DN.networkNodes().get(sourceID));
                const std::shared_ptr<NetworkNodeType> isSink(DN.networkNodes().get(sinkID));
                const auto isLink(DN.networkLinks().get(std::make_pair(sourceID,sinkID)));
                
                const auto& crosSlipSystem(DN.poly.grain(grainID).slipSystems()[slipID]); // last element in map has highest pkGlide
                
                if(isSource && isSink && isLink)
                {
                    
                    // Align source and sink to perfect screw orientation
                    const VectorDim midPoint(0.5*(isSource->get_P()+isSink->get_P()));
                    const long int height(crosSlipSystem->n.closestPlaneIndexOfPoint(midPoint));
                    //
                    //                        const int height=LatticePlane::computeHeight(crosSlipSystem->n,midPoint).second;
                    const VectorDim planePoint(height*crosSlipSystem->n.planeSpacing()*crosSlipSystem->unitNormal);
                    
                    //const VectorDim planePoint2=midPoint-(midPoint-planePoint).dot(crosSlipSystem->unitNormal)*crosSlipSystem->unitNormal; // closest point to midPoint on the crossSlip plane
                    
                    //                        PlanePlaneIntersection<dim> ppi(midPoint,isLink.second->glidePlaneNormal(),
                    //                                                        planePoint2,crosSlipSystem->unitNormal);
                    
                    
                    
                    
                    PlanePlaneIntersection<dim> ppi((*isLink->loopLinks().begin())->loop->glidePlane->P,
                                                    (*isLink->loopLinks().begin())->loop->glidePlane->unitNormal,
                                                    planePoint,
                                                    crosSlipSystem->unitNormal);
                    
                    
                    const VectorDim newSourceP(ppi.P+(isSource->get_P()-ppi.P).dot(ppi.d)*ppi.d);
                    const VectorDim newSinkP(ppi.P+(isSink->get_P()-ppi.P).dot(ppi.d)*ppi.d);

                        //Only one loopNode
                        //New position do not belong to the patch boundary based on the glide plane intersection
                    bool sourcePositiononBoundary(false);
                    bool sinkPositiononBoundary(false);

                    for (const auto& gp : isSource->glidePlanes())
                    {
                        if (gp->meshIntersections.contains(newSourceP))
                        {
                            sourcePositiononBoundary=true;
                            break;
                        }
                    }
                        
                    for (const auto& gp : isSink->glidePlanes())
                    {
                        if (gp->meshIntersections.contains(newSourceP))
                        {
                            sinkPositiononBoundary=true;
                            break;
                        }
                    }

                    if (!sourcePositiononBoundary && !sinkPositiononBoundary)
                    {
                        if (isSource->isMovableTo(newSourceP) && isSink->isMovableTo(newSinkP))
                        {

                            VerboseCrossSlip(1, "cross-slip " << sourceID << "->" << sinkID << std::endl;);

                            // Re-align source and sink
                            isSource->trySet_P(newSourceP);
                            isSink->trySet_P(newSinkP);

                            DN.updateBoundaryNodes();
                            if ((isSource->get_P() - newSourceP).norm() < FLT_EPSILON && (isSink->get_P() - newSinkP).norm() < FLT_EPSILON)
                            {

                                // Check if source and sink are already part of loops on the conjugate plane

                                // Construct and insert new loop in conjugate plane
                                const VectorDim newNodeP(0.5 * (isSource->get_P() + isSink->get_P()));
                                //                                const size_t newNodeID=DN.insertDanglingNode(newNodeP,VectorDim::Zero(),1.0).first->first;
                                std::shared_ptr<NetworkNodeType> newNode(DN.networkNodes().create(newNodeP, VectorDim::Zero(), 1.0));

                                //                                std::vector<size_t> nodeIDs;
                                std::vector<std::shared_ptr<NetworkNodeType>> networkNodes;

                                networkNodes.push_back(isSink);   // insert in reverse order, sink first, source second
                                networkNodes.push_back(isSource); // insert in reverse order, sink first, source second
                                networkNodes.push_back(newNode);

                                GlidePlaneKey<dim> loopPlaneKey(newNodeP, DN.poly.grain(grainID).slipSystems()[slipID]->n);
                                const auto glidePlane(DN.glidePlaneFactory.getFromKey(loopPlaneKey));
                                auto glissileLoop(DN.loops().create(DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(), glidePlane));

                                std::vector<std::shared_ptr<LoopNodeType>> loopNodes;

                                const auto periodicGlidePlane(DN.periodicGlidePlaneFactory->get(glidePlane->key));
                                const auto periodicPatch(periodicGlidePlane->getPatch(VectorDim::Zero()));

                                if (isSink->isBoundaryNode())
                                {
                                    assert(!isSource->isBoundaryNode() && "Cross-slip cannot happen at the boundary");
                                    //Get the loopNodes of Sink
                                    std::set<short int> edgeIDs;
                                    for (const auto &edge : periodicPatch->edges())
                                    {
                                        if (((isSink->get_P() - edge->meshIntersection->P0).cross(isSink->get_P() - edge->meshIntersection->P1)).squaredNorm() < FLT_EPSILON)
                                        {
                                            edgeIDs.insert(edge->edgeID);
                                        }
                                    }
                                    assert(edgeIDs.size() == 1 && "Cross-Slip  at corner");
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSink, isSink->get_P(), periodicPatch, periodicPatch->edges()[*edgeIDs.begin()]));
                                }
                                else
                                {
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSink, isSink->get_P(), periodicPatch, nullptr));
                                }

                                if (isSource->isBoundaryNode())
                                {
                                    assert(!isSink->isBoundaryNode() && "Cross-slip cannot happen at the boundary");
                                    //Get the loopNodes of Sink
                                    std::set<short int> edgeIDs;
                                    for (const auto &edge : periodicPatch->edges())
                                    {
                                        if (((isSource->get_P() - edge->meshIntersection->P0).cross(isSource->get_P() - edge->meshIntersection->P1)).squaredNorm() < FLT_EPSILON)
                                        {
                                            edgeIDs.insert(edge->edgeID);
                                        }
                                    }
                                    assert(edgeIDs.size() == 1 && "Glissile Junction Intersection at corner");
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSource, isSource->get_P(), periodicPatch, periodicPatch->edges()[*edgeIDs.begin()]));
                                }
                                else
                                {
                                    loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, isSource, isSource->get_P(), periodicPatch, nullptr));
                                }

                                //New node cannot be a boundary node
                                loopNodes.emplace_back(DN.loopNodes().create(glissileLoop, newNode, newNode->get_P(), periodicPatch, nullptr));
                                //                                nodeIDs.push_back(sinkID);      // insert in reverse order, sink first, source second
                                //                                nodeIDs.push_back(sourceID);    // insert in reverse order, sink first, source second
                                //                                nodeIDs.push_back(newNodeID);

                                //                                LatticePlane loopPlane(newNodeP,DN.poly.grain(grainID).slipSystems()[slipID]->n);
                                //                                GlidePlaneKey<dim> loopPlaneKey(grainID,loopPlane);

                                //                                DN.insertLoop(nodeIDs,
                                //                                              DN.poly.grain(grainID).slipSystems()[slipID]->s.cartesian(),
                                //                                              DN.glidePlaneFactory.get(loopPlaneKey));

                                DN.insertLoop(glissileLoop, loopNodes);

                                executed++;
                            }
                        }
                    }
                }
            }
            model::cout<<executed<<" executed"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
        }
        

        
    };
    
//    template <typename DislocationNetworkType>
//    int DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip=0;
//
//    template <typename DislocationNetworkType>
//    double DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg=2.0;
    
}
#endif
