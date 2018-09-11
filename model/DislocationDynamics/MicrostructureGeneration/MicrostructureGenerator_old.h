/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_H_
#define model_MicrostructureGenerator_H_

#include <chrono>
#include <random>
#include <cmath>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>

#include <model/DislocationDynamics/Polycrystals/Polycrystal.h> // defines mode::cout
//#include <model/IO/EigenDataReader.h>
#include <model/Mesh/SimplicialMesh.h> // defines mode::cout
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/Mesh/PlaneMeshIntersection.h>
#include <model/DislocationDynamics/IO/DislocationNodeIO.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>
#include <model/DislocationDynamics/IO/DislocationEdgeIO.h>
#include <model/DislocationDynamics/IO/EVLio.h>
#include <model/DislocationDynamics/IO/DislocationLinkingNumber.h>


namespace model
{
    
    
    class MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<long int,dim,dim>	MatrixDimI;
        typedef Material<dim,Isotropic> MaterialType;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        
        /**********************************************************************/
        static double min(const double& a,const double& b)
        {
            return a<b? a : b;
        }
        
        /**********************************************************************/
        static double max(const double& a,const double& b)
        {
            return a>b? a : b;
        }
        
        /**********************************************************************/
        std::deque<VectorDimD> straightLineBoundaryClosure(const VectorDimD& P0,
                                                           const VectorDimD& d,
                                                           const VectorDimD& n,
                                                           const int& grainID)
        {
            // Define line AB containing dislocaiton and piercing the mesh
            const VectorDimD A=P0+3.0*maxSize*d;
            const VectorDimD B=P0-3.0*maxSize*d;
            
            // Compute interseciton between mesh and glide plane
            PlaneMeshIntersection<dim> pmi(mesh,P0,n,grainID);
            std::deque<std::pair<VectorDimD,VectorDimD>> segDeq;
            
            for(size_t k=0;k<pmi.size();++k)
            {
                const size_t k1=(k+1)<pmi.size()? k+1 :0;
                segDeq.emplace_back(pmi[k].second,pmi[k1].second);
            }
            
            std::deque<VectorDimD> nodePos;
            int nIntersections=0;
            for(const auto& pair : segDeq)
            {
                SegmentSegmentDistance<dim> ssi(A,B,pair.first,pair.second);
                
                if(nIntersections==0)
                {
                    if(ssi.dMin>FLT_EPSILON) // no intersection
                    {
                        
                    }
                    else //if(ssi.size==1)
                    {
                        nIntersections++;
                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
                    }
                }
                else if(nIntersections==1)
                {
                    const VectorDimD X((ssi.x0+ssi.x1)*0.5);
                    if((X-nodePos.back()).norm()>FLT_EPSILON)
                    {
                        if(ssi.dMin>FLT_EPSILON) // no intersection
                        {
                            nodePos.push_back(pair.first);
                        }
                        else //if(ssi.size==1)
                        {
                            nIntersections++;
                            nodePos.push_back(pair.first);
                            nodePos.push_back(X);
                        }
                    }
                }
                else
                {
                    
                }
            }
            return nodePos;
        }
        
        /**********************************************************************/
        void addStraightDislocations()
        {
            if(targetStraightDislocationDensity>0.0)
            {
                // init counters
                double density=0.0;
                double sessileDensity=0.0;
                
                std::cout<<greenBoldColor<<"Generating straight dislocations"<<defaultColor<<std::endl;
                while(density<targetStraightDislocationDensity)
                {
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const int& grainID=rp.second;   // random grain ID
                    const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                    const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem ID
                    const SlipSystem& slipSystem=poly.grain(grainID).slipSystems()[rSS];
                    const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
                    
                    
                    VectorDimD n=slipSystem.unitNormal; // slip plane normal
                    std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
                    const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
                    VectorDimD d=Eigen::AngleAxisd(theta, n)*b.normalized();
                    
                    bool isSessile=false;
                    if(sessileDensity/targetStraightDislocationDensity<fractionSessileStraightDislocationDensity)
                    {
                        n=b.normalized();
                        isSessile=true;
                        d=Eigen::AngleAxisd(theta, n)*n.cross(VectorDimD::Random()).normalized();
                    }

                    
                    std::deque<VectorDimD> nodePos=straightLineBoundaryClosure(P0,d,n,grainID);

                    
                    const double lineLength=(nodePos[nodePos.size()-1]-nodePos[0]).norm();
                    //                nodePos.push_back(nodePos[nodePos.size()-1]+1.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                    //                nodePos.push_back(nodePos[nodePos.size()-1]+2.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                    
                    // Write files
                    if(nodePos.size()>=3)
                    {
                        // write node and edge file
                        for(size_t k=0;k<nodePos.size();++k)
                        {
                            nodesIO.emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            //                        DislocationNodeIO<dim> dlIO(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            //                        vertexFile<<dlIO<<"\n";
                            //                        vertexFile.write(dlIO);
                            //                        vertexFile << nodeID+k<<"\t" << std::setprecision(15)<<std::scientific<<nodePos[k].transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                            
                            const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                            edgesIO.emplace_back(loopID,nodeID+k,nextNodeID,0);
                            //                        edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nextNodeID<<"\n";
                            
                        }
                        nodeID+=nodePos.size();
                        
                        // write loop file
                        loopsIO.emplace_back(loopID+0, b,n,P0,grainID);
                        //                    DislocationLoopIO<dim> dlIO(loopID+0, b,n,P0,grainID);
                        //                    loopFile<< dlIO<<"\n";
                        
                        loopID+=1;
                        snID+=1;
                        density += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                        if(isSessile)
                        {
                            sessileDensity += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                        }
                        std::cout<<"theta="<<theta*180.0/M_PI<<", density="<<density<<" (sessileDensity="<<sessileDensity<<")"<<std::endl;
                    }
                }
            }
        }
        
        /**********************************************************************/
        void addFrankReadSources()
        {
            if(targetFrankReadDislocationDensity>0.0)
            {
                std::cout<<greenBoldColor<<"Generating Frank-Read sources"<<defaultColor<<std::endl;
                
                double density=0.0;
                double edgeDensity=0.0;
                
                const double fractionEdge=1.0; // TEMPORARY
                
                
                while(density<targetFrankReadDislocationDensity)
                {
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    
                    const int rSS=distribution(generator); // a random SlipSystem
                    
                    const auto& slipSystem=poly.grain(grainID).slipSystems()[rSS];
                    
                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));
                    
                    bool isEdge=true;
                    
                    
                    
                    LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)*randomSign()));
                    double d1cNorm(d1.cartesian().norm());
                    int a1=randomSize()/d1cNorm;
                    LatticeVector<dim> L1=L0+d1*a1;
                    
                    
                    if(edgeDensity>fractionEdge*density) // overwrite with screw dislocaiton
                    {
                        isEdge=false;
                        d1cNorm=slipSystem.s.cartesian().norm();
                        a1=randomSize()/d1cNorm;
                        L1=L0+slipSystem.s*a1;
                    }
                    
                    // Compute the LatticeDireciton corresponding to -n
                    LatticeDirection<3> d2(poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*randomSign()));
                    double d2cNorm(d2.cartesian().norm());
                    
                    const int a2=2*a1; // aspect ratio of double FR source
                    LatticeVector<dim> L2=L1+d2*a2;
                    LatticeVector<dim> L3=L0+d2*a2;
                    
                    
                    const auto search1(mesh.search(L1.cartesian()));
                    const auto search2(mesh.search(L2.cartesian()));
                    const auto search3(mesh.search(L3.cartesian()));
                    
                    if(   search1.first && search1.second->region->regionID==grainID
                       && search2.first && search2.second->region->regionID==grainID
                       && search3.first && search3.second->region->regionID==grainID)
                    {
                        density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                        if(isEdge)
                        {
                            edgeDensity+=2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                        }
                        std::cout<<"density="<<density<<" (edge density="<<edgeDensity<<")"<<std::endl;
                        
                        const VectorDimD P0=L0.cartesian();
                        const VectorDimD P1=L1.cartesian();
                        const VectorDimD P2=L2.cartesian();
                        const VectorDimD P3=L3.cartesian();
                        const VectorDimD P4=0.5*(P0+P1);
                        const VectorDimD P5=0.5*(P2+P3);
                        
                        const VectorDimD n1=slipSystem.unitNormal;
                        const VectorDimD n2=d2.cross(d1).cartesian().normalized();
                        
                        nodesIO.emplace_back(nodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        nodesIO.emplace_back(nodeID+1,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        nodesIO.emplace_back(nodeID+2,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        nodesIO.emplace_back(nodeID+3,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        nodesIO.emplace_back(nodeID+4,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        nodesIO.emplace_back(nodeID+5,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        
                        loopsIO.emplace_back(loopID+0,slipSystem.s.cartesian(),n1,P0,grainID);
                        loopsIO.emplace_back(loopID+1,slipSystem.s.cartesian(),n2,P0,grainID);
                        loopsIO.emplace_back(loopID+2,slipSystem.s.cartesian(),n1,P3,grainID);
                        
                        edgesIO.emplace_back(loopID+0,nodeID+0,nodeID+1,0);
                        edgesIO.emplace_back(loopID+0,nodeID+1,nodeID+4,0);
                        edgesIO.emplace_back(loopID+0,nodeID+4,nodeID+0,0);
                        
                        edgesIO.emplace_back(loopID+1,nodeID+0,nodeID+3,0);
                        edgesIO.emplace_back(loopID+1,nodeID+3,nodeID+2,0);
                        edgesIO.emplace_back(loopID+1,nodeID+2,nodeID+1,0);
                        edgesIO.emplace_back(loopID+1,nodeID+1,nodeID+0,0);
                        
                        edgesIO.emplace_back(loopID+2,nodeID+3,nodeID+5,0);
                        edgesIO.emplace_back(loopID+2,nodeID+5,nodeID+2,0);
                        edgesIO.emplace_back(loopID+2,nodeID+2,nodeID+3,0);
                        
                        nodeID+=6;
                        loopID+=3;
                        snID++;
                        
                    }
                }
                
            }
        }
        
        /**********************************************************************/
        void addSingleArmDislocations()
        {
            if(targetFrankReadDislocationDensity>0.0)
            {
                std::cout<<greenBoldColor<<"Generating single-arm sources"<<defaultColor<<std::endl;
                
                double fractionEdge=1.0; // TEMPORARY
                
                double density=0.0;
                double edgeDensity=0.0;
                
                while(density<targetFrankReadDislocationDensity)
                {
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    
                    const int rSS=distribution(generator); // a random SlipSystem
                    
                    const auto& slipSystem=poly.grain(grainID).slipSystems()[rSS];
                    
                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));
                    
                    //                    bool isEdge=true;
                    
                    LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)*randomSign()));  // this is the projection direction
                    //     if(edgeDensity>fractionEdge*density) // overwrite with screw dislocaiton
                    //      {
                    //           d1=slipSystem.s;
                    //      }
                    // Compute the LatticeDireciton corresponding to -n
                    LatticeDirection<3> d2(poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*randomSign()));
                    double d2cNorm(d2.cartesian().norm());
                    const int a2=randomSize()/d2cNorm;
                    LatticeVector<dim> L3=L0+d2*a2;
                    
                    
                    const auto search2(mesh.search(L3.cartesian()));
                    
                    if(  search2.first && search2.second->region->regionID==grainID)
                    {
                        const VectorDimD P0=L0.cartesian();
                        const VectorDimD P3=L3.cartesian();
                        
                        
                        const VectorDimD n1=slipSystem.unitNormal;
                        const VectorDimD n2=d2.cross(d1).cartesian().normalized();
                        
                        
                        PlaneMeshIntersectionContainerType pmi01=PlaneMeshIntersection<dim>(mesh,P0,n2,grainID);
                        const VectorDimD P1=boundaryProjection(P0,d1.cartesian(),pmi01).second;
                        const VectorDimD P2=boundaryProjection(P3,d1.cartesian(),pmi01).second;
                        const std::map<double,VectorDimD> P12=boundaryProjection(P0,P3,d1.cartesian(),pmi01);
                        
                        // const VectorDimD P1=P12.begin().second;
                        // const VectorDimD P2=P12.rbegin().second;
                        const VectorDimD P4=(P0+P1)/2.0;
                        const VectorDimD P5=(P3+P2)/2.0;
                        
                        if ((P1-P0).norm()>a2*0.5 && (P3-P2).norm()>a2*0.5)  //not too small arm
                        {
                            density+=((P1-P0).norm()+(P0-P3).norm()+(P3-P2).norm())/mesh.volume()/std::pow(poly.b_SI,2);
                            /*! Vertex file format is:
                             * ID Px Py Pz Vx Vy Vz velReducCoeff snID meshLocation grainID
                             */
                            const size_t refNodeID=nodeID;
                            nodesIO.emplace_back(refNodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            nodesIO.emplace_back(refNodeID+1,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            nodesIO.emplace_back(refNodeID+2,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            nodesIO.emplace_back(refNodeID+3,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            nodesIO.emplace_back(refNodeID+4,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            nodesIO.emplace_back(refNodeID+5,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            nodeID+=6;
                            
                            if (P12.size()==0)
                            {
                                edgesIO.emplace_back(loopID,refNodeID+4,refNodeID+5,0);
                            }
                            else if (P12.size()==1)
                            {
                                nodesIO.emplace_back(nodeID,P12.begin()->second,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                edgesIO.emplace_back(loopID,refNodeID+4,nodeID,0);
                                edgesIO.emplace_back(loopID,nodeID,refNodeID+5,0);
                                nodeID++;
                            }
                            else
                            {
                                for(const auto pair : P12)
                                {
                                    nodesIO.emplace_back(nodeID,pair.second,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                    nodeID++;
                                    if (pair.first==P12.begin()->first)
                                    {
                                        edgesIO.emplace_back(loopID,refNodeID+4,nodeID,0);
                                        edgesIO.emplace_back(loopID,nodeID,nodeID+1,0);
                                    }
                                    else if (pair.first==P12.rbegin()->first )
                                    {
                                        edgesIO.emplace_back(loopID,nodeID,refNodeID+5,0);
                                    }
                                    else
                                    {
                                        edgesIO.emplace_back(loopID,nodeID,nodeID+1,0);
                                    }
                                }
                            }
                            loopsIO.emplace_back(loopID,slipSystem.s.cartesian(),n2,P0,grainID);
                            edgesIO.emplace_back(loopID,refNodeID+5,refNodeID+1,0);
                            edgesIO.emplace_back(loopID,refNodeID+1,refNodeID+0,0);
                            edgesIO.emplace_back(loopID,refNodeID+0,refNodeID+4,0);
                            
                            loopsIO.emplace_back(loopID+1,slipSystem.s.cartesian(),n1,P0,grainID);
                            edgesIO.emplace_back(loopID+1,refNodeID+0,refNodeID+2,0);
                            edgesIO.emplace_back(loopID+1,refNodeID+2,refNodeID+4,0);
                            edgesIO.emplace_back(loopID+1,refNodeID+4,refNodeID+0,0);
                            
                            loopsIO.emplace_back(loopID+2,slipSystem.s.cartesian(),n1,P3,grainID);
                            edgesIO.emplace_back(loopID+2,refNodeID+1,refNodeID+5,0);
                            edgesIO.emplace_back(loopID+2,refNodeID+5,refNodeID+3,0);
                            edgesIO.emplace_back(loopID+2,refNodeID+3,refNodeID+1,0);
                            
                            loopID+=3;
                            snID++;
                            
                            std::cout<<"density="<<density<<std::endl;
                        }
                    }
                }
            }
        }
        
        /**********************************************************************/
        void addPrismaticLoops()
        {
            if(targetPrismaticLoopDensity>0.0)
            {
                std::cout<<greenBoldColor<<"Generating prismatic loops"<<defaultColor<<std::endl;
                
                double density=0.0;
                
                while(density<targetPrismaticLoopDensity)
                {
                    
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem
                    const auto& slipSystem=poly.grain(grainID).slipSystems()[rSS];
                    
                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));
                    
                    // find prismatic planes
                    std::set<int> normalIDs;
                    for(size_t k=0;k<poly.grain(grainID).planeNormals().size();++k)
                    {
                        if(slipSystem.s.dot(poly.grain(grainID).planeNormals()[k])==0)
                        {
                            normalIDs.insert(k);
                        }
                    }
                    if(normalIDs.size()<2)
                    {
                        std::cout<<"normalIDs.size()="<<normalIDs.size()<<std::endl;
                        std::cout<<"Cannot generate a prismatic loop with less than 2 planes. EXITING."<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    assert(0 && "FINISH DIPOLE GENERATOR");
                }
                
            }
        }
        
        /**********************************************************************/
        void addIndividualStraightDislocations()
        {
            if(straightDislocationsSlipSystemIDs.size())
            {
                std::cout<<greenBoldColor<<"Generating individual straight dislocations"<<defaultColor<<std::endl;
                if(straightDislocationsSlipSystemIDs.size()!=straightDislocationsAngleFromScrewOrientation.size())
                {
                    std::cout<<"straightDislocationsSlipSystemIDs.size()="<<straightDislocationsSlipSystemIDs.size()<<std::endl;
                    std::cout<<"straightDislocationsAngleFromScrewOrientation.size()="<<straightDislocationsAngleFromScrewOrientation.size()<<std::endl;
                    std::cout<<"You must provide one angle for each dislocation. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                if(int(straightDislocationsSlipSystemIDs.size())!=pointsAlongStraightDislocations.rows())
                {
                    std::cout<<"straightDislocationsSlipSystemIDs.size()="<<straightDislocationsSlipSystemIDs.size()<<std::endl;
                    std::cout<<"pointsAlongStraightDislocations.rows()="<<pointsAlongStraightDislocations.rows()<<std::endl;
                    std::cout<<"You must provide one point for each dislocation. Each point is a row of the matrix pointsAlongStraightDislocations. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                for(size_t k=0;k<straightDislocationsSlipSystemIDs.size();++k)
                {
                    
                    std::pair<bool,const Simplex<dim,dim>*> found=mesh.search(pointsAlongStraightDislocations.row(k));
                    if(!found.first)
                    {
                        std::cout<<"Point "<<pointsAlongStraightDislocations.row(k)<<" is outside mesh. EXITING."<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    int grainID=found.second->region->regionID;
                    
                    
                    const int& rSS(straightDislocationsSlipSystemIDs[k]);

                    
                    if(rSS>=0)
                    {
                        
                        if(rSS>=int(poly.grain(grainID).slipSystems().size()))
                        {
                            std::cout<<"requested slip system ID="<<rSS<<std::endl;
                            std::cout<<"# of slip systems ="<<poly.grain(grainID).slipSystems().size()<<std::endl;
                            std::cout<<"Requested slip system does not exist. EXITING."<<std::endl;
                            exit(EXIT_FAILURE);
                        }
                        
                        const auto& slipSystem=poly.grain(grainID).slipSystems()[rSS];
                        const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(slipSystem.n,pointsAlongStraightDislocations.row(k));
                        
                        const VectorDimD P0=pointsAlongStraightDislocations.row(k).transpose()-pointsAlongStraightDislocations.row(k).dot(slipSystem.unitNormal)*slipSystem.unitNormal+slipSystem.unitNormal*slipSystem.n.planeSpacing()*heightPair.second;
                        
                        const double theta(straightDislocationsAngleFromScrewOrientation[k]*M_PI/180.0);
                        const VectorDimD& n(slipSystem.unitNormal);
                        const VectorDimD b(slipSystem.s.cartesian());
                        
                        const VectorDimD d=Eigen::AngleAxisd(theta, n)*b.normalized();
                        const std::deque<VectorDimD> nodePos=straightLineBoundaryClosure(P0,d,n,grainID);

                        
                        const double lineLength=(nodePos[nodePos.size()-1]-nodePos[0]).norm();
                        
                        // Write files
                        if(nodePos.size()>=3)
                        {
                            // write node and edge file
                            for(size_t k=0;k<nodePos.size();++k)
                            {
                                nodesIO.emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                //                        DislocationNodeIO<dim> dlIO(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                //                        vertexFile<<dlIO<<"\n";
                                //                        vertexFile.write(dlIO);
                                //                        vertexFile << nodeID+k<<"\t" << std::setprecision(15)<<std::scientific<<nodePos[k].transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                                
                                const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                                edgesIO.emplace_back(loopID,nodeID+k,nextNodeID,0);
                                //                        edgeFile << loopID<<"\t" <<    nodeID+k<<"\t"<< nextNodeID<<"\n";
                                
                            }
                            nodeID+=nodePos.size();
                            
                            // write loop file
                            loopsIO.emplace_back(loopID+0, b,n,P0,grainID);
                            //                    DislocationLoopIO<dim> dlIO(loopID+0, b,n,P0,grainID);
                            //                    loopFile<< dlIO<<"\n";
                            
                            loopID+=1;
                            snID+=1;
                            std::cout<<"["<<b.transpose()<<"]("<<slipSystem.unitNormal.transpose()<<") dislocation. Line dir="<<d.transpose()<<". Length="<<lineLength<<std::endl;
                            if(lineLength<FLT_EPSILON)
                            {
                                std::cout<<"Line too short. EXITING."<<std::endl;
                                exit(EXIT_FAILURE);
                            }
                            //                        density += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                            //                        if(isSessile)
                            //                        {
                            //                            sessileDensity += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                            //                        }
                            //                        std::cout<<"theta="<<theta*180.0/M_PI<<", density="<<density<<" (sessileDensity="<<sessileDensity<<")"<<std::endl;
                        }
                    
                    }
                    else
                    {
                        std::cout<<"negative slip system ID. Skipping entries."<<std::endl;
                    }
                    

                    
                    
                }
                
                
            }
            
        }
        
        /**********************************************************************/
        void addEshelbyInclusions()
        {
            
            if(targetInclusionDensities.size())
            {
                
                assert(targetInclusionDensities.size()==inclusionsDistribution_alpha.size());
                assert(targetInclusionDensities.size()==inclusionsDistribution_beta.size());
                assert(targetInclusionDensities.size()==inclusionsDistribution_lambda.size());
                assert(int(targetInclusionDensities.size())==inclusionsTransformationStrains.rows());
                assert(int(targetInclusionDensities.size())==inclusionsPatterns.rows());
                
                
                std::cout<<greenBoldColor<<"Generating Inclusions"<<defaultColor<<std::endl;
                
                std::ofstream inclusionsfile("E/E_0.txt");
                std::deque<std::pair<double,VectorDimD>> existingPrecipitates;
                
                size_t inclusionID=0;
                for(size_t f=0;f<targetInclusionDensities.size();++f)
                {
                    if(   targetInclusionDensities[f]>0.0)
                    {
                        assert(inclusionsDistribution_alpha[f]>0.0);
                        assert(inclusionsDistribution_beta[f]>0.0);
                        assert(inclusionsDistribution_lambda[f]>0.0);
                        
                        const VectorDimD currentPattern(inclusionsPatterns.row(f)/poly.b_SI); // normalize to length units
                        const double patternHeigth(currentPattern.norm());
                        const bool applyPattern(patternHeigth>0.0);
                        const VectorDimD patternDir(applyPattern? (currentPattern/patternHeigth).eval() : VectorDimD::Zero());
                        
                        double numberDensity=0.0;
                        
                        while(numberDensity<targetInclusionDensities[f])
                        {
                            
                            //                        std::default_random_engine generator;
                            std::gamma_distribution<double> distribution(inclusionsDistribution_alpha[f],inclusionsDistribution_beta[f]);
                            
                            const double size = distribution(generator)*inclusionsDistribution_lambda[f]/poly.b_SI;
                            //poly.grain(grainID)
                            std::pair<LatticeVector<dim>,int> pointPair=randomPointInMesh();
                            VectorDimD P=pointPair.first.cartesian();
                            const int& grainID(pointPair.second);
                            
                            if(applyPattern)
                            {
                                const VectorDimD globalVector(poly.grain(grainID).get_C2G()*currentPattern);
                                const VectorDimD globalDir(poly.grain(grainID).get_C2G()*patternDir);
                                
                                
                                const long long pointHeigth=std::round(P.dot(globalDir)/patternHeigth);
                                const VectorDimD O(pointHeigth*globalVector);
                                P-=(P-O).dot(globalDir)*globalDir;
                                
                            }
                            
                            bool isGoodPopision=mesh.searchRegion(grainID,P).first;
                            for(const auto& pair : existingPrecipitates)
                            {
                                isGoodPopision *= (P-pair.second).norm()>pair.first+size;
                                if(!isGoodPopision)
                                {
                                    break;
                                }
                            }
                            
                            if(isGoodPopision)
                            {
                                inclusionsfile<<inclusionID
                                /*          */<<" "<<P.transpose()
                                /*          */<<" "<<size
                                /*          */<<" "<<inclusionsTransformationStrains.row(f)
                                /*          */<<" "<<f
                                /*          */<<"\n";
                                
                                numberDensity+=1.0/mesh.volume()/std::pow(poly.b_SI,3);
                                inclusionID++;
                                existingPrecipitates.emplace_back(size,P);
                                
                                std::cout<<"inclusions density="<<numberDensity<<std::endl;
                            }
                        }
                    }
                }
                inclusionsfile.close();
            }
        }
        
        std::mt19937 generator;
        size_t nodeID;
        size_t snID;
        size_t loopID;
        std::vector<DislocationNodeIO<dim>> nodesIO;
        std::vector<DislocationLoopIO<dim>> loopsIO;
        std::vector<DislocationEdgeIO<dim>> edgesIO;
        std::deque<std::deque<VectorDimD>> loopPoints;
        std::deque<VectorDimD> loopBurgers;
        
    public:
        
        const bool outputBinary;
        const int meshID;
        const SimplicialMesh<dim> mesh;
        const double minSize;
        const double maxSize;
        GlidePlaneObserver<dim> gpo;
        Polycrystal<dim> poly;
        
        // Straight Dislocations
        const double targetStraightDislocationDensity;
        const double fractionSessileStraightDislocationDensity;
        
        // Frank-Read sources
        const double targetFrankReadDislocationDensity;
        
        // Single-arm sources
        const double targetSingleArmDislocationDensity;
        
        // Prismatic loops
        const double targetPrismaticLoopDensity;
        
        // Individual dislocations
        const std::vector<int> straightDislocationsSlipSystemIDs;
        const std::vector<double> straightDislocationsAngleFromScrewOrientation;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> pointsAlongStraightDislocations;
        
        
        // Defects
        const double targetIrradiationLoopDensity;
        const double averageLoopSize;
        const std::vector<double> targetInclusionDensities;
        const std::vector<double> inclusionsDistribution_alpha;
        const std::vector<double> inclusionsDistribution_beta;
        const std::vector<double> inclusionsDistribution_lambda;
        const Eigen::Matrix<double,Eigen::Dynamic,dim*dim> inclusionsTransformationStrains;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> inclusionsPatterns;
        
        /**********************************************************************/
        MicrostructureGenerator(int argc, char* argv[]) :
        /* init*/ generator(std::chrono::system_clock::now().time_since_epoch().count())
        /* init*/,nodeID(0)
        /* init*/,snID(0)
        /* init*/,loopID(0)
        /* init*/,outputBinary(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputBinary",true))
        /* init*/,meshID(TextFileParser("./inputFiles/DD.txt").readScalar<int>("meshID",true))
        /* init*/,mesh(meshID)
        /* init*/,minSize(0.1*min(mesh.xMax(0)-mesh.xMin(0),min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
        /* init*/,maxSize(max(mesh.xMax(0)-mesh.xMin(0),max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
        /* init*/,poly("./inputFiles/polycrystal.txt",mesh,gpo)
        /* Straight Dislocations */
        /* init*/,targetStraightDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetStraightDislocationDensity",true))
        /* init*/,fractionSessileStraightDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("fractionSessileStraightDislocationDensity",true))
        /* Frank-Read sources */
        /* init*/,targetFrankReadDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetFrankReadDislocationDensity",true))
        /* Single-arm sources */
        /* init*/,targetSingleArmDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetSingleArmDislocationDensity",true))
        /* Prismatic loops */
        /* init*/,targetPrismaticLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetPrismaticLoopDensity",true))
        /* Indivial straight dislocations */
        /* init*/,straightDislocationsSlipSystemIDs(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<int>("straightDislocationsSlipSystemIDs",true))
        /* init*/,straightDislocationsAngleFromScrewOrientation(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("straightDislocationsAngleFromScrewOrientation",true))
        /* init*/,pointsAlongStraightDislocations(TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("pointsAlongStraightDislocations",straightDislocationsSlipSystemIDs.size(),dim,true))
        /* Irradiation Loops */
        /* init*/,targetIrradiationLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetIrradiationLoopDensity",true))
        /* init*/,averageLoopSize(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("averageLoopSize",true))
        /* Inclusions */
        /* init*/,targetInclusionDensities(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("targetInclusionDensities",true))
        /* init*/,inclusionsDistribution_alpha(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDistribution_alpha",true))
        /* init*/,inclusionsDistribution_beta(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDistribution_beta",true))
        /* init*/,inclusionsDistribution_lambda(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDistribution_lambda",true))
        /* init*/,inclusionsTransformationStrains(TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsTransformationStrains",targetInclusionDensities.size(),dim*dim,true))
        /* init*/,inclusionsPatterns(TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsPatterns",targetInclusionDensities.size(),dim,true))
        {
            
            // Some sanity checks
            if(mesh.volume()==0.0)
            {
                std::cout<<"mesh "<<meshID<<" is empty. MicrostructureGenerator cannot run. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            // Call individual generators
            addStraightDislocations();
            addFrankReadSources();
            addSingleArmDislocations();
            addPrismaticLoops();
            addIndividualStraightDislocations();
            addIrradiationLoops();
            addEshelbyInclusions();
            
            // Output to evl/evl_0
            if(outputBinary)
            {
                EVLio<dim>::writeBin(0,nodesIO,loopsIO,edgesIO);
            }
            else
            {
                EVLio<dim>::writeTxt(0,nodesIO,loopsIO,edgesIO);
            }
            
        }
        
        /**********************************************************************/
        void addIrradiationLoops()
        {
            if(targetIrradiationLoopDensity>0.0)
            {
                
                size_t ndefects=0;
                double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                while(defectsDensity<targetIrradiationLoopDensity)
                {
                    
                    
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const int& grainID=rp.second;   // random grain ID
                    const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                    const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem ID
                    const SlipSystem& slipSystem=poly.grain(grainID).slipSystems()[rSS];
                    const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
                    const VectorDimD a(b.normalized());
                    //const VectorDimD c(slipSystem.s.cartesian());
                    
                    std::vector<VectorDimD> points;
                    
                    const int NP=6;
                    for(int k=0;k<NP;++k)
                    {
                        points.push_back(P0+Eigen::AngleAxis<double>(k*2.0*M_PI/NP,a)*slipSystem.unitNormal*0.5*averageLoopSize/poly.b_SI);
                    }
                    
                    bool pointsIncluded=true;
                    for(const auto& point : points)
                    {
                        pointsIncluded*=mesh.searchRegion(grainID,point).first;
                    }
                    
                    if(pointsIncluded)
                    {
                        
                        for(int k=0;k<NP;++k)
                        {
                            nodesIO.emplace_back(nodeID+k,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            
                            const int nextNodeID=(k+1)<NP? nodeID+k+1 : nodeID;
                            edgesIO.emplace_back(loopID,nodeID+k,nextNodeID,0);
                            
                        }
                        
                        loopsIO.emplace_back(loopID+0, b,a,P0,grainID);
                        
                        
                        snID++;
                        loopID++;
                        nodeID+=NP;
                        ndefects++;
                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"irradiation defects density="<<defectsDensity<<std::endl;
                    }
                }
            }
        }
        
        /**********************************************************************/
        double deltaHelicity(const std::deque<VectorDimD>& newPoints,
                             const VectorDimD& newBurgers) const
        {
            
            double h(0.0);
            assert(loopPoints.size()==loopBurgers.size());
            for(size_t k=0;k<loopPoints.size(); ++k)
            {
                h+=LinkingNumber<dim>::loopPairHelicity(loopPoints[k],loopBurgers[k],newPoints,newBurgers);
            }
            
            return h;
        }
        
        /**********************************************************************/
        std::pair<LatticeVector<dim>,int> randomPointInMesh() const
        {
            return poly.randomLatticePointInMesh();
        }
        
        /**********************************************************************/
        double randomSize()
        {
            std::uniform_real_distribution<double> dist(minSize,maxSize);
            return dist(generator);
        }
        
        /**********************************************************************/
        int randomSign()
        {
            std::uniform_int_distribution<> dis(0,1);
            return  dis(generator)*2-1;
        }
        
        /**********************************************************************/
        static std::map<double,VectorDimD> boundaryProjection(const VectorDimD& P0,
                                                              const VectorDimD& P1,
                                                              const VectorDimD& D,
                                                              const PlaneMeshIntersectionContainerType& pp)
        {
            
            
            
            const double dNorm(D.norm());
            assert(dNorm>FLT_EPSILON);
            const VectorDimD dir=D/dNorm;
            // Let a point v on the boundary be written as v=P0+u1*(P1-P0)+u2*d
            // then we have [P1-P0 d]*[u1 u2]^T=v-P0
            
            Eigen::Matrix<double,3,2> A;
            A.col(0)=P1-P0;
            A.col(1)=dir;
            const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(A.transpose()*A);
            assert(llt.info()==Eigen::Success);
            
            std::map<double,VectorDimD> temp; // keep points sorted by parameter u1
            for(size_t m=0;m<pp.size();++m)
            {
                const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*(pp[m].second-P0));
                if(x(0)>FLT_EPSILON && x(0)<1.0-FLT_EPSILON && x(1)>FLT_EPSILON)
                {
                    temp.emplace(x(0),pp[m].second);
                }
            }
            
            return temp;
            
            
        }
        
        /**********************************************************************/
        static std::pair<int,VectorDimD> boundaryProjection(const VectorDimD& P,
                                                            const VectorDimD& D,
                                                            const PlaneMeshIntersectionContainerType& pp)
        {
            const double dNorm(D.norm());
            assert(dNorm>FLT_EPSILON);
            const VectorDimD dir=D/dNorm;
            // line1 is P+u1*dir, u>0
            
            bool success=false;
            std::pair<int,VectorDimD> temp=std::make_pair(-1,VectorDimD::Zero());
            
            for(size_t k=0;k<pp.size();++k)
            {
                const size_t k1 = ((k==(pp.size()-1))? 0 : k+1);
                const VectorDimD& v0=pp[k].second;
                const VectorDimD& v1=pp[k1].second;
                // line2 is v0+u2*(v1-v0), 0<=u2<=1
                
                //P+u1*dir=v0+u2*(v1-v0)
                // [dir -(v1-v0)] [u1 u2] = [v0-P]
                // In least square sense
                // [dir -(v1-v0)]^T*[dir -(v1-v0)]* [u1 u2] = [v0-P]
                
                Eigen::Matrix<double,3,2> A;
                A.col(0)=dir;
                A.col(1)=-(v1-v0);
                
                const Eigen::Matrix<double,3,1> b=v0-P;
                
                const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(A.transpose()*A);
                //                std::cout<<"DO NOT USE LLT TO SEE IF SYSTEM HAS SOLUTION. See https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html#a858dc77b65dd48248299bb6a6a758abf"<<std::endl;
                
                
                if(llt.info()==Eigen::Success)
                {
                    
                    const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*b);
                    
                    if(x(0)>=0.0 && x(1)>=0.0 && x(1)<=1.0)
                    {
                        success=true;
                        temp=std::make_pair(k,v0+x(1)*(v1-v0));
                        break;
                    }
                }
                
            }
            
            assert(success);
            return temp;
        }
        
        
        
    };
    
}
#endif