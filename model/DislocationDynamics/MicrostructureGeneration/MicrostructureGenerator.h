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
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <DislocatedMaterial.h>
#include <LatticeMath.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationEdgeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>

//#include <PeriodicGlidePlane.h>


namespace model
{
    
    
    class MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef LatticeDirection<dim> LatticeDirectionType;
        
        typedef Eigen::Matrix<double,dim,dim>    MatrixDimD;
        typedef Eigen::Matrix<long int,dim,dim>    MatrixDimI;
        typedef DislocatedMaterial<dim,Isotropic> MaterialType;
        
        typedef BoundingMeshSegments<dim> MeshBoundaryContainerType;
        
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
        static VectorDimD randomOrthogonalUnitVector(VectorDimD v)
        {
            const double vNorm(v.norm());
            assert(vNorm>FLT_EPSILON);
            v/=vNorm;
            
            VectorDimD temp(v.cross(VectorDimD::Random()));
            double tempNorm(temp.norm());
            while(tempNorm<FLT_EPSILON)
            {
                temp=v.cross(VectorDimD::Random());
                tempNorm=temp.norm();
            }
            return temp/tempNorm;
        }
        
        /**********************************************************************/
        bool addSingleLoop(const bool randomizeBurgersSense,
                           const std::vector<VectorDimD>& nodePos,
                           VectorDimD b,
                           const VectorDimD& unitNormal,
                           const VectorDimD& P0,
                           const int& grainID,
                           const int& loopType)
        {
            if(allPointsInGrain(nodePos,grainID))
            {
                if(randomizeBurgersSense)
                {
                    b*=randomSign();
                }
                
                for(size_t k=0;k<nodePos.size();++k)
                {
                    const size_t nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                    configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                }
//                configIO.loops().emplace_back(loopID, b,unitNormal,P0,grainID,loopType,periodicLoopID,periodicShift);
                configIO.loops().emplace_back(loopID, b,unitNormal,P0,grainID,loopType);
                nodeID+=nodePos.size();
                loopID+=1;
                snID+=1;
                return true;
            }
            else
            {
                return false;
            }
        }
        
        /**********************************************************************/
        void addJoggedDislocations()
        {
            if(targetJoggedtDislocationDensity>0.0)
            {
                if(poly.crystalStructure=="HEX")
                {
                    
                    const double jogLengthM(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("jogLength_M",true) );
                    const double jogLengthS(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("jogLength_S",true) );

                    const double jogHeightM(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("jogHeight_M",true) );
                    const double jogHeightS(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("jogHeight_S",true) );

                    const double fractionSessile(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("jogFractionSessile",true));
                    
                    std::lognormal_distribution<double> jogLengthDistribution(jogLengthM,jogLengthS);
                    std::lognormal_distribution<double> jogHeightDistribution(jogHeightM,jogHeightS);

                    
                    std::cout<<magentaBoldColor<<"Generating jogged dislocations"<<defaultColor<<std::endl;
                    double density=0.0;
                    double sessileDensity(0.0);
                    while(density<targetJoggedtDislocationDensity)
                    {
                        const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                        const int& grainID=rp.second;   // random grain ID
                        std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                        const int rSS=distribution(generator); // a random SlipSystem ID
                        const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                        
                        if(fabs(slipSystem.n.planeSpacing()-sqrt(3.0)/2.0)<FLT_EPSILON)
                        {// prismatic planes
                            const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                            const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                            const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
                            const VectorDimD n=slipSystem.unitNormal; // slip plane normal
                            std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
                            //                    const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
                            const double theta(0.0); // screw dislocations
                            const VectorDimD d=Eigen::AngleAxisd(theta, n)*b.normalized();
                            
                            MeshPlane<dim> plane(mesh,grainID,P0,n);
                            std::vector<VectorDimD> nodePos(DislocationInjectorBase<dim>::straightLineBoundaryClosure(P0,d,plane,mesh));
                            
                            const auto FS(nodePos.front()-nodePos.back()); // vector goinf from end of loop to start
                            const double lineLength(FS.norm());
                            const double jogLength_SI(jogLengthDistribution(generator));
                            const double jogLength(jogLength_SI/poly.b_SI);

                            const auto midpoint(0.5*(nodePos.front()+nodePos.back())); // vector going from end of loop to start
                            const VectorDimD jogPt0(midpoint-0.5*jogLength/lineLength*FS);
                            const VectorDimD jogPt1(midpoint+0.5*jogLength/lineLength*FS);
                            nodePos.push_back(jogPt0);
                            nodePos.push_back(jogPt1);
                            
                            const double jogHeight_SI(jogHeightDistribution(generator));
                            const double jogHeight(std::round(jogHeight_SI/poly.b_SI/slipSystem.n.planeSpacing())*slipSystem.n.planeSpacing());

                            if(lineLength>jogLength)
                            {
                                std::vector<VectorDimD> jogPos;
//                                jogPos.push_back(jogPt1);
//                                jogPos.push_back(jogPt0);
                                jogPos.push_back(jogPt0+jogHeight*n);
                                jogPos.push_back(jogPt1+jogHeight*n);
                                
                                if(allPointsInGrain(nodePos,grainID) && allPointsInGrain(jogPos,grainID))
                                {
                                 
                                    if(addSingleLoop(false,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP))
                                    {
                                        configIO.nodes().emplace_back(nodeID+0,jogPt1+jogHeight*n,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                        configIO.nodes().emplace_back(nodeID+1,jogPt0+jogHeight*n,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                        configIO.nodes().emplace_back(nodeID+2,0.5*(jogPt0+jogPt1)+jogHeight*n,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                        
                                        configIO.links().emplace_back(loopID+0,nodeID-1,nodeID-2,0);
                                        configIO.links().emplace_back(loopID+0,nodeID-2,nodeID+1,0);
                                        configIO.links().emplace_back(loopID+0,nodeID+1,nodeID+0,0);
                                        configIO.links().emplace_back(loopID+0,nodeID+0,nodeID-1,0);
                                        
                                        configIO.links().emplace_back(loopID+1,nodeID+0,nodeID+1,0);
                                        configIO.links().emplace_back(loopID+1,nodeID+2,nodeID+0,0);
                                        configIO.links().emplace_back(loopID+1,nodeID+1,nodeID+2,0);

                                        

                                        
                                        if(sessileDensity>fractionSessile*density)
                                        {
                                            configIO.loops().emplace_back(loopID, b,slipSystem.unitNormal.cross(b).normalized(),jogPt0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                            density += (lineLength-jogLength)/mesh.volume()/std::pow(poly.b_SI,2);
                                            std::cout<<"density="<<density<<std::endl;
                                        }
                                        else
                                        {
                                            configIO.loops().emplace_back(loopID, b,slipSystem.unitNormal.cross(b).normalized(),jogPt0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                                            density += (lineLength-jogLength)/mesh.volume()/std::pow(poly.b_SI,2);
                                            sessileDensity += (lineLength-jogLength)/mesh.volume()/std::pow(poly.b_SI,2);
                                            std::cout<<"density="<<density<<std::endl;
                                            std::cout<<"sessileDensity="<<sessileDensity<<std::endl;
                                        }
                                        
                                        configIO.loops().emplace_back(loopID+1, b,slipSystem.unitNormal,jogPt0+jogHeight*n,grainID,DislocationLoopIO<dim>::GLISSILELOOP);

                                        
                                        nodeID+=3;
                                        loopID+=2;

                                        
                                        
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    std::cout<<magentaBoldColor<<"JoggedDislocations implemented on HEX only"<<defaultColor<<std::endl;
                }

            }
        }
        
        
        /**********************************************************************/
        void addStraightDislocations()
        {
            if(targetStraightDislocationDensity>0.0)
            {
                // init counters
                double density=0.0;
                double sessileDensity=0.0;
                
                std::cout<<magentaBoldColor<<"Generating straight dislocations"<<defaultColor<<std::endl;
                while(density<targetStraightDislocationDensity)
                {
                    
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const int& grainID=rp.second;   // random grain ID
                    const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                    const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem ID
                    const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
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
                    
                    MeshPlane<dim> plane(mesh,grainID,P0,n);
                    const std::vector<VectorDimD> nodePos(DislocationInjectorBase<dim>::straightLineBoundaryClosure(P0,d,plane,mesh));
                    
                    const double lineLength=(nodePos.back()-nodePos.front()).norm();
                    //                nodePos.push_back(nodePos[nodePos.size()-1]+1.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                    //                nodePos.push_back(nodePos[nodePos.size()-1]+2.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                    
                    double dh=0.0;
                    if(enforceMonotonicHelicity)
                    {
                        dh=deltaHelicity(nodePos,b);
                    }
                    
                    // Write files
                    if(   nodePos.size()>=3
                       && ((fabs(helicity+dh)>fabs(helicity) && fabs(dh)>FLT_EPSILON ) || helicity==0.0 || !enforceMonotonicHelicity))
                    {
                        // write node and edge file
                        for(size_t k=0;k<nodePos.size();++k)
                        {
                            configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                            configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                        }
                        nodeID+=nodePos.size();
                        
                        // write loop file
                        configIO.loops().emplace_back(loopID+0, b,n,P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                        
                        if(enforceMonotonicHelicity)
                        {
                            loopPoints.push_back(nodePos);
                            loopBurgers.push_back(b);
                        }
                        
                        loopID+=1;
                        snID+=1;
                        density += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                        if(isSessile)
                        {
                            sessileDensity += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                        }
                        std::cout<<"theta="<<theta*180.0/M_PI<<", density="<<density<<" (sessileDensity="<<sessileDensity<<")"<<std::endl;
                        if(enforceMonotonicHelicity)
                        {
                            loopPoints.push_back(nodePos);
                            loopBurgers.push_back(b);
                            helicity+=dh;
                            std::cout<<"helicity="<<helicity<<std::endl;
                        }
                        
                    }
                }
            }
        }
        
        /**********************************************************************/
        void addFrankReadSources()
        {
            if(targetFrankReadDislocationDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Frank-Read sources"<<defaultColor<<std::endl;
                
                double density=0.0;
                double edgeDensity=0.0;
                
                const double fractionEdge=1.0; // TEMPORARY
                
                std::normal_distribution<double> sizeDistribution(FrankReadSizeMean/poly.b_SI,FrankReadSizeStd/poly.b_SI);
                std::normal_distribution<double> aspectRatioDistribution(FrankReadAspectRatioMean,FrankReadAspectRatioStd);
                
                while(density<targetFrankReadDislocationDensity)
                {
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    
                    const int rSS=distribution(generator); // a random SlipSystem
                    
                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(slipSystem.s.cartesian());
                    
                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(b));
                    
                    bool isEdge=true;
                    
                    
                    
                    LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(slipSystem.n)*randomSign()));
                    double d1cNorm(d1.cartesian().norm());
                    //                    const double size = distribution(generator)*inclusionsDiameterLognormalDistribution_A[f]/poly.b_SI;
                    
                    int a1=sizeDistribution(generator)/d1cNorm;
                    if(a1>0)
                    {
                        LatticeVector<dim> L1=L0+d1*a1;
                        
                        
                        if(edgeDensity>fractionEdge*density) // overwrite with screw dislocaiton
                        {
                            isEdge=false;
                            d1cNorm=b.norm();
                            a1=randomSize()/d1cNorm;
                            L1=L0+slipSystem.s.dir*a1;
                        }
                        
                        // Compute the LatticeDireciton corresponding to -n
                        LatticeDirection<3> d2(poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*randomSign()));
                        double d2cNorm(d2.cartesian().norm());
                        
                        const int a2=aspectRatioDistribution(generator)*a1; // aspect ratio of double FR source
                        if(a2>0)
                        {
                            LatticeVector<dim> L2=L1+d2*a2;
                            LatticeVector<dim> L3=L0+d2*a2;
                            
                            const VectorDimD P0=L0.cartesian();
                            const VectorDimD P1=L1.cartesian();
                            const VectorDimD P2=L2.cartesian();
                            const VectorDimD P3=L3.cartesian();
                            const auto search1(mesh.search(P1));
                            const auto search2(mesh.search(P2));
                            const auto search3(mesh.search(P3));
                            
                            double dh=0.0;
                            std::vector<VectorDimD> nodePos;
                            nodePos.push_back(P0);
                            nodePos.push_back(P1);
                            nodePos.push_back(P2);
                            nodePos.push_back(P3);
                            if(enforceMonotonicHelicity)
                            {
                                dh=deltaHelicity(nodePos,-b); // central loop is opposite direction
                            }
                            
                            if(allPointsInGrain(nodePos,grainID)
                               && ((fabs(helicity+dh)>fabs(helicity) && fabs(dh)>FLT_EPSILON ) || helicity==0.0 || !enforceMonotonicHelicity)
                               )
                            {
                                density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                                if(isEdge)
                                {
                                    edgeDensity+=2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                                }
                                std::cout<<"density="<<density<<" (edge density="<<edgeDensity<<")"<<std::endl;
                                
                                const VectorDimD P4=0.5*(P0+P1);
                                const VectorDimD P5=0.5*(P2+P3);
                                
                                const VectorDimD n1=slipSystem.unitNormal;
                                const VectorDimD n2=d2.cross(d1).cartesian().normalized();
                                
                                configIO.nodes().emplace_back(nodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+1,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+2,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+3,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+4,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.nodes().emplace_back(nodeID+5,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                
                                configIO.loops().emplace_back(loopID+0,b,n1,P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                configIO.loops().emplace_back(loopID+1,b,n2,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                                configIO.loops().emplace_back(loopID+2,b,n1,P3,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                
                                configIO.links().emplace_back(loopID+0,nodeID+0,nodeID+1,0);
                                configIO.links().emplace_back(loopID+0,nodeID+1,nodeID+4,0);
                                configIO.links().emplace_back(loopID+0,nodeID+4,nodeID+0,0);
                                
                                configIO.links().emplace_back(loopID+1,nodeID+0,nodeID+3,0);
                                configIO.links().emplace_back(loopID+1,nodeID+3,nodeID+2,0);
                                configIO.links().emplace_back(loopID+1,nodeID+2,nodeID+1,0);
                                configIO.links().emplace_back(loopID+1,nodeID+1,nodeID+0,0);
                                
                                configIO.links().emplace_back(loopID+2,nodeID+3,nodeID+5,0);
                                configIO.links().emplace_back(loopID+2,nodeID+5,nodeID+2,0);
                                configIO.links().emplace_back(loopID+2,nodeID+2,nodeID+3,0);
                                
                                nodeID+=6;
                                loopID+=3;
                                snID++;
                                
                                if(enforceMonotonicHelicity)
                                {
                                    loopPoints.push_back(nodePos);
                                    loopBurgers.push_back(-b);
                                    helicity+=dh;
                                    std::cout<<"helicity="<<helicity<<std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /**********************************************************************/
        void addSingleArmDislocations()
        {
            if(targetSingleArmDislocationDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating single-arm sources"<<defaultColor<<std::endl;
                
                //                double fractionEdge=1.0; // TEMPORARY
                
                double density=0.0;
                //                double edgeDensity=0.0;
                
                while(density<targetSingleArmDislocationDensity)
                {
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    
                    const int rSS=distribution(generator); // a random SlipSystem
                    
                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    
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
                    
                    if(enforceMonotonicHelicity)
                    {
                        assert(false && "enforceMonotonicHelicity NOT YET IMPLEMENTED FOR SINGLE-ARM SOURCES");
                    }
                    
                    if(  search2.first && int(search2.second->region->regionID)==grainID)
                    {
                        const VectorDimD P0=L0.cartesian();
                        const VectorDimD P3=L3.cartesian();
                        
                        
                        const VectorDimD n1=slipSystem.unitNormal;
                        const VectorDimD n2=d2.cross(d1).cartesian().normalized();
                        
                        MeshPlane<dim> plane01(mesh,grainID,P0,n2);
                        //                        PlaneMeshIntersectionContainerType pmi01=PlaneMeshIntersection<dim>(mesh,P0,n2,grainID);
                        const VectorDimD P1=boundaryProjection(P0,d1.cartesian(),plane01.meshIntersections).second;
                        const VectorDimD P2=boundaryProjection(P3,d1.cartesian(),plane01.meshIntersections).second;
                        const std::map<double,VectorDimD> P12=boundaryProjection(P0,P3,d1.cartesian(),plane01.meshIntersections);
                        
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
                            configIO.nodes().emplace_back(refNodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            configIO.nodes().emplace_back(refNodeID+1,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            configIO.nodes().emplace_back(refNodeID+2,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            configIO.nodes().emplace_back(refNodeID+3,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            configIO.nodes().emplace_back(refNodeID+4,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            configIO.nodes().emplace_back(refNodeID+5,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            nodeID+=6;
                            
                            if (P12.size()==0)
                            {
                                configIO.links().emplace_back(loopID,refNodeID+4,refNodeID+5,0);
                            }
                            else if (P12.size()==1)
                            {
                                configIO.nodes().emplace_back(nodeID,P12.begin()->second,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                configIO.links().emplace_back(loopID,refNodeID+4,nodeID,0);
                                configIO.links().emplace_back(loopID,nodeID,refNodeID+5,0);
                                nodeID++;
                            }
                            else
                            {
                                for(const auto pair : P12)
                                {
                                    configIO.nodes().emplace_back(nodeID,pair.second,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                    nodeID++;
                                    if (pair.first==P12.begin()->first)
                                    {
                                        configIO.links().emplace_back(loopID,refNodeID+4,nodeID,0);
                                        configIO.links().emplace_back(loopID,nodeID,nodeID+1,0);
                                    }
                                    else if (pair.first==P12.rbegin()->first )
                                    {
                                        configIO.links().emplace_back(loopID,nodeID,refNodeID+5,0);
                                    }
                                    else
                                    {
                                        configIO.links().emplace_back(loopID,nodeID,nodeID+1,0);
                                    }
                                }
                            }
                            configIO.loops().emplace_back(loopID,slipSystem.s.cartesian(),n2,P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                            configIO.links().emplace_back(loopID,refNodeID+5,refNodeID+1,0);
                            configIO.links().emplace_back(loopID,refNodeID+1,refNodeID+0,0);
                            configIO.links().emplace_back(loopID,refNodeID+0,refNodeID+4,0);
                            
                            configIO.loops().emplace_back(loopID+1,slipSystem.s.cartesian(),n1,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                            configIO.links().emplace_back(loopID+1,refNodeID+0,refNodeID+2,0);
                            configIO.links().emplace_back(loopID+1,refNodeID+2,refNodeID+4,0);
                            configIO.links().emplace_back(loopID+1,refNodeID+4,refNodeID+0,0);
                            
                            configIO.loops().emplace_back(loopID+2,slipSystem.s.cartesian(),n1,P3,grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                            configIO.links().emplace_back(loopID+2,refNodeID+1,refNodeID+5,0);
                            configIO.links().emplace_back(loopID+2,refNodeID+5,refNodeID+3,0);
                            configIO.links().emplace_back(loopID+2,refNodeID+3,refNodeID+1,0);
                            
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
                std::cout<<magentaBoldColor<<"Generating prismatic loops"<<defaultColor<<std::endl;
                
                double density=0.0;
                
                while(density<targetPrismaticLoopDensity)
                {
                    
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem
                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(slipSystem.s.cartesian());
                    
                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(b));
                    
                    // find prismatic planes
                    std::vector<int> normalIDs;
                    for(size_t k=0;k<poly.grain(grainID).planeNormals().size();++k)
                    {
                        if(slipSystem.s.dot(poly.grain(grainID).planeNormals()[k])==0)
                        {
                            normalIDs.push_back(k);
                        }
                    }
                    if(normalIDs.size()<2)
                    {
                        std::cout<<"normalIDs.size()="<<normalIDs.size()<<std::endl;
                        std::cout<<"Cannot generate a prismatic loop with less than 2 planes. EXITING."<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                    
                    std::vector<LatticeDirection<dim>> dirVector;
                    std::vector<int> sizeVector;
                    for(const int& normalID : normalIDs)
                    {
                        dirVector.emplace_back(sr,poly.grain(grainID).planeNormals()[normalID]);
                        sizeVector.emplace_back(randomSize()*randomSign()/dirVector.back().cartesian().norm());
                        //                        sizeVector.emplace_back(10000);
                    }
                    
                    std::vector<VectorDimD> posVector;
                    std::vector<VectorDimD> normalsVector;
                    
                    posVector.push_back(L0.cartesian());
                    for(size_t k=0;k<dirVector.size();++k)
                    {
                        posVector.push_back(posVector[k]+dirVector[k].cartesian()*sizeVector[k]);
                        normalsVector.push_back(poly.grain(grainID).planeNormals()[normalIDs[k]].cartesian().normalized());
                    }
                    for(size_t k=0;k<dirVector.size();++k)
                    {
                        posVector.push_back(posVector[k+dirVector.size()]-dirVector[k].cartesian()*sizeVector[k]);
                        normalsVector.push_back(-poly.grain(grainID).planeNormals()[normalIDs[k]].cartesian().normalized());
                    }
                    posVector.pop_back();
                    //                    std::cout<<"posVector.size()="<<posVector.size()<<std::endl;
                    //                    for(const auto& pos : posVector)
                    //                    {
                    //                        std::cout<<pos.cartesian().transpose()<<std::endl;
                    //                    }
                    
                    bool allInside=true;
                    for(const auto& pos : posVector)
                    {
                        allInside*=mesh.searchRegion(grainID,pos).first;
                        if(!allInside)
                        {
                            break;
                        }
                    }
                    
                    double dh=0.0;
                    if(enforceMonotonicHelicity)
                    {
                        dh=deltaHelicity(posVector,-b);
                    }
                    
                    if(allInside
                       && ((fabs(helicity+dh)>fabs(helicity) && fabs(dh)>FLT_EPSILON ) || helicity==0.0 || !enforceMonotonicHelicity))
                    {
                        
                        // Add nodes (two for-loops are needed)
                        for(size_t k=0;k<posVector.size();++k)
                        {
                            configIO.nodes().emplace_back(nodeID+k,posVector[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        }
                        for(size_t k=0;k<posVector.size();++k)
                        {
                            configIO.nodes().emplace_back(nodeID+k+posVector.size(),posVector[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                        }
                        
                        
                        // Add lateral loops
                        for(size_t k=0;k<posVector.size();++k)
                        {
                            const int nextNodeID=(k+1)<posVector.size()? nodeID+k+1 : nodeID;
                            configIO.links().emplace_back(loopID+k,nodeID+k,nextNodeID,0);
                            configIO.links().emplace_back(loopID+k,nextNodeID,nextNodeID+posVector.size(),0);
                            configIO.links().emplace_back(loopID+k,nextNodeID+posVector.size(),nodeID+k+posVector.size(),0);
                            configIO.links().emplace_back(loopID+k,nodeID+k+posVector.size(),nodeID+k,0);
                            
                            configIO.loops().emplace_back(loopID+k,b,normalsVector[k],posVector[k],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                            
                        }
                        
                        
                        // Add back loop (sessile)
                        for(size_t k=0;k<posVector.size();++k)
                        {
                            const size_t nextNodeID=(k+1)<posVector.size()? nodeID+k+1 : nodeID;
                            configIO.links().emplace_back(loopID+posVector.size(),nodeID+k,nextNodeID,0);
                        }
                        configIO.loops().emplace_back(loopID+posVector.size(),-b,b.normalized(),posVector[0],grainID,DislocationLoopIO<dim>::SESSILELOOP);
                        
                        
                        
                        for(size_t k=0;k<dirVector.size();++k)
                        {
                            density+=2.0*(dirVector[k]*sizeVector[k]).cartesian().norm()/mesh.volume()/std::pow(poly.b_SI,2);
                        }
                        
                        std::cout<<"density="<<density<<std::endl;
                        
                        nodeID+=2*posVector.size();
                        snID+=1;
                        loopID+=posVector.size()+1;
                        
                        if(enforceMonotonicHelicity)
                        {
                            loopPoints.push_back(posVector);
                            loopBurgers.push_back(-b);
                            helicity+=dh;
                            std::cout<<"helicity="<<helicity<<std::endl;
                        }
                        
                    }
                }
                
            }
        }
        
        /**********************************************************************/
        void addIndividualStraightDislocations()
        {
            if(straightDislocationsSlipSystemIDs.size())
            {
                std::cout<<magentaBoldColor<<"Generating individual straight dislocations"<<defaultColor<<std::endl;
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
                    const int& rSS(straightDislocationsSlipSystemIDs[k]);
                    
                    
                    if(rSS>=0)
                    {
                        std::pair<bool,const Simplex<dim,dim>*> found=mesh.search(pointsAlongStraightDislocations.row(k));
                        if(!found.first)
                        {
                            std::cout<<"Point "<<pointsAlongStraightDislocations.row(k)<<" is outside mesh. EXITING."<<std::endl;
                            exit(EXIT_FAILURE);
                        }
                        
                        int grainID=found.second->region->regionID;
                        
                        
                        
                        std::cout<<"generating individual straight dislocation "<<k<<defaultColor<<std::endl;
                        
                        
                        if(rSS>=int(poly.grain(grainID).slipSystems().size()))
                        {
                            std::cout<<"requested slip system ID="<<rSS<<std::endl;
                            std::cout<<"# of slip systems ="<<poly.grain(grainID).slipSystems().size()<<std::endl;
                            std::cout<<"Requested slip system does not exist. EXITING."<<std::endl;
                            exit(EXIT_FAILURE);
                        }
                        
                        const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                        const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(slipSystem.n,pointsAlongStraightDislocations.row(k));
                        
                        const VectorDimD P0=pointsAlongStraightDislocations.row(k).transpose()-pointsAlongStraightDislocations.row(k).dot(slipSystem.unitNormal)*slipSystem.unitNormal+slipSystem.unitNormal*slipSystem.n.planeSpacing()*heightPair.second;
                        
                        const double theta(straightDislocationsAngleFromScrewOrientation[k]*M_PI/180.0);
                        //const VectorDimD& n(slipSystem.unitNormal);
                        const VectorDimD b(slipSystem.s.cartesian());
                        
                        const VectorDimD d=Eigen::AngleAxisd(theta,slipSystem.unitNormal)*b.normalized();
                        const std::vector<VectorDimD> nodePos(DislocationInjectorBase<dim>::straightLineBoundaryClosure(P0,d,slipSystem.unitNormal,grainID,mesh));
                        
                        const double lineLength=(nodePos.back()-nodePos.front()).norm();
                        
                        
                        
                        
                        if(nodePos.size()>=3)
                        {// Write files
                            if(lineLength<FLT_EPSILON)
                            {
                                std::cout<<"Line too short. EXITING."<<std::endl;
                                exit(EXIT_FAILURE);
                            }
                            
                            if(addSingleLoop(false,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<dim>::GLISSILELOOP))
                            {
                                std::cout<<"["<<b.transpose()<<"]("<<slipSystem.unitNormal.transpose()<<") dislocation. Line dir="<<d.transpose()<<". Length="<<lineLength<<std::endl;
                            }
                            
                        }
                        else
                        {
                            std::cout<<"nodePos.size="<<nodePos.size()<<std::endl;
                            assert(false && "LOOP DOES NOT HAVE ENOUGH POINTS");
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
        void addNonPlanarLoops()
        {
            if(targetNonPlanarLoopDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating non-planar sessile loops"<<defaultColor<<std::endl;
                double density=0.0;
                while(density<targetNonPlanarLoopDensity)
                {
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const VectorDimD P0(L0.cartesian());
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem
                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(poly.grain(grainID).latticeDirection(slipSystem.n.cartesian()).cartesian()); // Frank loop
                    const VectorDimD sessileAxis(b.normalized());
                    
                    std::normal_distribution<double> sizeDistribution(nonPlanarLoopRadiusMean/poly.b_SI,nonPlanarLoopRadiusStd/poly.b_SI);
                    const double radius(sizeDistribution(generator));
                    
                    std::normal_distribution<double> heightDistribution(nonPlanarLoopRadiusMean/poly.b_SI,0.5*nonPlanarLoopRadiusMean/poly.b_SI);
                    
                    
                    std::vector<VectorDimD> nodePos;
                    for(int k=0;k<nonPlanarLoopSides;++k)
                    {
                        const double height(heightDistribution(generator));
                        nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/nonPlanarLoopSides, sessileAxis)*slipSystem.s.cartesian().normalized()*radius+height*sessileAxis);
                    }
                    
                    if(addSingleLoop(true,nodePos,b,VectorDimD::Zero(),P0,grainID,DislocationLoopIO<dim>::SESSILELOOP))
                    {
                        density += 2.0*radius*sin(M_PI/nonPlanarLoopSides)/mesh.volume()/std::pow(poly.b_SI,2);
                        std::cout<<"non-planar loop density="<<density<<std::endl;
                    }
                    
                    //                    if(allPointsInGrain(nodePos,grainID))
                    //                    {
                    //                        for(size_t k=0;k<nonPlanarLoopSides;++k)
                    //                        {
                    //                            const size_t nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                    //                            configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    //                            configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                    //                        }
                    //                        configIO.loops().emplace_back(loopID+0, b,VectorDimD::Zero(),P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);  // write loop file
                    //                        nodeID+=nodePos.size();
                    //                        loopID+=1;
                    //                        snID+=1;
                    //                        density += 2.0*radius*sin(M_PI/nonPlanarLoopSides)/mesh.volume()/std::pow(poly.b_SI,2);
                    //                        std::cout<<"non-planar loop density="<<density<<std::endl;
                    //                    }
                    
                }
                
            }
            
        }
        
//        /**********************************************************************/
//        void addPeriodicLoops()
//        {
//
//
//            if(targetPeriodicLoopDensity>0.0)
//            {
//                std::cout<<magentaBoldColor<<"Generating periodic loops"<<defaultColor<<std::endl;
//                assert(poly.grains().size()==1 && "PeriodicLoops only supported in single crystals.");
//                double density=0.0;
//                size_t periodicNodeID(0);
//                //                size_t periodicLoopID(0);
//                while(density<targetPeriodicLoopDensity)
//                {
//                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
//                    const LatticeVector<dim> L0=rp.first;
//                    const VectorDimD P0(L0.cartesian());
//                    const int grainID=rp.second;
//                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
//                    const int rSS=distribution(generator); // a random SlipSystem
//                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
//                    const VectorDimD b(slipSystem.s.cartesian());
//
//                    std::normal_distribution<double> sizeDistribution(periodicLoopRadiusMean/poly.b_SI,periodicLoopRadiusStd/poly.b_SI);
//                    const double radius(sizeDistribution(generator));
//
//                    //                    std::vector<VectorDimD> nodePos;
//                    std::vector<DislocationNodeIO<dim>> periodicNodesIO;
//
//                    for(int k=0;k<periodicLoopSides;++k)
//                    {
//                        const VectorDimD P(P0+Eigen::AngleAxisd(k*2.0*M_PI/periodicLoopSides, slipSystem.unitNormal)*slipSystem.s.cartesian().normalized()*radius);
//                        periodicNodesIO.emplace_back(periodicNodeID,P,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
//                        periodicNodeID++;
//                    }
//
//
//                    //                    GlidePlaneKey<3> referencePlaneKey(grainID,P0,slipSystem.n);
//
//                    PeriodicGlidePlaneFactory<dim> pgpf(poly,glidePlaneFactory);
//
//
//                    //                    GlidePlaneKey<3> pointsPlaneKey(grainID,P0,slipSystem.n);
//                    //                    std::shared_ptr<GlidePlane<3>> pointsPlane(glidePlaneFactory.get(pointsPlaneKey));
//                    ////                    PeriodicGlidePlane<dim>& pgp(periodicGlidePlaneContainer.emplace_back(glidePlaneFactory,referencePlane->key,b));
//                    //                    auto periodicGlidePlane(pgpf.get(*pointsPlane));
//                    //                    periodicGlidePlane->addPatchesContainingPolygon(periodicNodesIO);
//                    //
//                    //                    // Yash: Here decompose periodic loop into individual loops
//                    //
//                    ////                    configIO.periodicLoops().emplace_back(periodicLoopID);
//                    //
//                    //                    auxIO.setGlidePlaneBoundaries(glidePlaneFactory);
//                    //                    auxIO.addPeriodicGlidePlane(*periodicGlidePlane);
//
//                    //                    periodicLoopID++;
//                    density += 2.0*radius*sin(M_PI/periodicLoopSides)/mesh.volume()/std::pow(poly.b_SI,2);
//                    std::cout<<"periodicLoop density="<<density<<std::endl;
//                }
//            }
//
//        }
        
        /**********************************************************************/
        void addFrankLoops()
        {
            if(targetFrankLoopDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Frank loops"<<defaultColor<<std::endl;
                double density=0.0;
                while(density<targetFrankLoopDensity)
                {
                    const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                    const LatticeVector<dim> L0=rp.first;
                    const VectorDimD P0(L0.cartesian());
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> planesDistribution(0,poly.grain(grainID).planeNormals().size()-1);
                    const int rSS=planesDistribution(generator); // a random SlipSystem
                    const auto& n(poly.grain(grainID).planeNormals()[rSS]);
                    const VectorDimD unitNormal(n.cartesian().normalized());
                    const VectorDimD b(density<fractionFrankVacancyLoops*targetFrankLoopDensity? n.planeSpacing()*unitNormal : -n.planeSpacing()*unitNormal); // if b*n > 0 then the loop is a vacancy loop
                    
                    std::normal_distribution<double> sizeDistribution(frankLoopRadiusMean/poly.b_SI,frankLoopRadiusStd/poly.b_SI);
                    const double radius(sizeDistribution(generator));
                    const VectorDimD R(randomOrthogonalUnitVector(b)*radius);
                    
                    std::vector<VectorDimD> nodePos;
                    for(int k=0;k<frankLoopSides;++k)
                    {
                        nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/frankLoopSides, unitNormal)*R);
                    }
                    
//                    if(density>fractionFrankVacancyLoops*targetFrankLoopDensity)
//                    {// make an interstitial loop, that is b*n<0
//                        b*=-1.0;
//                    }
                    
                    if(addSingleLoop(false,nodePos,b,unitNormal,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP))
                    {
                        density += 2.0*radius*sin(M_PI/frankLoopSides)/mesh.volume()/std::pow(poly.b_SI,2);
                        std::cout<<"Frank loop density="<<density<<std::endl;
                    }
                    

                }
            }
            
        }
        
        /**********************************************************************/
        void addEshelbyInclusions()
        {
            
            if(targetInclusionDensities.size())
            {
                double totalDensity(0.0);
                for(const double& val : targetInclusionDensities)
                {
                    if(val>0.0)
                    {
                        totalDensity+=val;
                    }
                }
                
                if(totalDensity)
                {
                    assert(targetInclusionDensities.size()==inclusionsDiameterLognormalDistribution_M.size());
                    assert(targetInclusionDensities.size()==inclusionsDiameterLognormalDistribution_S.size());
                    assert(targetInclusionDensities.size()==inclusionsDiameterLognormalDistribution_A.size());
                    assert(int(targetInclusionDensities.size())==inclusionsTransformationStrains.rows());
                    assert(int(targetInclusionDensities.size())==inclusionsPatterns.rows());
                    
                    
                    std::cout<<magentaBoldColor<<"Generating Inclusions"<<defaultColor<<std::endl;
                    
                    std::ofstream inclusionsfile("E/E_0.txt");
                    std::deque<std::pair<double,VectorDimD>> existingPrecipitates;
                    
                    size_t inclusionID=0;
                    for(size_t f=0;f<targetInclusionDensities.size();++f)
                    {
                        if(   targetInclusionDensities[f]>0.0)
                        {
                            assert(inclusionsDiameterLognormalDistribution_M[f]>0.0);
                            assert(inclusionsDiameterLognormalDistribution_S[f]>0.0);
                            assert(inclusionsDiameterLognormalDistribution_A[f]>0.0);
                            
                            const VectorDimD currentPattern(inclusionsPatterns.row(f)/poly.b_SI); // normalize to length units
                            const double patternHeight(currentPattern.norm());
                            const bool applyPattern(patternHeight>0.0);
                            const VectorDimD patternDir(applyPattern? (currentPattern/patternHeight).eval() : VectorDimD::Zero());
                            
                            double numberDensity=0.0;
                            
                            std::lognormal_distribution<double> distribution(log(inclusionsDiameterLognormalDistribution_M[f]/inclusionsDiameterLognormalDistribution_A[f]),inclusionsDiameterLognormalDistribution_S[f]);
                            
                            while(numberDensity<targetInclusionDensities[f])
                            {
                                const double diameter = distribution(generator)*inclusionsDiameterLognormalDistribution_A[f]/poly.b_SI;
                                const double radius(0.5*diameter);
                                std::pair<LatticeVector<dim>,int> pointPair=randomPointInMesh();
                                VectorDimD P=pointPair.first.cartesian();
                                const int& grainID(pointPair.second);
                                
                                if(applyPattern)
                                {
                                    const VectorDimD globalVector(poly.grain(grainID).C2G*currentPattern);
                                    const VectorDimD globalDir(poly.grain(grainID).C2G*patternDir);
                                    const long long pointHeight=std::round(P.dot(globalDir)/patternHeight);
                                    const VectorDimD O(pointHeight*globalVector);
                                    P-=(P-O).dot(globalDir)*globalDir;
                                }
                                
                                bool isGoodPosition=mesh.searchRegion(grainID,P).first;
                                for(const auto& pair : existingPrecipitates)
                                {
                                    isGoodPosition *= (P-pair.second).norm()>pair.first+radius;
                                    if(!isGoodPosition)
                                    {
                                        break;
                                    }
                                }
                                
                                if(isGoodPosition)
                                {
                                    inclusionsfile<<inclusionID
                                    /*          */<<" "<<P.transpose()
                                    /*          */<<" "<<radius
                                    /*          */<<" "<<inclusionsTransformationStrains.row(f)
                                    /*          */<<" "<<f
                                    /*          */<<"\n";
                                    
                                    numberDensity+=1.0/mesh.volume()/std::pow(poly.b_SI,3);
                                    inclusionID++;
                                    existingPrecipitates.emplace_back(radius,P);
                                    
                                    std::cout<<"inclusions density="<<numberDensity<<std::endl;
                                }
                            }
                        }
                    }
                    inclusionsfile.close();
                }
            }
        }
        
        std::mt19937 generator;
        size_t nodeID;
        size_t snID;
        size_t loopID;
//        std::list<PeriodicGlidePlane<dim>> periodicGlidePlaneContainer;
        
        std::deque<std::vector<VectorDimD>> loopPoints;
        std::deque<VectorDimD> loopBurgers;
        const bool enforceMonotonicHelicity;
        double helicity;
        
    public:
        
        DDconfigIO<3> configIO;
        DDauxIO<dim> auxIO;
        
        
        const bool outputBinary;
        //        const int meshID;
        const std::string meshFilename;
        const SimplicialMesh<dim> mesh;
        const double minSize;
        const double maxSize;
        Polycrystal<dim> poly;
        GlidePlaneFactory<dim> glidePlaneFactory;
        
        // Straight Dislocations
        const double targetStraightDislocationDensity;
        const double fractionSessileStraightDislocationDensity;
        
        const double targetJoggedtDislocationDensity;
        
        // Frank-Read sources
        const double targetFrankReadDislocationDensity;
        const double FrankReadSizeMean;
        const double FrankReadSizeStd;
        const double FrankReadAspectRatioMean;
        const double FrankReadAspectRatioStd;
        
        // Single-arm sources
        const double targetSingleArmDislocationDensity;
        
        // Prismatic loops
        const double targetPrismaticLoopDensity;
        
        // Individual dislocations
        const std::vector<int> straightDislocationsSlipSystemIDs;
        const std::vector<double> straightDislocationsAngleFromScrewOrientation;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> pointsAlongStraightDislocations;
        
        // Frank  Loops
        const double targetFrankLoopDensity;
        const double frankLoopRadiusMean;
        const double frankLoopRadiusStd;
        const double frankLoopSides;
        const double fractionFrankVacancyLoops;
        
        // NonPlanar  Loops
        const double targetNonPlanarLoopDensity;
        const double nonPlanarLoopRadiusMean;
        const double nonPlanarLoopRadiusStd;
        const double nonPlanarLoopSides;
        
//        // Periodic  Loops
//        const double targetPeriodicLoopDensity;
//        const double periodicLoopRadiusMean;
//        const double periodicLoopRadiusStd;
//        const double periodicLoopSides;
        
        // Irradiation Loops
        const double targetIrradiationLoopDensity;
        //        const double averageLoopSize;
        const double irradiationLoopsDiameterLognormalDistribution_M;
        const double irradiationLoopsDiameterLognormalDistribution_S;
        const double irradiationLoopsDiameterLognormalDistribution_A;
        const double fraction111Loops;  // fraction of [111] glissile loop;
        const bool mobile111Loops;
        
        // SFTs
        const double targetSFTdensity;
        const double sftSizeMean;
        const double sftSizeStd;
        
        const std::vector<int> sftPlaneIDs;
        const std::vector<int> sftIsInverted;
        const std::vector<double> sftSizes;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> sftBasePoints;
        
        
        // Inclusions
        const std::vector<double> targetInclusionDensities;
        const std::vector<double> inclusionsDiameterLognormalDistribution_M;
        const std::vector<double> inclusionsDiameterLognormalDistribution_S;
        const std::vector<double> inclusionsDiameterLognormalDistribution_A;
        const Eigen::Matrix<double,Eigen::Dynamic,dim*dim> inclusionsTransformationStrains;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> inclusionsPatterns;
        //        const std::vector<double> inclusionsMobilityReduction;
        
        bool isInclusionsUsed() const
        {
            bool temp(false);
            for(const auto& density : targetInclusionDensities)
            {
                if(density>0.0)
                {
                    temp=true;
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        MicrostructureGenerator(int argc, char* argv[]) :
        /* init*/ generator(std::chrono::system_clock::now().time_since_epoch().count())
        /* init*/,nodeID(0)
        /* init*/,snID(0)
        /* init*/,loopID(0)
        /* init*/,enforceMonotonicHelicity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("enforceMonotonicHelicity",true))
        /* init*/,helicity(0.0)
        /* init*/,outputBinary(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputBinary",true))
        //        /* init*/,meshID(TextFileParser("./inputFiles/DD.txt").readScalar<int>("meshID",true))
        /* init*/,meshFilename(TextFileParser("./inputFiles/polycrystal.txt").readString("meshFile",true))
        /* init*/,mesh(meshFilename,TextFileParser("./inputFiles/polycrystal.txt").readMatrix<double>("A",3,3,true),TextFileParser("./inputFiles/polycrystal.txt").readMatrix<double>("x0",1,3,true).transpose())
        /* init*/,minSize(0.1*min(mesh.xMax(0)-mesh.xMin(0),min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
        /* init*/,maxSize(max(mesh.xMax(0)-mesh.xMin(0),max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
        /* init*/,poly("./inputFiles/polycrystal.txt",mesh)
        /* init*/,glidePlaneFactory(poly)
        /* Straight Dislocations */
        /* init*/,targetStraightDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetStraightDislocationDensity",true))
        /* init*/,fractionSessileStraightDislocationDensity(targetStraightDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("fractionSessileStraightDislocationDensity",true) : 0.0)
        /* Jogged Dislocations */
        /* init*/,targetJoggedtDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetJoggedtDislocationDensity",true))
        /* Frank-Read sources */
        /* init*/,targetFrankReadDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetFrankReadDislocationDensity",true))
        /* init*/,FrankReadSizeMean(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadSizeMean",true) : 0.0)
        /* init*/,FrankReadSizeStd(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadSizeStd",true) : 0.0)
        /* init*/,FrankReadAspectRatioMean(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadAspectRatioMean",true) : 0.0)
        /* init*/,FrankReadAspectRatioStd(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadAspectRatioStd",true) : 0.0)
        /* Single-arm sources */
        /* init*/,targetSingleArmDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetSingleArmDislocationDensity",true))
        /* Prismatic loops */
        /* init*/,targetPrismaticLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetPrismaticLoopDensity",true))
        /* Indivial straight dislocations */
        /* init*/,straightDislocationsSlipSystemIDs(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<int>("straightDislocationsSlipSystemIDs",true))
        /* init*/,straightDislocationsAngleFromScrewOrientation(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("straightDislocationsAngleFromScrewOrientation",true))
        /* init*/,pointsAlongStraightDislocations(TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("pointsAlongStraightDislocations",straightDislocationsSlipSystemIDs.size(),dim,true))
        /* Frank Loops */
        /* init*/,targetFrankLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetFrankLoopDensity",true))
        /* init*/,frankLoopRadiusMean(targetFrankLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("frankLoopRadiusMean",true) : 0.0)
        /* init*/,frankLoopRadiusStd(targetFrankLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("frankLoopRadiusStd",true) : 0.0)
        /* init*/,frankLoopSides(targetFrankLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("frankLoopSides",true) : 0.0)
        /* init*/,fractionFrankVacancyLoops(targetFrankLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("fractionFrankVacancyLoops",true) : 1.0)
        /* NonPlanar Loops */
        /* init*/,targetNonPlanarLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetNonPlanarLoopDensity",true))
        /* init*/,nonPlanarLoopRadiusMean(targetNonPlanarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("nonPlanarLoopRadiusMean",true) : 0.0)
        /* init*/,nonPlanarLoopRadiusStd(targetNonPlanarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("nonPlanarLoopRadiusStd",true) : 0.0)
        /* init*/,nonPlanarLoopSides(targetNonPlanarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("nonPlanarLoopSides",true) : 0.0)
//        /* PeriodicLoops */
//        /* init*/,targetPeriodicLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetPeriodicLoopDensity",true))
//        /* init*/,periodicLoopRadiusMean(targetPeriodicLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("periodicLoopRadiusMean",true) : 0.0)
//        /* init*/,periodicLoopRadiusStd(targetPeriodicLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("periodicLoopRadiusStd",true) : 0.0)
//        /* init*/,periodicLoopSides(targetPeriodicLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("periodicLoopSides",true) : 0.0)
        /* Irradiation Loops */
        /* init*/,targetIrradiationLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetIrradiationLoopDensity",true))
        /* init*/,irradiationLoopsDiameterLognormalDistribution_M(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("irradiationLoopsDiameterLognormalDistribution_M",true) : 0.0)
        /* init*/,irradiationLoopsDiameterLognormalDistribution_S(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("irradiationLoopsDiameterLognormalDistribution_S",true) : 0.0)
        /* init*/,irradiationLoopsDiameterLognormalDistribution_A(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("irradiationLoopsDiameterLognormalDistribution_A",true) : 0.0)
        /* init*/,fraction111Loops(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("fraction111Loops",true) : 0.0)
        /* init*/,mobile111Loops(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("mobile111Loops",true) : 0.0)
        /* SFTs */
        /* init*/,targetSFTdensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetSFTdensity",true))
        /* init*/,sftSizeMean(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("sftSizeMean",true))
        /* init*/,sftSizeStd(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("sftSizeStd",true))
        /* individual SFTs */
        /* init*/,sftPlaneIDs(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<int>("sftPlaneIDs",true))
        /* init*/,sftIsInverted(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<int>("sftIsInverted",true))
        /* init*/,sftSizes(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("sftSizes",true))
        /* init*/,sftBasePoints(TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("sftBasePoints",sftPlaneIDs.size(),dim,true))
        /* Inclusions */
        /* init*/,targetInclusionDensities(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("targetInclusionDensities",true))
        /* init*/,inclusionsDiameterLognormalDistribution_M(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_M",true) : std::vector<double>())
        /* init*/,inclusionsDiameterLognormalDistribution_S(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_S",true) : std::vector<double>())
        /* init*/,inclusionsDiameterLognormalDistribution_A(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_A",true) : std::vector<double>())
        /* init*/,inclusionsTransformationStrains(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsTransformationStrains",targetInclusionDensities.size(),dim*dim,true) : Eigen::Matrix<double,Eigen::Dynamic,dim*dim>::Zero(1,dim*dim))
        /* init*/,inclusionsPatterns(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsPatterns",targetInclusionDensities.size(),dim,true) : Eigen::Matrix<double,Eigen::Dynamic,dim>::Zero(1,dim))
        //        /* init*/,inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true))
        {
            
            
            
            
            // Some sanity checks
            if(mesh.volume()==0.0)
            {
                std::cout<<"mesh "<<meshFilename<<" is empty. MicrostructureGenerator cannot run. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            std::cout<<greenBoldColor<<"Generating initial microstructure"<<defaultColor<<std::endl;
            // Call individual generators
            addStraightDislocations();
            addJoggedDislocations();
            addFrankReadSources();
            addSingleArmDislocations();
            addPrismaticLoops();
            addIndividualStraightDislocations();
            addFrankLoops();
            addNonPlanarLoops();
//            addPeriodicLoops();
            addIrradiationLoops();
            addStackingFaultTetrahedra();
            addIndividualStackingFaultTetrahedra();
            addEshelbyInclusions();
            writeConfigFiles(0);
            
        }
        
        /**********************************************************************/
        void writeConfigFiles(const size_t& fileID)
        {
            
            auxIO.setGlidePlaneBoundaries(glidePlaneFactory); // change this function to take a GlidePlaneFactory during write
            
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
        
        /**********************************************************************/
        void addIndividualStackingFaultTetrahedra()
        {
            if(sftPlaneIDs.size())
            {
                std::cout<<magentaBoldColor<<"Generating individual Stacking Fault Tetrahedra"<<defaultColor<<std::endl;
                if(int(sftPlaneIDs.size())!=sftBasePoints.rows())
                {
                    std::cout<<"sftPlaneIDs()="<<sftPlaneIDs.size()<<std::endl;
                    std::cout<<"sftBasePoints.rows()="<<sftBasePoints.rows()<<std::endl;
                    std::cout<<"You must provide one point for each SFT. Each point is a row of the matrix sftBasePoints. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                if(sftPlaneIDs.size()!=sftSizes.size())
                {
                    std::cout<<"sftPlaneIDs()="<<sftPlaneIDs.size()<<std::endl;
                    std::cout<<"sftSizes.size()="<<sftSizes.size()<<std::endl;
                    std::cout<<"You must provide one size for each SFT. Each size is an element of the vector sftSizes. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                if(sftPlaneIDs.size()!=sftIsInverted.size())
                {
                    std::cout<<"sftPlaneIDs()="<<sftPlaneIDs.size()<<std::endl;
                    std::cout<<"sftIsInverted()="<<sftSizes.size()<<std::endl;
                    std::cout<<"You must provide one boolean value of sftIsInverted for each SFT. Each value is an element of the vector sftIsInverted. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                int generated(0);
                for(size_t k=0;k<sftPlaneIDs.size();++k)
                {
                    const int& rSS(sftPlaneIDs[k]);
                    
                    
                    if(rSS>=0)
                    {
                        // Get grain
                        std::pair<bool,const Simplex<dim,dim>*> found=mesh.search(sftBasePoints.row(k));
                        if(!found.first)
                        {
                            std::cout<<"Point "<<sftBasePoints.row(k)<<" is outside mesh. EXITING."<<std::endl;
                            exit(EXIT_FAILURE);
                        }
                        const int grainID=found.second->region->regionID;
                        
                        // Get slip system
                        if(rSS>=int(poly.grain(grainID).planeNormals().size()))
                        {
                            std::cout<<"requested planeNormal ID="<<rSS<<std::endl;
                            std::cout<<"# of plane normals ="<<poly.grain(grainID).planeNormals().size()<<std::endl;
                            std::cout<<"Requested planeNormal does not exist. EXITING."<<std::endl;
                            exit(EXIT_FAILURE);
                        }
//                        const auto& planeNormal(poly.grain(grainID).planeNormals()[rSS]);
                        
                        const auto basePlane(sftIsInverted[k]? LatticePlaneBase(poly.grain(grainID).planeNormals()[rSS].primitiveVectors.second*(-1),poly.grain(grainID).planeNormals()[rSS].primitiveVectors.first*(-1)) : poly.grain(grainID).planeNormals()[rSS]);
                        const auto unitNormal(basePlane.cartesian().normalized());
                        const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(basePlane,sftBasePoints.row(k));
                        const VectorDimD basePoint(sftBasePoints.row(k).transpose()-sftBasePoints.row(k).transpose().dot(unitNormal)*unitNormal+unitNormal*basePlane.planeSpacing()*heightPair.second);
                        generated+=addSingleStackingFaultTetrahedron(grainID,basePlane,basePoint,std::round(sftSizes[k]/poly.b_SI));
                        std::cout<<"# individual SFTs="<<generated<<std::endl;
                    }
                    else
                    {
                        std::cout<<"negative slip system ID. Skipping entries."<<std::endl;
                    }
                }
            }
        }
        
        /**********************************************************************/
        bool addSingleStackingFaultTetrahedron(const int& grainID,const LatticePlaneBase& basePlane,const VectorDimD& tentativeBasePoint,const int& Li)
        {
            if(Li>1)
            {
                const VectorDimD unitNormal(basePlane.cartesian().normalized());
                std::vector<VectorDimD> nodePos;
                const auto faceBaryShift((basePlane.primitiveVectors.first+basePlane.primitiveVectors.second).cartesian()*Li/3.0);
                const VectorDimD corner(poly.grain(grainID).snapToLattice(tentativeBasePoint-faceBaryShift).cartesian());
                nodePos.push_back(corner);
                nodePos.push_back(corner+basePlane.primitiveVectors.first.cartesian()*Li);
                nodePos.push_back(corner+basePlane.primitiveVectors.second.cartesian()*Li);
                nodePos.push_back((nodePos[0]+nodePos[1]+nodePos[2])/3.0-unitNormal*Li*sqrt(2.0/3.0));
                
                if(allPointsInGrain(nodePos,grainID))
                {
                    
                    configIO.nodes().emplace_back(nodeID+0,nodePos[0],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0); // base point   0
                    configIO.nodes().emplace_back(nodeID+1,nodePos[1],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0); // base point   1
                    configIO.nodes().emplace_back(nodeID+2,nodePos[2],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0); // base point   2
                    configIO.nodes().emplace_back(nodeID+3,nodePos[3],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0); // vertex point 3
                    
                    // base loop 0->1->2 (the Frank loop)
                    const VectorDimD b012(basePlane.interplaneVector()); // Frank loops are platelets of vanancies
                    configIO.links().emplace_back(loopID+0,nodeID+0,nodeID+1,0);
                    configIO.links().emplace_back(loopID+0,nodeID+1,nodeID+2,0);
                    configIO.links().emplace_back(loopID+0,nodeID+2,nodeID+0,0);
                    configIO.loops().emplace_back(loopID+0,b012,unitNormal,nodePos[0],grainID,DislocationLoopIO<dim>::SESSILELOOP);
                    
                    // loop 1->3->2
                    const VectorDimD b132((0.5*(nodePos[1]+nodePos[2])-nodePos[3]).normalized()*sqrt(3.0)/3.0);
                    configIO.links().emplace_back(loopID+1,nodeID+1,nodeID+3,0);
                    configIO.links().emplace_back(loopID+1,nodeID+3,nodeID+2,0);
                    configIO.links().emplace_back(loopID+1,nodeID+2,nodeID+1,0);
                    configIO.loops().emplace_back(loopID+1,b132,(nodePos[3]-nodePos[1]).cross(nodePos[2]-nodePos[3]).normalized(),nodePos[1],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                    
                    // loop 2->3->0
                    const VectorDimD b230((0.5*(nodePos[2]+nodePos[0])-nodePos[3]).normalized()*sqrt(3.0)/3.0);
                    configIO.links().emplace_back(loopID+2,nodeID+2,nodeID+3,0);
                    configIO.links().emplace_back(loopID+2,nodeID+3,nodeID+0,0);
                    configIO.links().emplace_back(loopID+2,nodeID+0,nodeID+2,0);
                    configIO.loops().emplace_back(loopID+2,b230,(nodePos[3]-nodePos[2]).cross(nodePos[0]-nodePos[3]).normalized(),nodePos[2],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                    
                    // loop 0->3->1
                    const VectorDimD b031((0.5*(nodePos[0]+nodePos[1])-nodePos[3]).normalized()*sqrt(3.0)/3.0);
                    configIO.links().emplace_back(loopID+3,nodeID+0,nodeID+3,0);
                    configIO.links().emplace_back(loopID+3,nodeID+3,nodeID+1,0);
                    configIO.links().emplace_back(loopID+3,nodeID+1,nodeID+0,0);
                    configIO.loops().emplace_back(loopID+3,b031,(nodePos[3]-nodePos[0]).cross(nodePos[1]-nodePos[3]).normalized(),nodePos[0],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                    
                    
                    snID++;
                    loopID+=4;
                    nodeID+=4;
                    return true;
                }
            }
            return false;
            
        }
        
        
        /**********************************************************************/
        void addStackingFaultTetrahedra()
        {
            
            if(targetSFTdensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Stacking Fault Tetrahedra"<<defaultColor<<std::endl;
                
                if(poly.crystalStructure=="FCC")
                {
                    size_t ndefects=0;
                    double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                    
                    while(defectsDensity<targetSFTdensity)
                    {
                        const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                        const int& grainID=rp.second;   // random grain ID
                        const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                        
                        std::uniform_int_distribution<> planeNormalDistribution(0,poly.grain(grainID).planeNormals().size()-1);
                        const int rSS=planeNormalDistribution(generator); // a random SlipSystem ID
                        
                        std::uniform_int_distribution<> invertedThompsonTetrahedronDistribution(0,1);
                        const bool invertedSFT(invertedThompsonTetrahedronDistribution(generator));
                        const auto basePlane(invertedSFT? LatticePlaneBase(poly.grain(grainID).planeNormals()[rSS].primitiveVectors.second*(-1),poly.grain(grainID).planeNormals()[rSS].primitiveVectors.first*(-1)) : poly.grain(grainID).planeNormals()[rSS]);
//
//
//                        const LatticePlaneBase& basePlane(poly.grain(grainID).planeNormals()[rSS]);
                        
                        std::normal_distribution<double> sftSizeDistribution(sftSizeMean,sftSizeStd);
                        
                        const int Li(std::round(sftSizeDistribution(generator)/poly.b_SI));
                        
                        ndefects+=addSingleStackingFaultTetrahedron(grainID,basePlane,L0.cartesian(),Li);
                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"SFT density="<<defectsDensity<<std::endl;
                    }
                }
                else
                {
                    std::cout<<"SFTs can only be generated for FCC crystals"<<std::endl;
                }
            }
            
        }
        
        /**********************************************************************/
        bool allPointsInGrain(const std::vector<VectorDimD>& points,const int& grainID)
        {
            bool temp=true;
            for(const auto& point : points)
            {
                temp*=mesh.searchRegion(grainID,point).first;
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        void addIrradiationLoopsFCC()
        {// Irradiation loops in FCC are Frank loops
            
            const size_t irradiationLoopsNumberOfNodes(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("irradiationLoopsNumberOfNodes",true));
            std::lognormal_distribution<double> sizeDistribution(log(irradiationLoopsDiameterLognormalDistribution_M/irradiationLoopsDiameterLognormalDistribution_A),(irradiationLoopsDiameterLognormalDistribution_S));
            
            size_t ndefects=0;
            double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
            while(defectsDensity<targetIrradiationLoopDensity)
            {
                const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                
                std::uniform_int_distribution<> planesDistribution(0,poly.grain(grainID).planeNormals().size()-1);
                const int rSS=planesDistribution(generator);
                const auto& n(poly.grain(grainID).planeNormals()[rSS]); // a random {111} plane
                const VectorDimD unitNormal(n.cartesian().normalized());
                const VectorDimD b(n.planeSpacing()*unitNormal); // Frank loops are plateletslane of vanancies/interstitials
                
                const double diameter_SI = sizeDistribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
                const double radius(0.5*diameter_SI/poly.b_SI);
                const VectorDimD R(randomOrthogonalUnitVector(b)*radius);
                
                std::vector<VectorDimD> nodePos;
                for(size_t k=0;k<irradiationLoopsNumberOfNodes;++k)
                {
                    nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/irradiationLoopsNumberOfNodes, unitNormal)*R);
                }
                
                if(addSingleLoop(true,nodePos,b,unitNormal,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP))
                {
                    ndefects++;
                    defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                    std::cout<<"irradiation loops density="<<defectsDensity<<std::endl;
                }
            }
        }
        
        /**********************************************************************/
        void addIrradiationLoopsHCP()
        {// Irradiation loops in HCP are prismatic+basal loops
            
            const double PrismaticLoopSizeMean(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("PrismaticLoopSizeMean",true) : 0.0);
            const double PrismaticLoopSizeStd(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("PrismaticLoopSizeStd",true) : 0.0);
            const double BasalLoopSizeMean(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("BasalLoopSizeMean",true) : 0.0);
            const double BasalLoopSizeStd(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("BasalLoopSizeStd",true) : 0.0);
            const size_t irradiationLoopsNumberOfNodes(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("irradiationLoopsNumberOfNodes",true));
            
            
            // fraction111Loops represents the fracton of prismatic loops for HCP
            double targetPrismLoopDensity=targetIrradiationLoopDensity*fraction111Loops;
            double targetBasalLoopDesnity=targetIrradiationLoopDensity*(1.-fraction111Loops);
            size_t ndefects=0, nprism=0, nbasal=0;
            double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
            double prismDensity=nprism/mesh.volume()/std::pow(poly.b_SI,3);
            double basalDensity=nbasal/mesh.volume()/std::pow(poly.b_SI,3);
            
            std::cout<<"HCP loop generator with target prismatic loop density:"<<targetPrismLoopDensity<<" and target basal loop density:"<<targetBasalLoopDesnity<<std::endl;
            
            while(defectsDensity<targetIrradiationLoopDensity)
            {
                const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                
                std::uniform_int_distribution<> slipSystemDistribution(0,poly.grain(grainID).slipSystems().size()-1);
                const int rSS=slipSystemDistribution(generator); // a random SlipSystem ID
                const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                
                // For HCP type, the number of defects (without units) in loops follows the normal distribution.
                std::normal_distribution<double> prismsizeDistribution(PrismaticLoopSizeMean,PrismaticLoopSizeStd);
                std::normal_distribution<double> basalsizeDistribution(BasalLoopSizeMean,BasalLoopSizeStd);
                
                const double planeSpacing(slipSystem.n.planeSpacing());
                if(fabs(planeSpacing-sqrt(3.0)/2.0)<FLT_EPSILON && prismDensity<targetPrismLoopDensity)
                {// prismatic plane spacing
                    
                    // For square shape prismatic loops, radius is half of the edge length
                    // Omega/pi/ba in m^-2
                    const double size2radius=2.33e-29/(poly.b_SI*4.);
                    double xsize=prismsizeDistribution(generator);
                    // if((xsize<(PrismaticLoopSizeMean-3.0*PrismaticLoopSizeStd) || xsize<0) || xsize>(PrismaticLoopSizeMean+3.0*PrismaticLoopSizeStd)) continue;
                    if(xsize<0) continue;
                    const double radius=sqrt(xsize*size2radius)/poly.b_SI;
                    
                    const VectorDimD e1=slipSystem.unitNormal;
                    const VectorDimD e2=slipSystem.s.cartesian().cross(slipSystem.n.cartesian()).normalized();                 // "edge" directions
                    
                    
                    std::vector<VectorDimD> nodePos;
                    nodePos.push_back(P0);
                    nodePos.push_back(P0+e2*std::round(2.0*radius/sqrt(8.0/3.0))*sqrt(8.0/3.0));
                    nodePos.push_back(P0+e1*std::round(2.0*radius/sqrt(3.0/4.0))*sqrt(3.0/4.0)+e2*std::round(2.0*radius/sqrt(8.0/3.0))*sqrt(8.0/3.0));
                    nodePos.push_back(P0+e1*std::round(2.0*radius/sqrt(3.0/4.0))*sqrt(3.0/4.0));
                    VectorDimD centerofloop(0.,0.,0.);
                    for(size_t k=0;k<4;++k) centerofloop += nodePos[k];
                    
                    
                    VectorDimD b=slipSystem.s.cartesian();   // Prism axis
                    
                    // make sure it is interstitial type
                    VectorDimD Cycle_plane=(nodePos[1]-nodePos[0]).cross(nodePos[2]-nodePos[1]);
                    if (b.dot(Cycle_plane)>0)
                    {
                        b*=-1.0;
                    }
                    
                    /*
                     std::cout<<"Radius is:"<<radius<<std::endl;
                     std::cout<<"Nodepos 0:"<<nodePos[0]<<std::endl;
                     std::cout<<"Nodepos 1:"<<nodePos[1]<<std::endl;
                     std::cout<<"Nodepos 2:"<<nodePos[2]<<std::endl;
                     std::cout<<"Nodepos 3:"<<nodePos[3]<<std::endl;
                     */
                    
                    if(allPointsInGrain(nodePos,grainID))
                    {
                        
                        for(int k=0;k<4;++k)
                        {// inser the back loop
                            configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            const int nextNodeID=(k+1)<4? nodeID+k+1 : nodeID;
                            configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                        }
                        configIO.loops().emplace_back(loopID+0, b,slipSystem.s.cartesian().normalized(),P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                        loopID++;
                        
                        if(mobile111Loops)
                        {// insert loops on the four sides
                            for(int k=0;k<4;++k)
                            {// inser lateral loops
                                configIO.nodes().emplace_back(nodeID+k+4,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                const int nextNodeID=(k+1)<4? nodeID+k+1 : nodeID;
                                configIO.links().emplace_back(loopID,nodeID+k,nodeID+k+4,0);
                                configIO.links().emplace_back(loopID,nodeID+k+4,nextNodeID+4,0);
                                configIO.links().emplace_back(loopID,nextNodeID+4,nextNodeID,0);
                                configIO.links().emplace_back(loopID,nextNodeID,nodeID+k,0);
                                
                                configIO.loops().emplace_back(loopID+0, b,Eigen::AngleAxis<double>(k*2.0*M_PI/4,b.normalized())*slipSystem.unitNormal,nodePos[k],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                loopID++;
                            }
                            nodeID+=2*nodePos.size();
                        }
                        else
                        {
                            nodeID+=nodePos.size();
                        }
                        snID++;
                        nprism++;
                        prismDensity=nprism/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"Prismatci loops density="<<prismDensity<<std::endl;
                    }
                }
                
                else if(fabs(planeSpacing-sqrt(8.0/3.0))<FLT_EPSILON && basalDensity<targetBasalLoopDesnity)
                {// basal plane spacing
                    VectorDimD b(0.5*slipSystem.n.planeSpacing()*slipSystem.n.cartesian().normalized()); // 1/2 c-type loop
                    
                    // Omega/pi/bc in m^-2
                    const double size2radius=2.33e-29/(0.5*poly.b_SI*sqrt(8.0/3.0)*M_PI);
                    double xsize=basalsizeDistribution(generator);
                    // if((xsize<(BasalLoopSizeMean-3.0*BasalLoopSizeStd) || xsize<0) || xsize>(BasalLoopSizeMean+3.0*BasalLoopSizeStd)) continue;
                    if(xsize<0) continue;
                    const double radius=sqrt(xsize*size2radius)/poly.b_SI;
                    
                    const VectorDimD R(slipSystem.s.cartesian().normalized()*radius);
                    std::vector<VectorDimD> nodePos;
                    for(size_t k=0;k<irradiationLoopsNumberOfNodes;++k)
                    {
                        nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/irradiationLoopsNumberOfNodes, slipSystem.unitNormal)*R);
                    }
                    
                    // make sure it is vacancy type
                    VectorDimD Cycle_plane=(nodePos[1]-nodePos[0]).cross(nodePos[2]-nodePos[1]);
                    if (b.dot(Cycle_plane)<0)
                    {
                        b*=-1.0;
                    }
                    
                    if(addSingleLoop(false,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP))
                    {
                        nbasal++;
                        basalDensity=nbasal/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"Basal loops density="<<basalDensity<<std::endl;
                    }
                    
                    /*
                     std::cout<<"Basal radius is:"<<radius<<std::endl;
                     std::cout<<"Nodepos 0:"<<nodePos[0]<<std::endl;
                     std::cout<<"Nodepos 1:"<<nodePos[1]<<std::endl;
                     */
                    
                }
                else
                {
                    
                }
                ndefects = nprism + nbasal;
                defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                std::cout<<"Total loop density="<<defectsDensity<<std::endl;
            }
        }
        
        //        /**********************************************************************/
        //        void addIrradiationLoopsHCP()
        //        {// Irradiation loops in FCC are Frank loops.
        //         // TO DO: actaully Frank loops can be unfaulted and become prismatic (exagonal) loops bounded by pairs of {111} planes and {001} planes. See P. B. Hirsch , J. Silcox , R. E. Smallman & K. H. Westmacott
        //
        //            const size_t irradiationLoopsNumberOfNodes(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("irradiationLoopsNumberOfNodes",true));
        //            std::lognormal_distribution<double> sizeDistribution(log(irradiationLoopsDiameterLognormalDistribution_M/irradiationLoopsDiameterLognormalDistribution_A),(irradiationLoopsDiameterLognormalDistribution_S));
        //
        //            size_t ndefects=0;
        //            double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
        //            while(defectsDensity<targetIrradiationLoopDensity)
        //            {
        //                const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
        //                const int& grainID=rp.second;   // random grain ID
        //                const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
        //                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
        //
        //                std::uniform_int_distribution<> slipSystemDistribution(0,poly.grain(grainID).slipSystems().size()-1);
        //                const int rSS=slipSystemDistribution(generator); // a random SlipSystem ID
        //                const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
        //
        //                const double diameter_SI = sizeDistribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
        //                const double radius(0.5*diameter_SI/poly.b_SI);
        //                const VectorDimD R(slipSystem.s.cartesian().normalized()*radius);
        //
        //
        //                std::vector<VectorDimD> nodePos;
        //                for(size_t k=0;k<irradiationLoopsNumberOfNodes;++k)
        //                {
        //                    nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/irradiationLoopsNumberOfNodes, slipSystem.unitNormal)*R);
        //                }
        //
        //                const double planeSpacing(slipSystem.n.planeSpacing());
        //                if(fabs(planeSpacing-sqrt(3.0)/2.0)<FLT_EPSILON)
        //                {// prismatic plane spacing
        //                    const VectorDimD e=slipSystem.s.cartesian().cross(slipSystem.n.cartesian()).normalized();                 // "edge" direction, along prism axis
        //                    const VectorDimD b=Eigen::AngleAxisd(randomSign()*M_PI/3.0, e)*(slipSystem.s.cartesian());                  // rotate slip direction out of the plane by 60 deg
        //                    if(addSingleLoop(true,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP)
        //                    {
        //                        ndefects++;
        //                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
        //                        std::cout<<"irradiation loops density="<<defectsDensity<<std::endl;
        //                    }
        //                }
        //                else if(fabs(planeSpacing-sqrt(8.0/3.0))<FLT_EPSILON)
        //                {// basal plane spacing
        //                    const VectorDimD b(0.5*slipSystem.n.planeSpacing()*slipSystem.n.cartesian().normalized()); // 1/2 c-type loop
        //                    if(addSingleLoop(true,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP)
        //                    {
        //                        ndefects++;
        //                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
        //                        std::cout<<"irradiation loops density="<<defectsDensity<<std::endl;
        //                    }
        //                }
        //                else
        //                {
        //
        //                }
        //            }
        //        }
        
        
        /**********************************************************************/
        void addIrradiationLoopsBCC()
        {
            size_t ndefects=0;
            double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
            //                    int NP=6;
            std::lognormal_distribution<double> distribution(log(irradiationLoopsDiameterLognormalDistribution_M/irradiationLoopsDiameterLognormalDistribution_A),(irradiationLoopsDiameterLognormalDistribution_S));
            
            while(defectsDensity<targetIrradiationLoopDensity)
            {
                const std::pair<LatticeVector<dim>,int> rp=randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                
                if (defectsDensity<targetIrradiationLoopDensity*fraction111Loops)
                {// add [111] loops
                    std::uniform_int_distribution<> slipSystemDistribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=slipSystemDistribution(generator); // a random SlipSystem ID
                    const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
                    const VectorDimD a=b.normalized();
                    
                    const double diameter_SI = distribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
                    const double radius_SI(0.5*diameter_SI);
                    const double radius(std::round(radius_SI/poly.b_SI/(2.0*sqrt(2.0)/3.0))*(2.0*sqrt(2.0)/3.0)); // radius must be an integer multiple of 2.0*sqrt(2.0)/3.0, or lateral sides will not be on glide planes
                    
                    if(radius>0.0)
                    {
                        std::vector<VectorDimD> points;
                        for(int k=0;k<6;++k)
                        {
                            points.push_back(P0+Eigen::AngleAxis<double>(k*2.0*M_PI/6-M_PI/6,a)*slipSystem.unitNormal*radius);
                        }
                        
                        if(allPointsInGrain(points,grainID))
                        {
                            
                            for(int k=0;k<6;++k)
                            {// inser the back loop
                                configIO.nodes().emplace_back(nodeID+k,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                const int nextNodeID=(k+1)<6? nodeID+k+1 : nodeID;
                                configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                            }
                            configIO.loops().emplace_back(loopID+0, b,a,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                            loopID++;
                            
                            
                            if(mobile111Loops)
                            {// insert loops on the six sides
                                for(int k=0;k<6;++k)
                                {// inser lateral loops
                                    configIO.nodes().emplace_back(nodeID+k+6,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                                    const int nextNodeID=(k+1)<6? nodeID+k+1 : nodeID;
                                    configIO.links().emplace_back(loopID,nodeID+k,nodeID+k+6,0);
                                    configIO.links().emplace_back(loopID,nodeID+k+6,nextNodeID+6,0);
                                    configIO.links().emplace_back(loopID,nextNodeID+6,nextNodeID,0);
                                    configIO.links().emplace_back(loopID,nextNodeID,nodeID+k,0);
                                    
                                    configIO.loops().emplace_back(loopID+0, b,Eigen::AngleAxis<double>(k*2.0*M_PI/6,a)*slipSystem.unitNormal,points[k],grainID,DislocationLoopIO<dim>::GLISSILELOOP);
                                    loopID++;
                                }
                                nodeID+=2*points.size();
                            }
                            else
                            {
                                nodeID+=points.size();
                            }
                            snID++;
                            ndefects++;
                            defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                            std::cout<<"irradiation defects density="<<defectsDensity<<std::endl;
                        }
                    }
                }
                else
                {// add [100] loop (by Yinan Cui)
                    std::vector<LatticeDirectionType> sessileb;
                    sessileb.emplace_back(LatticeVector<dim>(VectorDimI(0,1,1),poly.grain(grainID).lattice())); // is ( 1, 0, 0) in cartesian
                    sessileb.emplace_back(LatticeVector<dim>(VectorDimI(1,0,1),poly.grain(grainID).lattice())); // is ( 0, 1, 0) in cartesian
                    sessileb.emplace_back(LatticeVector<dim>(VectorDimI(1,1,0),poly.grain(grainID).lattice())); // is ( 0, 0, 1) in cartesian
                    
                    std::uniform_int_distribution<> dist(0,2);
                    typedef Eigen::Matrix<long int,2,1> Vector2I;
                    std::vector<Vector2I> sslinedirection;
                    sslinedirection.emplace_back(Vector2I(1,2));
                    sslinedirection.emplace_back(Vector2I(0,2));
                    sslinedirection.emplace_back(Vector2I(0,1));
                    
                    //NP=4;
                    const int rSS_sessile=dist(generator); // a random sessile plane
                    LatticeDirection<3> d1(sessileb[sslinedirection[rSS_sessile][0]]);
                    LatticeDirection<3> d2(sessileb[sslinedirection[rSS_sessile][1]]);
                    
                    const double diameter_SI = distribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
                    double a1=diameter_SI/poly.b_SI;
                    double a2=diameter_SI/poly.b_SI;
                    
                    std::vector<VectorDimD> points;
                    points.push_back(P0);
                    points.push_back(P0+d1.cartesian().normalized()*a1);
                    points.push_back(P0+d1.cartesian().normalized()*a1+d2.cartesian().normalized()*a2);
                    points.push_back(P0+d2.cartesian().normalized()*a2);
                    
                    // make sure it is interstitial
                    VectorDimD Cycle_plane=d1.cartesian().cross(d2.cartesian());
                    VectorDimD b=sessileb[rSS_sessile].cartesian();
                    if (b.dot(Cycle_plane)<0)
                    {
                        b*=-1.0;
                    }
                    const VectorDimD a=(b.normalized());
                    
                    if(allPointsInGrain(points,grainID))
                    {
                        
                        for(int k=0;k<4;++k)
                        {
                            configIO.nodes().emplace_back(nodeID+k,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                            
                            const int nextNodeID=(k+1)<4? nodeID+k+1 : nodeID;
                            configIO.links().emplace_back(loopID,nodeID+k,nextNodeID,0);
                            
                        }
                        
                        configIO.loops().emplace_back(loopID+0, b,a,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                        
                        
                        snID++;
                        loopID++;
                        nodeID+=4;
                        ndefects++;
                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"irradiation defects density="<<defectsDensity<<std::endl;
                    }
                }
            }
        }
        
        /**********************************************************************/
        void addIrradiationLoops()
        {
            if(targetIrradiationLoopDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Irradiation Loops"<<defaultColor<<std::endl;
                
                if(poly.crystalStructure=="BCC")
                {
                    addIrradiationLoopsBCC();
                }
                else if(poly.crystalStructure=="FCC")
                {
                    addIrradiationLoopsFCC();
                }
                else if(poly.crystalStructure=="HEX")
                {
                    addIrradiationLoopsHCP();
                }
                else
                {
                    std::cout<<"irradiation loops supported for FCC/BCC/HEX crystals only"<<std::endl;
                }
                
            }
        }
        
        /**********************************************************************/
        double deltaHelicity(const std::vector<VectorDimD>& newPoints,
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
                                                              //                                                              const PlaneMeshIntersectionContainerType& pp
                                                              const MeshBoundaryContainerType& pp)
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
                const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*(pp[m]->P0-P0));
                if(x(0)>FLT_EPSILON && x(0)<1.0-FLT_EPSILON && x(1)>FLT_EPSILON)
                {
                    temp.emplace(x(0),pp[m]->P0);
                }
                
                const Eigen::Matrix<double,2,1> y=llt.solve(A.transpose()*(pp[m]->P1-P0));
                if(y(0)>FLT_EPSILON && y(0)<1.0-FLT_EPSILON && y(1)>FLT_EPSILON)
                {
                    temp.emplace(y(0),pp[m]->P1);
                }
                
            }
            
            return temp;
            
            
        }
        
        /**********************************************************************/
        static std::pair<int,VectorDimD> boundaryProjection(const VectorDimD& P,
                                                            const VectorDimD& D,
                                                            //                                                            const PlaneMeshIntersectionContainerType& pp,
                                                            const MeshBoundaryContainerType& pp)
        {
            const double dNorm(D.norm());
            assert(dNorm>FLT_EPSILON);
            const VectorDimD dir=D/dNorm;
            // line1 is P+u1*dir, u>0
            
            bool success=false;
            std::pair<int,VectorDimD> temp=std::make_pair(-1,VectorDimD::Zero());
            
            for(size_t k=0;k<pp.size();++k)
            {
                //                const size_t k1 = ((k==(pp.size()-1))? 0 : k+1);
                const VectorDimD& v0=pp[k]->P0;
                const VectorDimD& v1=pp[k]->P1;
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
