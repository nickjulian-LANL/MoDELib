/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DipolarMicrostructureGenerator_H_
#define model_DipolarMicrostructureGenerator_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <model/DislocationDynamics/MicrostructureGeneration/MicrostructureGenerator.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/Mesh/PlaneMeshIntersection.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/IO/SequentialOutputFile.h>
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>


namespace model
{
    
    class DipolarMicrostructureGenerator : public MicrostructureGenerator
    {
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename PlaneMeshIntersection<dim>::PlaneMeshIntersectionContainerType PlaneMeshIntersectionContainerType;
        
        std::random_device rd;
        std::mt19937 generator;
        
//        SequentialOutputFile<'E',1> edgeFile;
//        SequentialOutputFile<'V',1> vertexFile;
//        SequentialOutputFile<'L',1> loopFile;
        
    public:
        DipolarMicrostructureGenerator(int argc, char* argv[]) :
        MicrostructureGenerator(argc,argv),
        /* init list */ generator(rd())
        {
            
            EigenDataReader EDR;
            
            double targetDensity=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","targetDensity",targetDensity);
            
            double fractionSessile=0.0;
            EDR.readScalarInFile("./microstructureInput.txt","fractionSessile",fractionSessile);
            
            std::cout<<"Generating dipoles..."<<std::endl;
            double density=0.0;
            size_t nodeID=0;
            size_t loopID=0;
            size_t snID=0;
            
            double helicity=0.0;
            
            while(density<targetDensity)
            {
                std::pair<LatticeVector<dim>,int> rp=this->randomPointInMesh();
                const int grainID=rp.second;
                
                const LatticeVector<dim> L0=rp.first;
                
                std::uniform_int_distribution<> distribution(0,this->poly.grain(grainID).slipSystems().size()-1);
                
                
                const int rSS=distribution(generator); // a random SlipSystem
                
                const SlipSystem& slipSystem=this->poly.grain(grainID).slipSystems()[rSS];
                
                std::set<int> planeIDs;
                for (unsigned int k=0;k<this->poly.grain(grainID).planeNormals().size();++k)
                {
                    if(slipSystem.s.dot(this->poly.grain(grainID).planeNormals()[k])==0)
                    {
                        planeIDs.insert(k);
                    }
                }
                assert(planeIDs.size()==2 && "ONLY FCC IS SUPPORTED AT THE MOMENT.");
                
                
                
                ReciprocalLatticeDirection<3> sr(this->poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));
                
                LatticeDirection<3> d1(LatticeVector<dim>(sr.cross(this->poly.grain(grainID).planeNormals()[*planeIDs.begin()]))*this->randomSign());
                LatticeDirection<3> d2(LatticeVector<dim>(sr.cross(this->poly.grain(grainID).planeNormals()[*planeIDs.rbegin()])*this->randomSign()));
                LatticeDirection<3> d3(slipSystem.s*this->randomSign());
                
                if(density/targetDensity<fractionSessile)
                { // overwrite d2
                    assert(0 && "SESSILE LOOPS NOT SUPPORTED YET.");
                }
                
                
                const double d1cNorm(d1.cartesian().norm());
                const double d2cNorm(d2.cartesian().norm());
                const double d3cNorm(d3.cartesian().norm());
                
                assert(d1cNorm>0.0);
                assert(d2cNorm>0.0);
                assert(d3cNorm>0.0);

                int a1=this->randomSize()/d1cNorm;
                int a2=this->randomSize()/d2cNorm;
                int a3=this->randomSize()/d3cNorm;
                
                assert(a1!=0);
                assert(a2!=0);
                assert(a3!=0);
                
                LatticeVector<dim> L1=L0+d1*a1;
                LatticeVector<dim> L2=L1+d2*a2;
                LatticeVector<dim> L3=L2-d1*a1;
                
                const VectorDimD P0=L0.cartesian();
                const VectorDimD P1=L1.cartesian();
                const VectorDimD P2=L2.cartesian();
                const VectorDimD P3=L3.cartesian();
                
                bool useHelicity=false;
                double deltaHelicity=0.0;
                std::deque<VectorDimD> tempPoints;
                if(useHelicity)
                {
                    
                    tempPoints.push_back(P0);
                    tempPoints.push_back(P1);
                    tempPoints.push_back(P2);
                    tempPoints.push_back(P3);
//                    this->loopPoints.emplace_back({P0,P1,P2,P3});
//                    this->loopBurgers.push_back(slipSystem.s.cartesian());
                    
                     deltaHelicity=this->deltaHelicity(tempPoints,slipSystem.s.cartesian());
                    
                }
                
                
                if(   mesh.searchRegion(grainID,P1).first
                   && mesh.searchRegion(grainID,P2).first
                   && mesh.searchRegion(grainID,P3).first
                   && fabs(helicity+deltaHelicity)>=fabs(helicity)
                   )
                {
                    
                    if(useHelicity)
                    {

                                            this->loopPoints.emplace_back(tempPoints);
                                            this->loopBurgers.push_back(slipSystem.s.cartesian());
                        helicity+=deltaHelicity;
                        
                    }
                    
                    const VectorDimD n1=d1.cross(slipSystem.s).cartesian().normalized();
                    const VectorDimD n2=d2.cross(slipSystem.s).cartesian().normalized();
                    
                    
                    
                    PlaneMeshIntersectionContainerType pmi01=PlaneMeshIntersection<dim>(this->mesh,P0,n1,grainID);
                    const VectorDimD P4=this->boundaryProjection(P0,d3.cartesian(),pmi01).second;
                    PlaneMeshIntersectionContainerType pmi12=PlaneMeshIntersection<dim>(this->mesh,P1,n2,grainID);
                    const VectorDimD P5=this->boundaryProjection(P1,d3.cartesian(),pmi12).second;
                    PlaneMeshIntersectionContainerType pmi23=PlaneMeshIntersection<dim>(this->mesh,P2,n1,grainID);
                    const VectorDimD P6=this->boundaryProjection(P2,d3.cartesian(),pmi23).second;
                    PlaneMeshIntersectionContainerType pmi30=PlaneMeshIntersection<dim>(this->mesh,P3,n2,grainID);
                    const VectorDimD P7=this->boundaryProjection(P3,d3.cartesian(),pmi30).second;
                    
                    const auto v54=this->boundaryProjection(P1,P0,d3.cartesian(),pmi01);
                    const auto v65=this->boundaryProjection(P2,P1,d3.cartesian(),pmi12);
                    const auto v76=this->boundaryProjection(P3,P2,d3.cartesian(),pmi23);
                    const auto v47=this->boundaryProjection(P0,P3,d3.cartesian(),pmi30);
                    
                    /*! Vertex file format is:
                     * ID Px Py Pz Vx Vy Vz velReducCoeff snID meshLocation grainID
                     */
                    const size_t refNodeID=nodeID;
                    this->nodesIO.emplace_back(refNodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    this->nodesIO.emplace_back(refNodeID+1,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    this->nodesIO.emplace_back(refNodeID+2,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    this->nodesIO.emplace_back(refNodeID+3,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    this->nodesIO.emplace_back(refNodeID+4,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    this->nodesIO.emplace_back(refNodeID+5,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    this->nodesIO.emplace_back(refNodeID+6,P6,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
                    this->nodesIO.emplace_back(refNodeID+7,P7,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
//                    vertexFile << refNodeID+0<<"\t" << std::setprecision(15)<<std::scientific<<P0.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+1<<"\t" << std::setprecision(15)<<std::scientific<<P1.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+2<<"\t" << std::setprecision(15)<<std::scientific<<P2.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+3<<"\t" << std::setprecision(15)<<std::scientific<<P3.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+4<<"\t" << std::setprecision(15)<<std::scientific<<P4.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+5<<"\t" << std::setprecision(15)<<std::scientific<<P5.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+6<<"\t" << std::setprecision(15)<<std::scientific<<P6.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                    vertexFile << refNodeID+7<<"\t" << std::setprecision(15)<<std::scientific<<P7.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
                    nodeID+=8;
                    
                    
                    writeLoop(P0,0,1,5,v54,4,slipSystem, n1,grainID,refNodeID,nodeID,loopID,snID);
                    writeLoop(P1,1,2,6,v65,5,slipSystem, n2,grainID,refNodeID,nodeID,loopID,snID);
                    writeLoop(P2,2,3,7,v76,6,slipSystem,-n1,grainID,refNodeID,nodeID,loopID,snID);
                    writeLoop(P3,3,0,4,v47,7,slipSystem,-n2,grainID,refNodeID,nodeID,loopID,snID);
                    
                    snID+=1;
                    density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/this->mesh.volume()/std::pow(Material<Isotropic>::b_real,2);
                    std::cout<<"density="<<density<<std::endl;
                    if(useHelicity)
                    {
                        std::cout<<"helicity="<<helicity<<std::endl;

                    }
                }
            }
            this->write();
        }
        
        
        /**********************************************************************/
        void writeLoop(const VectorDimD& P0,
                       const int& id0,
                       const int& id1,
                       const int& id2,
//                       const std::deque<std::pair<int,VectorDimD>>& bndVtx,
                       const std::map<double,VectorDimD>& bndVtx,
                       const int& id3,
                       const SlipSystem& slipSystem,
                       const VectorDimD& n,
                       const int& grainID,
                       const int& refNodeID,
                       size_t& nodeID,
                       size_t& loopID,
                       const size_t& snID)
        {
            // write Loop file
//            DislocationLoopIO<dim> dlIO(loopID+0, slipSystem.s.cartesian(),n,P0,grainID);
            //loopFile<< dlIO<<"\n";
            this->loopsIO.emplace_back(loopID+0, slipSystem.s.cartesian(),n,P0,grainID);

            
            // write to edge and node files
            this->edgesIO.emplace_back(loopID,refNodeID+id0,refNodeID+id1,0);
            this->edgesIO.emplace_back(loopID,refNodeID+id1,refNodeID+id2,0);

//            edgeFile << loopID<<"\t" << refNodeID+id0<<"\t"<< refNodeID+id1<<"\n";
//            edgeFile << loopID<<"\t" << refNodeID+id1<<"\t"<< refNodeID+id2<<"\n";
            
            size_t oldID=refNodeID+id2;
            for(const auto& pair : bndVtx)
            {
                this->nodesIO.emplace_back(nodeID,pair.second,Eigen::Matrix<double,1,3>::Zero(),1.0,snID,0);
//                vertexFile << nodeID<<"\t" << std::setprecision(15)<<std::scientific<<pair.second.transpose()<<"\t"<<Eigen::Matrix<double,1,3>::Zero()<<"\t"<<1.0<<"\t"<< snID <<"\t"<< 0<<"\n";
//                edgeFile << loopID<<"\t" << oldID<<"\t"<< nodeID<<"\n"; // CHANGE THIS
                this->edgesIO.emplace_back(loopID,oldID,nodeID,0);
                oldID=nodeID;
                nodeID++;
            }
            this->edgesIO.emplace_back(loopID,oldID,refNodeID+id3,0);
            this->edgesIO.emplace_back(loopID,refNodeID+id3,refNodeID+id0,0);

//            edgeFile << loopID<<"\t" <<    oldID<<"\t"<< refNodeID+id3<<"\n"; // CHANGE THIS
//            edgeFile << loopID<<"\t" << refNodeID+id3<<"\t"<< refNodeID+id0<<"\n";
            
            loopID+=1;
            
        }
        
        
    };
    
}
#endif