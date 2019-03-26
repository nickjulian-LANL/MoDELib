/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationInjector_H_
#define model_DislocationInjector_H_




namespace model
{
    template <int dim>
    class DislocationInjector
    {
//        constexpr static int dim=DislocationNetworkType::dim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
  //      DislocationNetworkType& DN;
        
    public:
        
//        DislocationInjector(DislocationNetworkType& DN_in) :
//        DN(DN_in)
//        {
//
//        }
        
        /**********************************************************************/
        static std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> straightLineBoundaryClosure(const VectorDim& P0,
                                                                                                       const VectorDim& d,
                                                                                                       const VectorDim& n,
                                                                                                       const int& grainID,
                                                                                                       const SimplicialMesh<dim>& mesh)
        {
            
            
            const double maxSize=std::max(mesh.xMax(0)-mesh.xMin(0),std::max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));

            
            // Define line AB containing dislocaiton and piercing the mesh
            const VectorDim A=P0+3.0*maxSize*d;
            const VectorDim B=P0-3.0*maxSize*d;
            
            MeshPlane<dim> plane(mesh,grainID,P0,n);
//            const auto& segDeq(plane.meshIntersections);
            
            std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> nodePos;
            int nIntersections=0;
            for(const auto& pair : plane.meshIntersections)
            {
                
                SegmentSegmentDistance<dim> ssi(A,B,pair.P0,pair.P1);
                
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
                    if(ssi.dMin>FLT_EPSILON)
                    {// no intersection
                        if((pair.P0-nodePos.back()).norm()>FLT_EPSILON)
                        {
                            nodePos.push_back(pair.P0);
                        }
                    }
                    else
                    {// intersection
                        const VectorDim X((ssi.x0+ssi.x1)*0.5);
                        if((X-nodePos.back()).norm()>FLT_EPSILON)
                        {
                            nIntersections++;
                            nodePos.push_back(pair.P0);
                            if((X-pair.P0).norm()>FLT_EPSILON)
                            {
                                nodePos.push_back(X);
                            }
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
        template <typename DislocationNetworkType>
        void insertRandomStraightDislocation(DislocationNetworkType& DN)
        {
            const std::pair<LatticeVector<dim>,int> rp=DN.poly.randomLatticePointInMesh();
            const int& grainID=rp.second;   // random grain ID
            const LatticeVector<dim>& L0=rp.first; // random lattice position in the grain
            const VectorDim P0(L0.cartesian());   // cartesian position of L0
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            std::uniform_int_distribution<> distribution(0,DN.poly.grain(grainID).slipSystems().size()-1);
            const int rSS=distribution(generator); // a random SlipSystem ID
            const SlipSystem& slipSystem=DN.poly.grain(grainID).slipSystems()[rSS];
            const VectorDim b=slipSystem.s.cartesian();    // Burgers vector
            const VectorDim n=slipSystem.n.cartesian().normalized(); // slip plane normal

            std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
            const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
            const VectorDim d=Eigen::AngleAxisd(theta, n)*b.normalized();

//            const double maxSize=std::max(DN.mesh.xMax(0)-DN.mesh.xMin(0),std::max(DN.mesh.xMax(1)-DN.mesh.xMin(1),DN.mesh.xMax(2)-DN.mesh.xMin(2)));
//
//
//            const VectorDim A=P0+3.0*maxSize*d;
//            const VectorDim B=P0-3.0*maxSize*d;
//
//            // Compute interseciton between mesh and glide plane
//            //PlaneMeshIntersection<dim> pmi(DN.mesh,P0,n,grainID);
//            MeshPlane<dim> pmi(DN.mesh,grainID,P0,n);
//
//            std::deque<std::pair<VectorDim,VectorDim>> segDeq;
//
//            for(size_t k=0;k<pmi.size();++k)
//            {
//                const int k1=(k+1)<pmi.size()? k+1 :0;
//                segDeq.emplace_back(pmi[k].second,pmi[k1].second);
//            }
//            int nIntersections=0;
//            for(const auto& pair : segDeq)
//            {
//                SegmentSegmentDistance<dim> ssi(A,B,pair.first,pair.second);
//
//                if(nIntersections==0)
//                {
//                    if(ssi.dMin>FLT_EPSILON) // no intersection
//                    {
//
//                    }
//                    else //if(ssi.size==1)
//                    {
//                        nIntersections++;
//                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
//                    }
//                }
//                else if(nIntersections==1)
//                {
//                    if(ssi.dMin>FLT_EPSILON) // no intersection
//                    {
//                        nodePos.push_back(pair.first);
//                    }
//                    else //if(ssi.size==1)
//                    {
//                        nIntersections++;
//                        nodePos.push_back(pair.first);
//                        nodePos.push_back((ssi.x0+ssi.x1)*0.5);
//                    }
//                }
//                else
//                {
//
//                }
//            }
            
            std::deque<VectorDim> nodePos(straightLineBoundaryClosure(P0,d,n,grainID,DN.mesh));
            std::vector<size_t> nodeIDs;
            for (const auto& node : nodePos)
            {
                nodeIDs.push_back(DN.insertDanglingNode(node,VectorDim::Zero(),1.0).first->first);
            }
            DN.insertLoop(nodeIDs,b,n,P0,grainID);
            DN.clearDanglingNodes();
            
            
        }
        
    };
    
}
#endif
