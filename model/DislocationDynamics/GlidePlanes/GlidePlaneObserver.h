/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GLIDEPLANEOBSERVER_H_
#define model_GLIDEPLANEOBSERVER_H_

#include <algorithm>
#include <map>
#include <tuple>
#include <memory> // std::shared_ptr (c++11)
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/Utilities/CompareVectorsByComponent.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>
#include <model/Geometry/PlanePlaneIntersection.h>


namespace model
{
    
    // Class Predeclaration
    template <typename LoopType>
    class GlidePlane;
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename LoopType>
    struct GlidePlaneObserver : private std::map<Eigen::Matrix<long int,LoopType::dim+2,1>,
    /*                                       */ const GlidePlane<LoopType>* const,
    /*                                       */ CompareVectorsByComponent<long int,LoopType::dim+2,long int>,
    /*                                       */ Eigen::aligned_allocator<std::pair<const Eigen::Matrix<long int,LoopType::dim+2,1>,const GlidePlane<LoopType>* const> > >,
    /*                       */ private std::map<std::pair<size_t,size_t>,PlanePlaneIntersection<LoopType::dim>>
    {
        
        
        static constexpr int dim=LoopType::dim;
        typedef GlidePlaneObserver<LoopType> GlidePlaneObserverType;
        typedef GlidePlane<LoopType> GlidePlaneType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim+2,1> GlidePlaneKeyType;
        typedef std::map<Eigen::Matrix<long int,LoopType::dim+2,1>,
        /*            */ const GlidePlane<LoopType>* const,
        /*            */ CompareVectorsByComponent<long int,LoopType::dim+2,long int>,
        /*            */ Eigen::aligned_allocator<std::pair<const Eigen::Matrix<long int,LoopType::dim+2,1>,const GlidePlane<LoopType>* const> > > GlidePlaneMapType;
        typedef std::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;
        typedef PlanePlaneIntersection<dim> PlanePlaneIntersectionType;
        typedef std::map<std::pair<size_t,size_t>,PlanePlaneIntersectionType> GlidePlaneIntersectionContainerType;
        
        
        
    public:
        
        /**********************************************************************/
        GlidePlaneIntersectionContainerType& glidePlaneIntersections()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GlidePlaneIntersectionContainerType& glidePlaneIntersections() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const PlanePlaneIntersectionType& glidePlaneIntersection(const GlidePlaneType* const p1,
                                                                 const GlidePlaneType* const p2)
        {/*@param[in] p1 first plane
          *@param[in] p2 second plane
          *\returns the infinite line of interseciton between the two planes, 
          * if already computed. Otherwise it computes the interseciton, stores it,
          * and returns it.
          */
            const auto key=std::make_pair(std::max(p1->sID,p2->sID),std::min(p1->sID,p2->sID));
            const auto iter=glidePlaneIntersections().find(key);
            if(iter==glidePlaneIntersections().end())
            {
                glidePlaneIntersections().emplace(key,PlanePlaneIntersectionType(p1->P.cartesian(),p1->n.cartesian(),p2->P.cartesian(),p2->n.cartesian()));
            }
            return glidePlaneIntersections().at(key);
        }
        
        /**********************************************************************/
        GlidePlaneMapType& glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GlidePlaneMapType& glidePlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        static GlidePlaneKeyType getGlidePlaneKey(const Grain<dim>& grain,
                                                  const VectorDimD& P,
                                                  const VectorDimD& N)
        {
            const LatticeVector<dim> p=grain.latticeVector(P);
            const ReciprocalLatticeDirection<dim> n=grain.reciprocalLatticeDirection(N);
            assert(n.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane key");
            return (GlidePlaneKeyType()<<grain.grainID,n,p.dot(n)).finished();
        }
        
        /**********************************************************************/
        std::shared_ptr<GlidePlaneType> sharedGlidePlane(const SimplicialMesh<dim>& mesh,
                                                         const Grain<dim>& grain,
                                                         const VectorDimD& P,
                                                         const VectorDimD& N)
        {
            const GlidePlaneKeyType key=getGlidePlaneKey(grain,P,N);
            const auto planeIter=glidePlanes().find(key);
            return (planeIter!=glidePlanes().end())? planeIter->second->loops().begin()->second->_glidePlane :
            /*                            */ std::make_shared<GlidePlaneType>(this,mesh,grain,P,N);
        }
        
        /**********************************************************************/
        void addGlidePlane(const GlidePlaneType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Adds pS to *this GLidePlane
          */
            const bool success=glidePlanes().emplace(pL->glidePlaneKey,pL).second;
            assert( success && "COULD NOT INSERT GLIDE PLANE POINTER IN GLIDE PLANE OBSERVER.");
        }
        
        /**********************************************************************/
        void removeGlidePlane(const GlidePlaneType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Removes pS from *this GLidePlane
          */
            const int success=glidePlanes().erase(pL->glidePlaneKey);
            assert(success==1 && "COULD NOT ERASE GLIDE PLANE POINTER FROM GLIDE PLANE OBSERVER.");
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const GlidePlaneObserverType& gpo)
        {
            for (const auto& glidePlane : gpo.glidePlanes())
            {
                os << (*glidePlane.second);
            }
            return os;
        }
        
        
        
        /**********************************************************************/
        static std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> planeSegmentIntersection(const VectorDimD& P0,
                                                                                                                                                const VectorDimD& N,
                                                                                                                                                const VectorDimD& v0,
                                                                                                                                                const VectorDimD& v1)
        {
            //            std::deque<VectorDim> temp;
            
            const double nNorm(N.norm());
            assert(nNorm>FLT_EPSILON);
            const VectorDimD n(N/nNorm);
            
            std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> temp;
            
            
            // check intersection of v0->v1 with plane
            // x=v0+u(v1-v0)
            // (x-P0).n=0
            // (v0+u(v1-v0)-P0).n=0
            // u=(P0-v0).n/(v1-v0).n;
            const double edgeNorm=(v1-v0).norm();
            if(edgeNorm<FLT_EPSILON)
            {
                
            }
            else
            {
                const double den=(v1-v0).dot(n);
                const double num=(P0-v0).dot(n);
                const double P0v0norm=(P0-v0).norm();
                
                const double numCheck= (P0v0norm<FLT_EPSILON)? 0.0 : num/P0v0norm;
                
                if (fabs(den/edgeNorm)>FLT_EPSILON)
                {
                    // edge intersects plane
                    const double u=num/den;
                    
                    if(fabs(u)<FLT_EPSILON)
                    {
                        temp.emplace_back(v0,v0);
                    }
                    else if (u>=FLT_EPSILON && u<=1.0-FLT_EPSILON)
                    {
                        const VectorDimD x((1.0-u)*v0 + u*v1);
                        temp.emplace_back(x,x);
                    }
                    else if (fabs(1.0-u)<FLT_EPSILON)
                    {
                        temp.emplace_back(v1,v1);
                    }
                    else
                    {// no roots
                        
                    }
                    
                }
                else
                {
                    if (fabs(numCheck)>FLT_EPSILON)
                    {// edge is parallel to plane, no intersection
                        
                    }
                    else
                    {// edge is coplanar
                        temp.emplace_back(v0,v1);
                        //                        temp.emplace_back(v1,&edge.child(1));
                    }
                }
            }
            
            
            
            
            return temp;
        }
        
        
    };
    
}	// close namespace
#endif










