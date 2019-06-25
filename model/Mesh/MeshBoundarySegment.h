/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshBoundarySegment_H_
#define model_MeshBoundarySegment_H_

#include <cfloat>
#include <set>
#include <Eigen/Dense>
#include <FiniteLineSegment.h>
#include <PlanarMeshFace.h>
#include <LineLineIntersection.h>

namespace model
{
    /*!\brief Class template representing a straight line at the intersection
     * of PlanarMeshFaces
     */
    template <int dim>
    struct MeshBoundarySegment : public FiniteLineSegment<dim>
    //    /*                        */,public std::set<const PlanarMeshFace<dim>*>
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef FiniteLineSegment<dim> FiniteLineSegmentType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFace<dim>*> FaceContainerType;
        
        
        static VectorDim getBoundaryNormal(const FaceContainerType& fcs)
        {
            assert(fcs.size() && "EMPY FACE CONTAINER");
            VectorDim temp(VectorDim::Zero());
            for(const auto& face : fcs)
            {
                temp+=face->outNormal();
            }
            const double tempNorm(temp.norm());
            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
        }
        
        FaceContainerType faces;
        const VectorDim boundaryNormal;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        //        const PlanarMeshFaceType* const face;
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const PlanarMeshFaceType* const face) :
        /* init */ FiniteLineSegmentType(p0,p1)
        /* init */,faces(FaceContainerType({face}))
        /* init */,boundaryNormal(getBoundaryNormal(faces))
        //        /* init */,face(face_in)
        {
            //            assert(face.size() && "EMPY FACE CONTAINER");
        }
        
        MeshBoundarySegment(const VectorDim& p0,
                            const VectorDim& p1,
                            const FaceContainerType& faces_in) :
        /* init */ FiniteLineSegmentType(p0,p1)
        /* init */,faces(faces_in)
        /* init */,boundaryNormal(getBoundaryNormal(faces))
        {
            //            assert(faces().size() && "EMPY FACE CONTAINER");
        }
        
        //        FaceContainerType& faces()
        //        {
        //            return *this;
        //
        //        }
        //
        //        const FaceContainerType& faces() const
        //        {
        //            return *this;
        //
        //        }
        
        //        /**********************************************************************/
        //        VectorDim boundaryNormal() const
        //        {
        //            VectorDim temp(VectorDim::Zero());
        //            for(const auto& face : faces())
        //            {
        //                temp+=face->outNormal();
        //            }
        //            const double tempNorm(temp.norm());
        //            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
        //        }
        
        
        
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const MeshBoundarySegment<dim>& seg)
        {
            os<<seg.faces.size()<<", "<<std::setprecision(15)<<std::scientific<<seg.P0.transpose()<<","<<seg.P1.transpose();
            return os;
        }
        
        
    };
    
    template <int dim>
    struct BoundingMeshSegments : public std::vector<MeshBoundarySegment<dim>, Eigen::aligned_allocator<MeshBoundarySegment<dim>>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::set<Eigen::Matrix<double,dim,1>,CompareVectorsByComponent<double,dim,float> > UniquePointContainer;

        
        /**********************************************************************/
        void emplace_unique(const VectorDim& P0,const VectorDim& P1,const PlanarMeshFace<dim>* const face)
        {
            if((P0-P1).norm()>FLT_EPSILON)
            {// ignore degenerate segments
                bool found(false);
                for(auto& seg : *this)
                {
                    if(   ((seg.P0-P0).norm()<FLT_EPSILON && (seg.P1-P1).norm()<FLT_EPSILON)
                       || ((seg.P0-P1).norm()<FLT_EPSILON && (seg.P1-P0).norm()<FLT_EPSILON)
                       )
                    {// coincident segments on different faces, need to merge faces
                        seg.faces.insert(face);
                        found=true;
                    }
                }
                if (!found)
                {
                    this->emplace_back(P0,P1,face);
                }
            }
        }
        
        /**********************************************************************/
        BoundingMeshSegments()
        {
            
        }
        
        /**********************************************************************/
        BoundingMeshSegments(const SimplicialMesh<dim>& mesh,
                             const int& rID,
                             const Plane<dim>& plane)
        {
            //std::cout<<"I'm here A"<<std::endl;
            const MatrixDim R(plane.localRotationMatrix());
            
            for(const auto& face : mesh.region(rID)->faces())
            {
                PlanePlaneIntersection<dim> ppi(plane,face.second->asPlane());
                //std::cout<<"I'm here A1"<<std::endl;

                if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                {// plane and mesh-face are coincident
                    
                    UniquePointContainer roots;
                    for(size_t k=0;k<face.second->convexHull().size();++k)
                    {
                        //std::cout<<"I'm here A2"<<std::endl;

                        // Pick face boundary segment
                        const size_t k1(k==face.second->convexHull().size()-1? 0 : k+1);
                        const VectorDim& P0(face.second->convexHull()[k]->P0);
                        const VectorDim& P1(face.second->convexHull()[k1]->P0);
                        const double segLength((P1-P0).norm());
                        const VectorDim D0((P1-P0)/segLength);
                        // Compute intersection between bonudary segment and line of intersection of the two planes
                        LineLineIntersection<dim> lli(P0,D0,ppi.P,ppi.d);
                        if(lli.type==LineLineIntersection<dim>::INCIDENT)
                        {
                            //std::cout<<"I'm here A3"<<std::endl;

                            const double u0((lli.x0-P0).dot(D0));
                            if(u0>0.0-FLT_EPSILON && u0<segLength+FLT_EPSILON)
                            {// intersection within segment
                                roots.insert(lli.x0);
                            }
                        }
                        else if(lli.type==LineLineIntersection<dim>::COINCIDENT)
                        {// a coincident line was found, which means that the glide planes intersec on a boundary face
                            //std::cout<<"I'm here A4"<<std::endl;

                            roots.insert(P0);
                            roots.insert(P1);
                        }
                    }
                    //std::cout<<"I'm here A5"<<std::endl;

                    switch (roots.size())
                    {
                        case 0:
                        {// no intersaction between plane and face
                            break;
                        }
                            
                        case 2:
                        {// an intersection segment
                            const VectorDim& P0(*roots.begin());
                            const VectorDim& P1(*roots.rbegin());
                            emplace_unique(P0,P1,face.second.get());
                            break;
                        }
                            
                        default:
                        {
                            std::cout<<"FAILED TO FIND A BOUNDARY PERIMETER FOR PLANE"<<std::endl;
                            std::cout<<"plane.P="<<plane.P.transpose()<<std::endl;
                            std::cout<<"plane.unitNormal="<<plane.unitNormal.transpose()<<std::endl;
                            std::cout<<"IN INTERSECTING PLANE AND FACE"<<std::endl;
                            std::cout<<face.second->outNormal()<<std::endl;
                            std::cout<<"ROOTS ARE"<<std::endl;
                            for(const auto& root : roots)
                            {
                                std::cout<<root.transpose()<<std::endl;
                            }
                            break;
                        }
                    }
                    //std::cout<<"I'm here A6"<<std::endl;

                }
                else if (ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                {
//                    this->clear();
                    //std::cout<<"I'm here A7"<<std::endl;

                    for(size_t k=0;k<face.second->convexHull().size();++k)
                    {
                        //std::cout<<"I'm here A8"<<std::endl;

                        const size_t k1(k==face.second->convexHull().size()-1? 0 : k+1);
                        const VectorDim& P0(face.second->convexHull()[k]->P0);
                        const VectorDim& P1(face.second->convexHull()[k1]->P0);
                        emplace_unique(P0,P1,face.second.get());
                    }
//                    break; // loop over faces
                }
            }
            
                        //std::cout<<"I'm here B"<<std::endl;
            
            // Now sort segments
            ConvexHull<2,MeshBoundarySegment<dim>> finalHull;
            //std::cout<<"unsorted hull"<<std::endl;
            for(const auto& pt : *this)
            {
                
                VectorDim x(R*(0.5*(pt.P0+pt.P1)-plane.P));
                finalHull.emplace(std::array<double,2>{x[0],x[1]},&pt);
                // THE PROBLEM HERE IS THAT IF COINCIDENT POINTS FROM DIFFERENCE FACES EXIST, THEN ONLY ONE OF THEM IS KEPT. E.G. A PLANE CUTTING AT THE INTERSECTION OF TWO FACES. IF WE HAD UNIQUE FACE EDGES WITH POINTERS TO THE ADJECENT FACES WE COULD SOLVE THIS
                //std::cout<<pt.P0.transpose()<<","<<pt.P1.transpose()<<std::endl;
            }
            const auto hullPts=finalHull.getPoints();
            
//            if(hullPts.size()!=temp.size())
//            {
//
//                std::cout<<"plane.P="<<plane.P.transpose()<<std::endl;
//                std::cout<<"plane.n="<<plane.unitNormal.transpose()<<std::endl;
//
//                std::cout<<"temp:\n";
//                for(const auto& pt : temp)
//                {
//                    VectorDim x(R*(0.5*(pt.P0+pt.P1)-plane.P));
//                    std::cout<<x.transpose()<<" "<<pt.P0.transpose()<<" "<<pt.P1.transpose()<<"\n";
//                }
//
//                std::cout<<"finalHull:\n";
//                for(const auto& pt : finalHull)
//                {
//                    std::cout<<pt[0]<<" "<<pt[1]<<"\n";
//                }
//
//                std::cout<<"hullPts:\n";
//                for(const auto& pt : hullPts)
//                {
//                    std::cout<<pt[0]<<" "<<pt[1]<<"\n";
//                }
//
//                assert(hullPts.size()==temp.size());
//            }
            
                        //std::cout<<"I'm here C"<<std::endl;
            BoundingMeshSegments<dim> sortedTemp;
            //            for(const auto& seg : hullPts)
            //            {
            //                sortedTemp.push_back(*seg.t);
            //            }
            //std::cout<<"sorted hull"<<std::endl;
            for(size_t k=0;k<hullPts.size();++k)
            {
                const auto& seg(*hullPts[k].t);
                if(k==0)
                {
                    const auto& seg1(*hullPts[1].t);
                    
                    if((seg.P1-seg1.P0).norm()<FLT_EPSILON || (seg.P1-seg1.P1).norm()<FLT_EPSILON)
                    {// first segment is oriented in the right direction
                        sortedTemp.push_back(seg);
                    }
                    else
                    {
                        sortedTemp.emplace_back(seg.P1,seg.P0,seg.faces);
                    }
                }
                else
                {
                    if((sortedTemp.back().P1-seg.P0).norm()<FLT_EPSILON)
                    {
                        sortedTemp.push_back(seg);
                    }
                    else if((sortedTemp.back().P1-seg.P1).norm()<FLT_EPSILON)
                    {
                        sortedTemp.emplace_back(seg.P1,seg.P0,seg.faces);
                    }
                    else
                    {
                        //std::cout<<"k="<<k<<" of "<<hullPts.size()<<std::endl;
                        assert(false && "DISCONNECTED FACE BOUNDARY");
                    }
                }
                //std::cout<<sortedTemp.back().P0.transpose()<<","<<sortedTemp.back().P1.transpose()<<std::endl;
                
            }
            assert((sortedTemp.back().P1-sortedTemp.front().P0).norm()<FLT_EPSILON && "OPEN FACE BOUNDARY");
            
            assert(sortedTemp.size()==this->size());
            
                        //std::cout<<"I'm here D"<<std::endl;
            this->swap(sortedTemp);
//            return sortedTemp;
            
        }
        
        
        //        BoundingMeshSegments(const BoundingMeshSegments<dim>& bb1,const BoundingMeshSegments<dim>& bb2)
        //        {
        //            assert(0 && "FINISH HERE");
        //        }
        
        /**********************************************************************/
        std::set<const MeshBoundarySegment<dim>*> containingSegments(const VectorDim& P) const
        {
            std::set<const MeshBoundarySegment<dim>*> temp;
            
            for(const auto& seg : *this)
            {
                //                std::cout<<pair.first<<", d="<<pair.second.distanceTo(P)<<std::endl;
                if(seg.contains(P))
                {
                    temp.insert(&seg);
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool contains(const VectorDim& P) const
        {
            return containingSegments(P).size();
        }
        
        /**********************************************************************/
        const BoundingMeshSegments<dim>& boundingBoxSegments() const
        {
            return *this;
        }
        
        BoundingMeshSegments<dim>& boundingBoxSegments()
        {
            return *this;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const BoundingMeshSegments<dim>& bls)
        {
            for(const auto& seg : bls)
            {
                os<<seg<<std::endl;
            }
            return os;
        }
        
    };
    
    
    //    template<int dim>
    //    struct MeshBoundarySegmentIntersection : public std::unique_ptr<MeshBoundarySegment<dim>>
    //    {
    //
    //
    //        MeshBoundarySegmentIntersection(const MeshBoundarySegment<dim>& s1,const MeshBoundarySegment<dim>& s2)
    //        {// Constructs the intersection of two MeshBoundarySegment(s)
    //            const SegmentSegmentDistance<dim> ssd(s1.P0,s1.P1,s2.P0,s2.P1);
    //            const auto iSeg(ssd.intersectionSegment());
    //
    //            switch (iSeg.size())
    //            {
    //                    std::set<const PlanarMeshFace<dim>*> allFaces;
    //                    for(const auto& tup : iSeg)
    //                    {// Make sure that each face cointais end points, and add to allFaces
    //                        const auto& x(std::get<0>(tup));
    //                        for(const auto& face : s1.faces)
    //                        {
    //                            assert(face->asPlane().contains(x) && "FACE DOES NOT CONTAIN INTERSECTION POINT");
    //                            allFaces.insert(face);
    //                        }
    //                        for(const auto& face : s2.faces)
    //                        {
    //                            assert(face->asPlane().contains(x) && "FACE DOES NOT CONTAIN INTERSECTION POINT");
    //                            allFaces.insert(face);
    //                        }
    //                    }
    //
    //                case 1:
    //                {// Single intersection point. This point belongs must be common to all faces of the origina segments
    //                    const auto& x(std::get<0>(iSeg[0]));
    //                    this->reset(new MeshBoundarySegment<dim>(x,x,allFaces));
    //                    break;
    //                }
    //
    //                case 2:
    //                {// extended intersection segment
    //                    const auto& x0(std::get<0>(iSeg[0]));
    //                    const auto& x1(std::get<0>(iSeg[1]));
    //                    this->reset(new MeshBoundarySegment<dim>(x0,x1,allFaces));
    //                    break;
    //                }
    //
    //                default:
    //                    break;
    //            }
    //
    //
    //
    //        }
    //
    //        MeshBoundarySegmentIntersection(const MeshBoundarySegment<dim>& seg,const PlanarMeshFace<dim>& face)
    //        {// Constructs the intersection of a MeshBoundarySegment and a PlanarMeshFace
    //
    //            const PlaneSegmentIntersection<dim> psi(face.asPlane(),seg);
    //            if(psi.type==PlaneSegmentIntersection<dim>::COINCIDENT || psi.type==PlaneSegmentIntersection<dim>::INCIDENT)
    //            {
    //                std::set<const PlanarMeshFace<dim>*> allFaces(seg.faces);
    //                allFaces.insert(face);
    //                this->reset(new MeshBoundarySegment<dim>(psi.x0,psi.x1,allFaces));
    //            }
    ////
    ////            switch (psi.type)
    ////            {
    ////                case PlaneSegmentIntersection::COINCIDENT:
    ////                {
    ////                    std::set<const PlanarMeshFace<dim>*> allFaces(seg.faces);
    ////                    allFaces.insert(face);
    ////                    this->reset(new MeshBoundarySegment<dim>(psi.x0,psi.x1,allFaces));
    ////                    break;
    ////                }
    ////
    ////                case PlaneSegmentIntersection::INCIDENT:
    ////                {
    ////                    this->reset(new MeshBoundarySegment<dim>(psi.x0,psi.x1,allFaces));
    ////                    break;
    ////                }
    ////
    ////                default:
    ////                    break;
    ////            }
    //        }
    //
    //    };
    
    
}
#endif

