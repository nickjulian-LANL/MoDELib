/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SegmentSegmentDistance_H_
#define model_SegmentSegmentDistance_H_

#include <cfloat>
#include <deque>
#include <Eigen/Dense>

namespace model
{
    /*! Class template which computes the distance between two lines segments
     * in dimension dim. Algorithm adapted from
     * V. LUMELSKY. ON FAST COMPUTATION OF DISTANCE BETWEEN LINE SEGMENTS.
     * Information Processing Letters 21 (1985) 55-61.
     */
    template <int dim>
    struct SegmentSegmentDistance
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        typedef std::pair<double,VectorDim> PointPairType;
        typedef std::deque<std::pair<PointPairType,PointPairType>,Eigen::aligned_allocator<std::pair<PointPairType,PointPairType>>> IntersectionContainerType;
        
        static constexpr double tol=FLT_EPSILON;
        
        /**********************************************************************/
        static double limitRange(const double& v)
        {
            if(v<0.0)
            {
                return 0.0;
            }
            else if(v>1.0)
            {
                return 1.0;
            }
            else
            {
                return v;
            }
        }
        
        /**********************************************************************/
        std::pair<double,double> step4U(const double& T) const
        {/*! Compute t from eq. (10) and limit by eq. (12)
          */
            return std::make_pair(T,limitRange((T*R-S2)/D2));
        }
        
        /**********************************************************************/
        std::pair<double,double> step4(const double& U) const
        {/*! Compute t from eq. (10) and limit by eq. (12)
          */
            return std::make_pair(limitRange((U*R+S1)/D1),U);
        }
        
        /**********************************************************************/
        std::pair<double,double> step3(const double& T) const
        {
            const double U=(T*R-S2)/D2;       // compute U from eq. (10)
            if(U<0.0 || U>1.0)
            {// u is not in [0,1]
                return step4(limitRange(U));
            }
            else
            {// u is in [0,1]
                return std::make_pair(T,U);
            }
        }
        
        /**********************************************************************/
        std::pair<double,double> step2() const
        {
            const double T=limitRange((S1*D2-S2*R)/den); // compute t from eq. (11) and limit by eq. (12)
            return step3(T);
        }
        
        /**********************************************************************/
        std::pair<double,double> getTU() const
        {
            if(D1<tol && D2<tol)
            {// Step 1b: both segments are degenerate
                return std::make_pair(0.0,0.0);
            }
            else if(D1<tol && D2>=tol)
            {// Step 1a: first segment is degenerate
                return step4U(0.0);
            }
            else if(D1>=tol && D2<tol)
            {// Step 1a: second segment is degenerate
                return step4(0.0);
            }
            else
            {// both segments are not degenerate
                if(fabs(den)>tol)
                {// Step 1d: skew segments
                    return step2();
                }
                else
                {//  Step 1c: parallel or coincident segments
                    return step3(0.0);
                }
            }
        }
        
        
        const VectorDim d1;
        const VectorDim d2;
        const VectorDim d12;
        
        const double R;
        const double S1;
        const double S2;
        const double D1;
        const double D2;
        const double den;
        
        const std::pair<double,double> tu;
        const double& t;
        const double& u;
        const double dMin;
        const VectorDim x0;
        const VectorDim x1;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /**********************************************************************/
        SegmentSegmentDistance(const VectorDim& A,
                                   const VectorDim& B,
                                   const VectorDim& C,
                                   const VectorDim& D) :
        /* init */ d1(B-A),
        /* init */ d2(D-C),
        /* init */ d12(C-A),
        /* init */ R(d1.dot(d2)),
        /* init */ S1(d1.dot(d12)),
        /* init */ S2(d2.dot(d12)),
        /* init */ D1(d1.squaredNorm()),
        /* init */ D2(d2.squaredNorm()),
        /* init */ den(D1*D2-R*R),
        /* init */ tu(getTU()),
        /* init */ t(tu.first),
        /* init */ u(tu.second),
        /* init */ dMin((d1*t-d2*u-d12).norm()),
        /* init */ x0(A+t*d1),
        /* init */ x1(C+u*d2)
        {/*!\param[in] A start point of first segment
          * \param[in] B   end point of first segment
          * \param[in] C start point of second segment
          * \param[in] D   end point of second segment
          */
            
        }
        
//        /**********************************************************************/
//        IntersectionContainerType intersectionPoints(const double& dMax) const
//        {
//        
//            IntersectionContainerType temp;
//            
//            if(dMin<dMax)
//            {
//                
//                if(D1<tol && D2<tol)
//                {// Step 1b: both segments are degenerate
//                    const PointPairType i0=std::make_pair(t,x0);
//                    const PointPairType i1=std::make_pair(u,x1);
//                    temp.emplace_back(i0,i1);
//                }
//                else if(D1<tol && D2>=tol)
//                {// Step 1a: first segment is degenerate
//                    const PointPairType i0=std::make_pair(t,x0);
//                    const PointPairType i1=std::make_pair(u,x1);
//                    temp.emplace_back(i0,i1);
//                }
//                else if(D1>=tol && D2<tol)
//                {// Step 1a: second segment is degenerate
//                    const PointPairType i0=std::make_pair(t,x0);
//                    const PointPairType i1=std::make_pair(u,x1);
//                    temp.emplace_back(i0,i1);
//
//                }
//                else
//                {// both segments are not degenerate
//                    if(fabs(den)>tol)
//                    {// Step 1d: skew segments
//                        const PointPairType i0=std::make_pair(t,x0);
//                        const PointPairType i1=std::make_pair(u,x1);
//                        temp.emplace_back(i0,i1);
//                    }
//                    else
//                    {//  Step 1c: parallel or coincident segments
//                        return step3(0.0);
//                    }
//                }
//            
//            }
//            
//            return temp;
//        }
        
        
    };
    
}
#endif
