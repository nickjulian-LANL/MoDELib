/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EvaluationPoint_H_
#define model_EvaluationPoint_H_

#include <Eigen/Dense>

namespace model
{
    
 
    /******************************************************************************/
    template<int rows,int cols>
    struct EvaluationPoint : public Eigen::Matrix<double,dim,cols>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        const size_t pointID;
        const VectorDim P;
        
        /**********************************************************************/
        DisplacementPoint(const size_t& _pointID,const VectorDim& _P) :
        /* init */ Eigen::Matrix<double,dim,cols>(Eigen::Matrix<double,dim,cols>::Zero()),
        /* init */ pointID(_pointID),
        /* init */ P(_P)
        {
        }
        
    };

}
#endif
