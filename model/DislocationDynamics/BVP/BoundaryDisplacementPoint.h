/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_BoundaryDisplacementPoint_h
#define _model_BoundaryDisplacementPoint_h

#include <assert.h>
#include <model/ParticleInteraction/FieldPoint.h>


namespace model
{
    
    /******************************************************************************/
    template<typename DislocationNetworkType>
    struct BoundaryDisplacementPoint :
    /* inheritance */ public FieldPoint<BoundaryDisplacementPoint<DislocationNetworkType>,DislocationNetworkType::DisplacementField::dim,typename DislocationNetworkType::DisplacementField>
    {
        typedef typename DislocationNetworkType::DisplacementField DisplacementField;
        typedef typename DislocationNetworkType::ElementType::NodeType NodeType;

        constexpr static int dim=DisplacementField::dim;
        
        const size_t gID;
        const Eigen::Matrix<double,dim,1> P;
        const Eigen::Matrix<double,dim,1> S;
        
        BoundaryDisplacementPoint(const NodeType& node) :
        /*   */ gID(node.gID),
        /*   */ P(node.P0),
        /*   */ S(node.outNormal())
        {
            if(std::abs(S.norm()-1.0)>=10.0*DBL_EPSILON)
            {
                model::cout<<"FEM node "<<gID<<"\n"<<"P="<<P.transpose()<<"\n"<<"S="<<S.transpose()<<"\n"<<"norm(S)="<<S.norm()<<std::endl;
                assert(0 && "S-vector has non-unit norm");
            }
        }
        
//        BoundaryDisplacementPoint(const size_t gID_in,
//                                  const Eigen::Matrix<double,dim,1>& P0,
//                                  const Eigen::Matrix<double,dim,1>& s_in) :
//        /*   */ gID(gID_in),
//        /*   */ P(P0),
//        /*   */ S(s_in)
//        {
//            assert(std::abs(S.norm()-1.0)>10.0*DBL_EPSILON && "S-vector has non-unit norm");
//        }
        
    };
    
} // end namespace
#endif
