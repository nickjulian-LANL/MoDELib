/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialDef_H_
#define model_TrialDef_H_

#include <utility> // for std::move
#include <type_traits> //static_assert
#include <model/Utilities/TypeTraits.h>
#include <model/Mesh/Simplex.h>
#include <model/FEM/TrialOperators/TrialBase.h>
#include <model/FEM/TrialOperators/TrialExpressionBase.h>
#include <model/FEM/TrialOperators/ExpressionRef.h>

namespace model
{
    
	
    /**************************************************************************/
	/**************************************************************************/
    template <typename T>
    struct TrialDef : //public TrialBase<typename T::TrialFunctionType>,
    /*            */ public TrialExpressionBase<TrialDef<T> >
    {
        
        typedef typename T::TrialFunctionType TrialFunctionType;
        typedef TrialBase<TrialFunctionType> TrialBaseType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;


        static_assert(TypeTraits<TrialFunctionType>::dim==TypeTraits<TrialFunctionType>::nComponents,"SYMMETRIC GRADIENT (DEF) CAN ONLY BE COMPUTED IF nComponents==dim.");
        
//        const T& op;
        ExpressionRef<T> op;

        
    //public:
        
        constexpr static int rows=(TrialFunctionType::dim*(TrialFunctionType::dim+1))/2;

        /**********************************************************************/
        TrialDef(const T& x) :
//        /* base init */ TrialBaseType(x.derived().trial()),
        /* init list */ op(x)
        {
            
        }
        
        /**********************************************************************/
        TrialDef(T&& x) :
//        /* base init */ TrialBaseType(x.derived().trial()),
        /* init list */ op(std::move(x))
        {
            
        }
        
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,dofPerElement> sfm(const ElementType& ele,
                                                     const BaryType& bary) const
        {
            return op().sfmDef(ele,bary);
        }
        
        //        /**********************************************************************/
        //        void sfmGrad(const typename TF::ElementType& ele, const Eigen::Matrix<double,TF::dim+1,1>& bary) const
        //        {
        //            static_assert(0,"SECOND GRADIENTS ARE NOT SUPPORTED YET.");
        //        }
        
        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const ElementType& ele,
//                                                const BaryType& bary) const
//        {/*!@param[in] ele the element
//          * @param[in] bary the vector of barycentric coordinates
//          * \returns the value of the Derived expression at bary.
//          *
//          * \todo: in order to be optimized, this function should be Derived-specific
//          */
//            return sfm(ele,bary)*this->trial().dofs(ele);
//        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P,
//                                                const Simplex<dim,dim>* guess) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            const std::pair<bool,const ElementType*> temp=this->trial().fe.searchWithGuess(P,guess);
//            //            return (temp.first? sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second)) : Eigen::Matrix<double,rows,1>::Zero());
//            Eigen::Matrix<double,rows,1> val(Eigen::Matrix<double,rows,1>::Zero());
//            if(temp.first)
//            {
//                val=sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second));
//            }
//            return val;
//        }

    };
    
    
    /**************************************************************************/
    template <typename T>
    TrialDef<T> def(const TrialExpressionBase<T>& x)
    {
        return TrialDef<T>(x);
    }
    
    /**************************************************************************/
    template <typename T>
    TrialDef<T> def(TrialExpressionBase<T>&& x)
    {
        return TrialDef<T>(std::move(x));
    }
    
}	// close namespace
#endif


//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            return operator()(P,&(op.trial.fe.mesh.begin()->second));
//
//        }
