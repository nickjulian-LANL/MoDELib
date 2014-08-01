/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakList_H_
#define model_LinearWeakList_H_

#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/TerminalColors.h>
#include <model/FEM/TrialOperators/TestExpression.h>
#include <model/FEM/WeakForms/LinearWeakExpression.h>
#include <model/FEM/WeakForms/JGNselector.h>
#include <model/FEM/Domains/IntegrationList.h>


namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename DerivedTest, typename _IntegrationListType>
	class LinearWeakList : public LinearWeakExpression<LinearWeakList<DerivedTest,_IntegrationListType> >
    {
        
    public:
        
        typedef typename DerivedTest::TrialFunctionType TrialFunctionType;
        
        constexpr static int dim=TrialFunctionType::dim;
        
    private:
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        typedef Eigen::Matrix<double,dofPerElement,1> ElementVectorType;
        
        
        
    public:
        
        const DerivedTest& testExp;
        const _IntegrationListType& list;
        
        /**********************************************************************/
        LinearWeakList(const TestExpression<DerivedTest>& _testExp,
                       const _IntegrationListType& _list) :
        /*init list */ testExp(_testExp.trial()),
        /*init list */ list(_list)
        {
            model::cout<<greenColor<<"Creating LinearWeakList "<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        size_t gSize() const
        {
            return testExp.gSize();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            
            model::cout<<"Assembling LinearWeakList (list size="<<list.size()<<") ..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            Eigen::Matrix<double,Eigen::Dynamic,1> _globalVector(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(gSize()));
            
            
            for (int k=0;k<list.size();++k)
            {
                const ElementType& ele(list[k].ele);  // element
                const int f(list[k].boundaryFace); //    face ID
                const Eigen::Matrix<double,dim+1,1>& bary(list[k].domainBary);
                
                const ElementVectorType ve(testExp.sfm(ele,bary).transpose()
                                           *list[k]()
                                           *JGNselector<3>::jGN(ele.jGN(bary,f))
                                           *list[k].weight
                                           );
                
                for (int i=0;i<dofPerElement;++i)
                {
                    const size_t  nodeID_I(i/dofPerNode);
                    const size_t nodeDof_I(i%dofPerNode);
                    const size_t gI= ele.node(nodeID_I).gID*dofPerNode+nodeDof_I;
                    _globalVector(gI) += ve(i);
                }
            }
            model::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            
            return _globalVector;
        }
        
    };
    
    
    //Operator ,
    template<typename T,int dimMinusDomainDim, typename PointType>
    LinearWeakList<T,IntegrationList<dimMinusDomainDim,PointType> > operator,(const TestExpression<T>& testExpr,const IntegrationList<dimMinusDomainDim,PointType>& list)
    {
        return LinearWeakList<T,IntegrationList<dimMinusDomainDim,PointType> >(testExpr,list);
    }
    
}	// close namespace
#endif
