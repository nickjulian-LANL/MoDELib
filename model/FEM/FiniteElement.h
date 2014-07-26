/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FiniteElement_H_
#define model_FiniteElement_H_

#include <deque>
#include <map>
#include <list>

#include <Eigen/Dense>

#include <model/Mesh/SimplicialMesh.h>
#include <model/Utilities/TerminalColors.h>
#include <model/Utilities/CompareVectorsByComponent.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/FEM/Elements/LagrangeElement.h>
#include <model/FEM/TrialFunction.h>
#include <model/FEM/WeakForms/LinearWeakForm.h>
#include <model/FEM/WeakForms/BilinearWeakForm.h>
#include <model/FEM/Domains/IntegrationDomain.h>
#include <model/FEM/Domains/EntireDomain.h>
#include <model/FEM/Boundaries/ExternalBoundary.h>
#include <model/FEM/Boundaries/NodeList.h>
#include <model/MPI/MPIcout.h>


namespace model
{

    template<typename _ElementType>
    class FiniteElement :
    /* inherits        */ public std::map<Eigen::Matrix<size_t,_ElementType::dim+1,1>, // key
    /*                                  */ _ElementType, // value
    /*                                  */ CompareVectorsByComponent<size_t,_ElementType::dim+1> // key compare
    /*                                  */ >,
    /* inherits        */ public std::deque<typename _ElementType::NodeType>,
    /* inherits        */ public std::map<Eigen::Matrix<double,_ElementType::dim,1>, // key
    /*                                  */ const typename _ElementType::NodeType* const, // value
    /*                                  */ CompareVectorsByComponent<double,_ElementType::dim,float> // key compare
    /*                                  */ >
    {
        
    private:
        Eigen::Matrix<double,_ElementType::dim,1> _xMin;
        Eigen::Matrix<double,_ElementType::dim,1> _xMax;
        
    public:
        
        typedef _ElementType ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef typename ElementType::NodeType NodeType;
        typedef std::list<const NodeType*> NodeListType;
        typedef std::map<size_t,NodeListType> NodeMapType;
        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        typedef SimplicialMesh<dim> MeshType;
        typedef std::map<Eigen::Matrix<size_t,_ElementType::dim+1,1>, // key
        /*                                  */ _ElementType, // value
        /*                                  */ CompareVectorsByComponent<size_t,_ElementType::dim+1> // key compare
        /*                                  */ > ElementContainerType;
        typedef std::deque<typename _ElementType::NodeType> NodeContainerType;
        typedef std::map<Eigen::Matrix<double,_ElementType::dim,1>, // key
        /*                                  */ const typename _ElementType::NodeType* const, // value
        /*                                  */ CompareVectorsByComponent<double,_ElementType::dim,float> // key compare
        /*                                  */ > NodeFinderType;
        
        
        const MeshType& mesh;
        
        /**********************************************************************/
        FiniteElement(const SimplicialMesh<dim>& m) :
        /* init list */ _xMin(Eigen::Matrix<double,ElementType::dim,1>::Constant( DBL_MAX)),
        /* init list */ _xMax(Eigen::Matrix<double,ElementType::dim,1>::Constant(-DBL_MAX)),
        /* init list */ mesh(m)
        {/*!@param[in] s A const reference to a SimplicialMesh on which *this 
          * FiniteElement is constructed.
          */
            
             model::cout<<greenColor<<"Creating FiniteElement:\n"<<defaultColor<<std::flush;
            
            // THIS IS NECESSARY TO AVOID "STATIC INITIALIZATION FIASCO"
             model::cout<<"Element barycentric coordinates:\n"<<ElementType::baryNodalCoordinates<<std::endl;
            
            
            // Insert elements
            for (typename SimplicialMesh<dim>::const_iterator eIter=mesh.begin();eIter!=mesh.end();++eIter)
            {
                auto temp=ElementContainerType::emplace(eIter->first,ElementType(eIter->second,*this,*this));
                assert(temp.second && "UNABLE TO INSERT ELEMENT.");
            }
            
            // Compute _xMin and _xMax
            for (int n=0;n<nodeSize();++n)
            {
                for(int d=0;d<dim;++d)
                {
                    if (node(n).p0(d)<_xMin(d))
                    {
                        _xMin(d)=node(n).p0(d);
                    }
                    if (node(n).p0(d)>_xMax(d))
                    {
                        _xMax(d)=node(n).p0(d);
                    }
                }
                
            }
            
             model::cout<<"   # elements: "<<elementSize()    <<"\n";
             model::cout<<"   # nodes: "   <<nodeSize()       <<"\n";
             model::cout<<"   xMin= "    <<_xMin.transpose()<<"\n";
             model::cout<<"   xMax= "    <<_xMax.transpose()<<std::endl;
        }
        
        /**********************************************************************/
        template <int nComponents>
        TrialFunction<nComponents,FiniteElementType> trial() const
        {
            return TrialFunction<nComponents,FiniteElementType>(*this);
        }
        
        /**********************************************************************/
        typename ElementContainerType::const_iterator elementBegin() const
        {
            return ElementContainerType::begin();
        }
        
        /**********************************************************************/
        typename ElementContainerType::const_iterator elementEnd() const
        {
            return ElementContainerType::end();
        }
        
        /**********************************************************************/
        size_t elementSize() const
        {/*!\returns the number of elements in the FiniteElement
          */
            return ElementContainerType::size();
        }
        
        /**********************************************************************/
        typename NodeContainerType::const_iterator nodeBegin() const
        {
            return NodeContainerType::begin();
        }
        
        /**********************************************************************/
        typename NodeContainerType::const_iterator nodeEnd() const
        {
            return NodeContainerType::end();
        }
        
        /**********************************************************************/
        const NodeContainerType& nodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        size_t nodeSize() const
        {
            return NodeContainerType::size();
        }
        
        /**********************************************************************/
        const NodeType& node(const size_t& n) const
        {/*!@param[in] n the n-th node stored in *this FiniteElement
          * \returns a reference to the n-th node in *this FiniteElement
          */
            return NodeContainerType::operator[](n);
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMin() const
        {
            return _xMin;
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMax() const
        {
            return _xMax;
        }
        
        /**********************************************************************/
        template <typename BndType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> boundary() const
        {
            return BndType::template boundary<FiniteElementType,qOrder,QuadratureRule>(*this);
        }
        
        /**********************************************************************/
        template <typename DomainType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        IntegrationDomain<FiniteElementType,0,qOrder,QuadratureRule> domain() const
        {
            return DomainType::template domain<FiniteElementType,qOrder,QuadratureRule>(*this);
        }
        
        /**********************************************************************/
        template <typename NodeSelectorType>
        NodeList<FiniteElementType> getNodeList() const
        {
            NodeList<FiniteElementType> temp(*this);
            const NodeSelectorType nodeSelector(*this);
            for(typename NodeContainerType::const_iterator nIter=nodeBegin(); nIter!=nodeEnd();++nIter)
            {
                if(nodeSelector(*nIter))
                {
                    temp.emplace_back(&(*nIter));
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        std::pair<bool,const ElementType*> search(const Eigen::Matrix<double,dim,1>& P) const
        {/*!@param[in] P position to search for
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          *
          * By default the search starts at this->begin()->second
          */
            return searchWithGuess(P,&(elementBegin()->second.simplex));
        }
        
        /**********************************************************************/
        std::pair<bool,const ElementType*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess) const
        {/*!@param[in] P position to search for
          * @param[in] guess Simplex* where the search starts
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(P,guess));
            const typename ElementContainerType::const_iterator eIter(ElementContainerType::find(temp.second->xID));
            assert(eIter!=elementEnd() && "ELEMENT NOT FOUND");
            return std::pair<bool,const ElementType*>(temp.first,&(eIter->second));
        }
        
    };
    
}	// close namespace
#endif