/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNode_h_
#define model_DislocationNode_h_

//#include <PlanarDislocationNode.h>

#include <DislocationDynamicsModule.h>

#ifndef NDEBUG
#define VerboseDislocationNode(N,x) if(this->network().verboseDislocationNode>=N){std::cout<<x;}
#else
#define VerboseDislocationNode(N,x)
#endif

namespace model
{
    
    template <int dim, short unsigned int corder>
    class DislocationNode : public NetworkNode<DislocationNode<dim,corder>>
    /*                   */,public SplineNode<DislocationNode<dim,corder>,dim,corder,Hermite>
    /*                   */,public ConfinedDislocationObject<dim>
    {
        
        public:
        
        typedef TypeTraits<DislocationLoop<dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef SplineNode<DislocationNode<dim,corder>,dim,corder,Hermite> SplineNodeType;
        typedef ConfinedDislocationObject<dim> ConfinedDislocationObjectType;
        typedef typename TypeTraits<NetworkNodeType>::MeshLocation MeshLocation;
        typedef std::vector<VectorDim> VectorOfNormalsType;

        // static bool use_velocityFilter;
        // static double velocityReductionFactor;
        // static int verboseDislocationNode;

        const Simplex<dim,dim>* p_Simplex;
        VectorDim velocity;
        VectorDim vOld;         //! The previous velocity vector of *this PlanarDislocationNode
        double velocityReductionCoeff;
        std::shared_ptr<NetworkNodeType> virtualNode;
        NetworkNodeType* const masterNode;

        
        DislocationNode(LoopNetworkType* const,const VectorDim&,const VectorDim&,const double&);
        ~DislocationNode();
        std::shared_ptr<DislocationNode> clone() const;
        const Simplex<dim,dim>* get_includingSimplex(const VectorDim&,const Simplex<dim,dim>* const) const;
        const Simplex<dim,dim>* includingSimplex() const;
        bool isMovableTo(const VectorDim&) const;
        const double& velocityReduction() const;
        bool isVirtualBoundaryNode() const __attribute__ ((deprecated));
        bool isBoundaryNode() const;
        bool isGrainBoundaryNode() const;
        void updateGeometry();
        bool set_P(const VectorDim&);
        bool trySet_P(const VectorDim&);
        const VectorDim& get_V() const;
        // const VectorDim& get_P() const;
        MeshLocation meshLocation() const;
        void set_V(const VectorDim& vNew);
        void projectVelocity();
        // size_t uniqueLoopNodes() const;
        void addLoopNode(LoopNodeType* const);
        void removeLoopNode(LoopNodeType* const);
        // bool isRemovable(const double& Lmin,const double& relAreaThIn);
        VectorDim invariantDirectionOfMotion() const;
        const std::shared_ptr<NetworkNodeType>& virtualBoundaryNode() const;
        // static void initFromFile(const std::string&);
        // bool isZeroBurgersNode() const;

    };
    
}
#endif