/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SplineSegment_H_
#define model_SplineSegment_H_

#include <math.h>
#include <Eigen/Dense>



//#include <model/Geometry/Splines/SplineConsts.h>

//#include <model/Network/NetworkLink.h>
#include <model/LoopNetwork/NetworkLink.h>
#include <model/Geometry/ParametricCurve.h>
//#include <model/Geometry/Splines/Intersection/SplineIntersection.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/Geometry/Splines/SplineShapeFunction.h>
#include <model/Network/Operations/EdgeFinder.h>
//#include <Eigen/Sparse>
#include <model/Math/MatrixCompanion.h>

namespace model
{
    
    /*!
     A SplineSegment is a special kind of parametric curve where the position vector can be written
     as:
     \verbatim
     r(u)=N(u)*q
     \endverbatim
     where N(u) is a (dim x NdofXsegment) matrix of shape functions and q is a (NdofXsegment x 1) column
     vector of degrees of freedom, where NdofXsegment= dim x n and n=porder+1 is the number of shape functions:
     \verbatim
     | N1  0  0		Nn  0  0 |
     N= |  0 N1  0	...  0 Nn  0 |
     |  0  0 N1		 0  0 Nn |
     \endverbatim
     
     Each shape function is a porder:
     \verbatim
     [N1 ... Nn]= [1 u ... n^porder] * SFC
     \endverbatim
     where SFC is a n x n matrix of shape function coefficients that in general dependes on the knot vector. It's posible
     to take advantage of this particular form of shape functions to calculate ruu and ru therefore these functions are
     redefined here (virtual in ParametricCurve) and the flow chart of their implementation is shown below:
     \verbatim
     make_UPOWuu()	  make_UPOWu()	make_UPOW()
     ^				  ^				^
     |				  |				|
     make_SFuu()	  make_SFu()	make_SF()
     ^				  ^				^
     |				  |				|
     make_Nuu()		  make_Nu()		make_N()
     ^				  ^				^
     |				  |				|
     make_ruu()		  make_ru()		make_r()	(this level is virtual in ParametricCurve)
     \endverbatim
     */
    
    
    
    
    /************************************************************************/
    /* SplineSegment, general case **************************************/
    /************************************************************************/
    template <typename Derived, short unsigned int dim,short unsigned int corder>
    class SplineSegment : public NetworkLink<Derived>,
    /*                 */ public ParametricCurve<Derived,dim>
    {
        
        static_assert(dim>=1 && dim <=3,"DIMENSION MUST BE 1, 2 or 3."); // requires c++11
        static_assert(corder>=0 && corder <=2,"CONTINUITY ORDER MUST BE 0, 1 or 2."); // requires c++11
        
        
    public:
        //enum  {dim=_dim};
        static constexpr int Ncoeff= 2*(corder+1);
        static constexpr int pOrder= 2*corder+1;
        static constexpr int Ndof  = dim*Ncoeff;
        static constexpr int eigenSize=pOrder*pOrder;
        
        typedef Eigen::Matrix<double, 1, Ncoeff> RowNcoeff;
        typedef Eigen::Matrix<double, 1, Ncoeff-1> RowNcoeffu;
        typedef Eigen::Matrix<double, 1, Ncoeff-2> RowNcoeffuu;
        typedef Eigen::Matrix<double, Ncoeff, Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double, Ncoeff, dim> MatrixNcoeffDim;
        typedef Eigen::Matrix<double, Ndof,1> VectorNdof;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        typedef Eigen::Matrix<double, dim, Ndof> MatrixDimNdof;
        
        typedef Eigen::Matrix<double,Ndof,Ndof>	MatrixNdof;
        
        
        typedef typename NetworkLink<Derived>::NodeType NodeType;
        typedef typename NetworkLink<Derived>::LinkType LinkType;
        typedef typename NetworkLink<Derived>::FlowType FlowType;
        
        //typedef SplineSegment<Derived,dim,corder,alpha,qOrder,QuadratureRule> SplineSegmentType;
        //typedef SplineSegment<Derived,dim,corder,alpha> SplineSegmentType;
        //typedef SplineSegment<Derived,dim,corder> SplineSegmentType;
        
        typedef SplineShapeFunction<dim,corder> SplineShapeFunctionType;
        typedef ParametricCurve<Derived,dim> ParametricCurveType;
        
        //#include<model/Geometry/Splines/SplineEnums.h>
        
        
        
        //protected:
        //
        //Eigen::Matrix<double,qOrder,Ncoeff> SFgauss;
        
        //        int sourceTfactor; // LEAVE THIS UNINITIALIZED: this is calculated in TopologyChangeActions of source node, which happens before constructor of this
        //        int sinkTfactor;   // LEAVE THIS UNINITIALIZED: this is calculated in TopologyChangeActions of   sink node, which happens before constructor of this
        //
        //
        /**********************************************************************/
        VectorDim sourceT() const
        {
            VectorDim temp=VectorDim::Zero();
            int n=0;
            for(const auto& link : this->loopLinks())
            {
                if(this->source->sID==link->source()->sID && this->sink->sID==link->sink()->sID)
                {
                    temp+=this->source->tangents().at(link->loop()->sID);
                    n++;
                }
                else if(this->source->sID==link->sink()->sID && this->sink->sID==link->source()->sID)
                {
                    temp-=this->source->tangents().at(link->loop()->sID);
                    n++;
                }
                else
                {
                    assert(0);
                }
            }
            
            return n==0? temp : temp/n;
            
//            return this->source->get_T()*sourceTfactor;
        }
        
        /**********************************************************************/
        VectorDim sinkT() const
        {
            VectorDim temp=VectorDim::Zero();
            int n=0;
            for(const auto& link : this->loopLinks())
            {
                if(this->source->sID==link->source()->sID && this->sink->sID==link->sink()->sID)
                {
                    temp+=this->sink->tangents().at(link->loop()->sID);
                    n++;
                }
                else if(this->source->sID==link->sink()->sID && this->sink->sID==link->source()->sID)
                {
                    temp-=this->sink->tangents().at(link->loop()->sID);
                    n++;
                }
                else
                {
                    assert(0);
                }
            }
            
            return n==0? temp : temp/n;
            //AFTER INTRODUCING THE ENERGY CRITERION CHANGE THIS IN A -1
            //return -this->sink->get_T()*sinkTfactor;
            //return this->sink->get_T()*sinkTfactor;
        }
        
    public:
        
        static double alpha;
        
        
        /******************************************************************************/
        SplineSegment(const std::shared_ptr<NodeType>& nI,
                      const std::shared_ptr<NodeType>& nJ) :
        NetworkLink<Derived>(nI,nJ)
//        sourceTfactor(1),
//        sinkTfactor(-1)
        {/*! Constructor with Nodes and flow
          */
            
        }
        
        /**********************************************************************/
        MatrixNcoeffDim get_qH() const
        {/*!\returns The matrix of Hermite dof of this spline segment.
          *  [P0x P0y P0z;T0x T0y T0z;P1x P1y P1z;T1x T1y T1z]
          */
            return (MatrixNcoeffDim()<< this->source->get_P().transpose(),
                    /*            */	sourceT().transpose(),
                    /*            */	this->  sink->get_P().transpose(),
                    /*            */	sinkT().transpose()).finished();
        }
        
        
        /******************************************************************************/
        VectorDim chord() const
        {/*!\returns the chord vector (source -> sink)
          */
            return this->sink->get_P()-this->source->get_P();
        }
        
        /******************************************************************************/
        double chordLength() const
        {/*!\returns the length of the chord vector
          */
            return chord().norm();
        }
        
        /**********************************************************************/
        double parametricChordLength() const
        {//!\returns  the length of the chord vector to the power alpha
            return std::pow(chordLength(),alpha);
        }
        
        
        
        /******************************************************************************/
        VectorDim get_r(const double & u) const
        {/*!\returns The position vector at parameter u
          *  @param[in] u the parametrization variable in [0:1]
          *	\f[
          *		\mathbf{r} = \mathbf{q}\mathbf{H}\mathbf{u} \\
          *		r_i = q_{ik}H_{km}u^{m} = q_{ik} N_k
          *	\f]
          * with i=0..dim-1, k,m = 0... Ncoeff
          * ACTUALLY IN THE CODE WE HAVE THE TRANSPOSE OF THIS !!!!
          */
            return SplineShapeFunctionType::sf(u,parametricChordLength())*get_qH();
        }
        
        /******************************************************************************/
        VectorDim get_ru(const double & uin) const
        {
            return SplineShapeFunctionType::sfDiff1(uin,parametricChordLength())*get_qH();
        }
        
        /******************************************************************************/
        VectorDim get_rmu(const double & uin) const
        {
            return this->get_ru(uin)/parametricChordLength();
        }
        
        /******************************************************************************/
        VectorDim get_ruu(const double & uin) const
        {
            return SplineShapeFunctionType::sfDiff2(uin,parametricChordLength())*get_qH();
        }
        
        /******************************************************************************/
        VectorDim get_rmumu(const double & uin) const
        {
            return this->get_ruu(uin)/std::pow(parametricChordLength(),2);
        }
        
        //        /******************************************************************************/
        //        Eigen::Matrix<double, Ndof, Eigen::Dynamic>  get_G2H() const
        //        {
        //            //make_G2H();
        //
        //            size_t gDof(this->pSN()->nodeOrder()*this->source->NdofXnode); // CHANGE HERE, NdofXnode should be available directly
        //
        //            Eigen::Matrix<double, Ndof, Eigen::Dynamic> G2H(Eigen::Matrix<double, Ndof, Eigen::Dynamic>::Zero(Ndof,gDof));
        //
        //            //G2H.setZero(Ndof,gDof);
        //
        //            Eigen::VectorXi dofid(this->source->dofID());
        //            Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> M(this->source->W2H());
        //            M.block(dim,0,dim,M.cols())*=sourceTfactor;
        //            //	std::cout<<"M source=\n"<<M<<std::endl;
        //
        //            for (int k=0;k<dofid.size();++k)
        //            {
        //                G2H.template block<Ndof/2,1>(0,dofid(k))=M.col(k);
        //            }
        //
        //            dofid=this->sink->dofID();
        //            M=this->sink->W2H();
        //            M.block(dim,0,dim,M.cols())*=(-sinkTfactor);
        //            //	std::cout<<"M sink=\n"<<M<<std::endl;
        //
        //            for (int k=0;k<dofid.size();++k)
        //            {
        //                G2H.template block<Ndof/2,1>(Ndof/2,dofid(k))=M.col(k);
        //            }
        //
        //            return G2H;
        //        }
        
        /******************************************************************************/
        std::pair<double,std::pair<double,VectorDim> > closestPoint(const VectorDim& P0) const
        {/*!@param[in] P0 reference point
          * \returns The closesest point to P0 along this segment. The return value is a
          * pair, where pair.first is the parameter value u, and pair.second is the position P
          * of the closest point, so that get_r(u)=P.
          *
          *
          */
            
            // solve (P-P0)*dP/du=0
            
            // The polynomial coefficients of this spine segment
            Eigen::Matrix<double,dim,Ncoeff> coeffs(this->polynomialCoeff());
            coeffs.col(0)-=P0;
            
            
            // The derivative of the polynomial coefficients
            Eigen::Matrix<double,dim,Ncoeff-1> dcoeffs(Eigen::Matrix<double,dim,Ncoeff-1>::Zero());
            for (int i=0;i<Ncoeff-1;++i)
            {
                dcoeffs.col(i)=(i+1)*coeffs.col(i+1);
            }
            
            Eigen::Matrix<double,1,2*Ncoeff-2> pcoeffs(Eigen::Matrix<double,1,2*Ncoeff-2>::Zero()); // degree of product = pOrder+(pOrder-1)=2*pOrder-1. nCoeffs of product = 2*pOrder-1+1= 2*pOrder = 2*Ncoeff-2
            
            // The polynomial coefficients of (P-P0)*dP/du, in reverse order
            for (int i=0;i<Ncoeff;++i)
            {
                for (int j=0;j<Ncoeff-1;++j)
                {
                    pcoeffs(2*Ncoeff-3-i-j) += coeffs.col(i).dot(dcoeffs.col(j));
                }
            }
            
            // Compute roots using the eigenvalue method
            MatrixCompanion<2*Ncoeff-3> mc(pcoeffs);
            
            // sort roots according to distance to P0
            std::map<double,std::pair<double,VectorDim> > rootMap;
            
            //    for (int k=0;k<2*Ncoeff-3;++k)
            for (size_t k=0;k<mc.rootSize;++k)
            {
                if (std::fabs(mc.root(k).imag())<FLT_EPSILON && mc.root(k).real()>0.0 && mc.root(k).real()<1.0 )
                {
                    
                    VectorDim P(this->get_r(mc.root(k).real()));
                    rootMap.insert(std::make_pair((P-P0).norm(), std::make_pair(mc.root(k).real(),P) ));
                    
                }
                
            }
            
            // check distance to limits of interval
            rootMap.insert(std::make_pair((this->source->get_P()-P0).norm(), std::make_pair(0.0,this->source->get_P()) ));
            rootMap.insert(std::make_pair((this->  sink->get_P()-P0).norm(), std::make_pair(1.0,this->  sink->get_P()) ));
            
            return *rootMap.begin();
            
        }
        
    };
    
    
    //static data
    template <typename Derived, short unsigned int dim,short unsigned int corder>
    double SplineSegment<Derived,dim,corder>::alpha=0.5;
    
    
    
}
#endif

