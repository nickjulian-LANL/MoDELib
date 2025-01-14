/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CSL_h_
#define model_CSL_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <SmithDecomposition.h>
#include <Lattice.h>
#include <LatticeTransitionMatrix.h>
#include <LLL.h>
#include <RLLL.h>


namespace model
{
    /*!Class template that computes the Coincident Site Lattice (CSL) of two
     * parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class CSL : public LatticeTransitionMatrix<dim>
    /*       */,public Lattice<dim>
    {
        static_assert(dim>0,"dim must be > 0.");
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Lattice<dim> LatticeType;
        typedef long long int IntValueType;
        typedef Eigen::Matrix<IntValueType,dim,dim> MatrixInt;
        
        /**********************************************************************/
        static MatrixDimD getLatticeBasis(const LatticeTransitionMatrix<dim>& rm,
                                          const LatticeType& A,
                                          const LatticeType& B,
                                          const bool& useRLLL)
        {
            // The transition matrix is T=P/sigma, where P=rm.integerMatrix is
            // an integer matrix and sigma=rm.sigma is an integer
            // The integer matrix P can be decomposed as P=X*D*Y using the Smith decomposition.
            // X and Y are unimodular, and D is diagonal with D(k,k) dividing D(k+1,k+1)
            // The decomposition also computes the matices U and V such that D=U*P*V
            SmithDecomposition<dim> sd(rm.integerMatrix);
            
            // From T=inv(A)*B=P/sigma=X*D*Y/sigma=X*D*inv(V)/sigma, we have
            // B1*(sigma*I)=A1*D
            // where
            // B1=B*V
            // A1=A*X
            // Since V and X are unimodular matrices, B1 and A1 are new bases
            // of the lattices B and A, respectively. Moreover, since
            // (sigma*I) and D are diagonal, the columns of B1 and A1 are
            // proportional, with rational proportionality factors different for each column.
            // For example, columns "i" read
            // b1_i*sigma=a1_i*D(i,i)
            // Therefore, the i-th primitive vectors of the CSL is
            // c_i=b1_i*sigma/gcd(sigma,D(i,i))=a1_i*D(i,i)/gcd(sigma,D(i,i))
            // or, in matrix form
            // C=B1*N=A1*M, that is
            // C=B*V*N=A*X*M
            // where M=diag(D(i,i)/gcd(sigma,D(i,i))) and
            //       N=diag(sigma/gcd(sigma,D(i,i))) and
            
            MatrixInt M(MatrixInt::Identity());
            MatrixInt N(MatrixInt::Identity());
            for(int i=0;i<dim;++i)
            {
                const IntValueType& dii=sd.matrixD()(i,i);
                M(i,i)=dii/rm.gcd(rm.sigma,dii);
                N(i,i)=rm.sigma/rm.gcd(rm.sigma,dii);
            }
            
            const MatrixDimD C1=A.latticeBasis*(sd.matrixX()*M).template cast<double>();
            const MatrixDimD C2=B.latticeBasis*(sd.matrixV()*N).template cast<double>();
            assert((C1-C2).norm()<FLT_EPSILON && "CSL calculation failed.");
            
            if(useRLLL)
            {
                return RLLL(0.5*(C1+C2),0.75).reducedBasis();
            }
            else
            {
                return 0.5*(C1+C2);
            }
        }
        
    public:
        
//        const LatticeType& A;
//        const LatticeType& B;
        
        /**********************************************************************/
        CSL(const LatticeType& A_in,
            const LatticeType& B_in,
            const bool& useRLLL=true) :
        /* init */ LatticeTransitionMatrix<dim>(A_in,B_in)
        /* init */,Lattice<dim>(getLatticeBasis(*this,A_in,B_in,useRLLL),MatrixDimD::Identity())
//        /* init */,A(A_in)
//        /* init */,B(B_in)
        {
            
            if(true)
            {//verify that CSL can be obtained as multiple of A and B
                
                const MatrixDimD tempA(A_in.reciprocalBasis.transpose()*this->latticeBasis);
                assert((tempA-tempA.array().round().matrix()).norm()<FLT_EPSILON && "CSL is not a multiple of A");
                
                const MatrixDimD tempB(B_in.reciprocalBasis.transpose()*this->latticeBasis);
                assert((tempB-tempB.array().round().matrix()).norm()<FLT_EPSILON && "CSL is not a multiple of B");
            }
            //            update(useRLLL);
        }
        
        
    };
    
    
} // end namespace
#endif


//
//
//        /**********************************************************************/
//        const IntValueType& sigma() const
//        {
////            return _sigma;
//            return rationalMatrix().den;
//        }
//
//        const RationalMatrix<dim>& transitionMatrix() const
//        {
//            return *this;
//        }

//        IntValueType _sigma;

//        /**********************************************************************/
//        static IntValueType gcd(const IntValueType& a,const IntValueType& b)
//        {
//            return b>0? gcd(b, a % b) : a;
//        }


//        /**********************************************************************/
//        static RationalMatrix<dim> getTransitionMatrix(const LatticeType& A,
//                                                     const LatticeType& B)
//        {
//            const MatrixDimD R(B.covBasis()*A.contraBasis().transpose());
//
//            // Check that R is a proper rotation
//            const MatrixDimD RRT=R*R.transpose();
//            const double RRTmInorm=(RRT-Eigen::Matrix<double,dim,dim>::Identity()).norm()/Eigen::Matrix<double,dim,dim>::Identity().norm();
//            if(RRTmInorm>FLT_EPSILON)
//            {
//                std::cout<<"R="<<std::endl<<R<<std::endl;
//                std::cout<<"R*R^T="<<std::endl<<RRT<<std::endl;
//                std::cout<<"norm(R*R^T-I)/norm(I)="<<RRTmInorm<<", tol="<<FLT_EPSILON<<std::endl;
//                assert(0 && "R IS NOT ORTHOGONAL.");
//            }
//            // make sure that C2G is proper
//            assert(std::fabs(R.determinant()-1.0) < FLT_EPSILON && "R IS NOT PROPER.");
//
//
//            // Compute the transition matrix T=inv(A)*B
//            const MatrixDimD T=A.contraBasis().transpose()*B.covBasis();
//
//            // For the two lattices to have coincident sites, R must be a rational matrix
//            // Compute the integer matrix P and the integer sigma such that T=P/sigma
//            return RationalMatrix<dim> rm(T);
//        }
