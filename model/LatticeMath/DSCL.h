/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DSCL_h_
#define model_DSCL_h_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <SmithDecomposition.h>
#include <Lattice.h>
#include <LatticeTransitionMatrix.h>
#include <LLL.h>
#include <RLLL.h>


namespace model
{
    /*!Class template that computes the Displacement Shift Complete Lattice (DSCL)
     * of two parent lattices using the Smith decomposition method [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class DSCL : public LatticeTransitionMatrix<dim>
    /*       */ ,public Lattice<dim>
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
            // Therefore, the i-th primitive vectors of the DSCL is
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
            
            const MatrixDimD D1=A.latticeBasis*sd.matrixX().template cast<double>()*N.template cast<double>().inverse();
            const MatrixDimD D2=B.latticeBasis*sd.matrixV().template cast<double>()*M.template cast<double>().inverse();
            assert((D1-D2).norm()<FLT_EPSILON && "DSCL calculation failed.");
            
            if(useRLLL)
            {
                return RLLL(0.5*(D1+D2),0.75).reducedBasis();
            }
            else
            {
                return 0.5*(D1+D2);
            }
        }
        
    public:
        
//        const LatticeType& A;
//        const LatticeType& B;
        
        /**********************************************************************/
        DSCL(const LatticeType& A_in,
             const LatticeType& B_in,
             const bool& useRLLL=true) :
        /* init */ LatticeTransitionMatrix<dim>(A_in,B_in)
        /* init */,Lattice<dim>(getLatticeBasis(*this,A_in,B_in,useRLLL),MatrixDimD::Identity())
//        /* init */,A(A_in)
//        /* init */,B(B_in)
        {
            //            update(useRLLL);
            
            if(true)
            {//verify that A and B are multiples of DSCL
                
                const MatrixDimD tempA(this->reciprocalBasis.transpose()*A_in.latticeBasis);
                assert((tempA-tempA.array().round().matrix()).norm()<FLT_EPSILON && "A is not a multiple of DSCL");
                
                const MatrixDimD tempB(this->reciprocalBasis.transpose()*B_in.latticeBasis);
                assert((tempB-tempB.array().round().matrix()).norm()<FLT_EPSILON && "B is not a multiple of DSCL");
            }
            
        }
        
        
    };
    
    
}
#endif




//        /**********************************************************************/
//        const IntValueType& sigma() const
//        {
//            return _sigma;
//        }

//        IntValueType _sigma;

//        /**********************************************************************/
//        static IntValueType gcd(const IntValueType& a,const IntValueType& b)
//        {
//            return b>0? gcd(b, a % b) : a;
//        }
