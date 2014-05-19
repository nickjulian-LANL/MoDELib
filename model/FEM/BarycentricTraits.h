/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_BarycentricTraits_H_
#define  model_BarycentricTraits_H_

#include <set>
#include <Eigen/Dense>

namespace model {
	
    
    /**************************************************************************/
	/**************************************************************************/
    /*!Class template that performs auxiliary barycentric-to-standard
     * coordinate transformation in a Simplex.
     */
	template<short int dim>
	struct BarycentricTraits
    {
        static_assert(dim>0,"dim must be > 0");

        
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef Eigen::Matrix<double,dim+1,dim+1> MatrixHigherDimD;
        typedef Eigen::Matrix<double,dim+1,dim> MatrixHigherDimDimD;
        typedef Eigen::Matrix<double,dim,dim+1> MatrixDimHigherDimD;

        typedef Eigen::Matrix<double,dim  ,1> VectorDimD;
        typedef Eigen::Matrix<double,dim+1,1> VectorHigherDimD;

        static const MatrixHigherDimD L2X;
        static const MatrixHigherDimD X2L;
        static const MatrixHigherDimDimD dLdX;
        static const MatrixDimHigherDimD NdA;
        
        /**********************************************************************/
        static VectorHigherDimD x2l(const VectorDimD& X)
        {/*!\param[in] X vector of standard coordinates
          *\returns the barycentric coordinates corresponding to X
          */
            return X2L* (VectorHigherDimD()<< X, 1.0).finished();
        }
        
        /**********************************************************************/
        static VectorHigherDimD x2l(const double& X)
        {/*!\param[in] X vector of standard coordinates
          *\returns the barycentric coordinates corresponding to X
          */
            assert(dim==1);
            return X2L* (VectorHigherDimD()<< X, 1.0).finished();
        }
        
        /**********************************************************************/
        static VectorDimD l2x(const VectorHigherDimD& L)
        {/*!\param[in] L vector of barycentric coordinates
             *\returns the vector of standard coordinates corresponding to L
             */
            return (L2X* L).template segment<dim>(0);
        }

    private:
        
        /**********************************************************************/
        static MatrixHigherDimD get_L2X()
        {/*!/returns the barycentric-to-standard coordinate transformation matrix.
          *
          * For example, in 2d, the transformation is:
          *	\f[
          *		\left[\begin{array}{lll}
          *     v_{11}&v_{21}&v_{31}\\
          *     v_{12}&v_{22}&v_{32}\\
          *     1&1&1\\
          *     \end{array}
          *     \right]
          *     \left[\begin{array}{l}
          *     \lambda_1\\
          *     \lambda_2\\
          *     \lambda_3\\
          *     \end{array}
          *     \right]=
          *     \left[\begin{array}{c}
          *     s_1\\
          *     s_2\\
          *     1\\
          *     \end{array}
          *     \right]
          *	\f]
          * 
          * In the sandard domain, we now adopt the convention that, for i<n, 
          * the i-th vertex is on the i-th axis at a unit distance from the 
          * origin, while the n-th vertex is the origin. With this convention,
          * plugging in the coordinates of the standard vertices, we obtain:
          *	\f[
          *		\underbrace{\left[\begin{array}{lll}
          *     1&0&0\\
          *     0&1&0\\
          *     1&1&1\\
          *     \end{array}
          *     \right]}_{L2X}
          *     \left[\begin{array}{l}
          *     \lambda_1\\
          *     \lambda_2\\
          *     \lambda_3\\
          *     \end{array}
          *     \right]=
          *     \left[\begin{array}{c}
          *     s_1\\
          *     s_2\\
          *     1\\
          *     \end{array}
          *     \right]
          *	\f]
          */
            MatrixHigherDimD temp(MatrixHigherDimD::Identity());
            temp.template block<1,dim>(dim,0)=Eigen::Matrix<double,1,dim>::Ones();
            return temp;
        }
        
        /**********************************************************************/
        static MatrixHigherDimD get_X2L()
        {/*!/returns the standard-2-barycentric coordinate transformation matrix,
          * that is the inverse ofget_L2X().
          */
            MatrixHigherDimD temp(MatrixHigherDimD::Identity());
            temp.template block<1,dim>(dim,0)= -Eigen::Matrix<double,1,dim>::Ones();
            return temp;
        }
        
        /**********************************************************************/
        static MatrixDimHigherDimD get_NdA()
        {/*!/returns a matrix having in columns the non-unitary outward normals 
          * in the standard simplex. The norm of the normal is the ratio between 
          * the area of the face to the area of a dim-1 Simplex.
          */
            MatrixHigherDimDimD temp(-get_L2X().template block<dim+1,dim>(0,0));
            for (int r=0;r<dim;++r)
            {
                temp.row(r).normalize(); // norm of first dim-1 normals is unitary
            }
            temp.row(dim)*=-1.0; // flip direction of the inclined face normal
            return temp.transpose();
        }
        
	};
    
    // Declare static data member
    template<short int dim>
    const typename BarycentricTraits<dim>::MatrixHigherDimD BarycentricTraits<dim>::L2X=BarycentricTraits<dim>::get_L2X();
    
	template<short int dim>
    const typename BarycentricTraits<dim>::MatrixHigherDimD BarycentricTraits<dim>::X2L=BarycentricTraits<dim>::get_X2L();

    template<short int dim>
    const typename BarycentricTraits<dim>::MatrixHigherDimDimD BarycentricTraits<dim>::dLdX=BarycentricTraits<dim>::X2L. template block<dim+1,dim>(0,0);
    
    template<short int dim>
    const typename BarycentricTraits<dim>::MatrixDimHigherDimD BarycentricTraits<dim>::NdA=BarycentricTraits<dim>::get_NdA();

	
}
#endif
