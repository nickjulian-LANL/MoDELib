/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GAUSSLEGENDRE_H_
#define model_GAUSSLEGENDRE_H_

#include <assert.h>
#include <model/Quadrature/GaussLegendre/GaussLegendreRecursiveRelation.h>
#include <model/Quadrature/QuadratureRecursiveRelation.h>

namespace model
{
	
	/**************************************************/
	/* GaussLegendre: general case                    */
	/**************************************************/
	/*! \brief Class template defining the GaussLegendre rules for determination 
	 *	of quadrature abscissas and weights.
	 */
	template<short unsigned int dim, short unsigned int qOrder>
	struct GaussLegendre
    {
		
	};
    
    /* Template specialization for cases not included in includeGaussLegendre1D.h
     */
    template<short unsigned int qOrder>
    struct GaussLegendre<1,qOrder>
    {
        
        static Eigen::Matrix<double,2,qOrder> abcsissasAndWeights()
        {
            return QuadratureRecursiveRelation::template abcsissasAndWeights<GaussLegendreRecursiveRelation>(qOrder);
        }
        
    };
    
} // namespace model

//#include <model/Quadrature/GaussLegendre/includeGaussLegendre1D.h>
#include <model/Quadrature/GaussLegendre/includeGaussLegendre2D.h>
#include <model/Quadrature/GaussLegendre/includeGaussLegendre3D.h>

#endif

