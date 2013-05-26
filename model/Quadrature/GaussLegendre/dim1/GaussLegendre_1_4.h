 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_1_4_H_ 
#define model_GAUSSLEGENDRE_1_4_H_ 

namespace model{

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
   template<>
   struct GaussLegendre<1,4>{
       static Eigen::Matrix<double,2,4> abcsissasAndWeights(){
           Eigen::Matrix<double,4,2> aw;
           aw<<6.943184420297371e-02, 1.739274225687270e-01, 
               3.300094782075718e-01, 3.260725774312732e-01, 
               6.699905217924280e-01, 3.260725774312732e-01, 
               9.305681557970260e-01, 1.739274225687266e-01; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 
