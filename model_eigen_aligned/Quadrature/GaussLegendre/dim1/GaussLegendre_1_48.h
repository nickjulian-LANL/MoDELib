 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
#ifndef model_GAUSSLEGENDRE_1_48_H_ 
#define model_GAUSSLEGENDRE_1_48_H_ 

namespace model
{

   template<>
   struct GaussLegendre<1,48>
   {
       static Eigen::Matrix<double,2,48> abcsissasAndWeights()
       {
           Eigen::Matrix<double,48,2> aw;
           aw<<6.144963737866882e-04, 1.576673026150361e-03, 
               3.234913866824563e-03, 3.663776950668462e-03, 
               7.937708138586907e-03, 5.738617289018582e-03, 
               1.470420372687653e-02, 7.789657862048529e-03, 
               2.350614841978443e-02, 9.808080228628958e-03, 
               3.430665464672256e-02, 1.178538042090235e-02, 
               4.706043164221535e-02, 1.371325485320402e-02, 
               6.171398986287585e-02, 1.558361391618433e-02, 
               7.820586918780331e-02, 1.738861128235488e-02, 
               9.646689798527897e-02, 1.912067553355732e-02, 
               1.164204837421297e-01, 2.077254147148780e-02, 
               1.379829345380928e-01, 2.233728042810313e-02, 
               1.610638101836677e-01, 2.380832924623162e-02, 
               1.855663016117429e-01, 2.517951777699943e-02, 
               2.113876369580137e-01, 2.644509474245263e-02, 
               2.384195126388833e-01, 2.759975184999257e-02, 
               2.665485476245206e-01, 2.863864605020110e-02, 
               2.956567590046416e-01, 2.955741984926540e-02, 
               3.256220568539195e-01, 3.035221958285786e-02, 
               3.563187563222723e-01, 3.101971157992198e-02, 
               3.876181048026555e-01, 3.155709614312747e-02, 
               4.193888219655542e-01, 3.196211929232415e-02, 
               4.514976503952686e-01, 3.223308221797508e-02, 
               4.838099145185653e-01, 3.236884840634200e-02, 
               5.161900854814346e-01, 3.236884840634192e-02, 
               5.485023496047312e-01, 3.223308221797496e-02, 
               5.806111780344457e-01, 3.196211929232399e-02, 
               6.123818951973444e-01, 3.155709614312687e-02, 
               6.436812436777277e-01, 3.101971157994614e-02, 
               6.743779431460803e-01, 3.035221958294705e-02, 
               7.043432409953586e-01, 2.955741984919762e-02, 
               7.334514523754794e-01, 2.863864605020167e-02, 
               7.615804873611168e-01, 2.759975184999290e-02, 
               7.886123630419861e-01, 2.644509474259696e-02, 
               8.144336983882565e-01, 2.517951777692672e-02, 
               8.389361898163321e-01, 2.380832924624538e-02, 
               8.620170654619079e-01, 2.233728042834729e-02, 
               8.835795162578710e-01, 2.077254147173219e-02, 
               9.035331020147223e-01, 1.912067553291489e-02, 
               9.217941308121976e-01, 1.738861128238572e-02, 
               9.382860101371249e-01, 1.558361391640075e-02, 
               9.529395683577853e-01, 1.371325485417946e-02, 
               9.656933453532768e-01, 1.178538041966096e-02, 
               9.764938515802153e-01, 9.808080228676221e-03, 
               9.852957962731234e-01, 7.789657861471936e-03, 
               9.920622918614134e-01, 5.738617289617233e-03, 
               9.967650861331753e-01, 3.663776950638680e-03, 
               9.993855036262138e-01, 1.576673026153069e-03; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 
