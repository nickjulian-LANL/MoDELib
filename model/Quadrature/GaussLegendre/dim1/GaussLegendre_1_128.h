 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_1_128_H_ 
#define model_GAUSSLEGENDRE_1_128_H_ 

namespace model{

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
   template<>
   struct GaussLegendre<1,128>{
       static Eigen::Matrix<double,2,128> abcsissasAndWeights(){
           Eigen::Matrix<double,128,2> aw;
           aw<<8.755602643395477e-05, 2.246904801448384e-04, 
               4.612700113116874e-04, 5.229063396697922e-04, 
               1.133375687242311e-03, 8.212515093328697e-04, 
               2.103620732509359e-03, 1.119144215481767e-03, 
               3.371443549892772e-03, 1.416375735730765e-03, 
               4.936090754132261e-03, 1.712763020455615e-03, 
               6.796628637706248e-03, 2.008127491869299e-03, 
               8.951945782139925e-03, 2.302292128351538e-03, 
               1.140075426804532e-02, 2.595080916337736e-03, 
               1.414159062643100e-02, 2.886318771432751e-03, 
               1.717281678401689e-02, 3.175831580853050e-03, 
               2.049262107314925e-02, 3.463446283451159e-03, 
               2.409901932936709e-02, 3.748990962817218e-03, 
               2.798985608488957e-02, 4.032294945247378e-03, 
               3.216280586104153e-02, 4.313188899304265e-03, 
               3.661537456052577e-02, 4.591504935830539e-03, 
               4.134490095951909e-02, 4.867076707502472e-03, 
               4.634855829912110e-02, 5.139739507836143e-03, 
               5.162335597542089e-02, 5.409330369853157e-03, 
               5.716614132730097e-02, 5.675688162035219e-03, 
               6.297360152098386e-02, 5.938653686370326e-03, 
               6.904226553022552e-02, 6.198069771976567e-03, 
               7.536850621101537e-02, 6.453781369634732e-03, 
               8.194854246954619e-02, 6.705635644309276e-03, 
               8.877844152217779e-02, 6.953482066477246e-03, 
               9.585412124604314e-02, 7.197172502085450e-03, 
               1.031713526189030e-01, 7.436561301076757e-03, 
               1.107257622467940e-01, 7.671505384435479e-03, 
               1.185128349779530e-01, 7.901864329702142e-03, 
               1.265279166014693e-01, 8.127500454896627e-03, 
               1.347662166290459e-01, 8.348278900798228e-03, 
               1.432228111582069e-01, 8.564067711559794e-03, 
               1.518926458152428e-01, 8.774737913563768e-03, 
               1.607705387761403e-01, 8.980163592507943e-03, 
               1.698511838636768e-01, 9.180221968664177e-03, 
               1.791291537188460e-01, 9.374793470328635e-03, 
               1.885989030447077e-01, 9.563761804941130e-03, 
               1.982547719207257e-01, 9.747014029406790e-03, 
               2.080909891856185e-01, 9.924440616404089e-03, 
               2.181016758866910e-01, 1.009593552085555e-02, 
               2.282808487935948e-01, 1.026139624338005e-02, 
               2.386224239744124e-01, 1.042072389037528e-02, 
               2.491202204319281e-01, 1.057382323411059e-02, 
               2.597679637979141e-01, 1.072060276960456e-02, 
               2.705592900832238e-01, 1.086097476908124e-02, 
               2.814877494814478e-01, 1.099485533428433e-02, 
               2.925468102238625e-01, 1.112216444692532e-02, 
               3.037298624833666e-01, 1.124282601640766e-02, 
               3.150302223250705e-01, 1.135676792514504e-02, 
               3.264411357011814e-01, 1.146392207184255e-02, 
               3.379557824877932e-01, 1.156422441219353e-02, 
               3.495672805611619e-01, 1.165761499703156e-02, 
               3.612686899110480e-01, 1.174403800826822e-02, 
               3.730530167886527e-01, 1.182344179222387e-02, 
               3.849132178866699e-01, 1.189577889050139e-02, 
               3.968422045489604e-01, 1.196100606835144e-02, 
               4.088328470073315e-01, 1.201908434051221e-02, 
               4.208779786428876e-01, 1.206997899450997e-02, 
               4.329704002694061e-01, 1.211365961140755e-02, 
               4.451028844361781e-01, 1.215010008398602e-02, 
               4.572681797477421e-01, 1.217927863234513e-02, 
               4.694590151979302e-01, 1.220117781692452e-02, 
               4.816681045156332e-01, 1.221578454892481e-02, 
               4.938881505196921e-01, 1.222309009813126e-02, 
               5.061118494803080e-01, 1.222309009813127e-02, 
               5.183318954843668e-01, 1.221578454892522e-02, 
               5.305409848020698e-01, 1.220117781692459e-02, 
               5.427318202522577e-01, 1.217927863234521e-02, 
               5.548971155638220e-01, 1.215010008398612e-02, 
               5.670295997305939e-01, 1.211365961140763e-02, 
               5.791220213571124e-01, 1.206997899450985e-02, 
               5.911671529926685e-01, 1.201908434051147e-02, 
               6.031577954510399e-01, 1.196100606835202e-02, 
               6.150867821133300e-01, 1.189577889050199e-02, 
               6.269469832113472e-01, 1.182344179222429e-02, 
               6.387313100889522e-01, 1.174403800826842e-02, 
               6.504327194388386e-01, 1.165761499703152e-02, 
               6.620442175122065e-01, 1.156422441219311e-02, 
               6.735588642988175e-01, 1.146392207184282e-02, 
               6.849697776749293e-01, 1.135676792511821e-02, 
               6.962701375166336e-01, 1.124282601637225e-02, 
               7.074531897761375e-01, 1.112216444689922e-02, 
               7.185122505185524e-01, 1.099485533423064e-02, 
               7.294407099167759e-01, 1.086097476902580e-02, 
               7.402320362020862e-01, 1.072060276960328e-02, 
               7.508797795680726e-01, 1.057382323411068e-02, 
               7.613775760255880e-01, 1.042072389037573e-02, 
               7.717191512064054e-01, 1.026139624348149e-02, 
               7.818983241133093e-01, 1.009593552106589e-02, 
               7.919090108143814e-01, 9.924440616416036e-03, 
               8.017452280792738e-01, 9.747014029353265e-03, 
               8.114010969552918e-01, 9.563761804975165e-03, 
               8.208708462811534e-01, 9.374793470272376e-03, 
               8.301488161363226e-01, 9.180221968665265e-03, 
               8.392294612238589e-01, 8.980163592504867e-03, 
               8.481073541847566e-01, 8.774737913558017e-03, 
               8.567771888417932e-01, 8.564067711554629e-03, 
               8.652337833709542e-01, 8.348278900794001e-03, 
               8.734720833985312e-01, 8.127500454892087e-03, 
               8.814871650220477e-01, 7.901864329699972e-03, 
               8.892742377532061e-01, 7.671505384434269e-03, 
               8.968286473810958e-01, 7.436561301073849e-03, 
               9.041458787539575e-01, 7.197172502083081e-03, 
               9.112215584778225e-01, 6.953482066477137e-03, 
               9.180514575304535e-01, 6.705635644308015e-03, 
               9.246314937889846e-01, 6.453781369631923e-03, 
               9.309577344697750e-01, 6.198069771975543e-03, 
               9.370263984790159e-01, 5.938653686370035e-03, 
               9.428338586726990e-01, 5.675688162041160e-03, 
               9.483766440245793e-01, 5.409330369752567e-03, 
               9.536514417008785e-01, 5.139739507915691e-03, 
               9.586550990404805e-01, 4.867076707504079e-03, 
               9.633846254394735e-01, 4.591504935830873e-03, 
               9.678371941389574e-01, 4.313188899308077e-03, 
               9.720101439151095e-01, 4.032294945242497e-03, 
               9.759009806706316e-01, 3.748990962815882e-03, 
               9.795073789268498e-01, 3.463446283449422e-03, 
               9.828271832159824e-01, 3.175831580853165e-03, 
               9.858584093735683e-01, 2.886318771432173e-03, 
               9.885992457319535e-01, 2.595080916337841e-03, 
               9.910480542178594e-01, 2.302292128352489e-03, 
               9.932033713622925e-01, 2.008127491868935e-03, 
               9.950639092458675e-01, 1.712763020454675e-03, 
               9.966285564501070e-01, 1.416375735729250e-03, 
               9.978963792674904e-01, 1.119144215482236e-03, 
               9.988666243127569e-01, 8.212515093330048e-04, 
               9.995387299886878e-01, 5.229063396712613e-04, 
               9.999124439735656e-01, 2.246904801468211e-04; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 

