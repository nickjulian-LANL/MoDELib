 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_1_256_H_ 
#define model_GAUSSLEGENDRE_1_256_H_ 

namespace model{

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
   template<>
   struct GaussLegendre<1,256>{
       static Eigen::Matrix<double,2,256> abcsissasAndWeights(){
           Eigen::Matrix<double,256,2> aw;
           aw<<2.197499050371476e-05, 5.639450891035658e-05, 
               1.157812953684334e-04, 1.312674721473316e-04, 
               2.845312668698918e-04, 2.062316272118132e-04, 
               5.282370782954682e-04, 2.811744770153659e-04, 
               8.468667634971005e-04, 3.560770817369917e-04, 
               1.240373621639534e-03, 4.309268507125937e-04, 
               1.708698988309110e-03, 5.057121966038366e-04, 
               2.251772759452042e-03, 5.804217787842603e-04, 
               2.869513538794721e-03, 6.550443409526445e-04, 
               3.561828695588742e-03, 7.295686666557487e-04, 
               4.328614396208685e-03, 8.039835653730389e-04, 
               5.169755627467854e-03, 8.782778681650507e-04, 
               6.085126217569903e-03, 9.524404267480332e-04, 
               7.074588856937347e-03, 1.026460113983829e-03, 
               8.137995119842312e-03, 1.100325824918307e-03, 
               9.275185487268522e-03, 1.174026478164924e-03, 
               1.048598937118916e-02, 1.247551017352516e-03, 
               1.177022514038317e-02, 1.320888412713640e-03, 
               1.312770014781500e-02, 1.394027662664674e-03, 
               1.455821075962843e-02, 1.466957795415547e-03, 
               1.606154238575502e-02, 1.539667870599339e-03, 
               1.763746951214712e-02, 1.612146980897435e-03, 
               1.928575573463398e-02, 1.684384253658535e-03, 
               2.100615379441079e-02, 1.756368852528540e-03, 
               2.279840561514174e-02, 1.828089979073074e-03, 
               2.466224234168612e-02, 1.899536874382374e-03, 
               2.659738438043646e-02, 1.970698820703570e-03, 
               2.860354144126870e-02, 2.041565143026066e-03, 
               3.068041258109255e-02, 2.112125210690663e-03, 
               3.282768624899857e-02, 2.182368438983998e-03, 
               3.504504033299705e-02, 2.252284290723613e-03, 
               3.733214220834219e-02, 2.321862277839315e-03, 
               3.968864878742667e-02, 2.391091962947081e-03, 
               4.211420657125442e-02, 2.459962960906266e-03, 
               4.460845170246786e-02, 2.528464940392603e-03, 
               4.717101001992746e-02, 2.596587625436331e-03, 
               4.980149711484783e-02, 2.664320796958686e-03, 
               5.249951838845762e-02, 2.731654294321260e-03, 
               5.526466911119421e-02, 2.798578016841772e-03, 
               5.809653448342056e-02, 2.865081925302048e-03, 
               6.099468969764638e-02, 2.931156043461653e-03, 
               6.395870000225595e-02, 2.996790459557045e-03, 
               6.698812076672350e-02, 3.061975327783651e-03, 
               7.008249754831297e-02, 3.126700869770436e-03, 
               7.324136616024957e-02, 3.190957376053376e-03, 
               7.646425274135282e-02, 3.254735207526530e-03, 
               7.975067382711937e-02, 3.318024796890696e-03, 
               8.310013642224834e-02, 3.380816650087656e-03, 
               8.651213807459468e-02, 3.443101347721859e-03, 
               8.998616695054146e-02, 3.504869546485727e-03, 
               9.352170191178411e-02, 3.566111980537782e-03, 
               9.711821259350129e-02, 3.626819462917817e-03, 
               1.007751594839154e-01, 3.686982886905223e-03, 
               1.044919940052274e-01, 3.746593227402832e-03, 
               1.082681585959086e-01, 3.805641542271530e-03, 
               1.121030867943373e-01, 3.864118973692293e-03, 
               1.159962033237769e-01, 3.922016749471005e-03, 
               1.199469241786731e-01, 3.979326184375935e-03, 
               1.239546567122545e-01, 4.036038681436356e-03, 
               1.280187997254445e-01, 4.092145733220095e-03, 
               1.321387435570409e-01, 4.147638923119145e-03, 
               1.363138701751743e-01, 4.202509926610564e-03, 
               1.405435532700151e-01, 4.256750512510763e-03, 
               1.448271583477289e-01, 4.310352544199753e-03, 
               1.491640428256574e-01, 4.363307980849739e-03, 
               1.535535561287113e-01, 4.415608878624466e-03, 
               1.579950397869614e-01, 4.467247391879129e-03, 
               1.624878275344185e-01, 4.518215774331182e-03, 
               1.670312454089757e-01, 4.568506380225944e-03, 
               1.716246118535134e-01, 4.618111665477913e-03, 
               1.762672378181374e-01, 4.667024188811612e-03, 
               1.809584268635448e-01, 4.715236612870716e-03, 
               1.856974752654929e-01, 4.762741705315049e-03, 
               1.904836721203695e-01, 4.809532339920095e-03, 
               1.953162994518330e-01, 4.855601497636521e-03, 
               2.001946323185157e-01, 4.900942267629975e-03, 
               2.051179389227727e-01, 4.945547848353762e-03, 
               2.100854807204585e-01, 4.989411548517069e-03, 
               2.150965125317157e-01, 5.032526788166977e-03, 
               2.201502826527594e-01, 5.074887099572664e-03, 
               2.252460329686406e-01, 5.116486128205294e-03, 
               2.303829990669702e-01, 5.157317633972599e-03, 
               2.355604103525888e-01, 5.197375491592900e-03, 
               2.407774901631627e-01, 5.236653692085751e-03, 
               2.460334558856926e-01, 5.275146343290796e-03, 
               2.513275190739098e-01, 5.312847670948115e-03, 
               2.566588855665549e-01, 5.349752019480217e-03, 
               2.620267556065082e-01, 5.385853852848918e-03, 
               2.674303239607601e-01, 5.421147755552826e-03, 
               2.728687800412051e-01, 5.455628433066478e-03, 
               2.783413080262362e-01, 5.489290712870853e-03, 
               2.838470869831293e-01, 5.522129545412701e-03, 
               2.893852909911881e-01, 5.554140004507183e-03, 
               2.949550892656417e-01, 5.585317288277954e-03, 
               3.005556462822693e-01, 5.615656720379298e-03, 
               3.061861219027431e-01, 5.645153747378694e-03, 
               3.118456715006416e-01, 5.673803947775960e-03, 
               3.175334460881728e-01, 5.701603021519007e-03, 
               3.232485924435158e-01, 5.728546799045363e-03, 
               3.289902532388149e-01, 5.754631238519795e-03, 
               3.347575671687917e-01, 5.779852427021105e-03, 
               3.405496690799472e-01, 5.804206581127550e-03, 
               3.463656901003408e-01, 5.827690047472395e-03, 
               3.522047577699327e-01, 5.850299303310289e-03, 
               3.580659961714596e-01, 5.872030957030034e-03, 
               3.639485260618319e-01, 5.892881748672312e-03, 
               3.698514650040287e-01, 5.912848550412427e-03, 
               3.757739274994715e-01, 5.931928367035445e-03, 
               3.817150251208578e-01, 5.950118336383097e-03, 
               3.876738666454337e-01, 5.967415729781587e-03, 
               3.936495581886868e-01, 5.983817952452843e-03, 
               3.996412033384366e-01, 5.999322543902522e-03, 
               4.056479032893056e-01, 6.013927178291120e-03, 
               4.116687569775490e-01, 6.027629664780221e-03, 
               4.177028612162230e-01, 6.040427947861979e-03, 
               4.237493108306718e-01, 6.052320107670632e-03, 
               4.298071987943122e-01, 6.063304360263855e-03, 
               4.358756163646965e-01, 6.073379057897033e-03, 
               4.419536532198337e-01, 6.082542689267532e-03, 
               4.480403975947455e-01, 6.090793879740993e-03, 
               4.541349364182403e-01, 6.098131391557536e-03, 
               4.602363554498835e-01, 6.104554124018636e-03, 
               4.663437394171419e-01, 6.110061113652012e-03, 
               4.724561721526830e-01, 6.114651534354820e-03, 
               4.785727367318106e-01, 6.118324697519909e-03, 
               4.846925156100106e-01, 6.121080052136477e-03, 
               4.908145907605932e-01, 6.122917184874232e-03, 
               4.969380438124053e-01, 6.123835820144786e-03, 
               5.030619561875949e-01, 6.123835820144888e-03, 
               5.091854092394069e-01, 6.122917184874180e-03, 
               5.153074843899894e-01, 6.121080052136605e-03, 
               5.214272632681896e-01, 6.118324697519877e-03, 
               5.275438278473170e-01, 6.114651534355389e-03, 
               5.336562605828581e-01, 6.110061113652108e-03, 
               5.397636445501165e-01, 6.104554124018432e-03, 
               5.458650635817598e-01, 6.098131391556726e-03, 
               5.519596024052549e-01, 6.090793879740700e-03, 
               5.580463467801665e-01, 6.082542689267967e-03, 
               5.641243836353036e-01, 6.073379057897740e-03, 
               5.701928012056880e-01, 6.063304360263787e-03, 
               5.762506891693283e-01, 6.052320107670379e-03, 
               5.822971387837770e-01, 6.040427947862376e-03, 
               5.883312430224508e-01, 6.027629664780021e-03, 
               5.943520967106944e-01, 6.013927178291010e-03, 
               6.003587966615638e-01, 5.999322543902921e-03, 
               6.063504418113130e-01, 5.983817952452899e-03, 
               6.123261333545662e-01, 5.967415729782046e-03, 
               6.182849748791421e-01, 5.950118336383842e-03, 
               6.242260725005284e-01, 5.931928367035524e-03, 
               6.301485349959713e-01, 5.912848550411948e-03, 
               6.360514739381684e-01, 5.892881748670986e-03, 
               6.419340038285410e-01, 5.872030957030514e-03, 
               6.477952422300678e-01, 5.850299303310433e-03, 
               6.536343098996595e-01, 5.827690047472671e-03, 
               6.594503309200531e-01, 5.804206581126029e-03, 
               6.652424328312085e-01, 5.779852427021521e-03, 
               6.710097467611860e-01, 5.754631238520093e-03, 
               6.767514075564850e-01, 5.728546799045556e-03, 
               6.824665539118273e-01, 5.701603021519189e-03, 
               6.881543284993582e-01, 5.673803947772599e-03, 
               6.938138780972579e-01, 5.645153747938251e-03, 
               6.994443537177296e-01, 5.615656719824861e-03, 
               7.050449107343584e-01, 5.585317288277089e-03, 
               7.106147090088120e-01, 5.554140004505284e-03, 
               7.161529130168711e-01, 5.522129545407256e-03, 
               7.216586919737640e-01, 5.489290712863832e-03, 
               7.271312199587955e-01, 5.455628433024554e-03, 
               7.325696760392395e-01, 5.421147755557103e-03, 
               7.379732443934913e-01, 5.385853852902375e-03, 
               7.433411144334447e-01, 5.349752019490587e-03, 
               7.486724809260903e-01, 5.312847670948376e-03, 
               7.539665441143082e-01, 5.275146343290163e-03, 
               7.592225098368373e-01, 5.236653692085601e-03, 
               7.644395896474112e-01, 5.197375491605836e-03, 
               7.696170009330294e-01, 5.157317633967467e-03, 
               7.747539670313590e-01, 5.116486128239363e-03, 
               7.798497173472405e-01, 5.074887099547152e-03, 
               7.849034874682842e-01, 5.032526788152444e-03, 
               7.899145192795407e-01, 4.989411548516805e-03, 
               7.948820610772269e-01, 4.945547848346794e-03, 
               7.998053676814837e-01, 4.900942267627266e-03, 
               8.046837005481671e-01, 4.855601497633372e-03, 
               8.095163278796305e-01, 4.809532339921271e-03, 
               8.143025247345074e-01, 4.762741705313811e-03, 
               8.190415731364556e-01, 4.715236612869424e-03, 
               8.237327621818626e-01, 4.667024188811324e-03, 
               8.283753881464867e-01, 4.618111665478691e-03, 
               8.329687545910245e-01, 4.568506380224819e-03, 
               8.375121724655819e-01, 4.518215774331824e-03, 
               8.420049602130384e-01, 4.467247391880710e-03, 
               8.464464438712891e-01, 4.415608878625174e-03, 
               8.508359571743433e-01, 4.363307980848721e-03, 
               8.551728416522717e-01, 4.310352544201167e-03, 
               8.594564467299852e-01, 4.256750512511149e-03, 
               8.636861298248250e-01, 4.202509926609729e-03, 
               8.678612564429582e-01, 4.147638923117275e-03, 
               8.719812002745557e-01, 4.092145733218499e-03, 
               8.760453432877457e-01, 4.036038681437195e-03, 
               8.800530758213275e-01, 3.979326184377731e-03, 
               8.840037966762231e-01, 3.922016749469367e-03, 
               8.878969132056631e-01, 3.864118973691497e-03, 
               8.917318414040920e-01, 3.805641542272890e-03, 
               8.955080059947731e-01, 3.746593227404092e-03, 
               8.992248405160850e-01, 3.686982886907498e-03, 
               9.028817874064990e-01, 3.626819462916511e-03, 
               9.064782980882160e-01, 3.566111980538239e-03, 
               9.100138330494582e-01, 3.504869546485126e-03, 
               9.134878619254059e-01, 3.443101347723610e-03, 
               9.168998635777517e-01, 3.380816650087917e-03, 
               9.202493261728799e-01, 3.318024796888678e-03, 
               9.235357472586467e-01, 3.254735207526021e-03, 
               9.267586338397510e-01, 3.190957376052337e-03, 
               9.299175024516870e-01, 3.126700869769986e-03, 
               9.330118792332769e-01, 3.061975327783936e-03, 
               9.360412999977438e-01, 2.996790459557663e-03, 
               9.390053103023527e-01, 2.931156043462104e-03, 
               9.419034655165781e-01, 2.865081925299684e-03, 
               9.447353308888043e-01, 2.798578016841509e-03, 
               9.475004816115409e-01, 2.731654294323007e-03, 
               9.501985028851510e-01, 2.664320796959142e-03, 
               9.528289899800717e-01, 2.596587625434701e-03, 
               9.553915482975321e-01, 2.528464940391648e-03, 
               9.578857934287441e-01, 2.459962960905523e-03, 
               9.603113512125718e-01, 2.391091962944836e-03, 
               9.626678577916576e-01, 2.321862277838815e-03, 
               9.649549596670027e-01, 2.252284290723800e-03, 
               9.671723137510013e-01, 2.182368438985840e-03, 
               9.693195874189069e-01, 2.112125210691704e-03, 
               9.713964585587314e-01, 2.041565143026326e-03, 
               9.734026156195630e-01, 1.970698820704664e-03, 
               9.753377576583134e-01, 1.899536874384073e-03, 
               9.772015943848577e-01, 1.828089979070901e-03, 
               9.789938462055878e-01, 1.756368852526528e-03, 
               9.807142442653649e-01, 1.684384253656905e-03, 
               9.823625304878528e-01, 1.612146980897064e-03, 
               9.839384576142440e-01, 1.539667870600918e-03, 
               9.854417892403712e-01, 1.466957795415074e-03, 
               9.868722998521846e-01, 1.394027662661028e-03, 
               9.882297748596167e-01, 1.320888412714057e-03, 
               9.895140106288107e-01, 1.247551017352175e-03, 
               9.907248145127325e-01, 1.174026478163490e-03, 
               9.918620048801583e-01, 1.100325824919919e-03, 
               9.929254111430637e-01, 1.026460113983701e-03, 
               9.939148737824314e-01, 9.524404267516015e-04, 
               9.948302443725330e-01, 8.782778681651860e-04, 
               9.956713856037914e-01, 8.039835653741221e-04, 
               9.964381713044106e-01, 7.295686666562332e-04, 
               9.971304864612037e-01, 6.550443409521282e-04, 
               9.977482272405481e-01, 5.804217787818576e-04, 
               9.982913010116907e-01, 5.057121966039054e-04, 
               9.987596263783614e-01, 4.309268507095519e-04, 
               9.991531332365047e-01, 3.560770817360450e-04, 
               9.994717629217046e-01, 2.811744770166419e-04, 
               9.997154687331317e-01, 2.062316272119270e-04, 
               9.998842187046328e-01, 1.312674721496317e-04, 
               9.999780250094974e-01, 5.639450891445142e-05; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 

