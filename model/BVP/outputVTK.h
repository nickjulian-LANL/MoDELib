//====================================================================================
// program to output the elastic field data in vtk format
//
// Mamdouh Mohamed 02/07/2012  FSU
//====================================================================================

std::vector<Eigen::Matrix<double,3,3> > tetGrad_U;                   // displacement gradient field inside tetrahedron
std::vector<Eigen::Matrix<double,3,3> > tetStrain;                  // strain field inside tetrahedron
std::vector<Eigen::Matrix<double,3,1> > tetRotation;                // lattice rotation field inside tetrahedron

//==================================================================================
// function to write the header part of the VTK file
//==================================================================================
void writeHeader (char* VTKfile) {
  FILE *fout =fopen(VTKfile, "w");
  
  //-------------- write header lines -------------
  fprintf (fout, "# vtk DataFile Version 1.0\n");
  fprintf (fout, "Indentation Elastic Fields\n");
  fprintf (fout, "ASCII\n");
  fprintf (fout, "\n");
  
  fclose(fout);
}

//==================================================================================
// function to write the nodes coordinates part of the VTK file
//==================================================================================
void writeNodes(char* VTKfile) {
  FILE *fout =fopen(VTKfile, "a");
  
  unsigned int nn = nodeContainer.size();
  
  fprintf (fout, "DATASET UNSTRUCTURED_GRID\n");
  fprintf (fout, "POINTS %u float\n",nn);
  
  for (unsigned int i= 0 ; i < nodeContainer.size() ; i++){
    fprintf (fout, "%f %f %f\n", nodeContainer[i].Pc(0) ,  nodeContainer[i].Pc(1)  , nodeContainer[i].Pc(2) );
  }
  
  fprintf (fout, "\n");
  
  fclose(fout);
}

//==================================================================================
// function to write the tetrahedron connections part of the VTK file
//==================================================================================
void writeTetrahedrons(char* VTKfile) {
    
  FILE *fout =fopen(VTKfile, "a");
  
  unsigned int nn = tetContainer.size();
  
  fprintf (fout, "CELLS %u %u\n",nn, nn*5 );
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
    fprintf (fout,"%u %lu %lu %lu %lu\n", 4 ,tetContainer[i].eleNodes[0]->sID,tetContainer[i].eleNodes[1]->sID,tetContainer[i].eleNodes[2]->sID,tetContainer[i].eleNodes[3]->sID);
  }
  
  fprintf (fout, "\n");
  
  
  
  fprintf (fout, "CELL_TYPES %u \n",nn);
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
    fprintf (fout,"%u\n", 10);
  }
  
  fprintf (fout, "\n");
  
  fclose(fout);
}

//==================================================================================
// function to write the strain field inside tetrahedrons for the VTK file
//==================================================================================
void writeTetStrainField(char* VTKfile) {
    
  FILE *fout =fopen(VTKfile, "a");
  
  fprintf (fout, "TENSORS Strain float\n");
  
  Eigen::Matrix<double,3,3> strain ;
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
       
    strain = tetStrain[i];
    
    for (unsigned int j = 0; j<3; j++) {
      fprintf (fout,"%f  %f %f \n", strain(j,0), strain(j,1), strain(j,2));
    }
    fprintf (fout, "\n");
  }
  
  fprintf (fout, "\n");
  fclose(fout);
}

//==================================================================================
// function to write the lattice rotation field inside tetrahedrons for the VTK file
//==================================================================================
void writeTetRotationField(char* VTKfile) {
    
  FILE *fout =fopen(VTKfile, "a");
  
  fprintf (fout, "VECTORS Rotation float\n");
  
  Eigen::Matrix<double,3,1> omega;
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
    omega = tetRotation[i];
    
    fprintf (fout,"%f %f %f \n", omega(0), omega(1), omega(2));
  }
  
  fprintf (fout, "\n");
  fclose(fout);
}


//==================================================================================
// function to calculate the total (infinite+image) displacement gradient field 
//==================================================================================
template<typename T>
void calDisplacementGrad(const T* const pT){
  Eigen::Matrix<double,3,3> uprim , rotation;
  Eigen::Matrix<double,3,1> omega;
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
    uprim = tetContainer[i].getTetInfiniteUprim<4,T>(pT) + tetContainer[i].getUprim();     // infinite medium + image
    
    tetGrad_U.push_back(uprim);
    tetStrain.push_back  (0.5e00*(uprim.transpose()+uprim));
    
    rotation = 0.5e00*(uprim.transpose()-uprim);
    omega(0) = rotation(1,2); omega(1) = rotation(2,0); omega(2) = rotation(0,1);
    tetRotation.push_back(omega);
  }
  
}

//==================================================================================
// function to calculate and write the dislocation density tensor data
//==================================================================================

void writeGND(Eigen::Vector3d box,unsigned int npnts , char* GNDfile){
  
   
  
  //Eigen::Matrix<double,3,1> dx = box/double(npnts-1);
  //---vectors hold the (i,j,k) position of sample points and their global index (iv)
  std::vector< std::vector< std::vector<unsigned int> > > vindx(npnts, std::vector< std::vector<unsigned int> > (npnts, std::vector<unsigned int>(npnts)) );
  
  std::vector< std::vector<unsigned int> > vpos (npnts*npnts*npnts, std::vector<unsigned int>(3));
  
  unsigned int iv = 0;
  
  for (unsigned int k=0; k<npnts ; k++){
    for (unsigned int j=0; j<npnts ; j++){
      for (unsigned int i=0; i<npnts ; i++) {
	vindx [i][j][k] = iv;
	vpos[iv][0]=i;   vpos[iv][1]=j;    vpos[iv][2]=k;
	iv++;
      }
    }
  }
  
  //------------ calculate the elastic strain and lattice rotation on sample points --------------
  
  FILE *fomega  =fopen("omega.txt", "w"); 
  FILE *fstrain =fopen("strain.txt", "w"); 
  
  std::vector<Eigen::Matrix<double,3,3> > pntStrain;                  // strain field at sample points
  std::vector<Eigen::Matrix<double,3,1> > pntRotation;                // lattice rotation field at sample points
  
  double dx = box(0)/double(npnts-1);
  Eigen::Vector3d P ;
  
  unsigned int itet;
  Eigen::Matrix<double,3,3> strain;
  Eigen::Matrix<double,3,1> omega;
  
  for (unsigned int k=0; k<npnts ; k++){
    P(2) = k*dx;
    for (unsigned int j=0; j<npnts ; j++){
      P(1) = j*dx;
      for (unsigned int i=0; i<npnts ; i++) {
	P(0) = i*dx;
	
	isTetrahedronType isT = findIncludingTet(P);
	
	if(!isT.first) {std::cout << P.transpose() << std::endl;   assert(0 &&"Error in GND calculations. Sample point outside domain");}
	else itet = isT.second->sID;
	
	strain = tetStrain[itet];
	omega  = tetRotation[itet];
	pntStrain.push_back(strain);
	pntRotation.push_back(omega);
	
	fprintf (fomega ,"%22.15e %22.15e %22.15e \n", omega(0), omega(1), omega(2));
	fprintf (fstrain,"%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", strain(0,0), strain(0,1), strain(0,2), strain(1,1), strain(1,2), strain(2,2));
      }
    }
  }
   
  fclose(fomega);    fclose(fstrain);     
  
  //----------------- calculate the GND at sample points ----------------------
   
  writeGNDFileHeader (GNDfile);
  writeGNDNodes (box,npnts ,GNDfile);
  
  FILE *fout =fopen(GNDfile, "a");
  fprintf (fout, "POINT_DATA %u \n",(npnts)*(npnts)*(npnts));

  
  fprintf (fout, "TENSORS STRAIN double\n");
  
  for (unsigned int i=0; i<pntStrain.size(); i++) {
    strain = pntStrain[i];
    for (unsigned int jj = 0; jj<3; jj++) {
      fprintf (fout,"%22.15e  %22.15e %22.15e \n", strain(jj,0), strain(jj,1), strain(jj,2));
    }
    fprintf (fout, "\n");
  }
  
  /*
  fprintf (fout, "VECTORS ROTATION double\n");
  
  for (unsigned int i=0; i<pntRotation.size(); i++) {
    omega = pntRotation[i];
    fprintf (fout,"%22.15e  %22.15e %22.15e \n", omega(0), omega(1), omega(2));
  }
  
  
  fprintf (fout, "VECTORS ROTATION_Cyl double\n");
  
  Eigen::Matrix<double,3,1> Px , Omega_C, cntr, R , phi;
  
  cntr = 0.5e00*box;
  
  iv = 0;
  
  for (unsigned int k=0; k<npnts ; k++){
    Px(2) = k*dx;
    cntr(2) = Px(2);
    for (unsigned int j=0; j<npnts ; j++){
      Px(1) = j*dx;
      for (unsigned int i=0; i<npnts ; i++) {
	Px(0) = i*dx;
	
	omega = pntRotation[iv];
	
	R = Px - cntr;
	R = R.normalized();
	
	phi(0) = -R(1);   phi(1) =  R(0);    phi(2) = R(2);
	
	Omega_C(0) = R.transpose()*omega;
	Omega_C(1) = phi.transpose()*omega;
	Omega_C(2) = omega(2);
	
	fprintf (fout,"%22.15e  %22.15e %22.15e \n", Omega_C(0), Omega_C(1), Omega_C(2));
	
	iv++;
      }
    }
  }
  
  
  
  
  Eigen::Matrix<double,3,3> omega_prim , alpha_Nye;
  
  Eigen::Matrix<double,3,3> I = Eigen::Matrix<double,3,3>::Zero();
  for (unsigned int i=0;i<3;i++) I(i,i) = 1.0e00;
  
  Eigen::Matrix<int,3,1> index, index_p , index_m;
  unsigned int ivp , ivm;
  
  double dx1;
  
  fprintf (fout, "TENSORS GND_Nye double\n");
  for (unsigned int k=0; k<npnts ; k++){
    index(2) = k;
    for (unsigned int j=0; j<npnts ; j++){
      index(1) = j;
      for (unsigned int i=0; i<npnts ; i++) {
	index(0) = i;
	
	//----- derivative in 3 directions ----------------
	for (unsigned int ix=0;ix<3;ix++) {
	  index_p = index;        index_m = index;
	  index_p(ix) += 1;       index_m(ix) -= 1;
	  
	  dx1 = dx;
	  
	  if(index_p(ix)> int(npnts-1)) {
	    index_p(ix)  = index(ix);
	    dx1 = dx/2.0e00;
	  } 
	  if(index_m(ix)< 0) {
	    index_m(ix)  = index(ix);
	    dx1 = dx/2.0e00;
	  } 
	  
	  ivp = vindx[index_p(0)][index_p(1)][index_p(2)];
	  ivm = vindx[index_m(0)][index_m(1)][index_m(2)];
	  
	  omega_prim.col(ix) = (pntRotation[ivp] - pntRotation[ivm])*1.0e6/(dx1*0.2556e00);           // millirad / micron
	}
	
	//---------The dislocation density tensor -----------
	
	alpha_Nye = omega_prim - (omega_prim.trace()*I);
	
	for (unsigned int jj = 0; jj<3; jj++) {
	  fprintf (fout,"%22.15e  %22.15e %22.15e \n", alpha_Nye(jj,0), alpha_Nye(jj,1), alpha_Nye(jj,2));
	}
	fprintf (fout, "\n");
	
      }
    }
  }
  */
  fclose(fout);
}




//==================================================================================
// function to write the header part of the VTK file
//==================================================================================
void writeGNDFileHeader (char* GNDfile) {
  FILE *fout =fopen(GNDfile, "w");
  
  //-------------- write header lines -------------
  fprintf (fout, "# vtk DataFile Version 1.0\n");
  fprintf (fout, "Indentation Dislocation Density Tensor\n");
  fprintf (fout, "ASCII\n");
  fprintf (fout, "\n");
  
  fclose(fout);
}

//==================================================================================
// function to write the nodes coordinates part of the VTK file
//==================================================================================
void writeGNDNodes(Eigen::Vector3d box,unsigned int npnts , char* GNDfile) {
  
  //unsigned int n = npnts-2;
  unsigned int n = npnts;
  
  FILE *fout =fopen(GNDfile, "a");
  
  fprintf (fout, "DATASET STRUCTURED_GRID\n");
  fprintf (fout, "DIMENSIONS %u %u %u\n",n , n , n);
  fprintf (fout, "POINTS %u double\n",(n)*(n)*(n));
  
  Eigen::Vector3d P;
  Eigen::Vector3d step = box/(npnts-1);
  
  
  for (unsigned int k=0; k<n ; k++){
    P(2) = k * step(2);
    for (unsigned int j=0; j<n ; j++){
      P(1) = j* step(1);
      for (unsigned int i=0; i<n ; i++) {
	P(0) = i * step(0);
	fprintf (fout, "%22.15e %22.15e %22.15e\n", P(0) ,  P(1)  , P(2) );
      }
    }
  }
  
  fprintf (fout, "\n");
  
  fclose(fout);
}



//==================================================================================
// function to write mesh field data in VTK format
//==================================================================================
template<typename T>
void writeVTK_file (const T* const pT) {
  
  //char VTKfile[] =  "data.vtk";
  char VTK_GND_file[] =  "GND_double.vtk";
  Eigen::Vector3d box(4000.0e00, 4000.0e00, 4000.0e00);
  unsigned int npnts = 50;
  
  //writeHeader(VTKfile);
  //writeNodes(VTKfile);
  //writeTetrahedrons(VTKfile);
  
  //FILE *fout =fopen(VTKfile, "a"); 
  //unsigned int nn = tetContainer.size();
  //fprintf (fout, "CELL_DATA %u \n",nn);
  //fclose(fout);
  
  //writeTetStressField(VTKfile,pT);
  
  calDisplacementGrad(pT);
  
  //writeTetStrainField(VTKfile);
  
  //writeTetRotationField(VTKfile);
  
  writeGND(box,npnts,VTK_GND_file);
  
  
}



//==================================================================================
// function to write the stress field inside tetrahedrons for the VTK file
//==================================================================================
template<typename T>
void writeTetStressField(char* VTKfile,const T* const pT) {
    
  FILE *fout =fopen(VTKfile, "a");
  
  fprintf (fout, "TENSORS Stress float\n");
  
  Eigen::Matrix<double,3,3> stress;
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
    //stress = tetContainer[i].getStress() + tetContainer[i].getTetInfiniteStress<4,T>(pT);
    //stress = tetContainer[i].getTetInfiniteStress<4,T>(pT);
    stress = tetContainer[i].getStress();
    
    for (unsigned int j = 0; j<3; j++) {
      fprintf (fout,"%f  %f %f \n", stress(j,0), stress(j,1), stress(j,2));
    }
    fprintf (fout, "\n");
  }
  
  fprintf (fout, "\n");
  fclose(fout);
}



/*
 
 //=================================================================================
// function to calculate the strain from stress by Hooke's law
//================================================================================

Eigen::Matrix<double,3,3> stress2strain(Eigen::Matrix<double,3,3> stress, double mu, double nu,double E){
  Eigen::Matrix<double,3,3> strain;
  
  strain(0,0) = (stress(0,0)-nu*(stress(1,1)+stress(2,2)))/E;
  strain(1,1) = (stress(1,1)-nu*(stress(0,0)+stress(2,2)))/E;
  strain(2,2) = (stress(2,2)-nu*(stress(1,1)+stress(0,0)))/E;
  
  strain(0,1) = stress(0,1)/(2.0e00*mu);     strain(1,0) =  strain(0,1); 
  strain(0,2) = stress(0,2)/(2.0e00*mu);     strain(2,0) =  strain(0,2);
  strain(1,2) = stress(1,2)/(2.0e00*mu);     strain(2,1) =  strain(1,2);
  
  return strain;
}
 
 
 

 
//==================================================================================
// function to write the strain field inside tetrahedrons for the VTK file
//==================================================================================
template<typename T>
void writeRotationDifference(char* VTKfile,const T* const pT) {
    
  FILE *fout =fopen(VTKfile, "a");
  
  fprintf (fout, "TENSORS Rotation_diff float\n");
  
  Eigen::Matrix<double,3,3> omega, omegau ,  uprim , diff;
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
    
    omega  = tetContainer[i].getTetInfiniteRotation<4,T>(pT);
    
    uprim = tetContainer[i].getTetInfiniteUprim<4,T>(pT);
    
    omegau  = 0.5e00*(uprim.transpose()-uprim);
    
    diff = omega - omegau;
        
    for (unsigned int j = 0; j<3; j++) {
      fprintf (fout,"%f  %f %f \n", diff(j,0), diff(j,1), diff(j,2));
    }
    fprintf (fout, "\n");
  }
  
  fprintf (fout, "\n");
  fclose(fout);
}



//==================================================================================
// function to write the strain field inside tetrahedrons for the VTK file
//==================================================================================
template<typename T>
void writeStrainDifference(char* VTKfile,const T* const pT) {
    
  FILE *fout =fopen(VTKfile, "a");
  
  fprintf (fout, "TENSORS Strain_diff float\n");
  
  double mu = tetContainer[0].mu;
  double nu = tetContainer[0].nu;
  double E = 2.0e00*mu*(1.0e00+nu);
  
  Eigen::Matrix<double,3,3> strain , stress , strainu, uprim , diff;
  
  for (unsigned int i= 0 ; i < tetContainer.size() ; i++){
    stress = tetContainer[i].getTetInfiniteStress<4,T>(pT);
    strain = stress2strain(stress,mu,nu,E);
    
    uprim = tetContainer[i].getTetInfiniteUprim<4,T>(pT);
    
    strainu = 0.5e00*(uprim.transpose()+uprim);
    
    diff = strainu - strainu;
        
    for (unsigned int j = 0; j<3; j++) {
      fprintf (fout,"%f  %f %f \n", diff(j,0), diff(j,1), diff(j,2));
    }
    fprintf (fout, "\n");
  }
  
  fprintf (fout, "\n");
  fclose(fout);
}*/
