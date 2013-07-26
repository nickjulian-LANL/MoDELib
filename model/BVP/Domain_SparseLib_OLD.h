/* This file is part of finite element solution of BVP attached with mmdl "the Mechanics of Material Defects Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_domain_H_
#define bvpfe_domain_H_

//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include <mmdl/BVP/UmfPackSupport.h>
//#include "mmdl/BVP/SuperLUSupport.h"

#include <map>

#include "mmdl/Quadrature/Quadrature.h"
#include "mmdl/Utilities/SequentialOutputFile.h"

#include "mmdl/BVP/Node.h"
#include "mmdl/BVP/Tetrahedron.h"
#include "mmdl/BVP/Triangle.h"
#include "mmdl/BVP/Face.h"
#include "mmdl/BVP/SearchData.h"

//		  #include "mmdl/BVP/UpdateBoundaryConditions.h"


#include <stdio.h>

#include "SparseLib/CG.h"

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>


namespace bvpfe{
	
	
	class Domain {
		
		enum {dim=3};
		
#include <mmdl/BVP/commonTypeDefs.h>
		
		
		
		typedef boost::ptr_vector<Node<dim> >  NodeContainerType;
		typedef boost::ptr_vector<Tetrahedron> TetContainerType;
		typedef boost::ptr_vector<Triangle>    TriContainerType;
		typedef boost::ptr_map<size_t,Face>    FaceContainerType;
		typedef std::pair<bool,Tetrahedron*>   isTetrahedronType;
		
		
		
	public:
		NodeContainerType nodeContainer;                  // container for the nodes' pointers
		TetContainerType tetContainer;		         // container for the tetrahedrons' pointers
		std::vector<Triangle*> triContainer;
	private:
		FaceContainerType faceContainer;                 // container for the faces pointers
		
		//TriContainerType triContainer;                 // container for triangle pointers
		//std::vector<Triangle*> triContainer;
		
		std::vector<double> Finf;                  // holds the r.h.s. traction (infinite + external) force vector for all the nodes
		
		std::map<bvpfe::Face*,std::vector<double> > facesWtraction;    // container for the faces that have traction, and their traction vector
		
		unsigned int sysDim;                             // dimension of the global linear system	
		
		public :
		Domain () {}

		//==================================================================================
		// function to read the mesh data from input files mesh.*
		//==================================================================================
		
		void readMesh(){
		  
			readNodes();
			
			Finf.resize(nodeContainer.size()*3, 0.0);
			
			readVolElements();
			
			readSurfElements();
			
			output();
		}
		
		//==================================================================================
		// ---------- function to read & store the nodes data in Node objects and container ---------------
		//==================================================================================
		
		void readNodes()
		{
			// int nNodes , d1, d2 , d3, in; // Giacomo 09-30-2011
			unsigned int nNodes , d1, d2 , d3, in;
			VectorDim tempP; 
			
			FILE *fp =fopen("mesh.node", "r");
			
			 assert(fscanf (fp, "%d%d%d%d", &nNodes , &d1, &d2, &d3)==4);
			
			if ((d1!=3)&&(d2!=0)&&(d3!=0)) assert(0 &&"Error in .node file format");
			
			float X[3];
  
			for (unsigned int i = 0; i< nNodes; i++){
			   assert(fscanf (fp, "%d%f%f%f", &in , &X[0], &X[1], &X[2])==4);
			  
			  for (unsigned int ii = 0; ii<3; ii++) {tempP(ii) = X[ii];}
			  			  
			  std::auto_ptr<Node<dim> > pNode (new Node<dim>(tempP) );
			  if(pNode->sID != in) assert(0&&"Error in .node file format. use -z option when creating the mesh to start numbering from 0");
			  nodeContainer.push_back(pNode);
			}
			
			fclose(fp);
			
		}
		
		//==================================================================================
		//----------------- function to read & distribute the volume elements data ------------
		//==================================================================================
		
		void readVolElements()
		{
		  
		  unsigned int nTets , d1, nn , ti , nt;
		  
		  FILE *fTet =fopen("mesh.ele", "r");
		  FILE *fneigh =fopen("mesh.neigh", "r");
		  
		   assert(fscanf (fTet, "%d%d%d", &nTets , &nn, &d1)==3);
		  if ((nn!=4)&&(d1!=0)) assert(0&&"Error in .ele file format");
		  
		   assert(fscanf (fneigh, "%d%d",  &nt, &nn)==2);
		  if ((nn!=4)&&(nt!=nTets)) assert(0&&"Error in .neigh file format");
		  
		  int neighbor[4];
		  int tetNodes[4];
			
		  for (unsigned int i=0; i<nTets; i++) {
				
			std::auto_ptr<Tetrahedron> pTet (new Tetrahedron ); 
			
			 assert(fscanf (fneigh, "%d%d%d%d%d",  &ti, &neighbor[0], &neighbor[1], &neighbor[2], &neighbor[3])==5);
			 assert(fscanf (fTet, "%d%d%d%d%d",  &ti, &tetNodes[0], &tetNodes[1], &tetNodes[2], &tetNodes[3])==5);	
			
			if(pTet->sID != ti) assert(0&&"Error in .ele or .neigh file format. use -z option when creating the mesh to start numbering from 0");

			for(unsigned int j = 0; j<nn; j++) {
				pTet->insertNode(&nodeContainer[tetNodes[j]]); // save the node ptr in the element's eleNodes array
			}
				
			//----------- set neighbor elements -----------------
			pTet->addNeighbor(neighbor);
				
			//--------- set the nodes neighborlist ---------
			pTet->setNodesNeighbors();
				
			//assert(tetContainer.insert(std::make_pair(pTet->sID,pTet)).second);
			tetContainer.push_back(pTet);
		  }
		  			
		  fclose(fTet);
		  fclose(fneigh);
		}
		
		//==================================================================================
		//------------ function to read & distribute the triangular surface elements data----------
		//==================================================================================
		
		void readSurfElements()
		{
		  int triNodes[3];			
		  unsigned int nTris , di , iFc, iTet , ti;
		  
		  FILE *fTri =fopen("mesh.face", "r");
		  
		   assert(fscanf (fTri, "%d%d", &nTris , &di)==2);
		  
		  if (di!=0) assert(0&&"Error in .face file format");
		  
		  Triangle* pTri;
		  
		  for (unsigned int iTri = 0; iTri < nTris; iTri++ ){
		     assert(fscanf (fTri, "%d%d%d%d%d%d", &ti , &triNodes[0], &triNodes[1], &triNodes[2], &iTet  , &iFc )==6);
		    
		    //----------- if the face does NOT exist, creat it ----------
		    if(faceContainer.find(iFc) == faceContainer.end()) {
		      std::auto_ptr<Face> pFace (new Face);
		      assert(faceContainer.insert(iFc,pFace).second);
		    }
		    
		    //---------- add triangles ---------
		    pTri = new Triangle;
		    
		    if(pTri->sID != ti) {
		      assert(0&&"Error in .face file format. use -z option when creating the mesh to start numbering from 0");
		    }

		    
		    for (unsigned int j = 0; j<3; j++ ) {      // insert triangle nodes
		      nodeContainer[triNodes[j]].isBoundaryNode=true;
		      pTri->insertNode(&nodeContainer[triNodes[j]]);
		    }   
		    
		    triContainer.push_back(pTri);
		    
		    pTri->outNormal = pTri->triNormal();                       // set the outward normal for the triangle
		    
		    faceContainer.at(iFc).insertTri(pTri);                     // insert triangle in the face
				
		    // ------------ set neighbor tetrahedron for the face triangle ----------------
		    
		    pTri->neighTetIndx = iTet;         // index of the triangle's neighbor tetrahedron 			
		    tetContainer[iTet].insertSurfTri(pTri);   // set tetrahedron's neighbor surface triangle
		  }
			
		  setNeighborTriangles();                 // set for each triangle an array of pointers to the 3 neighbor triangles
			
		  fclose(fTri);
		  }
		
		//===================================================================================
		// function to set the 3 neighbor triangles for each surface triangle
		//===================================================================================
		
		void setNeighborTriangles(){
			
			//Triangle* pTri1, pTri2;
			
			Eigen::Matrix<size_t,2,3>   edgi , edgj; 
			
			//unsigned int indxi, indxj;
			
			for (unsigned int i = 0; i<(triContainer.size()-1); i++){
				
				//----------- take here specific order, and take the opposite for the pTri2 to avoid the need to sort
				edgi(0,0)=triContainer[i]->eleNodes[1]->sID;     edgi(1,0)=triContainer[i]->eleNodes[2]->sID;
				edgi(0,1)=triContainer[i]->eleNodes[2]->sID;     edgi(1,1)=triContainer[i]->eleNodes[0]->sID;
				edgi(0,2)=triContainer[i]->eleNodes[0]->sID;     edgi(1,2)=triContainer[i]->eleNodes[1]->sID;
				
				for (unsigned int j = i+1; j<triContainer.size(); j++){
					
					edgj(0,0)=triContainer[j]->eleNodes[2]->sID;     edgj(1,0)=triContainer[j]->eleNodes[1]->sID;
					edgj(0,1)=triContainer[j]->eleNodes[1]->sID;     edgj(1,1)=triContainer[j]->eleNodes[0]->sID;
					edgj(0,2)=triContainer[j]->eleNodes[0]->sID;     edgj(1,2)=triContainer[j]->eleNodes[2]->sID;
					
					for (unsigned int ii = 0; ii<3; ii++){
						for (unsigned int jj = 0; jj<3; jj++){
							
							if(edgi.col(ii) == edgj.col(jj)){
								setTriCouples(i,j,ii,jj); break;
							}
							
						}
					}
					
				}
			}	
			
			
			/*for (unsigned int i = 0; i<triContainer.size(); i++){
			 std::cout<< "Triangle: "<< triContainer[i]->sID<< " : " << triContainer[i]->neighbor[0]->sID << " "<< triContainer[i]->neighbor[1]->sID << " "
			 << triContainer[i]->neighbor[2]->sID << std::endl;
			 }*/
			
			
		}
		
		//==================================================================================
		// function to set 2 neighbor triangles as neighbors
		//===================================================================================
		
		void setTriCouples(unsigned int ti,unsigned int tj,unsigned int ni,unsigned int nj){
			
			//-----temporary, and needs modification. Index modification is mainly becasue of different sequence convention for edgj 
			if(nj==2) {nj=1;}
			else if (nj==1) {nj=2;}
			
			//std::cout << ti << " " << tj << " " << ni << " " << nj << std::endl;
			
			triContainer[ti]->neighbor[ni] = triContainer[tj];
			triContainer[tj]->neighbor[nj] = triContainer[ti];
			
		}
		
		
		//==================================================================================
		// function to read & set the traction and displacement boundary conditions for the domain
		//===================================================================================
		
		void setBoundaryConditions() {
		  
		  Eigen::Matrix <float,dim,1>  u , tr;
		  Eigen::Matrix <unsigned int,dim,1> isBC; 
		  unsigned int ni;
		 		  
		  FILE *fp =fopen("BCs_0.txt", "r");
		 		  
		  while (!feof(fp)) {
		    if(fscanf(fp, "%d%d%d%d%f%f%f%f%f%f", &ni,&isBC(0), &isBC(1), &isBC(2),&u(0),&u(1),&u(2),&tr(0),&tr(1),&tr(2))==10){
		      if(nodeContainer[ni].isBoundaryNode) {
			for (unsigned int j = 0; j<3 ; j++) {
			  if(isBC(j))  nodeContainer[ni].setBC(j ,double(u(j)) );
			  nodeContainer[ni].traction(j) = double(tr(j));
			}
		      }
		      else  {assert(0&&"inputting boundary condition for non-boundary node");}
		    }
		  }	
		  
		  fclose(fp);
		}
		
		//================================================================================
		//  function to set the dof of a face (the whole face)
		//================================================================================
		//	template <typename T>
		void setFaceDof(int iface, int idof, double val)
		{
			bvpfe::Face* pFace;
			std::map<size_t,bvpfe::Triangle*>::iterator it;
			
			idof = idof -1;
			pFace = faceContainer.find(iface)->second;    // pointer to the targetted face
			
			//----- loop over the face triangles ---
			
			for (it= pFace->triContainer.begin() ; it != pFace->triContainer.end(); it++ )
			{
				//----- loop over the triangle nodes ---
				for(unsigned int i = 0; i<it->second->Nnodes; i++)
				{
					//std::cout<< it->second->eleNodes[i]->sID;
					//VectorDim uInf = pT->displacement(it->second->eleNodes[i]->P);
					it->second->eleNodes[i]->setBC (idof , val-it->second->eleNodes[i]->uInf(idof)*0);
				}
			}
			
			
		}
		
		//===================================================================================
		// function to assemble the global linear system, solve it, and distribute nodal displacement
		//===================================================================================
		
		template<typename T>
		void solveBVP(bool dislocations_Stress,const T* const pT)
		{

		  //-------------------- Calculate infinite medium stress field --------------------
		  for (unsigned int dN=0; dN<nodeContainer.size();++dN){
		    if(nodeContainer[dN].isBoundaryNode) nodeContainer[dN].uInf=pT->displacement(nodeContainer[dN].P);
		  }
		  
		  int* row;  int* col;   double* val;
		  std::vector<int> rr, cc;
		  int dum, ii, ir , nz;
		  
		  double tS1=clock();
		  
		  //------ set the equation number for each dof based on the boundary conditions
		  setEquationNumbers();
		  
		  //-------- initiate the sparse system ------------------
		  
		  initiateSparseSystem(nz, rr , cc);
		  
		  //------------- necessary to build the sparse system (when using SparseLib++ library) -----------------
		  col = new int [nz];
		  val = new double [nz];   for(int i =0; i<nz; i++) {val[i] = 0.0e+00;}
		  
		  row = new int [sysDim+1];  row[0] = 0;  row[sysDim] = nz;
		  
		  ir = 0;  ii = 0; dum = 0;
		  for (unsigned int i = 0; i<cc.size() ; i++)
		  {
			  col[i] = cc[i];
			  if(rr[i] != dum )  {ir ++;   row[ir] = ii; }
			  ii++;
			  dum = rr[i] ;
		  }
		  
		  //------- the global linear system K.U = F -------------
		  CompRow_Mat_double K(sysDim, sysDim, nz, val, row, col);                // The global stiffness matrix
		  
		  delete val;     delete row;   delete col;
		  VECTOR_double U(sysDim);                     // The displacement vector 
		  VECTOR_double F(sysDim);                     // The force vector 
		  for(unsigned int i=0; i<sysDim; i++) {F(i) = 0.0e+00; U(i) =  0.0e+00;}
		  
		  //-------- Assemble the the linear system ---------
		  
		  assemble (dislocations_Stress, K,F,pT);    // assmble stiffness matrix and force vector from displacment BCs
			  
		  //assembleTractionForce (F);   // assemble force vector from surface traction
		  
		  if (dislocations_Stress) assembleInfiniteTraction (F,pT);   // assemble force vector from infinite medium surface traction

		  //-----------solve the linear system  to get the displacement U ---------------
		  
		  int max_itr = 1000;     double tol = 1.0e-8;
		  
		  double tS0=clock();
		  //ICPreconditioner_double D(K);            // generate preconditioning matrix
		  CompRow_ILUPreconditioner_double D(K);
		  
		  ir = CG(K, U , F , D, max_itr, tol);
		  
		  std::cout<<"SparseLib time= "<<(clock()-tS0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
		  std::cout<<"total BVP time= "<<(clock()-tS1)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
		  //-------------- distribute displacement to nodes -----------			
		  for (unsigned int i = 0 ; i < nodeContainer.size(); i++ )
		  {	
			  for (int j=0; j<dim; j++)
			  {
				  if (!nodeContainer[i].isBC(j)) {nodeContainer[i].u(j) = U(nodeContainer[i].eqnNumber(j));}
			  }
			  nodeContainer[i].displaceNode();
			  //std::cout<< "node: "<< nodeContainer[i].sID<< ": " <<nodeContainer[i].u.value[0]<< " " <<nodeContainer[i].u.value[1]<< " " <<nodeContainer[i].u.value[2] << std::endl;
		  }
		  
		  //------ if the loading is contact type, check for overconstraint nodes
#ifdef Contact_Loading
		while(pT->overConstraintNodes()) {reSolveBVP(true,pT);}
#endif	
		  
		  output();
		}
		
		  //--------------get stress on tets ---------------
		  //std::auto_ptr<Tetrahedron> pTet;
		  /*Tetrahedron* pTet;
		    Eigen::Matrix<double,dim,dim> stress;
		    for (int i = 0 ; i < tetContainer.size(); i++ )
		    {
		    pTet = &tetContainer[i];
		    
		    stress = pTet->getStress();
		    
		    std::cout << "stress =" << stress(0,0)<<" "<<stress(0,1)<<" "<<stress(0,2)<<" "<<stress(1,1)<<" "<<stress(1,2)<<" "<<stress(2,2)<<" "<< std::endl;
		    }*/	
		
		
		
		//===================================================================================
		// function to assemble the global linear system, solve it, and distribute nodal displacement
		//===================================================================================
		
		template<typename T>
		void reSolveBVP(bool dislocations_Stress,const T* const pT)
		{
			int* row;  int* col;   double* val;
			std::vector<int> rr, cc;
			int dum, ii, ir , nz;
			
			double tS1=clock();
			
			//------ set the equation number for each dof based on the boundary conditions
			setEquationNumbers();
			
			//-------- initiate the sparse system ------------------
			
			initiateSparseSystem(nz, rr , cc);
			
			//------------- necessary to build the sparse system (when using SparseLib++ library) -----------------
			col = new int [nz];
			val = new double [nz];   for(int i =0; i<nz; i++) {val[i] = 0.0e+00;}
			
			row = new int [sysDim+1];  row[0] = 0;  row[sysDim] = nz;
			
			ir = 0;  ii = 0; dum = 0;
			for (unsigned int i = 0; i<cc.size() ; i++)
			{
				col[i] = cc[i];
				if(rr[i] != dum )  {ir ++;   row[ir] = ii; }
				ii++;
				dum = rr[i] ;
			}
			
			//------- the global linear system K.U = F -------------
			CompRow_Mat_double K(sysDim, sysDim, nz, val, row, col);                // The global stiffness matrix
			
			delete val;     delete row;   delete col;
			VECTOR_double U(sysDim);                     // The displacement vector 
			VECTOR_double F(sysDim);                     // The force vector 
			for(unsigned int i=0; i<sysDim; i++) {F(i) = 0.0e+00; U(i) =  0.0e+00;}
			
			//-------- Assemble the the linear system ---------
			
			assemble (dislocations_Stress, K,F,pT);    // assmble stiffness matrix and force vector from displacment BCs
							
			if (dislocations_Stress)  copyTractionVector(F);
			
			//-----------solve the linear system  to get the displacement U ---------------
			
			int max_itr = 1000;     double tol = 1.0e-8;
			
			double tS0=clock();
			//ICPreconditioner_double D(K);            // generate preconditioning matrix
			CompRow_ILUPreconditioner_double D(K);
			
			ir = CG(K, U , F , D, max_itr, tol);
			
			std::cout<<"SparseLib time= "<<(clock()-tS0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
			std::cout<<"total BVP time= "<<(clock()-tS1)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
			//-------------- distribute displacement to nodes -----------			
			for (unsigned int i = 0 ; i < nodeContainer.size(); i++ )
			{	
				for (int j=0; j<dim; j++)
				{
					if (!nodeContainer[i].isBC(j)) {nodeContainer[i].u(j) = U(nodeContainer[i].eqnNumber(j));}
				}
				nodeContainer[i].displaceNode();
				//std::cout<< "node: "<< nodeContainer[i].sID<< ": " <<nodeContainer[i].u.value[0]<< " " <<nodeContainer[i].u.value[1]<< " " <<nodeContainer[i].u.value[2] << std::endl;
			}
	
		}
		
		
		//=================================================================================================
		// Copy the r.h.s. traction (infinite and external) force vector that was calculated inside solveBVP() function
		//=================================================================================================
		
		void copyTractionVector(VECTOR_double& Fm) {
			  
		  int qi;
		  
		  for (unsigned int ni = 0; ni < nodeContainer.size(); ni ++ ) {
		    
		    for (unsigned int dofi = 0; dofi<3; dofi++){
		      
		      qi = nodeContainer[ni].eqnNumber(dofi);
		      
		      if ( qi != -1) Fm(qi) = Fm(qi) + Finf[(ni*3)+dofi];
		      
		    }
		      
		  }
		}

		//===================================================================================
		// function to set the equation numbers for all dofs
		//===================================================================================
		
		void setEquationNumbers()
		{
			//unsigned int dim = 3;       // 3 dofs
			int num;
			bvpfe::Node<dim>* pNode;
			
			num = -1;
			
			for (unsigned int i = 0 ; i < nodeContainer.size(); i++ )
			{	
				pNode = &nodeContainer[i];
				for (int j=0; j<dim; j++)
				{
					if (!pNode->isBC(j))
					{
						num = num + 1;
						pNode->eqnNumber(j) = num;
					}
				}
			}
			
			sysDim = num + 1;          // set linear system dimension
			
			std::cout<< "number of dof = " << sysDim << std::endl;
			
			/*for (int i = 0 ; i < nodeContainer.size(); i++ )
			 {	
			 pNode = nodeContainer[i];
			 std::cout<< "Node=   "<< pNode->sID << ":  " << pNode->u.isBC[0] << "  " << pNode->u.isBC[1] << "  "<< pNode->u.isBC[2] << "  "<< pNode->u.eqnNumber[0] << "  " << pNode->u.eqnNumber[1] << "  "<< pNode->u.eqnNumber[2] << "  "<< std::endl;
			 }*/
			
			
		}
		
		//===================================================================================
		// function to assemble the global stiffness matrix
		//===================================================================================
		
		template<typename T>
		void assemble (bool disStress, CompRow_Mat_double & Km, VECTOR_double& Fm, const T* const pT)
		{
			//bvpfe::Tetrahedron* pTet;
			//std::vector<Vector> tetMat;         // 12x12 stiffness matrix for each tetrahedron
			Eigen::Matrix<double,12,12> tetMat;         // 12x12 stiffness matrix for each tetrahedron
			
			//std::vector<Vector> Np;             // shape functions derivatives mapped to actual element
			Eigen::Matrix<double,dim,Tetrahedron::Nnodes> Np;   // shape functions derivatives mapped to actual element
			
			int qi, qj;
			int ie , je;
			
			double wv;
			
			for (unsigned int i=0; i<tetContainer.size() ; i++) {
				
				//pTet = &tetContainer[i];
				
				tetMat = tetContainer[i].getElementStiffness(Np, wv);
				
				//---------- assemble tetMat to the global matrix ---------------------------
				for ( unsigned int ai = 0 ; ai < 4 ; ai ++)            // no. of nodes per Tet
				{
					ie= ai*3;
					for ( unsigned int in = 0 ; in < 3 ; in ++)    //no. of DOFs per node
					{
						qi = tetContainer[i].eleNodes[ai]->eqnNumber[in];
						if (qi == -1) continue;
						
						for ( unsigned int bi = 0 ; bi < 4 ; bi ++)
						{
							je = bi*3 ; 
							for ( unsigned int jn = 0 ; jn < 3 ; jn ++)
							{
								qj = tetContainer[i].eleNodes[bi]->eqnNumber(jn);
								if (qj == -1)
								{
									Fm(qi)=Fm(qi)-(tetMat(ie+in,je+jn)*tetContainer[i].eleNodes[bi]->bcValue[jn]);
									//std::cout<< "out111=== " << tetMat(ie+in,je+jn) << "  " << tetContainer[i].eleNodes[bi]->u.bcValue[jn] << std::endl;
									continue;
								}
								
								Km.set(qi,qj) = Km(qi,qj) + (tetMat(ie+in,je+jn));
							}
						}
					}
				}
				
				
				//------------- assemble the dislocation force vector, if required ----------
				
				/*if(disStress)
				 {
				 
				 Eigen::Matrix<double,dim,Tetrahedron::Nnodes> temp=tetContainer[i].get_dislocationForceVector<4,T>(pT,Np);
				 
				 
				 //----------- assemble in the global force vector ------------------
				 for ( unsigned int in = 0 ; in < Tetrahedron::Nnodes ; in ++)            // no. of nodes per Tet
				 {						
				 for ( unsigned int id = 0 ; id < dim ; id ++)    //no. of DOFs per node
				 {
				 qi = tetContainer[i].eleNodes[in]->u.eqnNumber[id];
				 if (qi == -1) continue;
				 
				 Fm(qi)=Fm(qi)-temp(id,in);
				 
				 }
				 }
				 }*/
				
			}
		}
		
		
		//===================================================================================
		// function to initiate the sparse linear system
		//===================================================================================
		
		void initiateSparseSystem(int& nz, std::vector<int>& rr, std::vector<int>& cc )
		{
			int qi, qj;
			
			Node<dim>* pNodei;
			Node<dim>* pNodej;
			
			// ---- sort the neighbor vector for each node. This mkae sit easy to initiate sparse system ---
			
			for (unsigned int i=0; i<nodeContainer.size() ; i++) {
				nodeContainer[i].sortNeighbors();
				//				pNodei = &nodeContainer[i];
				//				pNodei->sortNeighbors();
				
				/*std::cout<< pNodei->sID<< " : " ;
				 for (int j=0; j<pNodei->neighbor.size() ; j++)
				 {
				 std::cout<< pNodei->neighbor[j]->sID<< "  " ;
				 }
				 std::cout<< std::endl;*/
			}
			
			//---- loop over all nodes to get the positions of nonzero elements------
			
			for (unsigned int i=0; i<nodeContainer.size() ; i++) {
				
				pNodei = &nodeContainer[i];
				
				for (int ii = 0; ii<dim; ii++)
				{
					qi = pNodei->eqnNumber(ii);
					if (qi == -1) continue;
					
					for (unsigned int j = 0; j<pNodei->neighbor.size(); j++)
					{
						pNodej = pNodei->neighbor[j];
						
						for (int jj = 0; jj<dim; jj++)
						{
							qj = pNodej->eqnNumber(jj);
							if (qj == -1) continue;
							
							rr.push_back(qi);     cc.push_back(qj);
						}
					}
				}
				
			}
			
			nz = rr.size();
			
			std::cout << "nz = " << nz << std::endl;
			
			/*for (int i = 0; i<rr.size() ; i++)
			 {
			 std::cout << rr[i] << "  "<< cc[i] << std::endl;
			 }*/	
			
			
		}
		
		//===================================================================================
		// function to assemble the force vector from infinite medium surface traction
		//===================================================================================
		template<typename T>
		void assembleInfiniteTraction (VECTOR_double& Fm , const T* const pT) {
			
			for (unsigned int i=0 ; i< Finf.size() ; i++ )  Finf[i]= 0.0;
		  
			Eigen::Matrix<double,dim,3> TriVec;       // 3 is the number of nodes per triangle
			
			unsigned int nID;
			Triangle* pTri;
			std::vector<Triangle*>::iterator itt;
				
			for (itt= triContainer.begin() ; itt != triContainer.end(); itt++ ) {
				
				pTri = *itt; 
				
				TriVec = pTri->getTriInfiniteForce<3,T>(pT);
				
				//------------- assemble the infinite traction vector --------------
				int  qi;
				for (int ai = 0 ; ai < 3; ai++ )             // loop over the 3 surface element's nodes
				{
					nID = pTri->eleNodes[ai]->sID;        // node global index
					
					for  (int in = 0 ; in < 3; in++ )   // loop over the 3 dofs
					{
						Finf[(nID*3)+in] = Finf[(nID*3)+in] + TriVec(in,ai);  // -ve sigen for the infinite field is considered in Triangle::dislocationStressKernel
					  
						qi = pTri->eleNodes[ai]->eqnNumber(in);
						if (qi != -1) Fm(qi) = Fm(qi) + TriVec(in,ai);  // -ve sigen for the infinite field is considered in Triangle::dislocationStressKernel
					}
				}
				
			}	  
			
			//boost::ptr_map<size_t,Face>::iterator itf;
// 			for (itf= faceContainer.begin() ; itf != faceContainer.end(); itf++ ) {
// 				
// 				for (itt= itf->second->triContainer.begin() ; itt != itf->second->triContainer.end(); itt++ ) {
// 					
// 					pTri = itt->second; 
// 					
// 					TriVec = pTri->getTriInfiniteForce<3,T>(pT);
// 					
// 					//------------- assemble the infinite traction vector --------------
// 					int ie, qi;
// 					for (int ai = 0 ; ai < 3; ai++ )             // loop over the 3 surface element's nodes
// 					{
// 						ie = ai*3;
// 						
// 						for  (int in = 0 ; in < 3; in++ )   // loop over the 3 dofs
// 						{
// 							qi = pTri->eleNodes[ai]->eqnNumber(in);
// 							Finf[qi] = Finf[qi] + TriVec(in,ai);  // -ve sigen for the infinite field is considered in Triangle::dislocationStressKernel
// 							if (qi == -1) continue;
// 							Fm(qi) = Fm(qi) + TriVec(in,ai);  // -ve sigen for the infinite field is considered in Triangle::dislocationStressKernel
// 						}
// 					}
// 					
// 				}	  
// 			}
			
		}
		
		//===================================================================================
		// function to assemble the force vector from surface traction
		// we assume here uniform traction over the surface, so one quadrature point is used
		//===================================================================================
		
		void assembleTractionForce (VECTOR_double& Fm )
		{
			bvpfe::Face* pFace;
			std::map<bvpfe::Face*,std::vector<double> >::iterator itf;
			std::vector<double> tr;
			
			bvpfe::Triangle* pTri;
			std::map<size_t,bvpfe::Triangle*>::iterator itt;
			
			std::vector<double> triVec;
			triVec.resize(9,0.0e+00);     //The force vector for each triangle
			
			int ie , qi;
			
			//--------- loop over all faces subjected to traction ----------
			for (itf= facesWtraction.begin() ; itf != facesWtraction.end(); itf++ )
			{
				pFace = itf->first;
				tr = itf->second;
				
				//std::cout<< "traction " << tr[0] << " " << tr[1] << " " << tr[2] << std::endl;
				
				// ---------- loop over all triangles of the face ------------
				for (itt= pFace->triContainer.begin() ; itt != pFace->triContainer.end(); itt++ )
				{
					pTri = itt->second;
					
					triVec = pTri->getForceVec(tr);
					
					for (int ai = 0 ; ai < 3; ai++ )             // loop over the 3 surface element's nodes
					{
						ie = ai*3;
						
						for  (int in = 0 ; in < 3; in++ )   // loop over the 3 dofs
						{
							qi = pTri->eleNodes[ai]->eqnNumber(in);
							if (qi == -1) continue;
							Fm(qi) = Fm(qi) +  triVec[ie+in]; 
						}
					}
				}
			}
			
			
		}
		
		//===================================================================================
		// function to calculate the stress field at any given point
		//==================================================================================
		Eigen::Matrix<double,dim,dim> stressAt(const VectorDim& P)
		{
			Eigen::Matrix<double,dim,dim> stress = Eigen::Matrix<double,dim,dim>::Zero(); 
			isTetrahedronType isT = findIncludingTet(P);
			
			/*std::cout << "tetrahedron number : " << pTet->sID<< std::endl;
			 std::cout<< "Point: "<< P[0]<< " "<< P[1]<< " "<< P[2]<< " "<<std::endl;
			 
			 std::cout<< pTet->eleNodes[0]->P[0]<< " "<< pTet->eleNodes[0]->P[1]<< " "<< pTet->eleNodes[0]->P[2]<< " "<<std::endl; 
			 std::cout<< pTet->eleNodes[1]->P[0]<< " "<< pTet->eleNodes[1]->P[1]<< " "<< pTet->eleNodes[1]->P[2]<< " "<<std::endl; 
			 std::cout<< pTet->eleNodes[2]->P[0]<< " "<< pTet->eleNodes[2]->P[1]<< " "<< pTet->eleNodes[2]->P[2]<< " "<<std::endl; 
			 std::cout<< pTet->eleNodes[3]->P[0]<< " "<< pTet->eleNodes[3]->P[1]<< " "<< pTet->eleNodes[3]->P[2]<< " "<<std::endl; */
			
			if (isT.first) {stress = isT.second->getStress();}
			
			return stress;
		}
		
		//===================================================================================
		// function to search for the including tetrahedron given the new position for the point, 
		// old point's tetrahedron index, and search direction vector
		//===================================================================================
		
		void SearchMovingNode (mmdl::SearchData<dim> & data){
			
		  // --------------NOT boundary node (coming from inside or outside) -----------------------
		  if(data.nodeMeshLocation != 2) { 
		    
		    data.newMeshID = data.currentMeshID;	    // initialize
		    int dI, sI;
			Eigen::Matrix<double,4,1>::Index ii;
			  // , jj;
		    VectorDim interP;                              // should save the velocity-surface intersection point
		    Triangle* pTri;
		    Tetrahedron* pTet;
			
		    Eigen::Matrix<double,4,1>  Bary;
		    
		    while (!data.found){
		      bool insideTet =  tetContainer[data.newMeshID].isInsideTet(data, sI, Bary);
		      
		      if(insideTet) break;
		      
		      //--------if node is found outside, find intersection with boundary, and set it as boundary node ----
		      if(data.nodeMeshLocation == 0){
			pTet = &tetContainer[data.newMeshID];
			
			while(!data.found){
			  
			  //---- find the vector-surface intersection point ---------
			  data.newMeshID = pTet->sID;
			  
			  interP = getSurfaceIntersection(pTet->TetSurfTris.find(sI)->second->outNormal,
							pTet->TetSurfTris.find(sI)->second->eleNodes[0]->P, -data.normalizedDir, data.P);
			  
			  //------ calculate the barycentric coordinates for this intersection point -------
			  double baryMin = pTet->getBarycentric(interP).minCoeff(&ii);
			  
			  //------ check if the point is inside, or just outside (but too close) the tetrahedron surface --> DONE ---------
			  //------ if (ii==sI) so the point is on the tetrahedron surface, with small error from calculations --- 
			  if((baryMin>=0.0)||(ii==sI)){
			    data.found=true;
			    data.projectedP = interP;
			    data.outwardFaceNormal= pTet->TetSurfTris.find(sI)->second->outNormal;
			    data.nodeMeshLocation = 2;      // boundary node
			    data.triIndex = pTet->TetSurfTris.find(ii)->second->sID;
			    break;
			  }
			  
			  //-------- if the intersection point is not in this tetrahedron, move to the next -----
			  
			  // --- find the index (from 0 to 2) for the point with lowest bary for the surface triangle ----
			  
			  for(unsigned int k=0; k<3; k++){
			    if(pTet->TetSurfTris.find(sI)->second->eleNodes[k]->sID == pTet->eleNodes[ii]->sID) {dI = k; break;}
			  }
			  
			  pTri = pTet->TetSurfTris.find(sI)->second->neighbor[dI];    // next neighbor triangle			  
			  pTet = &tetContainer[pTri->neighTetIndx];                   // new tetrahedron to move to
			  
			  for (std::map<size_t,Triangle*>::iterator it = pTet->TetSurfTris.begin(); it!= pTet->TetSurfTris.end(); it++){
			    if (it->second->sID == pTri->sID){sI = it->first; break;}
			  }
						
			}			 
		      }
		    }
		  }
		  
		  // ----------------- boundary node ---------------------
		  
		 else {   
		    
		    VectorDim org = data.P - data.Dir;      // old position of the node
		    unsigned int iTri;
		    
		    if(triContainer[data.triIndex]->intersectWithLine(org,data.P,data.projectedP,iTri)) {
		      data.outwardFaceNormal= triContainer[iTri]->outNormal;
		      data.newMeshID = triContainer[iTri]->neighTetIndx;
		      data.triIndex = iTri;
		      //std::cout<< "new mesh ID ================== "<<data.newMeshID  << std::endl;
		    }
		    else{
		      data.newMeshID = data.currentMeshID;
		      data.outwardFaceNormal= triContainer[data.triIndex]->outNormal;
		    }
		    data.found=true;
		 }
		  
		  
	      }
		//===================================================================================
		// function to calculate the intersection point between a line and a plane, given
		// the plane normal, point on the plane, line direction, point on the line
		//==================================================================================
		
		VectorDim getSurfaceIntersection (VectorDim pD, VectorDim pP,const VectorDim lD,const VectorDim lP){
			
			double dm = lD.dot(pD);
			
			if(dm == 0.0) assert(0&&"Dislocation Node is moving parallel to the boundary: unable to get intersection point");
			
			double d = ((pP-lP).dot(pD))/dm;
			
			return (lP+d*lD);
			
		}
		
		//===================================================================================
		// function to return a pointer to the including tetrahedron for any given point, 
		// or return null if the point is outside the domain
		//==================================================================================
		isTetrahedronType findIncludingTet(VectorDim P){
			bool found = false;
			int ci; // index of neighboor tet. -1 means outside.
			int ni; // 
			
			isTetrahedronType temp=std::make_pair(true, (Tetrahedron*) NULL);
			
			//----- starting Tet, random start point ------------
			ci = int(tetContainer.size()/2);     
			
			while ((!found)&&(ci>=0))
			{
				temp.second = &tetContainer[ci];
				ci = temp.second->nextNeighbor(P,found,ni);    // returns -10 if the point is inside the Tet
			}
			
			temp.first=found;
			
//			if (!found)   // this means it is outside the domain
//			{
//				temp.first=false;
//			//	temp.second = NULL;
//				//std::cout << "============= ERROR: Searching for point outside the domain ==============" << std ::endl; 
//				//assert(0);
//			}
			
			return temp;
		} 
		
		//===================================================================================
		// function to detect the intersection line between a slip plane, and the mesh surface
		//==================================================================================
		void get_planeMeshIntersection(const VectorDim& x0, const VectorDim& n, std::vector<std::pair<VectorDim,VectorDim> >& collisionContainer){
		
			
			double tol = 1.0e-10;
			
			for(unsigned int i=0; i<triContainer.size();i++){
				std::vector<VectorDim> intersectionPoints;
				
				for(unsigned int j = 0; j< 3; j++){
					VectorDim v0 = triContainer[i]->eleNodes[j]->P;
					VectorDim v1 = triContainer[i]->eleNodes[(j+1)%3]->P;
					double u = getSurfaceIntersectionPar (x0,n,v0,v1);
					
					if(u>=0 && u<=1.0){ 
						VectorDim P = v0 + u*(v1-v0);
						

						bool isDifferent=true;
						for (unsigned int k=0;k<intersectionPoints.size();k++){
							isDifferent*= (P-intersectionPoints[k]).squaredNorm()>tol;
						}
						
						if (isDifferent){
							intersectionPoints.push_back(P);
						}
					}
				}
				
				assert(intersectionPoints.size()<3);
				
				if (intersectionPoints.size()==2){
					collisionContainer.push_back(std::make_pair(intersectionPoints[0],intersectionPoints[1]));
				
				}
				
			}
			
			
			
		}
		
		
		//===================================================================================
		// function to calculate the intersection point between a line and a plane, given
		// the plane normal, point on the plane, line direction, point on the line
		//==================================================================================
		
		double getSurfaceIntersectionPar (VectorDim x0, VectorDim n,const VectorDim v0,const VectorDim v1){
			
			return ((x0-v0).dot(n)) / (v1-v0).dot(n);
		}		
		
		//===================================================================================
		// function to output Tetrahedron mesh
		//==================================================================================
		void output(){
			
			mmdl::SequentialOutputFile<'N',true> nodesFile;
			for (NodeContainerType::const_iterator iter=nodeContainer.begin();iter!=nodeContainer.end();++iter){
				nodesFile<< iter->sID<<"	" << iter->Pc.transpose()<<"	"<<iter->isBoundaryNode<<std::endl;
			}
			
			mmdl::SequentialOutputFile<'T',true> tetFile;
			for (TetContainerType::const_iterator iter=tetContainer.begin();iter!=tetContainer.end();++iter){
				for (unsigned int k=0;k<iter->eleNodes.size()-1;++k){
					for (unsigned int j=k+1;j<iter->eleNodes.size();++j){
					 	tetFile<< iter->eleNodes[k]->sID<<"	"<<iter->eleNodes[j]->sID<<" 0"<<std::endl;
					}
				}
				//tetFile<<std::endl;
			}
			
		}
		
		//==================================================================================
		// function to remove all previous traction and displacement boundary conditions
		// called before setting new ones
		//==================================================================================
		
		void removeBoundaryConditions(){
		  
		  for (unsigned int i =0; i<nodeContainer.size(); i++){
		    if(nodeContainer[i].isBoundaryNode){
		      nodeContainer[i].traction = VectorDim::Zero();         // set traction to be zero
		      nodeContainer[i].remove_BCs();                        // set displacement BC to be false
		    }
		  } 
		  
		}
		
		
		//===========================================================================================
		
		
		
	};
	
	
	
}  //  namespace bvpfe
#endif