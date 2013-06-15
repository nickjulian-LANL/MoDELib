/* This file is part of finite element solution of BVP attached with model "the Mechanics of Defects Evolution Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_dof_H_
#define bvpfe_dof_H_

namespace bvpfe{

	template<short unsigned int dim>
	class Dof {
		
		#include <model/BVP/commonTypeDefs.h>	
		public:
			VectorDim u;                                       // the dof value
			Eigen::Matrix<bool,dim,1> isBC;                    // =false if no BCs, =true if there is BCs
			VectorDim bcValue;                                 // BC value 
			Eigen::Matrix<int,dim,1> eqnNumber;                // the position of this dof in the global matrix. =-1 if there is BC
			VectorDim uInf;                                    // infinite medium displacement field at this node
			VectorDim uVir;                                    // infinite medium displacement field induced by Virtual dislocation segments
				

		public :

			Dof(){
			  remove_BCs();                                                      // initialize BCs to be false
			  uInf = VectorDim::Zero();
			  uVir = VectorDim::Zero();
			  //for(unsigned int i=0; i<dim;i++) {eqnNumber(i) = -1;}             // initialize the equation number to be -1             
			}
			

			//----------- set the dof value ---------------
			//inline void setDof (unsigned int i,  )

			// ---------- get the dof value --------------
			VectorDim & getDof()
			{
				return u;
			} 

			// --------- set the boundary condition for a DOF ------------
			void setBC (unsigned int i , double val)
			{
			  //if(!this->isBoundaryNode) assert(0&&"setting boundary condition for non-boundary node");
			  isBC(i) = true;
			  bcValue(i) = double(val-uInf(i));
			  u (i) = double(val-uInf(i));
			}
			
			// --------- remove the boundary condition for a DOF ------------
			void removeBC (unsigned int i)
			{
			  //if(!this->isBoundaryNode) assert(0&&"removing boundary condition for non-boundary node");
			  isBC(i) = false;
			  eqnNumber(i) = -1;
			}
			
			// --------- remove all boundary conditions from the node ------------
			void remove_BCs()
			{
				for(unsigned int i=0; i<dim;i++) {isBC(i) = false;  eqnNumber(i) = -1;}
			}
		
			

	};
		
	

}  //  namespace bvpfe
#endif