/* This file is part of model, the Mechanics of Defects Evolution  Library.
 *
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_SEARCHDATA_H_
#define model_SEARCHDATA_H_


#include <Eigen/Core>

namespace model {
		
		template <short unsigned int dim>
		struct SearchData {
			
			typedef Eigen::Matrix<double,dim,1> VectorDim;
			
			const VectorDim P;
			const VectorDim normalizedDir;
			const int currentMeshID;
			bool found;
			int nodeMeshLocation;      //=0 outside,   =1 inside,   =2 on mesh boundary

			
			VectorDim Dir;
			
			//bool isOutside;
			unsigned int triIndex;

			int newMeshID;
			VectorDim outwardFaceNormal;
			VectorDim projectedP;
			
			SearchData(const VectorDim& Pin, const VectorDim& dirin, const int& mID, const int& nodeMeshLocation_in, const int& triIndex_in, const VectorDim& boundaryNormal_in) : P(Pin), 
			normalizedDir(dirin.normalized()), 
			currentMeshID(mID),
			found(false),
			nodeMeshLocation(nodeMeshLocation_in),
			Dir(dirin),
			triIndex(triIndex_in),
			outwardFaceNormal(boundaryNormal_in)
			{
			  //outwardFaceNormal = VectorDim::Zero();
			  projectedP = VectorDim::Zero();
			}
			
			
			SearchData(const VectorDim& Pin) : P(Pin), 
			normalizedDir(VectorDim::Zero()), 
			currentMeshID(0),
			found(false),
			nodeMeshLocation(-1),
			Dir(VectorDim::Zero()),
			triIndex(0){
			  outwardFaceNormal = VectorDim::Zero();
			  projectedP = VectorDim::Zero();
			}
			
			
			/* some other face info*/
		};
		
} // namespace model
#endif
