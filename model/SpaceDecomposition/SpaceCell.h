/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPACECELL_H_
#define model_SPACECELL_H_

#include <assert.h>
#include <math.h>
#include <set>
#include <map>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>
#include <model/Math/CompileTimeMath/Pow.h>
#include <model/SpaceDecomposition/SpaceCellObserver.h>
#include <model/SpaceDecomposition/NeighborShift.h>
#include <model/Utilities/CompareVectorsByComponent.h>


namespace model {
	
	/********************************************************************************************/
	/********************************************************************************************/
	/*! \brief A dim-dimensional cell occupying the spatial region cellID<= x/cellSize < (cellID+1).
	 *  SpaceCell is aware off all ParticleType objects present inside it.
	 */
	template<typename ParticleType, short unsigned int dim, double & cellSize>
	struct SpaceCell : boost::noncopyable,
	/*              */ private SpaceCellObserver<SpaceCell<ParticleType,dim,cellSize>,dim,cellSize>{
		
		typedef SpaceCell<ParticleType,dim,cellSize> SpaceCellType;
		typedef SpaceCellObserver<SpaceCellType,dim,cellSize> SpaceCellObserverType;
		typedef typename SpaceCellObserverType::CellMapType  CellMapType;
		typedef typename SpaceCellObserverType::VectorDimD  VectorDimD;
		typedef typename SpaceCellObserverType::VectorDimI  VectorDimI;
        //        typedef typename SpaceCellObserverType::CellMapType  CellMapType;
        //        typedef std::map<VectorDimI,const SpaceCellType* const,CompareVectorsByComponent<int,dim> >  CellMapType;
        
		typedef std::set<const ParticleType*> ParticleContainerType; // PTR COMPARE IS NOT NECESSARY
        
        
        
        // typedef Eigen::Matrix<double,dim,dim>  MatrixDimD;	// remove this with Dislocation Stuff
        // MatrixDimD alpha;	// the dislocation density tensor !! remove this with Dislocation Stuff
        
        typedef NeighborShift<dim,1> NeighborShiftType;
        
        
        CellMapType nearCells;
        CellMapType  farCells;
        
        /* isNearCell *******************************************/
        bool isNearCell(const VectorDimI& otherCellID) const {
            bool temp(false);
            for (int k=0;k<neighborCellIDs.cols();++k){
                if(neighborCellIDs.col(k)==otherCellID){ // cell is a neighbor
                    temp=true;
                    break;
                }
            }
            return temp;
        }
        
        
	public:
        
        //static int nearestNeighborOrder;
		
		//! The container of pointers to particles in this cell
		ParticleContainerType particleContainer;
		
        //! The container of pointers to particles in this and first neighbor cells
		//ParticleContainerType neighborParticleContainer;
		
        //! The ID of this cell, defining the dim-dimensional spatial region cellID<= x/cellSize < (cellID+1).
		const VectorDimI cellID;
		//! The cellID(s) of the neighboring SpaceCell(s) (in column)
		const Eigen::Matrix<int,dim, NeighborShiftType::Nneighbors> neighborCellIDs;
        
		/* Constructor *******************************************/
		SpaceCell(const VectorDimI& cellID_in) : cellID(cellID_in),neighborCellIDs(NeighborShiftType::neighborIDs(cellID)){
			//! 1- Adds this to static SpaceCellObserver::cellMap
			assert(this->cellMap.insert(std::make_pair(cellID,this)).second && "CANNOT INSERT SPACE CELL IN STATIC cellMap.");
			
//            //! 2- Loops over neighboring SpaceCell(s) and copies the content of their particleContainer into this neighborParticleContainer
//			for (int n=0;n<NeighborShiftType::Nneighbors;++n){
//				typename CellMapType::const_iterator iter(this->cellMap.find(neighborCellIDs.col(n)));
//				if (iter!=this->cellMap.end()){ // neighbor cell exists
//					for (typename ParticleContainerType::iterator iterP=iter->second->particleContainer.begin(); iterP!=iter->second->particleContainer.end();++iterP){
//						assert(neighborParticleContainer.insert(*iterP).second && "CANNOT COPY NEIGHBOR PARTICLE");
//					}
//				}
//			}
            
            
            for (typename CellMapType::const_iterator cellIter=this->begin();cellIter!=this->end();++cellIter){
                if (isNearCell(cellIter->second->cellID)){
                    assert(                  nearCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN NEARCELLS");
                    if(cellID!=cellIter->second->cellID){
                        assert(cellIter->second->nearCells.insert(std::make_pair(         cellID,            this)).second && "CANNOT INSERT THIS IN NEARCELLS");
                    }
                }
                else{
                    assert(                  farCells.insert(std::make_pair(cellIter->first,cellIter->second)).second && "CANNOT INSERT CELL IN FARCELLS");
                    assert(cellIter->second->farCells.insert(std::make_pair(         cellID,            this)).second && "CANNOT INSERT THIS IN FARCELLS");
                }
            }
            
		}
		
		/* Destructor *******************************************/
		~SpaceCell(){
			//! Removes this from static SpaceCellObserver::cellMap
			this->cellMap.erase(cellID);
			assert(particleContainer.empty() && "DESTROYING NON-EMPTY SPACE CELL.");
		}
		
		/* addParticle *******************************************/
		void addParticle(const ParticleType* const pP){
			//! 1- Adds pP to the particleContainer
			assert(particleContainer.insert(pP).second && "CANNOT INSERT PARTICLE IN SPACECELL");
			
//            //! 2- Updates the first nieghboring cells calling addNeighborParticle(pP)
//			for (int n=0;n<NeighborShiftType::Nneighbors;++n){
//				typename CellMapType::iterator iter(this->cellMap.find(neighborCellIDs.col(n)));
//				if (iter!=this->cellMap.end()){ // neighbor cell exists
//					iter->second->addNeighborParticle(pP);
//				}
//			}
		}
		
//		/* addNeighborParticle ***********************************/
//		void addNeighborParticle(const ParticleType* const pP){
//			assert(neighborParticleContainer.insert(pP).second && "CANNOT INSERT NEIGHBOR PARTICLE IN SPACECELL");
//			
//		}
		
		/* removeParticle ****************************************/
		void removeParticle(const ParticleType* const pP){
			//! 1- Removes pP to the particleContainer
			assert(particleContainer.erase(pP)==1 && "CANNOT ERASE PARTICLE FROM particleContainer.");
			
//            //! 2- Updates the first nieghboring cells calling removeNeighborParticle(pP)
//			for (int n=0;n<NeighborShiftType::Nneighbors;++n){
//				typename CellMapType::iterator iter(this->cellMap.find(neighborCellIDs.col(n)));
//				if (iter!=this->cellMap.end()){ // neighbor cell exists
//					iter->second->removeNeighborParticle(pP);
//				}
//			}
		}
		
//		/* removeNeighborParticle ********************************/
//		void removeNeighborParticle(const ParticleType* const pP){
//			//! Removes pP to the neighborParticleContainer
//			assert(neighborParticleContainer.erase(pP)==1 && "CANNOT ERASE PARTICLE FROM neighborParticleContainer.");
//		}
		
//		/* check ************************************************/
//		void check() const{
//			unsigned int temp(0);
//			for (int n=0;n<NeighborShiftType::Nneighbors;++n){
//				typename CellMapType::iterator iter(this->cellMap.find(neighborCellIDs.col(n)));
//				if (iter!=this->cellMap.end()){ // neighbor cell exists
//					temp+=iter->second->particleContainer.size();
//				}
//			}
//			assert(temp==neighborParticleContainer.size());
//		}
        
        /* nearCellBegin ***************************************/
        typename CellMapType::const_iterator nearCellsBegin() const {
            return nearCells.begin();
        }
        
        /* nearCellEnd ***************************************/
        typename CellMapType::const_iterator nearCellsEnd() const {
            return nearCells.end();
        }
        
        /* nearCellBegin ***************************************/
        typename CellMapType::const_iterator farCellsBegin() const {
            return farCells.begin();
        }
        
        /* nearCellEnd ***************************************/
        typename CellMapType::const_iterator farCellsEnd() const {
            return farCells.end();
        }
        
        /* particleBegin ***************************************/
        typename ParticleContainerType::const_iterator particleBegin() const {
            return particleContainer.begin();
        }
        
        /* particleEnd *****************************************/
        typename ParticleContainerType::const_iterator particleEnd() const {
            return particleContainer.end();
        }
        
        
	};
    
    
    // Declare Static data
    //    template<typename ParticleType, short unsigned int dim, double & cellSize>
    //    int SpaceCell<ParticleType,dim,cellSize>::nearestNeighborOrder=1;
	
	////////////////////////////////////////////////////////////////////////////////
}	// close namespace model
#endif

