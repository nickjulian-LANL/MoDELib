/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPACECELLOBSERVER_H_
#define model_SPACECELLOBSERVER_H_

#include <map>
#include <memory> // std::shared_ptr (c++11)
//#include <Eigen/Dense>
//#include <model/Utilities/CompareVectorsByComponent.h>
//#include <model/SpaceDecomposition/SpatialCellTraits.h>
//#include <boost/ptr_container/ptr_map.hpp>
#include <model/SpaceDecomposition/SpatialCell.h>
#include <model/SpaceDecomposition/CellShift.h>


namespace model {
	
	
	/********************************************************************************************/
	/********************************************************************************************/
	template<typename ParticleType,short unsigned int dim>
	struct SpatialCellObserver
    {
		
//        typedef SpatialCell<ParticleType,dim> SpatialCellType;
//        typedef typename SpatialCellType::VectorDimD VectorDimD;
//        typedef typename SpatialCellType::CellIdType CellIdType;
//		typedef std::map<Eigen::Matrix<int,dim,1>, SpatialCell<ParticleType,dim>* const,CompareVectorsByComponent<int,dim> >  CellMapType;
//        typedef typename SpatialCellType::SharedPtrType SharedPtrType;
//        typedef typename SpatialCellType::isCellType isCellType;

        
//        SpatialCellObserver()
//        {
//
//        }
//        
//        SpatialCellObserver(const SpatialCellObserver& other)
//        {
//            
//        }
        
//    public:
        
        

        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<int,dim,1> CellIdType;
        typedef SpatialCell<ParticleType,dim> SpatialCellType;
		typedef std::map<CellIdType, SpatialCellType* const,CompareVectorsByComponent<int,dim> >  CellMapType;
		typedef std::shared_ptr<SpatialCellType> SharedPtrType;
        typedef std::pair<bool,SpatialCellType* const> isCellType;
            typedef CellShift<dim,1>    CellShiftType;
        

        
        //! The size of a SpatialCell
        static double cellSize;
        
		/* begin() ***************************************************/
		static typename CellMapType::const_iterator cellBegin()
        {/*! \returns A const_iterator to the first SpatialCellType cell.
          */
			return cellMap.begin();
		}
		
		/* end() *****************************************************/
		static typename CellMapType::const_iterator cellEnd()
        {/*! \returns A const_iterator to the past-the-last SpatialCellType cell.
          */
			return cellMap.end();
		}
        
        
        static size_t totalCells()
        {/*! \returns The number of observed SpatialCellType cells.
          */
            return cellMap.size();
        }
        
        /* getCellIDByPosition ************************************************/
		static CellIdType getCellIDByPosition(const VectorDimD& P)
        {/*! \returns The CellIdType ID of the cell that contains P. The ID
          *  satisfies cellID <= P/cellSize < (cellID+1).
          */
			return floorEigen<dim>(P/cellSize);
		}
        
		/* getCellByID ********************************************************/
		static SharedPtrType getCellByID(const CellIdType& cellID)
        {/*! @param[in] cellID The ID of the cell.
          *  \returns If a cell with ID=cellID exists, a shared-pointer to that
          *  cell is returned. Otherwise, a shared-pointer to a new cell with
          *  ID=cellID is returned.
          */
			typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? (*(iter->second->particleContainer.begin()))->pCell : SharedPtrType(new SpatialCellType(cellID));
//			return (iter!=cellMap.end())? (*((*(iter->second))->particleContainer.begin()))->pCell : SharedPtrType(new SpatialCellType(cellID));
		}
		
		/* getCellByPosition **************************************************/
		static SharedPtrType getCellByPosition(const VectorDimD& P)
        {/*! @param[in] P The position vector.
          *  \returns If a cell satisfying cellID <= P/cellSize < (cellID+1)
          *   exists, a shared-pointer to that cell is returned. Otherwise, a
          *   shared-pointer to a new cell satisfying cellID <= P/cellSize < (cellID+1)
          *   is returned.
          */
			return getCellByID(getCellIDByPosition(P));
		}
        
        /* isCell *************************************************************/
        static isCellType isCell(const CellIdType& cellID)
        {
            typename CellMapType::const_iterator iter(cellMap.find(cellID));
			return (iter!=cellMap.end())? std::make_pair(true,iter->second) : std::make_pair(false,( SpatialCellType*) NULL);
        }
        
        static const CellMapType& cells()
        {
            return cellMap;
        }
        
        /* neighborCells ******************************************************/
        static CellMapType neighborCells(const VectorDimD& P)
        {
            const CellIdType cellID(getCellIDByPosition(P));
            const Eigen::Matrix<int,dim, CellShiftType::Nneighbors> neighborCellIDs(CellShiftType::neighborIDs(cellID));
            
            CellMapType temp;
            
            for (unsigned int c=0;c<CellShiftType::Nneighbors;++c)
            {
                
                isCellType isC(isCell(neighborCellIDs.col(c)));
                if (isC.first)
                {
                    model_execAssert(temp.insert(std::make_pair(isC.second->cellID,isC.second)),.second,"CANNOT INSERT CELL IN NEIGHBORCELLS");

                }

            }
            
            return temp;
        }
    
        
        /*****************************************/
        template <class T>
        friend T& operator<< (T& os, const SpatialCellObserver<SpatialCellType,dim>& sCO)
        {/*! Operator << uses ParticleType-specific operator <<
          */
            for (typename CellMapType::const_iterator cIter=sCO.begin();cIter!=sCO.end();++cIter)
            {
                os << (*cIter->second) << std::endl;
            }
            return os;
        }
        
        
    private:
        
        friend class SpatialCell<ParticleType,dim>;
        
        static  CellMapType cellMap;
    };
    
    /////////////////////////////
    // Declare static data member
    template <typename ParticleType,short unsigned int dim>
    double SpatialCellObserver<ParticleType,dim>::cellSize=1.0;
    
    template <typename ParticleType,short unsigned int dim>
    std::map<Eigen::Matrix<int,dim,1>, SpatialCell<ParticleType,dim>* const,CompareVectorsByComponent<int,dim> > SpatialCellObserver<ParticleType,dim>::cellMap;

//                template <typename ParticleType,short unsigned int dim>
//                boost::ptr_map<Eigen::Matrix<int,dim,1>, SpatialCell<ParticleType,dim>* const,CompareVectorsByComponent<int,dim> > SpatialCellObserver<ParticleType,dim>::cellMap;

                
//                template <typename ParticleType,short unsigned int dim>
//                std::map<Eigen::Matrix<int,dim,1>,
//                /*    */ SpatialCell<ParticleType,dim>* const,
//                /*    */ CompareVectorsByComponent<int,dim>,
//                /*    */ Eigen::aligned_allocator<std::pair<const Eigen::Matrix<int,dim,1>, SpatialCell<ParticleType,dim>* const> > > SpatialCellObserver<ParticleType,dim>::cellMap;

                
}	// close namespace
#endif

                
                
                
                //        const SharedPtrType pCell;
                //
                //
                //        SpatialCellObserver(const VectorDimD& P) :
                //		/* init list */ pCell(getCellByPosition(P))
                //        //        /* init list */ mpiID(this->sID)
                //        {/*!\param[in] P the position of this SpatialCellObserver
                //          *
                //          * Creates and anchor to a SpatialCell at position P
                //          */
                //		}
                //
                //        /* neighborCellsBegin ***************************************/
                //        typename CellMapType::const_iterator neighborCellsBegin() const
                //        {/*!\returns a const iterator to the first neighbor SpatialCell
                //          */
                //            return pCell->neighborCellsBegin();
                //        }
                //
                //        /* neighborCellsEnd ***************************************/
                //        typename CellMapType::const_iterator neighborCellsEnd() const
                //        {/*!\returns a const iterator to the past-the-end neighbor SpatialCell
                //          */
                //            return pCell->neighborCellsEnd();
                //        }

