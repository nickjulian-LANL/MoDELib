/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MPIinit_h_
#define model_MPIinit_h_

#include <assert.h>
#include <mpi.h>
#include <model/Utilities/TerminalColors.h>

namespace model {
    
    
    
    class MPIinit
    {
    
        
        /*****************************************/
        static int getMPIrank()
        {
            //           assert(mpiInitialized && "MPI not initialized.");
            int mpiRank_temp;
            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank_temp);
            return mpiRank_temp;
        }
        
        /*****************************************/
        static int getMPInProcs()
        {
            int nProcs_temp;
            MPI_Comm_size(MPI_COMM_WORLD,&nProcs_temp);
            return nProcs_temp;
        }
        
        
        bool initMPI(int argc, char* argv[])
        {
        
            int temp(0);
            MPI_Initialized(&temp);
            if (!temp)
            {
                MPI_Init(&argc,&argv);
            }
            MPI_Initialized(&temp);
            assert(temp && "MPI not initialized");
            return temp;
        }
    
    public:
    
        const bool mpiInitialized;
        const int mpiRank;
        const int mpiProcs;

        
        MPIinit(int argc, char* argv[]) :
        /* init list */ mpiInitialized(initMPI(argc,argv)),
        /* init list */ mpiRank(getMPIrank()),
        /* init list */ mpiProcs(getMPInProcs())
        {
            std::cout<<greenBoldColor<<"MPI process "<<mpiRank<<" of "<<mpiProcs<<defaultColor<<std::endl;
        }
        
        ~MPIinit()
        {
            MPI_Finalize();
        }

        
    };
        
} // end namespace
#endif
