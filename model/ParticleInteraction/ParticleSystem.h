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

#ifndef _MODEL_ParticleSystem_h_
#define _MODEL_ParticleSystem_h_

#ifdef _MODEL_MPI_
#include <model/ParticleInteraction/ParticleSystemParallel.h>
#else
#include <model/ParticleInteraction/ParticleSystemSerial.h>
#endif

#include <model/ParticleInteraction/SystemProperties.h>

namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Class template that stores _ParticleType objects organized
     * in SpatialCell(s) and computes nearest-neighbor and far-field interactions
     * among them. 
     * 
     * ParticleSystem can be compiled in serial and parallel version, depending
     * whether _MODEL_MPI_ is defined or not.
     */
    template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    struct ParticleSystem :
#ifdef _MODEL_MPI_
    /* inheritance */  public ParticleSystemParallel<_ParticleType,UserSystemProperties>
#else
    /* inheritance */  public ParticleSystemSerial<_ParticleType,UserSystemProperties>
#endif
    {
        
        /*****************************************/
        ParticleSystem()
        {/*!
          */
        }
        
#ifdef _MODEL_MPI_
        /*****************************************/
        ParticleSystem(int argc, char* argv[]) :
        /* init list */  ParticleSystemParallel<_ParticleType,UserSystemProperties>(argc,argv)
        {/*!
          */
        }
#else
        /*****************************************/
        ParticleSystem(int argc, char* argv[])
        {/*!
          */
            argc+=0; // avoid unused warning
            argv=argv; // avoid unused warning

        }
#endif
        
    };
    
} // end namespace
#endif

