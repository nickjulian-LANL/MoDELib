/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicLoopGenerator_H_
#define model_PeriodicLoopGenerator_H_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterial.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGeneratorBase.h>


namespace model
{

    class PeriodicLoopGenerator : public MicrostructureGeneratorBase
    {
        
        
    public:
        
        PeriodicLoopGenerator(const std::string& fileName);

//        void generate(MicrostructureGenerator& mg) override;
        void generateIndividual(MicrostructureGenerator& mg) override;
        void generateDensity(MicrostructureGenerator& mg) override;

    };

}
#endif