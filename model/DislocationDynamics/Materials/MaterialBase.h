/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_MaterialBase_H_
#define model_MaterialBase_H_

#include <string>


namespace model
{
    

    
    struct MaterialBase
    {
        
        virtual std::string name() const =0;
        
        
    };
}
#endif
