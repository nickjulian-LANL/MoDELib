/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GL2tga_H_
#define model_GL2tga_H_

#include <fstream>
#include <string>
#include <sstream>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#endif



namespace model {
	
	class GL2tga  {
		

        
    public:
        

        static bool saveTGA;
        
        static void saveAs(const std::string& filename)
        {/*! Saves to file
          */
            
            if (saveTGA)
            {
                int x=glutGet(GLUT_WINDOW_WIDTH);
                int y=glutGet(GLUT_WINDOW_HEIGHT);
                
                // get the image data
                long imageSize = x * y * 3;
                unsigned char* data = new unsigned char[imageSize];
                glReadPixels(0,0,x,y, GL_BGR,GL_UNSIGNED_BYTE,data);
                
                // split x and y sizes into bytes
                int xa= x % 256;
                int xb= (x-xa)/256;
                
                int ya= y % 256;
                int yb= (y-ya)/256;
                
                //assemble the header
                unsigned char header[18]={0,0,2,0,0,0,0,0,0,0,0,0,(char)xa,(char)xb,(char)ya,(char)yb,24,0};
                
                
                std::stringstream filenameStream;
                filenameStream << filename << ".tga";
                std::string filenameWithExtension=filenameStream.str();
                std::cout<<"Saving file"<<filenameWithExtension<<std::endl;
                
                
                // write header and data to file
                std::fstream file(filenameWithExtension.c_str(), std::ios::out | std::ios::binary);
                file.write (reinterpret_cast<const char*>(header), sizeof (char)*18);
                file.write (reinterpret_cast<const char*>(data), sizeof (char)*imageSize);
                file.close();
                
                delete[] data;
                data=NULL;
            }
            
        }
		
	};
	
    // Declare Static data
    bool GL2tga::saveTGA=true;
    
}
#endif
/*********************************************************************/
/*********************************************************************/





