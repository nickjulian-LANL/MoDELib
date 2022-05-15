/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigVtk_h_
#define model_DDconfigVtk_h_

#include <QGroupBox>
#include <QGridLayout>
#include <QCheckBox>
#include <QLabel>
#include <QTabWidget>
#include <QTextEdit>
#include <QLineEdit>
#include <QPushButton>

#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <DDconfigIO.h>
#include <DDauxIO.h>
//#include <DDconfigVtkBase.h>
#include <SimplicialMesh.h>
#include <NetworkNodeActor.h>
#include <NetworkLinkActor.h>
#include <NetworkLoopActor.h>
#include <InclusionActor.h>
//#include <DislocationLoopActor.h>

namespace model
{
    struct DDconfigVtk : public QWidget,
                         public DDconfigIO<3>,
                         public DDauxIO<3>
    {
        
        Q_OBJECT
        
        vtkGenericOpenGLRenderWindow* const renderWindow;
        NetworkNodeActor* nodes;
        NetworkLinkActor* segments;
        NetworkLoopActor* loops;
        InclusionActor* inclusions;
//        DislocationLoopActor loops;
        private slots:
        bool updateConfiguration();
        void nextConfiguration();
        void prevConfiguration();

        public:
        
        QGridLayout* mainLayout;
        QLineEdit* frameIDedit;
        QPushButton* plusFrameButton;
        QPushButton* minusFrameButton;
        QLineEdit* frameIncrementEdit;

        QTabWidget* tabWidget;


        
        /**********************************************************************/
        DDconfigVtk(const std::string& evlFolder,const std::string& auxFolder,vtkGenericOpenGLRenderWindow* const, vtkRenderer* const ren,const SimplicialMesh<3>& mesh);
        bool updateConfiguration(const size_t& frameID);

        void modify();
        const DDconfigIO<3>& configIO() const;
        DDconfigIO<3>& configIO();
        const DDauxIO<3>& auxIO() const;
        DDauxIO<3>& auxIO();
    };
    
}
#endif