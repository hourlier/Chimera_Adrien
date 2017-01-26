/**
 * \file FindProtonTracks.h
 *
 * \ingroup ChimeraProtonTracks
 *
 * \brief Class def header for a class FindProtonTracks
 *
 * @author hourlier
 */

/** \addtogroup ChimeraProtonTracks

 @{*/

#ifndef LARLITE_FINDPROTONTRACKS_H
#define LARLITE_FINDPROTONTRACKS_H

#include "Analysis/ana_base.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TPad.h"
#include <vector>

namespace larlite {
    class FindProtonTracks : public ana_base{

        public:

        FindProtonTracks(){ _name="FindProtonTracks"; _fout=0;}
        virtual ~FindProtonTracks(){}
        virtual bool initialize();
        virtual bool analyze(storage_manager* storage);
        virtual bool finalize();

        protected:

        TCanvas *cWireSignal;
        TH1D    *hWireSignal;
        TH2D    *hROI[3];
        std::vector<TCanvas*> canvasContainer;
        
    };
}
#endif
/** @} */ // end of doxygen group
