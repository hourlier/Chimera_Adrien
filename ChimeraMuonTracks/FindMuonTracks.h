/**
 * \file FindMuonTracks.h
 *
 * \ingroup ChimeraMuonTracks
 *
 * \brief Class def header for a class FindMuonTracks
 *
 * @author hourlier
 */

/** \addtogroup ChimeraMuonTracks

 @{*/

#ifndef LARLITE_FINDMUONTRACKS_H
#define LARLITE_FINDMUONTRACKS_H

#include "Analysis/ana_base.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TPad.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>

#include "DataFormat/track.h"

namespace larlite {
    class FindMuonTracks : public ana_base{
        public :

        FindMuonTracks(){ _name="FindMuonTracks"; _fout=0;}
        virtual ~FindMuonTracks(){}
        virtual bool initialize();
        virtual bool analyze(storage_manager* storage);
        virtual bool finalize();
        void ReadCSVFile();

        protected:

        double Run;

    };
}

#endif
