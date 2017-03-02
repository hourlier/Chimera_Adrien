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
        void DrawTrack();

        protected:

        int _Run, _SubRun, _Event, _decayIndex;
        std::vector< std::vector<int> > TrackListInfo;
        TCanvas                             *cWireSignal;
        TTree                               *T;
        TGraph                         *XtremPoints[3][2];
        TGraph                         *gGausHits[3];
        TGraph                         *gTrackHits[3];
        //TH2D                           *hROI[3];

        double                              Phi;
        double                              Theta;
        double                              Length;
        double                              TrackID;
        size_t                              Npoints;
        size_t                              TrackPoint;
        TVector3                            Vertex;
        TVector3                            End;
        bool                                isOtherTrack;

    };
}

#endif
