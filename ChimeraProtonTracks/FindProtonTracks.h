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
    class FindProtonTracks : public ana_base{

        public:

        FindProtonTracks(){ _name="FindProtonTracks"; _fout=0;}
        virtual ~FindProtonTracks(){}
        virtual bool initialize();
        virtual bool analyze(storage_manager* storage);
        virtual bool finalize();

        bool ReadListFile();
        bool FindCorrespondingTrackInList(int run, int subrun, int evt);
        void DrawTrack();

        protected:

        TCanvas                             *cWireSignal;
        TTree                               *T;
        TGraph                              *XtremPoints[3][2];
        TGraph                              *gGausHits[3];
        TGraph                              *gTrackHits[3];
        std::vector< std::vector<float> >   TrackHitCoordinates[3];
        std::vector< std::vector<int> >     TrackListInfo;
        std::vector<int>                    TrackIDs;
        bool                                isProtonTrack;
        double                              Phi;
        double                              Theta;
        double                              Length;
        double                              Run, SubRun, Event, TrackID;
        size_t                              Npoints;
        size_t                              TrackPoint;
        TVector3                            Vertex;
        TVector3                            End;
        TVector3                            TrackPos;
        TVector3                            TrackDir;
        int                                 WireStart[3];
        int                                 WireEnd[3];
        int                                 TimeStart[3];
        int                                 TimeEnd[3];
        int                                 XtremLoc[3][2][2];
        int                                 iXtrem;
        int                                 jXtrem;
        bool                                isOtherTrack;
        std::string                         wrongTracks;
        std::string                         goodTracks;

        
    };
}
#endif
/** @} */ // end of doxygen group
