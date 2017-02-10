#ifndef LARLITE_FINDPROTONTRACKS_CXX
#define LARLITE_FINDPROTONTRACKS_CXX

#include "FindProtonTracks.h"
#include "DataFormat/wire.h"
#include "DataFormat/track.h" // Include track header
#include "DataFormat/hit.h"   // Include hit header
#include "DataFormat/partid.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include <cmath>
#include <algorithm>
#include <fstream>

#include "TStyle.h"

namespace larlite {

    bool FindProtonTracks::initialize() {
        cWireSignal = new TCanvas("cWireSignal","cWireSignal",900,300);
        cWireSignal->Divide(3,1);
        isProtonTrack = false;

        hROI[0] = new TH2D("hROI_U","U plane; wire;time tick",10000,0,10000,10000,0,10000);
        hROI[1] = new TH2D("hROI_V","V plane; wire;time tick",10000,0,10000,10000,0,10000);
        hROI[2] = new TH2D("hROI_Y","Y plane; wire;time tick",10000,0,10000,10000,0,10000);

        ReadListFile();
        return true;
    }

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

    bool FindProtonTracks::ReadListFile(){
        TrackListInfo.clear();
        std::ifstream ListFile("passedGBDT_extBNB_AnalysisTrees_cosmic_trained_only_on_mc_score_0.99.txt");
        std::vector<int> eventInfo;
        int Run, subrun, evt, trackID, wire_min_U, time_min_U, wire_max_U, time_max_U,  wire_min_V, time_min_V, wire_max_V, time_max_V,  wire_min_Y, time_min_Y, wire_max_Y, time_max_Y;
        double BDTscore;
        bool goOn = true;
        std::string line;
        getline(ListFile,line);
        while(goOn == true){
            ListFile >> Run >> subrun >> evt >> trackID >> wire_min_U >> time_min_U >> wire_max_U >> time_max_U >>  wire_min_V >> time_min_V >> wire_max_V >> time_max_V >> wire_min_Y >> time_min_Y >> wire_max_Y >> time_max_Y >> BDTscore;
            eventInfo.clear();
            eventInfo.push_back(Run);
            eventInfo.push_back(subrun);
            eventInfo.push_back(evt);
            eventInfo.push_back(trackID);
            /*eventInfo.push_back(wire_min_U);
            eventInfo.push_back(time_min_U);
            eventInfo.push_back(wire_max_U);
            eventInfo.push_back(time_max_U);
            eventInfo.push_back(wire_min_V);
            eventInfo.push_back(time_min_V);
            eventInfo.push_back(wire_max_V);
            eventInfo.push_back(time_max_V);
            eventInfo.push_back(wire_min_Y);
            eventInfo.push_back(time_min_Y);
            eventInfo.push_back(wire_max_Y);
            eventInfo.push_back(time_max_Y);*/
            TrackListInfo.push_back(eventInfo);
            //std::cout << eventInfo[0] << "\t" << eventInfo[1] << "\t" << eventInfo[2] << std::endl;
            if(ListFile.eof())goOn = false;
        }
        return true;
    }

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

    bool FindProtonTracks::FindCorrespondingTrackInList(int run, int subrun, int evt){
        TrackIDs.clear();

        for(size_t iEvtInfo = 0;iEvtInfo<TrackListInfo.size();iEvtInfo++){
            if(TrackListInfo[iEvtInfo][0] == run && TrackListInfo[iEvtInfo][1] == subrun && TrackListInfo[iEvtInfo][2] == evt){
                TrackIDs.push_back(TrackListInfo[iEvtInfo][3]);
            }
        }
        if(TrackIDs.size() == 0){std::cout << "ERROR, no correspondance found in list" << std::endl;return false;}
        else return true;
    }

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

    bool FindProtonTracks::analyze(storage_manager* storage) {
        gStyle->SetOptStat(0);
        double margin = 0.25;
        auto track_v       = storage->get_data<larlite::event_track>("pandoraNu"); // A
        if(!FindCorrespondingTrackInList((int)storage->run_id(),(int)storage->subrun_id(),(int)storage->event_id()))return false;

        if(!track_v){
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate tracks " << std::endl;
            return false;
        }

        larlite::event_hit* hit_v = nullptr;
        auto const& ass_info = storage->find_one_ass(track_v->id(),         // provide A's product ID info
                                                     hit_v,                 // provide B's receiver (pointer)
                                                     track_v->id().second); // provide B's producer name;

        if(!hit_v) {
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate hits associated with " << track_v->id().second
            << " track producer!" << std::endl;
            return false;
        }

        for(size_t track_index=0; track_index < ass_info.size(); ++track_index) {
            isProtonTrack = false;

            auto const& one_track   = (*track_v)[track_index];
            auto const& hit_index_v = ass_info[track_index];

            unsigned int WireMin[3] = {100000,100000,100000};
            unsigned int WireMax[3] = {0,0,0};
            double TimeMin[3] = {100000,100000,100000};
            double TimeMax[3] = {0,0,0};


            for(size_t iTrack = 0;iTrack<TrackIDs.size();iTrack++){
                if(one_track.ID() == TrackIDs[iTrack])isProtonTrack=true;
            }
            if(!isProtonTrack)continue;

            hROI[0]->Reset();
            hROI[1]->Reset();
            hROI[2]->Reset();

            int hitNum = 0;
            for(auto const& hit_index : hit_index_v) {
                hROI[(*hit_v)[hit_index].WireID().Plane]->SetBinContent( hROI[(*hit_v)[hit_index].WireID().Plane]->FindBin((*hit_v)[hit_index].WireID().Wire, (*hit_v)[hit_index].PeakTime()), (*hit_v)[hit_index].Integral() );

                WireMin[(*hit_v)[hit_index].WireID().Plane] = std::min(WireMin[(*hit_v)[hit_index].WireID().Plane], (*hit_v)[hit_index].WireID().Wire);
                WireMax[(*hit_v)[hit_index].WireID().Plane] = std::max(WireMax[(*hit_v)[hit_index].WireID().Plane], (*hit_v)[hit_index].WireID().Wire);

                TimeMin[(*hit_v)[hit_index].WireID().Plane] = std::min( (double)(TimeMin[(*hit_v)[hit_index].WireID().Plane]),  (double)((*hit_v)[hit_index].PeakTime()));
                TimeMax[(*hit_v)[hit_index].WireID().Plane] = std::max( (double)(TimeMax[(*hit_v)[hit_index].WireID().Plane]),  (double)((*hit_v)[hit_index].PeakTime()));

                hitNum++;
            }

            for(int iPlane = 0;iPlane < 3;iPlane++){
                //hROI[iPlane]->RebinY(8);
                hROI[iPlane]->GetXaxis()->SetRangeUser(WireMin[iPlane]-margin*(WireMax[iPlane]-WireMin[iPlane]),WireMax[iPlane]+margin*(WireMax[iPlane]-WireMin[iPlane]));
                hROI[iPlane]->GetYaxis()->SetRangeUser(TimeMin[iPlane]-margin*(TimeMax[iPlane]-TimeMin[iPlane]),TimeMax[iPlane]+margin*(TimeMax[iPlane]-TimeMin[iPlane]));
                cWireSignal->cd(iPlane+1);
                hROI[iPlane]->Draw("colz");
            }

            cWireSignal->Modified();
            cWireSignal->Update();
            cWireSignal->SaveAs(Form("plot/Tracks_%06d_%03d_%05d_%02d.png",storage->run_id(),storage->subrun_id(),storage->event_id(),one_track.ID()));
        }

        return true;
    }

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

    bool FindProtonTracks::finalize() {

        return true;
    }
    
}
#endif
