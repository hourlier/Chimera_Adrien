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

        hROI[0] = new TH2D("hROI_U","hROI_U",10000,0,10000,10000,0,10000);
        hROI[1] = new TH2D("hROI_V","hROI_V",10000,0,10000,10000,0,10000);
        hROI[2] = new TH2D("hROI_Y","hROI_Y",10000,0,10000,10000,0,10000);

        ReadListFile();
        return true;
    }

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

    bool FindProtonTracks::analyze(storage_manager* storage) {
        //return true;
        gStyle->SetOptStat(0);
        double margin = 0.25;

        //
        // This is an example program to demonstrate how to access association
        // information between data product A and B. For simplicity let us
        // assume A=track, B=hit. An association itself is a data product (so
        // you can think it's C, if you want to). So, in principle, what you
        // can do is to retrieve A, B, and C, then ask C about how elements of
        // A and B are correlated.
        //
        // But this means you have to know the producer name of A, B, and C
        // apriori, and often this is cumbersome to look up. It is also counter
        // intuitive because, when you want an association, you rther want
        // "association to tell you which data products are correlated", instead
        // of requiring you some knowledge beforehand.
        //
        // One of several ways to retrieve association solves this problem, and
        // I show that method below.
        //
        // This method is relevant when you know the followings:
        // 0) You know the type (track) and producer of A
        // 1) You know the type (hit) of B made by a producer of A
        //
        // Code below shows how you can directly retrieve the association information
        // in such situation. So this approach allows you to "retrieve associated 2D
        // hits with kalmanhit tracks (where association is made by kalmanhit)"
        //

        // Step 0: retrieve A

        auto track_v       = storage->get_data<larlite::event_track>("pandoraNu"); // A

        if(!FindCorrespondingTrackInList((int)storage->run_id(),(int)storage->subrun_id(),(int)storage->event_id()))return false;

        if(!track_v){
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate tracks " << std::endl;
            return false;
        }

        // Step 1: retrieve association info + B
        //         We use a function of storage manager that searches for B based
        //         on type and producer name of A. Upon successful finding, it
        //         provides you a pointer to B + returns association information.

        larlite::event_hit* hit_v = nullptr;

        auto const& ass_info = storage->find_one_ass(track_v->id(),         // provide A's product ID info
                                                     hit_v,                 // provide B's receiver (pointer)
                                                     track_v->id().second); // provide B's producer name;

        // Check if B is found (if not found it stays as nullptr)
        if(!hit_v) {
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate hits associated with " << track_v->id().second
            << " track producer!" << std::endl;
            return false;
        }

        // Step 2: loop over association information
        //         "ass_info" is a vector-of-vector-of-int, where the outer vector represents the index
        //         number of track_v (vector of track) element, and the inner vector corresponds to
        //         the indecies of associated elements in hit_v.

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

            hROI[0]->Reset();hROI[0]->Clear();
            hROI[1]->Reset();hROI[1]->Clear();
            hROI[2]->Reset();hROI[2]->Clear();

            hROI[0]->SetNameTitle(Form("hROI_%d_%d_%d_%zu_U",storage->run_id(),storage->subrun_id(),storage->event_id(),track_index),
                                  Form("hROI_%d_%d_%d_%zu_U;wire;time tick",storage->run_id(),storage->subrun_id(),storage->event_id(),track_index));
            hROI[1]->SetNameTitle(Form("hROI_%d_%d_%d_%zu_V",storage->run_id(),storage->subrun_id(),storage->event_id(),track_index),
                                  Form("hROI_%d_%d_%d_%zu_V;wire;time tick",storage->run_id(),storage->subrun_id(),storage->event_id(),track_index));
            hROI[2]->SetNameTitle(Form("hROI_%d_%d_%d_%zu_Y",storage->run_id(),storage->subrun_id(),storage->event_id(),track_index),
                                  Form("hROI_%d_%d_%d_%zu_Y;wire;time tick",storage->run_id(),storage->subrun_id(),storage->event_id(),track_index));
            
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
                hROI[iPlane]->GetXaxis()->SetRangeUser(WireMin[iPlane]-margin*(WireMax[iPlane]-WireMin[iPlane]),WireMax[iPlane]+margin*(WireMax[iPlane]-WireMin[iPlane]));
                hROI[iPlane]->GetYaxis()->SetRangeUser(TimeMin[iPlane]-margin*(TimeMax[iPlane]-TimeMin[iPlane]),TimeMax[iPlane]+margin*(TimeMax[iPlane]-TimeMin[iPlane]));
                cWireSignal->cd(iPlane+1);
                hROI[iPlane]->Draw("colz");
            }

            cWireSignal->Modified();
            cWireSignal->Update();
            cWireSignal->SaveAs(Form("plot/Tracks_%06d_%03d_%05d_%02zu.png",storage->run_id(),storage->subrun_id(),storage->event_id(),track_index));
        }

        return true;
    }

    bool FindProtonTracks::finalize() {

        return true;
    }
    
}
#endif
