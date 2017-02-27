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
        wrongTracks = "";
        goodTracks  = "";

        for(int iPlane = 0;iPlane<3;iPlane++){
            gGausHits[iPlane] = new TGraph();
            gTrackHits[iPlane] = new TGraph();
            gTrackHits[iPlane]->SetMarkerStyle(7);
            gGausHits[iPlane]->SetMarkerStyle(20);
            gGausHits[iPlane]->SetMarkerColor(3);

            XtremPoints[iPlane][0] = new TGraph();
            XtremPoints[iPlane][1] = new TGraph();
            XtremPoints[iPlane][0]->SetMarkerColor(4);
            XtremPoints[iPlane][1]->SetMarkerColor(4);
            XtremPoints[iPlane][0]->SetMarkerStyle(4);
            XtremPoints[iPlane][1]->SetMarkerStyle(4);
            XtremPoints[iPlane][0]->SetMarkerSize(2);
            XtremPoints[iPlane][1]->SetMarkerSize(2);
        }

        //hROI[0] = new TH2D("hROI_U","U plane; wire;time tick",10000,0,10000,10000,0,10000);
        //hROI[1] = new TH2D("hROI_V","V plane; wire;time tick",10000,0,10000,10000,0,10000);
        //hROI[2] = new TH2D("hROI_Y","Y plane; wire;time tick",10000,0,10000,10000,0,10000);

        T = new TTree("T","T");
        T->Branch("Run",&Run);
        T->Branch("SubRun",&SubRun);
        T->Branch("Event",&Event);
        T->Branch("TrackID",&TrackID);
        T->Branch("Phi",&Phi);
        T->Branch("Theta",&Theta);
        T->Branch("Length",&Length);
        //T->Branch("Npoints","size_t",&Npoints);
        T->Branch("Vertex","TVector3",&Vertex);
        T->Branch("End","TVector3",&End);


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
            eventInfo.push_back(wire_min_U);
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
            eventInfo.push_back(time_max_Y);
            TrackListInfo.push_back(eventInfo);
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
        double margin = 0.125;
        auto track_v   = storage->get_data<larlite::event_track>("pandoraNu"); // A
        auto gaushit_v = storage->get_data<larlite::event_hit>("gaushit");

        if(!FindCorrespondingTrackInList((int)storage->run_id(),(int)storage->subrun_id(),(int)storage->event_id()))return false;

        if(!track_v){
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate tracks " << std::endl;
            return false;
        }

        larlite::event_hit* hit_v = nullptr; // B
        auto const& ass_info = storage->find_one_ass(track_v->id(),         // provide A's product ID info
                                                     hit_v,                 // provide B's receiver (pointer)
                                                     track_v->id().second); // provide B's producer name;

        if(!hit_v) {
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate hits associated with " << track_v->id().second
            << " track producer!" << std::endl;
            return false;
        }

        if(!gaushit_v){
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not find gaushits " << std::endl;
            return false;
        }

        for(size_t track_index=0; track_index < ass_info.size(); ++track_index) {
            isProtonTrack = false;

            unsigned int WireMin[3] = {100000,100000,100000};
            unsigned int WireMax[3] = {0,0,0};
            double TimeMin[3] = {100000,100000,100000};
            double TimeMax[3] = {0,0,0};

            auto const& one_track   = (*track_v)[track_index];
            auto const& hit_index_v = ass_info[track_index];

            for(size_t iTrack = 0;iTrack<TrackIDs.size();iTrack++){
                if(one_track.ID() == TrackIDs[iTrack])isProtonTrack=true;
            }
            if(!isProtonTrack)continue;

            for(int iPlane = 0;iPlane<3;iPlane++){
                TrackHitCoordinates[iPlane].clear();
                gGausHits[iPlane]->SetMarkerColor(3);
                while(gGausHits[ iPlane]->GetN()>0){gGausHits[ iPlane]->RemovePoint(0);}
                while(gTrackHits[iPlane]->GetN()>0){gTrackHits[iPlane]->RemovePoint(0);}
                //hROI[iPlane]->Reset();
            }

            Phi = one_track.Phi();
            Theta = one_track.Theta();
            Length = one_track.Length(0);
            Run = storage->run_id();
            SubRun = storage->subrun_id();
            Event = storage->event_id();
            TrackID = one_track.ID();
            Npoints = one_track.NumberTrajectoryPoints();
            Vertex = one_track.Vertex();
            End = one_track.End();

            TrackPoint = 1;
            TrackPos = one_track.LocationAtPoint(TrackPoint);
            TrackDir = TrackPos - Vertex;


            int hitNum[3] = {0,0,0};
            std::vector<float> hitCoordinates;
            for(auto const& hit_index : hit_index_v) {
                int iPlane = (*hit_v)[hit_index].WireID().Plane;
                hitCoordinates.clear();
                hitCoordinates.push_back((*hit_v)[hit_index].WireID().Wire);
                hitCoordinates.push_back((*hit_v)[hit_index].PeakTime());
                hitCoordinates.push_back((*hit_v)[hit_index].Integral());

                TrackHitCoordinates[iPlane].push_back(hitCoordinates);


                WireMin[iPlane] = std::min(WireMin[iPlane], (*hit_v)[hit_index].WireID().Wire);
                WireMax[iPlane] = std::max(WireMax[iPlane], (*hit_v)[hit_index].WireID().Wire);

                TimeMin[iPlane] = std::min( (double)(TimeMin[iPlane]),  (double)((*hit_v)[hit_index].PeakTime()));
                TimeMax[iPlane] = std::max( (double)(TimeMax[iPlane]),  (double)((*hit_v)[hit_index].PeakTime()));

                gTrackHits[iPlane]->SetPoint(hitNum[iPlane],(*hit_v)[hit_index].WireID().Wire,(*hit_v)[hit_index].PeakTime());
                //hROI[iPlane]->SetBinContent( hROI[iPlane]->FindBin((*hit_v)[hit_index].WireID().Wire, (*hit_v)[hit_index].PeakTime()), (*hit_v)[hit_index].Integral() );
                hitNum[iPlane]++;
            }

            /*for(int iPlane = 0;iPlane<3;iPlane++){
                hROI[iPlane]->GetXaxis()->SetRangeUser(WireMin[iPlane]-margin*(WireMax[iPlane]-WireMin[iPlane]),WireMax[iPlane]+margin*(WireMax[iPlane]-WireMin[iPlane]));
                hROI[iPlane]->GetYaxis()->SetRangeUser(TimeMin[iPlane]-margin*(TimeMax[iPlane]-TimeMin[iPlane]),TimeMax[iPlane]+margin*(TimeMax[iPlane]-TimeMin[iPlane]));
                for(int i = 0;i<hROI[iPlane]->GetNbinsX();i++){
                    if(i < WireMin[iPlane]-margin*(WireMax[iPlane]-WireMin[iPlane])-1 || i > WireMax[iPlane]+margin*(WireMax[iPlane]-WireMin[iPlane])+1) continue;
                    for(int j = 0;j<hROI[iPlane]->GetNbinsY();j++){
                        if(j < TimeMin[iPlane]-margin*(TimeMax[iPlane]-TimeMin[iPlane])-1 || j > TimeMax[iPlane]+margin*(TimeMax[iPlane]-TimeMin[iPlane])+1) continue;
                        if(hROI[iPlane]->GetBinContent(i,j) > 0) continue;
                        hROI[iPlane]->SetBinContent(i,j,0.01);
                    }
                }
            }*/

            int NgausHit[3]  = {0,0,0};
            for(size_t iGausHit = 0; iGausHit<gaushit_v->size(); iGausHit++){
                int iPlane = gaushit_v->at(iGausHit).WireID().Plane;
                if(   gaushit_v->at(iGausHit).PeakTime()    > TimeMin[iPlane]-margin*(TimeMax[iPlane]-TimeMin[iPlane])
                   && gaushit_v->at(iGausHit).PeakTime()    < TimeMax[iPlane]+margin*(TimeMax[iPlane]-TimeMin[iPlane])
                   && gaushit_v->at(iGausHit).WireID().Wire > WireMin[iPlane]-margin*(WireMax[iPlane]-WireMin[iPlane])
                   && gaushit_v->at(iGausHit).WireID().Wire < WireMax[iPlane]+margin*(WireMax[iPlane]-WireMin[iPlane])){
                    gGausHits[iPlane]->SetPoint(NgausHit[iPlane],gaushit_v->at(iGausHit).WireID().Wire,gaushit_v->at(iGausHit).PeakTime());
                    NgausHit[iPlane]++;
                }
            }
            isOtherTrack = false;
            for(int iPlane = 0;iPlane < 3;iPlane++){
                if(NgausHit[iPlane]-hitNum[iPlane] > 15 || (NgausHit[0]-hitNum[0])+(NgausHit[1]-hitNum[1])+(NgausHit[2]-hitNum[2]) > 20){
                    isOtherTrack = true;
                }
            }

            if(isOtherTrack)return true;

            double dmax[3] = {0,0,0};
            double dLoc[3] = {0,0,0};
            int iHit_indexMax[3] = {0,0,0};
            int jHit_indexMax[3] = {0,0,0};
            for(auto const& iHit_index : hit_index_v) {
                int iPlane = (*hit_v)[iHit_index].WireID().Plane;
                for(auto const& jHit_index : hit_index_v) {
                    int jPlane = (*hit_v)[jHit_index].WireID().Plane;
                    if(jPlane!= iPlane) continue;
                    if(iHit_index <= jHit_index)continue;
                    dLoc[iPlane] = sqrt( pow((*hit_v)[iHit_index].WireID().Wire-(*hit_v)[jHit_index].WireID().Wire,2)+pow((double)((*hit_v)[iHit_index].PeakTime())-(double)((*hit_v)[jHit_index].PeakTime()),2) );
                    if(dmax[iPlane]<=dLoc[iPlane]){
                        dmax[iPlane] = dLoc[iPlane];
                        XtremPoints[iPlane][0]->SetPoint(0,(*hit_v)[iHit_index].WireID().Wire,(double)((*hit_v)[iHit_index].PeakTime()));
                        XtremPoints[iPlane][1]->SetPoint(0,(*hit_v)[jHit_index].WireID().Wire,(double)((*hit_v)[jHit_index].PeakTime()));
                        iHit_indexMax[iPlane] = iHit_index;
                        jHit_indexMax[iPlane] = jHit_index;
                    }
                }
            }

            //try and figure out the end/begining of track
            if(Vertex.X() > End.X()){
                for(int iPlane = 0;iPlane < 3;iPlane++){
                    if((double)((*hit_v)[iHit_indexMax[iPlane]].PeakTime()) > (double)((*hit_v)[jHit_indexMax[iPlane]].PeakTime())){//then iHit_indexMax is associated with the vertex
                        XtremPoints[iPlane][0]->SetMarkerColor(8);
                        XtremPoints[iPlane][1]->SetMarkerColor(2);
                    }
                    else{
                        XtremPoints[iPlane][0]->SetMarkerColor(2);
                        XtremPoints[iPlane][1]->SetMarkerColor(8);
                    }
                }

            }
            else if(Vertex.X() < End.X()) {
                for(int iPlane = 0;iPlane < 3;iPlane++){
                    if((double)((*hit_v)[iHit_indexMax[iPlane]].PeakTime()) > (double)((*hit_v)[jHit_indexMax[iPlane]].PeakTime())){//then jHit_indexMax is associated with the vertex
                        XtremPoints[iPlane][0]->SetMarkerColor(2);
                        XtremPoints[iPlane][1]->SetMarkerColor(8);
                    }
                    else{
                        XtremPoints[iPlane][0]->SetMarkerColor(8);
                        XtremPoints[iPlane][1]->SetMarkerColor(2);
                    }
                }
            }

            else{std::cout << "WOW! exact same time! what are the odds?" << std::endl;}

            DrawTrack();
            T->Fill();
        }

        return true;
    }

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

    void FindProtonTracks::DrawTrack(){

        gGausHits[0]->SetTitle(Form("%05d_%03d_%05d_%02d_U",(int)Run,(int)SubRun,(int)Event,(int)TrackID));
        gGausHits[1]->SetTitle(Form("%05d_%03d_%05d_%02d_V",(int)Run,(int)SubRun,(int)Event,(int)TrackID));
        gGausHits[2]->SetTitle(Form("%05d_%03d_%05d_%02d_Y",(int)Run,(int)SubRun,(int)Event,(int)TrackID));

        double Tmin = 1e9;
        double Tmax = 0;
        for(int iPlane = 0;iPlane < 3;iPlane++){
            if( gGausHits[iPlane]->GetYaxis()->GetXmin() < Tmin ){Tmin = gGausHits[iPlane]->GetYaxis()->GetXmin();}
            if( gGausHits[iPlane]->GetYaxis()->GetXmax() > Tmax ){Tmax = gGausHits[iPlane]->GetYaxis()->GetXmax();}
        }

        for(int iPlane = 0;iPlane < 3;iPlane++){

            cWireSignal->cd(iPlane+1);
            if(isOtherTrack){gGausHits[iPlane]->SetMarkerColor(2);}
            if(gTrackHits[iPlane]->GetN() == 0 && gGausHits[iPlane]->GetN() == 0){gGausHits[iPlane]->SetPoint(0,1,1); gTrackHits[iPlane]->SetPoint(0,1,1);}
            gGausHits[iPlane]->GetYaxis()->SetRangeUser(Tmin-0.1*(Tmax-Tmin),Tmax+0.1*(Tmax-Tmin));

            //hROI[iPlane]->Draw("colz");

            gGausHits[iPlane]->Draw("AP");
            gTrackHits[iPlane]->Draw("same LP");

            XtremPoints[iPlane][0]->Draw("P same");
            XtremPoints[iPlane][1]->Draw("P same");
        }

        cWireSignal->Modified();
        cWireSignal->Update();
        cWireSignal->SaveAs(Form("plot/Tracks_%05d_%03d_%05d_%02d.pdf",(int)Run,(int)SubRun,(int)Event,(int)TrackID));
    }


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

    bool FindProtonTracks::finalize() {
        T->Write();
        return true;
    }
    
}
#endif
