#ifndef LARLITE_FINDMUONTRACKS_CXX
#define LARLITE_FINDMUONTRACKS_CXX

#include "FindMuonTracks.h"
#include "DataFormat/wire.h"
#include "DataFormat/track.h" // Include track header
#include "DataFormat/hit.h"   // Include hit header
#include "DataFormat/partid.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

#include "TStyle.h"


//Work on the end points. so far I'll be just showing the end points of the last track in the decay, I may want both of them

namespace larlite {

    /*****************************/

    bool FindMuonTracks::initialize() {

        cWireSignal = new TCanvas("cWireSignal","cWireSignal",900,300);
        cWireSignal->Divide(3,1);

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
        T->Branch("Run",&_Run);
        T->Branch("SubRun",&_SubRun);
        T->Branch("Event",&_Event);
        T->Branch("TrackID",&TrackID);
        T->Branch("Phi",&Phi);
        T->Branch("Theta",&Theta);
        T->Branch("Length",&Length);
        //T->Branch("Npoints","size_t",&Npoints);
        T->Branch("Vertex","TVector3",&Vertex);
        T->Branch("End","TVector3",&End);

        ReadCSVFile();
        return true;
    }

    /*****************************/

    bool FindMuonTracks::analyze(storage_manager* storage){
        gStyle->SetOptStat(0);
        double margin = 0.125;
        auto track_v   = storage->get_data<larlite::event_track>("trackkalmanhit");
        auto gaushit_v = storage->get_data<larlite::event_hit>("gaushit");

        _Run    = storage->run_id();
        _SubRun = storage->subrun_id();
        _Event  = storage->event_id();

        //std::cout << _Run << "\t" << _SubRun << "\t" << _Event << std::endl;

        if(!track_v){
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate tracks " << std::endl;
            return false;
        }

        int Ndecays = -1;
        for(size_t itrack = 0;itrack<TrackListInfo.size();itrack++){
            if(_Run == TrackListInfo[itrack][0] && _SubRun == TrackListInfo[itrack][1] && _Event == TrackListInfo[itrack][2] && TrackListInfo[itrack][4] > 0){
                Ndecays = std::max(Ndecays,TrackListInfo[itrack][3]);
            }
        }
        Ndecays++;
        if(Ndecays == 0){return true;}

        std::cout << _Run << "\t" << _SubRun << "\t" << _Event << "\t" << Ndecays << std::endl;

        larlite::event_hit* hit_v = nullptr;
        auto const& ass_info = storage->find_one_ass(track_v->id(),
                                                     hit_v,
                                                     track_v->id().second);
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

        for(int iDecay = 0;iDecay < Ndecays;iDecay++){
            _decayIndex = iDecay;

            for(int iPlane = 0;iPlane<3;iPlane++){
                gGausHits[iPlane]->SetMarkerColor(3);
                while(gGausHits[ iPlane]->GetN()>0){gGausHits[ iPlane]->RemovePoint(0);}
                while(gTrackHits[iPlane]->GetN()>0){gTrackHits[iPlane]->RemovePoint(0);}
                //hROI[iPlane]->Reset();
            }

            int hitNum[3] = {0,0,0};
            unsigned int WireMin[3] = {100000,100000,100000};
            unsigned int WireMax[3] = {0,0,0};
            double TimeMin[3] = {100000,100000,100000};
            double TimeMax[3] = {0,0,0};

            double dmax[3] = {0,0,0};
            double dLoc[3] = {0,0,0};
            int iHit_indexMax[3] = {0,0,0};
            int jHit_indexMax[3] = {0,0,0};

            for(size_t track_index=0; track_index < ass_info.size(); ++track_index) {
                auto const& one_track   = (*track_v)[track_index];
                auto const& hit_index_v = ass_info[track_index];

                bool foundCorrespondingTrack = false;
                //std::cout << "\t track# " << one_track.ID();

                for(size_t itrack = 0;itrack<TrackListInfo.size();itrack++){
                    if(_Run == TrackListInfo[itrack][0] && _SubRun == TrackListInfo[itrack][1] && _Event == TrackListInfo[itrack][2] && _decayIndex == TrackListInfo[itrack][3] && (int)one_track.ID() == (int)TrackListInfo[itrack][5]){
                        //std::cout << "\t OK" << std::endl;
                        foundCorrespondingTrack = true;
                    }
                    else{continue;}
                }

                if(!foundCorrespondingTrack){/*std::cout << "\t -- " << std::endl;*/ continue;}

                /*for(int iPlane = 0;iPlane<3;iPlane++){
                    gGausHits[iPlane]->SetMarkerColor(3);
                    while(gGausHits[ iPlane]->GetN()>0){gGausHits[ iPlane]->RemovePoint(0);}
                    while(gTrackHits[iPlane]->GetN()>0){gTrackHits[iPlane]->RemovePoint(0);}
                    //hROI[iPlane]->Reset();
                }*/

                Phi = one_track.Phi();
                Theta = one_track.Theta();
                Length = one_track.Length(0);
                TrackID = one_track.ID();
                Npoints = one_track.NumberTrajectoryPoints();
                Vertex = one_track.Vertex();
                End = one_track.End();


                /*int hitNum[3] = {0,0,0};
                unsigned int WireMin[3] = {100000,100000,100000};
                unsigned int WireMax[3] = {0,0,0};
                double TimeMin[3] = {100000,100000,100000};
                double TimeMax[3] = {0,0,0};*/

                for(auto const& hit_index : hit_index_v) {
                    int iPlane = (*hit_v)[hit_index].WireID().Plane;
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

                /*double dmax[3] = {0,0,0};
                double dLoc[3] = {0,0,0};
                int iHit_indexMax[3] = {0,0,0};
                int jHit_indexMax[3] = {0,0,0};*/
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
            }

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

            //if(isOtherTrack)return true;

            DrawTrack();
        }


        return true;
    }

    /*****************************/

    bool FindMuonTracks::finalize(){
        T->Write();
        return true;
    }

    /*****************************/

    void FindMuonTracks::DrawTrack(){

        gGausHits[0]->SetTitle(Form("%05d_%03d_%05d_%02d_U",(int)_Run,(int)_SubRun,(int)_Event,(int)_decayIndex));
        gGausHits[1]->SetTitle(Form("%05d_%03d_%05d_%02d_V",(int)_Run,(int)_SubRun,(int)_Event,(int)_decayIndex));
        gGausHits[2]->SetTitle(Form("%05d_%03d_%05d_%02d_Y",(int)_Run,(int)_SubRun,(int)_Event,(int)_decayIndex));

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
        cWireSignal->SaveAs(Form("plot/Tracks_%05d_%03d_%05d_%02d.pdf",(int)_Run,(int)_SubRun,(int)_Event,(int)_decayIndex));
    }

    /*****************************/

    void FindMuonTracks::ReadCSVFile(){
        std::ifstream  file("data/muon_track_list.csv");
        int Run,SubRun,Event,decayIdx,Ntracks,trackIndex,Ymin,Tmin,Ymax,Tmax;
        char coma;
        std::vector<int> eventInfo;
        TrackListInfo.clear();
        bool goOn = true;

        while(goOn == true){
            file >> Run >> coma >> SubRun >> coma >> Event >> coma >> decayIdx >> coma >> Ntracks;
            if(Ntracks == 0){continue;}
            else{
                for(int iTrack = 0;iTrack < Ntracks;iTrack++){
                    file >> coma >> trackIndex >> coma >> Ymin >> coma >> Tmin >> coma >> Ymax >> coma >> Tmax;
                    if(file.eof()){goOn = false;return;}
                    eventInfo.clear();
                    eventInfo.push_back(Run);
                    eventInfo.push_back(SubRun);
                    eventInfo.push_back(Event);
                    eventInfo.push_back(decayIdx);
                    eventInfo.push_back(Ntracks);
                    eventInfo.push_back(trackIndex);
                    eventInfo.push_back(Ymin);
                    eventInfo.push_back(Ymax);
                    eventInfo.push_back(Tmin);
                    eventInfo.push_back(Tmax);
                    TrackListInfo.push_back(eventInfo);
                }
            }
            std::cout << Run << "\t" << SubRun << "\t" << Event << "\t" << decayIdx << "\t" << Ntracks << std::endl;
            for(int iTrack = 0 ; iTrack < Ntracks;iTrack++){
                std::cout << "\t" << trackIndex << std::endl;
                std::cout << "\t" << Ymin << "=>" << Ymax << std::endl;
                std::cout << "\t" << Tmin << "=>" << Tmax << std::endl;
                std::cout << std::endl;
            }
        }
    }
    
    /*****************************/
    
    
}

#endif
