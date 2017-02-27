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

namespace larlite {

    /*****************************/

    bool FindMuonTracks::initialize() {
        ReadCSVFile();
        return true;
    }

    /*****************************/

    bool FindMuonTracks::analyze(storage_manager* storage){
        gStyle->SetOptStat(0);
        double margin = 0.125;
        auto track_v   = storage->get_data<larlite::event_track>("trackkalmanhit");
        auto gaushit_v = storage->get_data<larlite::event_hit>("gaushit");
        std::cout << storage->run_id() << "\t" << storage->subrun_id() << "\t" << storage->event_id() << std::endl;

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

        if(!gaushit_v){
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not find gaushits " << std::endl;
            return false;
        }


        return true;
    }

    /*****************************/

    bool FindMuonTracks::finalize(){
        return true;
    }

    /*****************************/

    void FindMuonTracks::ReadCSVFile(){
        std::ifstream  file("data/muon_track_list.csv");
        int Run,SubRun,Event,decayIdx,Ntracks,trackIndex,Ymin,Tmin,Ymax,Tmax;
        char coma;
        while(!file.eof()){

            file >> Run >> coma >> SubRun >> coma >> Event >> coma >> decayIdx >> coma >> Ntracks;
            std::cout << Run << "\t" << SubRun << "\t" << Event << "\t" << decayIdx << "\t" << Ntracks << std::endl;
            if(Ntracks == 0){std::cout << std::endl; continue;}
            else{
                for(int iTrack = 0;iTrack < Ntracks;iTrack++){
                    file >> coma >> trackIndex >> coma >> Ymin >> coma >> Tmin >> coma >> Ymax >> coma >> Tmax;
                }
            }
            for(int iTrack = 0 ; iTrack < Ntracks;iTrack++){
                std::cout << "\t" << Ymin << "=>" << Ymax << std::endl;
                std::cout << "\t" << Tmin << "=>" << Tmax << std::endl;
                std::cout << std::endl;
            }
        }
    }

    /*****************************/

}

#endif
