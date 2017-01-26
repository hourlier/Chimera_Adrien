#ifndef LARLITE_FINDPROTONTRACKS_CXX
#define LARLITE_FINDPROTONTRACKS_CXX

#include "FindProtonTracks.h"
#include "DataFormat/wire.h"
#include "DataFormat/track.h" // Include track header
#include "DataFormat/hit.h"   // Include hit header
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"


namespace larlite {

    bool FindProtonTracks::initialize() {
        return true;
    }

    bool FindProtonTracks::analyze(storage_manager* storage) {

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
        std::cout << "New EVT" << std::endl;
        auto track_v = storage->get_data<larlite::event_track>("pandoraNu");
        if(!track_v){
            std::cout << "\033[93m" << "ERROR" << "\033[00m"
            << " could not locate tracks " << std::endl;
            return false;
        }
        std::cout << "FOUND " << track_v->size() << " tracks in event" << std::endl;

        // Step 1: retrieve association info + B
        //         We use a function of storage manager that searches for B based
        //         on type and producer name of A. Upon successful finding, it
        //         provides you a pointer to B + returns association information.
        larlite::event_hit* hit_v = nullptr;
        auto const& ass_info = storage->find_one_ass(track_v->id(), // provide A's product ID info
                                                     hit_v,         // provide B's receiver (pointer)
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

            auto const& one_track   = (*track_v)[track_index];
            auto const& hit_index_v = ass_info[track_index];
            
            std::cout <<"\033[93m" << "Track ID " << one_track.ID() << "\033[00m"
            << " associated with " << hit_index_v.size()
            << " hits (hit charge dump below)" << std::endl;
            
            for(auto const& hit_index : hit_index_v) {
                std::cout << (int)((*hit_v)[hit_index].Integral()) << " ";
            }
            std::cout << std::endl;
        }

        return true;
    }

    bool FindProtonTracks::finalize() {

        return true;
    }
    
}
#endif
