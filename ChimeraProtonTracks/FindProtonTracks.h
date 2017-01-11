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

namespace larlite {
  /**
     \class FindProtonTracks
     User custom analysis class made by SHELL_USER_NAME
   */
  class FindProtonTracks : public ana_base{
  
  public:

    /// Default constructor
    FindProtonTracks(){ _name="FindProtonTracks"; _fout=0;}

    /// Default destructor
    virtual ~FindProtonTracks(){}

    /** IMPLEMENT in FindProtonTracks.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FindProtonTracks.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FindProtonTracks.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
