#ifndef __TA2Moeller_h__
#define __TA2Moeller_h__

#include "TA2Detector.h"            // Acqu-Root histogrammer
#include <vector>

class TA2Moeller : public TA2Detector
{
private:
  UShort_t fNLeftChannels;
  UShort_t fNPairsPerCh;
  UInt_t fNBins;
  UShort_t fCounter;
  std::vector<UInt_t> fVupromOffset;
  std::vector<UInt_t*> fTDCs;          // pointers to NBins long histogram arrays
  std::vector<const char*> fTDCDisplayLines; // pointers to histogram Display line
  Bool_t fReadoutStarted;
  
  UInt_t fNvaluesRead;
public:
  TA2Moeller( const char*, TA2System* );  // pass ptr to analyser
  virtual ~TA2Moeller();                  // destructor
  virtual void PostInit();                    // some setup after parms read in
  virtual void LoadVariable();              // display setup
  virtual void Decode();                // event by event analysis
  virtual void Cleanup() {}                     // reset at end of event
  virtual void SetConfig( Char_t*, Int_t );   // setup decode in implement
  virtual void SaveDecoded( ) {}
  virtual void ReadDecoded( ) {}
 
  // Root needs this line for incorporation in dictionary
  ClassDef(TA2Moeller,1)
};

#endif
