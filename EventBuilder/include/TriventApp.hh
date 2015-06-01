#ifndef _TriventProc_hh_
#define _TriventProc_hh_

#define  HISTOGRAM_PARSER true

#include <string>
#include <iostream>
#include <fstream>
//#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/RawCalorimeterHit.h>
#include <EVENT/LCParameters.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/RawCalorimeterHitImpl.h>
//#include <TTree.h>
//#include <TFile.h>
//#include <TH1.h>
//#include <TH1F.h>
#include "IO/LCWriter.h"
#include <map>
#include <algorithm>
#include "Mapping.h"
#include <bits/stream_iterator.h>


class TriventApp  //: public marlin::Processor
{
public:

//  Processor*  newProcessor() { return new TriventApp ; }
  TriventApp();

  ~TriventApp() {};

  void init(const EVENT::LCParameters *pPar);

  void    processEvent( EVENT::LCEvent * evtP );
  void    ClearVector();
//  void    fillTree();
  void    processRunHeader( EVENT::LCRunHeader * runH);// added by me
  void    XMLReader(std::string xmlfile);
  void    readDifGeomFile(std::string geomfile);
  void    printDifGeom();

  uint    getCellDif_id(int cell_id);
  uint    getCellAsic_id(int cell_id);
  uint    getCellChan_id(int cell_id);

  void    getMaxTime();
  std::vector<int> getTimeSpectrum();
  uint*   getPadIndex(uint dif_id, uint asic_id, uint chan_id);
  void    eventBuilder(EVENT::LCCollection* col_event,int time_peak, int prev_time_peak);
  bool    peakOrNot(std::vector<int> time_spectrum, int itime ,int threshold);
  void    end();

protected:
//  TH1F *noise_dist; // Replace with Vector / LCEvent
//  TH1F *gain_chan;
//  TH1F *mean_hit_dif;
//  TH1F *time_hit_dif;


//  // xml test
  std::map<std::string,std::string> m_parameters;

  std::vector<EVENT::RawCalorimeterHit*> _trigger_raw_hit;
  std::vector<EVENT::LCEvent*> evts; // Vector of evts reconstructed in the Trigger

  bool GAIN_CORRECTION_MODE;
  std::string _outFileName;
  std::string _rootFileName;
  std::string _noiseFileName;
  std::string _treeName;
  std::string _logrootName;
  std::string _colName;
  std::string _fileName;
  std::string _mappingfile;
  std::string _geomXML;
  std::ostream *_output;
  std::vector<std::string> _hcalCollections;
  int _overwrite;

//  TFile *_file;
//  TTree *_tree;
//  TClonesArray *_hitArray;
//  TClonesArray *_rootArray;// added by me
  unsigned int _triggerNr;
  unsigned int _runNr;
  unsigned int _eventType;
  int _evt_id;

  std::vector<int> x;
  std::vector<int> y;
  std::vector<int> z;
  std::vector<int> fThr;
  std::vector<int> fDifId;
  std::vector<int> fNevt;
  std::vector<int> fEvtReconstructed;
  std::vector<int> fPrevTime;
  std::vector<int> fTime;
  std::vector<int> fTimePeak;
  std::vector<int> fPrevTimePeak;
  std::vector<int> fDeltaTimePeak;
  std::vector<int> fTriggerNr;
  std::vector<int> fNhit;
  std::vector<int> fCerenkovTag1;
  std::vector<int> fCerenkovTag2;
  std::vector<int> fCerenkovCount1;
  std::vector<int> fCerenkovCount2;
  std::vector<int> fCerenkovCountDecember;
  std::vector<int> fPion;
  std::vector<int> fElectron;
  std::vector<int> fMuon;


  unsigned int _eventNr;
  int _Nhit;
  int _elec_noise_cut;

  std::map<int, LayerID  > _mapping;
  std::map<int, double  > _chamber_pos;//camber , pos

  int _noiseCut;
  int _timeWin;
  int _LayerCut;
  int _time2prev_event_cut;
  double _layer_gap;
  float _beamEnergy;
  int trig_count;
  int _maxtime;
  int evtnum;
  int previousEvtNbr;
  int _rejectedNum;
  int _selectedNum;
  uint _index[3];
  uintVec zcut;
  IO::LCWriter* _lcWriter;
  int _bcid1;
  int _bcid2;
  bool _cerenkov1;
  bool _cerenkov2;
  int _cerenkCount1;
  int _cerenkCount2;
  int _cerenkCountTot1;
  int _cerenkCountTot2;

  int _cerenkWindow;
  int _cerenkLength;

  std::string normal  ;
  std::string red     ;
  std::string green   ;
  std::string yellow  ;
  std::string blue    ;
  std::string magenta ;
  std::string white   ;

};



#endif
