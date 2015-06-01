#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <EVENT/RawCalorimeterHit.h>

#include <iostream>
#include <cmath>
#include "Mapping.h"

#include <TH1.h>
#include <TH2.h>

using namespace std;

struct GainHist{
  TH1F *noise_on_dif;
  TH1F *gain_chan;
};

class GainCorrect{
 public:
  GainCorrect();
  ~GainCorrect();
  // methods
  uint    getCellDif_id(int cell_id);
  uint    getCellAsic_id(int cell_id);
  uint    getCellChan_id(int cell_id);

  void    gainCorrectionInit(std::map<int, LayerID  > _mapping);
  void    gainCorrectionParser(EVENT::LCCollection* col);
  void    gainCorrector(const string gain_file);

 protected:
  std::map<int, LayerID  > _mapping;
  std::map<int, GainHist > _h_gain ;

  TH1F *noise_dist;
  TH1F *gain_chan;
  TH1F *mean_hit_dif;
  TH1F *time_hit_dif;
};
