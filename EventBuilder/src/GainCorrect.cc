#include "GainCorrect.h"
#include <math.h>
#include "UTIL/CellIDDecoder.h"
#include <string>
#include <algorithm>
#include <fstream>
using namespace std;

uint GainCorrect::getCellDif_id(int cell_id){
  return cell_id & 0xFF;
}
uint GainCorrect::getCellAsic_id(int cell_id){
  return (cell_id & 0xFF00)>>8;
}
uint GainCorrect::getCellChan_id(int cell_id){
  return (cell_id & 0x3F0000)>>16;
}
void GainCorrect::gainCorrectionInit(std::map<int, LayerID  > _mapping){
  mean_hit_dif = new TH1F("mean_hit_dif","mean_hit_dif",500,0.,500.);
  GainHist _gain_tmp;

  for(std::map<int,LayerID>::iterator itt = _mapping.begin();itt!=_mapping.end();itt++){
    int dif_id = itt->first;
    _gain_tmp.noise_on_dif = new TH1F(Form("noise_on_dif_%i",dif_id),"noise_on_dif",3072,0,3072);
    _gain_tmp.gain_chan    = new TH1F(Form("gain_chan_%i",dif_id),"gain_chan",3072,0,3072);
    _h_gain[dif_id] = _gain_tmp;
  }
}

void GainCorrect::gainCorrectionParser(LCCollection* col){
  try{
    for (int ihit(0); ihit < col->getNumberOfElements(); ++ihit) {
      RawCalorimeterHit *rawhit=dynamic_cast<RawCalorimeterHit*>(col->getElementAt(ihit));
      int Dif_id  =  getCellDif_id(rawhit->getCellID0());
      int Asic_id =  getCellAsic_id(rawhit->getCellID0());
      int Chan_id =  getCellChan_id(rawhit->getCellID0());

      for(std::map<int, GainHist >::iterator itt = _h_gain.begin();itt!=_h_gain.end();itt++){
	if( Dif_id == itt->first ){
	  itt->second.noise_on_dif->Fill(Chan_id + (Asic_id-1)*64);
	}
      }

    }
  }catch (std::exception ec){
    std::cout<<"No hits "<<std::endl;
  }
}

void GainCorrect::gainCorrector(const string gain_file="gain_correction.txt"){
  ofstream gainFlux(gain_file.c_str());
  std::map<int,double> mean_gain_dif;
  TH1F * h_mean_gain_full =new  TH1F("h_mean_gain_full","mean_gain_full",100,0,10);
  for(std::map<int, GainHist >::iterator itt = _h_gain.begin();itt!=_h_gain.end();itt++){
    mean_hit_dif->Clear();

    for(int i=0; i < itt->second.noise_on_dif->GetNbinsX();i++){
      mean_hit_dif->Fill(itt->second.noise_on_dif->GetBinContent(i));
    }
    double mean =  mean_hit_dif->GetMean();
    TH1F * h_mean_gain =new TH1F("h_mean_gain","mean_gain",100,0,10);
    //double sigma = mean_hit_dif->GetRMS();
    //TH1F *noise_on_dif_copy = new TH1F("noise_on_dif_copy","noise_on_dif_copy",3072,0,3072);
    for(int i=0; i < itt->second.noise_on_dif->GetNbinsX();i++){ // loop over the pad
      //itt->second.noise_on_dif->GetBinContent(i);
      if(itt->second.noise_on_dif->GetBinContent(i) > 0.0){
	//double gain = 1./( 1 + TMath::Exp((itt->second.noise_on_dif->GetBinContent(i)-mean)/mean));
	//noise_on_dif_copy->Fill(i,itt->second.noise_on_dif->GetBinContent(i) * gain);
	double gain = mean / (itt->second.noise_on_dif->GetBinContent(i));
	itt->second.gain_chan->Fill(i,gain*itt->second.noise_on_dif->GetBinContent(i));
	if(gain > 2. )gain = 1.99;

	h_mean_gain->Fill(gain);
	//gainFlux<<"oa.ChangeGain("
	//	  << itt->first <<","
	//	  << (int(i/64) + 1)<<","
	//	  << (i - 64*int(i/64)) <<","
	//	  << gain <<")"
	//	  <<endl;
      }
    }//loop over pad

    double mean_gain = h_mean_gain->GetMean();
    mean_gain_dif.insert(std::pair<int,double>(itt->first,mean_gain));
    h_mean_gain_full->Fill(mean_gain);
    //for(int i=0; i < itt->second.noise_on_dif->GetNbinsX();i++){
    //if(itt->second.noise_on_dif->GetBinContent(i) > 0.0){
    //gain_chan->Fill(i,(mean + sigma)/(1.0*itt->second.noise_on_dif->GetBinContent(i)));
    //double factor_rms = (itt->second.noise_on_dif->GetBinContent(i)- mean)/sigma;
    //double factor_rms = (itt->second.noise_on_dif->GetBinContent(i)-mean)/(1.0*mean);
    //double gain = ( mean + TMath::Exp( - factor_rms)*mean )/(mean + factor_rms * mean);
    //  double gain_new = 1./( 1 + TMath::Exp((itt->second.noise_on_dif->GetBinContent(i)-mean_new)/mean_new));
    //itt->second.gain_chan->Fill(i,gain_new*itt->second.noise_on_dif->GetBinContent(i));
    //	  gainFlux<<"oa.ChangeGain("
    //	  << itt->first <<","
    //	  << (int(i/64) + 1)<<","
    //	  << (i - 64*int(i/64)) <<","
    //	  << gain_new <<")"
    //	  <<endl;
    //}
    //}
    delete h_mean_gain;
  }// loop over the diff

  double mean_gain_full = h_mean_gain_full->GetMean();
  for(std::map<int, GainHist >::iterator itt = _h_gain.begin();itt!=_h_gain.end();itt++){
    mean_hit_dif->Clear();
    for(int i=0; i < itt->second.noise_on_dif->GetNbinsX();i++){
      mean_hit_dif->Fill(itt->second.noise_on_dif->GetBinContent(i));
    }

    std::cout << "write file gain_correction.txt ..." << std::endl;

    double mean =  mean_hit_dif->GetMean();
    for(int i=0; i < itt->second.noise_on_dif->GetNbinsX();i++){
      if(itt->second.noise_on_dif->GetBinContent(i) > 0.0){
	double gain = mean / (itt->second.noise_on_dif->GetBinContent(i));
	double gain_new =  (mean_gain_full/(mean_gain_dif.find(itt->first)->second))*gain;
	if(gain_new > 2. )gain_new = 1.99;
	gainFlux<<"oa.ChangeGain("
		<< itt->first <<","
		<< (int(i/64) + 1)<<","
		<< (i - 64*int(i/64)) <<","
		<< gain_new <<")"
		<<endl;
      }
    }

  }

  delete h_mean_gain_full;
}
