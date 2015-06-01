#include "TriventApp.hh"

#include <EVENT/LCCollection.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

#include <limits.h>
#include <Rtypes.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <fstream>

#include "Mapping.h"
#include "tinyxml.h"
#include "streamlog/streamlog.h"
//#include <lcio.h>

//=========================================================
TriventApp::TriventApp():
    _output(0),
    _runNr(0),
    _eventType(0),
    _evt_id(0)
{

    streamlog_out( MESSAGE )<< "Trivent ... begin " << endl;
    _rejectedNum = 0;
    _selectedNum = 0;
}

void TriventApp::XMLReader(std::string xmlfile){
    TiXmlDocument doc(xmlfile.c_str());
    bool load_key = doc.LoadFile();
    if(load_key){
        streamlog_out( MESSAGE ) << green << "File : " << xmlfile.c_str() << normal <<std::endl;
        // tout ici
        TiXmlHandle hDoc(&doc);
        TiXmlElement* pElem;
        TiXmlHandle hRoot(0);
        // name block
        {
            pElem=hDoc.FirstChildElement().Element();
            // should always have a valid root but handle gracefully if it does
            if (!pElem) streamlog_out( WARNING ) << red << "error elem" << normal << std::endl;
            streamlog_out( MESSAGE ) << green << pElem->Value() << normal << std::endl;

            // save this for later
            hRoot=TiXmlHandle(pElem);
        }
        // parameters block
        {
            m_parameters.clear();
            pElem=hRoot.FirstChild("parameter").Element();
            std::string key = pElem->Attribute("name");
            streamlog_out( MESSAGE ) << green << key.c_str() << normal << std::endl;
            streamlog_out( DEBUG1 ) << green
                                    <<"parameter : "
                                   << pElem->Attribute("name")
                                   << normal
                                   << std::endl;

            std::vector<std::string> lines;
            {
                std::string value = pElem->GetText() ;
                std::vector<std::string> lines;
                istringstream iss(value);
                copy(istream_iterator<string>(iss),
                     istream_iterator<string>(),
                     back_inserter<vector<string> >(lines));
                for(unsigned int iline = 0; iline < lines.size(); iline++){
                    std::string line = lines.at(iline);
                    streamlog_out( MESSAGE ) << red << line << normal << std::endl;

                    stringstream ss( line.c_str() );
                    vector<string> result;

                    LayerID mapp;
                    int Dif_id;
                    while( ss.good() )
                    {
                        string substr;
                        getline( ss, substr, ',' );
                        result.push_back( substr );
                    }
                    istringstream ( result.at(0) ) >> Dif_id;
                    istringstream ( result.at(1) ) >> mapp.K;
                    istringstream ( result.at(2) ) >> mapp.DifX;
                    istringstream ( result.at(3) ) >> mapp.DifY;
                    istringstream ( result.at(4) ) >> mapp.IncX;
                    istringstream ( result.at(5) ) >> mapp.IncY;
                    _mapping[Dif_id] = mapp;
                }
            }
            pElem = pElem->NextSiblingElement();
            // ChamberGeom  Node.
            {
                streamlog_out( DEBUG1 ) << green
                                        <<"parameter : "
                                       << pElem->Attribute("name")
                                       << normal
                                       << std::endl;
                std::vector<std::string> lines;
                {
                    std::string value = pElem->GetText() ;
                    std::vector<std::string> lines;
                    istringstream iss(value);
                    copy(istream_iterator<string>(iss),
                         istream_iterator<string>(),
                         back_inserter<vector<string> >(lines));
                    for(unsigned int iline = 0; iline < lines.size(); iline++){
                        std::string line = lines.at(iline);
                        streamlog_out( MESSAGE ) << red << line << normal << std::endl;

                        stringstream ss( line.c_str() );
                        vector<string> result;

                        double position;
                        int Dif_id;
                        while( ss.good() )
                        {
                            string substr;
                            getline( ss, substr, ',' );
                            result.push_back( substr );
                        }
                        istringstream ( result.at(0) ) >> Dif_id;
                        istringstream ( result.at(3) ) >> position;

                        _chamber_pos[Dif_id] = position;
                    }
                }
            }
        }
    }else{
        streamlog_out( WARNING ) << red << "Failed to load file : " << xmlfile.c_str() << normal <<std::endl;
    }
}

void TriventApp::readDifGeomFile(std::string geomfile){

    cout << "read the mapping file .."<< endl;

    LayerID contenu;
    ifstream file(geomfile.c_str(), ios::in);
    if(file){
        while(!file.eof()){
            int Dif_id;
            char co;
            file >> Dif_id >> co
                    >> contenu.K >> co
                    >> contenu.DifX >> co
                    >> contenu.DifY >> co
                    >> contenu.IncX >> co
                    >> contenu.IncY ;
            _mapping [Dif_id] = contenu;
        }
        file.close();
    }
    else
        cerr << "ERROR ... maping file not correct !" << endl;
}

void TriventApp::printDifGeom(){

    for(std::map<int,LayerID>::iterator itt = _mapping.begin();itt!=_mapping.end();itt++)     {
        streamlog_out( MESSAGE ) << itt->first << "\t" << itt->second.K
                                 <<"\t"<<itt->second.DifX
                                <<"\t"<<itt->second.DifY
                               <<"\t"<<itt->second.IncX
                              <<"\t"<<itt->second.IncY
                             << std::endl;
    }
}

uint TriventApp::getCellDif_id(int cell_id){
    return cell_id & 0xFF; // Applique le masque 1111 1111
}
uint TriventApp::getCellAsic_id(int cell_id){
    return (cell_id & 0xFF00)>>8; //  Applique le masque 1111 1111 0000 0000 puis tronque les 8 derniers bits
}
uint TriventApp::getCellChan_id(int cell_id){
    return (cell_id & 0x3F0000)>>16; //  Applique le masque 0011 1111 0000 0000 0000 0000 puis tronque les 16 derniers bits
}

uint* TriventApp::getPadIndex(uint dif_id, uint asic_id, uint chan_id){
    _index[0]=_index[1]=_index[2]=0;
    double DifY = -1.,DifZ = -1.;
    DifZ = _mapping.find(dif_id)->second.K;
    DifY = _mapping.find(dif_id)->second.DifY;
    _index[0] = (1+MapILargeHR2[chan_id]+AsicShiftI[asic_id]);
    _index[1] = (32-(MapJLargeHR2[chan_id]+AsicShiftJ[asic_id]))+int(DifY);
    _index[2] = abs(int(DifZ));
    streamlog_out( DEBUG0 ) << " Dif_id == " << dif_id
                            << " Asic_id ==" << asic_id
                            << " Chan_id ==" << chan_id
                            << " I == " << _index[0]
                            << " J == " << _index[1]
                            << " K == " << _index[2]
                            << std::endl;
    return _index;
}

//===============================================
void TriventApp::getMaxTime()
{
    _maxtime = 0;
    try{
        for(std::vector<EVENT::RawCalorimeterHit*>::iterator raw_hit=_trigger_raw_hit.begin();raw_hit!= _trigger_raw_hit.end();raw_hit++){
            int time =  int((*raw_hit)->getTimeStamp());
            if(time >= 0) _maxtime = max(_maxtime, time);
        }
    }catch (std::exception ec){
        streamlog_out( WARNING )<<"No hits "<<std::endl;
    }
    streamlog_out( DEBUG1 ) << " maxtime before == " << _maxtime << std::endl;
}

//===============================================
std::vector<int> TriventApp::getTimeSpectrum() //__attribute__((optimize(0)))
{
    std::vector<int> time_spectrum(_maxtime+1);
    try{
        for(std::vector<EVENT::RawCalorimeterHit*>::iterator raw_hit=_trigger_raw_hit.begin();raw_hit!= _trigger_raw_hit.end();raw_hit++){
            int time =  int((*raw_hit)->getTimeStamp());
            if(time >= 0) time_spectrum[time]++;
        }
    }catch (std::exception ec){
        streamlog_out( WARNING )<<"No hits "<<std::endl;
    }
    return time_spectrum;
}

//===============================================
bool TriventApp::peakOrNot(std::vector<int> time_spectrum ,int itime ,int threshold){

#if HISTOGRAM_PARSER
    //    noise_dist->Fill(time_spectrum[itime]);
#endif

    if(time_spectrum[itime] >= threshold
            && time_spectrum[itime] >  time_spectrum[itime+1]
            && time_spectrum[itime] >= time_spectrum[itime+1]){
        return true;
    }else{
        return false;
    }
}

//===============================================
int IJKToKey(const int i,const int j,const int k){return 100*100*k+100*j+i;}


//===============================================
//===============================================
void TriventApp::eventBuilder(EVENT::LCCollection* col_event,int time_peak, int prev_time_peak){
    zcut.clear();

    col_event->setFlag(col_event->getFlag()|( 1 << EVENT::LCIO::RCHBIT_LONG));
    col_event->setFlag(col_event->getFlag()|( 1 << EVENT::LCIO::RCHBIT_TIME));
    UTIL::CellIDEncoder<IMPL::CalorimeterHitImpl> cd( "M:3,S-1:3,I:9,J:9,K-1:6" ,col_event) ;
    try{
        std::vector<int> hitKeys;
        int prevTime = 0;
        for(std::vector<EVENT::RawCalorimeterHit*>::iterator rawhit=_trigger_raw_hit.begin();rawhit!= _trigger_raw_hit.end();rawhit++){
            float pos[3];
            int time = (*rawhit)->getTimeStamp();

            int Dif_id  =  getCellDif_id ((*rawhit)->getCellID0());
            if (/*fabs(time-time_peak) <=20 && // Cerenkov signal arrives between 7 to 12 bins before event signal = Decembre
                                    fabs(time-time_peak) >=0  &&*/ // = 7-8 ou 17-18 en avril
                    Dif_id==3)
            {

                std::cout << (*rawhit)->getAmplitude()<< std::endl; // 3rd Threshold

                _cerenkCount1++;
                _cerenkCount2++;
                _cerenkCountTot1++;
                _cerenkCountTot2++;
                // if (_cerenkCount==4) ( Cerenkov de decembre dure 5clocks, Avril=1clock -> ajouter une condition sur le run# pour prendre en compte ces deux cas
                //   {
                _cerenkov1 = 1;
                _cerenkov2 = 1;
                streamlog_out ( MESSAGE ) << green << "Event #: " << evtnum << "\tTrigger #: " << trig_count <<  \
                                             "\tDelta t:" << time-time_peak << "\tAt Time:" << time << "\tMaTime: " << _maxtime << "\t MaxTime - CeremkovTime: " << _maxtime - time << normal << std::endl;
                // }

                x.push_back( -1 );
                y.push_back( -1 );
                z.push_back( -1 );
                fThr.push_back( -1 );
                fDifId.push_back( Dif_id );
                fTime.push_back( time );
                fPrevTime.push_back( prevTime );
                fTimePeak.push_back( time_peak );
                fPrevTimePeak.push_back( prev_time_peak );
                fDeltaTimePeak.push_back( time-time_peak );
                fNevt.push_back( evtnum );

                if(evtnum==previousEvtNbr) fEvtReconstructed.push_back( 0 );     // Only one entry per event
                else                       fEvtReconstructed.push_back( evtnum );

                fTriggerNr.push_back( trig_count );
                fNhit.push_back( -1 );
                fCerenkovTag1.push_back( _cerenkov1 );
                fCerenkovTag2.push_back( _cerenkov2 );

                if (_cerenkCount1>1) // Because in december CTag is 4-5TimeSlot long, only One Cerenkov was used
                {                    // With this, CTag is not counted multiple time
                    fCerenkovCountDecember.pop_back();
                    fCerenkovCountDecember.push_back( 0 );
                }

                fCerenkovCount1.push_back( _cerenkCount1 );
                fCerenkovCount1.push_back( _cerenkCount2 );
                fCerenkovCountDecember.push_back( _cerenkCount1 );
                fCerenkovCountDecember.push_back( _cerenkCount2 );


                // In the caloHit collection
                // Register the Cerenkov tag as a 3rd Threshold hit in the top left corner of the 1st layer -> Easily spotable in the visualisation
                IMPL::CalorimeterHitImpl* caloHit = new IMPL::CalorimeterHitImpl();

                caloHit->setTime(float((*rawhit)->getTimeStamp()));
                caloHit->setEnergy(float((*rawhit)->getAmplitude()&3)+2); // 3rd Threshold

                // set the cell id
                cd["I"] = 96 ;
                cd["J"] = 96 ;
                cd["K-1"] = 0 + (_cerenkCount1 - 1); // Add a hit on the following layer depending on the number of Cerenkov hit (if 4Tag -> 1 hit in the first 4 layers)
                cd["M"] = 0 ;
                cd["S-1"] = 3 ;

                int aHitKey=IJKToKey(cd["I"],cd["J"],cd["K-1"]);
                pos[0] = cd["I"]*10.*1.04125;
                pos[1] = cd["J"]*10.*1.04125;
                pos[2] = _chamber_pos.find(cd["K-1"])->second*10;
                streamlog_out( DEBUG0 ) << " I == " << cd["I"]
                                        << " J == " << cd["J"]
                                        << " K == " << cd["K-1"]
                                        << std::endl;
                cd.setCellID( caloHit ) ;
                caloHit->setPosition(pos);
                col_event->addElement(caloHit);
                hitKeys.push_back(aHitKey);

                continue;
            }
            else
            {
                if(fabs(time-time_peak) <=_timeWin &&
                        (time > prev_time_peak + _timeWin ))
                {
                    _cerenkov1 = 0;
                    _cerenkov2 = 0;
                    int Asic_id =  getCellAsic_id((*rawhit)->getCellID0());
                    int Chan_id =  getCellChan_id((*rawhit)->getCellID0());

                    uint I = getPadIndex(Dif_id, Asic_id, Chan_id)[0];
                    uint J = getPadIndex(Dif_id, Asic_id, Chan_id)[1];
                    uint K = getPadIndex(Dif_id, Asic_id, Chan_id)[2];

                    int aHitKey=IJKToKey(I,J,K);
                    pos[0] = I*10.*1.04125;
                    pos[1] = J*10.*1.04125;
                    pos[2] = _chamber_pos.find(K)->second*10;

                    if(K<=0||K>64)
                    {
                        streamlog_out( DEBUG0 ) << " Dif_id  == " << Dif_id
                                                << " Asic_id == " << Asic_id
                                                << " Chan_id == " << Chan_id
                                                << " I == " << I
                                                << " J == " << J
                                                << " K == " << K
                                                << std::endl;
                        continue;
                    }

                    x.push_back( pos[0] );
                    y.push_back( pos[1] );
                    z.push_back( pos[2] );
                    fThr.push_back( int((*rawhit)->getAmplitude()) );
                    fDifId.push_back( Dif_id );
                    fTime.push_back( time );
                    fPrevTime.push_back( prevTime );
                    fTimePeak.push_back( time_peak );
                    fPrevTimePeak.push_back( prev_time_peak );
                    fDeltaTimePeak.push_back( time-time_peak );
                    fNevt.push_back( evtnum );
                    streamlog_out ( DEBUG0 ) << green << "Evt: " << evtnum << "\t Previous: " << previousEvtNbr << normal << std::endl;

                    if(evtnum==previousEvtNbr) fEvtReconstructed.push_back( 0 ); // To have only 1 entry per evt reconstructed
                    else                       fEvtReconstructed.push_back( evtnum );

                    fTriggerNr.push_back( trig_count );
                    fNhit.push_back( col_event->getNumberOfElements() );
                    fCerenkovTag1.push_back( _cerenkov1 );
                    fCerenkovTag2.push_back( _cerenkov2 );
                    fCerenkovCount1.push_back ( 0 );
                    fCerenkovCount2.push_back ( 0 );
                    fCerenkovCountDecember.push_back( 0 );

                    IMPL::CalorimeterHitImpl* caloHit = new IMPL::CalorimeterHitImpl();
                    caloHit->setTime(float((*rawhit)->getTimeStamp())); //

                    if(float((*rawhit)->getAmplitude()&3)>2.5) caloHit->setEnergy(float((*rawhit)->getAmplitude()&3));         // 3rd treshold
                    else if(float((*rawhit)->getAmplitude()&3)>1.5) caloHit->setEnergy(float((*rawhit)->getAmplitude()&3)-1);  // 2nd treshold -1 to shift color to green
                    else caloHit->setEnergy(float((*rawhit)->getAmplitude()&3)+1);                                             // 1st treshold +1 to shift color to blue

                    //avoid two hit in the same cell
                    if(std::find(hitKeys.begin(),hitKeys.end(),aHitKey)!=hitKeys.end()){
                        IMPL::CalorimeterHitImpl* hit=
                                dynamic_cast<IMPL::CalorimeterHitImpl*>(col_event->getElementAt(std::distance(hitKeys.begin(),std::find(hitKeys.begin(),hitKeys.end(),aHitKey))));
                        float hitTime=hit->getTime();
                        if( fabs(time_peak-hitTime)>fabs(time_peak-time) ){
                            hit->setEnergy(caloHit->getEnergy());
                        }
                        continue;
                    }
                    // set the cell id
                    cd["I"] = I ;
                    cd["J"] = J ;
                    cd["K-1"] = K-1 ;
                    cd["M"] = 0 ;
                    cd["S-1"] = 3 ;
                    streamlog_out( DEBUG0 ) << " I == " << I
                                            << " J == " << J
                                            << " K == " << K
                                            << std::endl;
                    cd.setCellID( caloHit ) ;
                    if(std::find(zcut.begin(), zcut.end(), K)==zcut.end())
                        zcut.push_back(K);
                    caloHit->setPosition(pos);
                    col_event->addElement(caloHit);
                    hitKeys.push_back(aHitKey);

                    streamlog_out ( DEBUG0 ) << red << "Evt: " << evtnum << "\t Previous: " << previousEvtNbr << normal << std::endl;
                    previousEvtNbr=evtnum;
                    streamlog_out ( DEBUG0 ) << yellow << "Evt: " << evtnum << "\t Previous: " << previousEvtNbr << normal << std::endl;
                } // End if (timeWin)
            }  // End Dif_id != 1

            prevTime = time;
        }//loop over the hit
        hitKeys.clear();
    }catch(EVENT::DataNotAvailableException &e){
        streamlog_out(WARNING) << " collection not available" << std::endl;
    }
}


//===============================================
//===============================================
void TriventApp::init(const EVENT::LCParameters *pParameters) {
    trig_count = 0;
    _cerenkov1 = 0;
    _cerenkov2 = 0;
    _cerenkCount1=0;
    _cerenkCount2=0;
    _cerenkCountTot1=0;
    _cerenkCountTot2=0;
    //========================
    //readDifGeomFile(_mappingfile.c_str());

    char cnormal[8] =  {0x1b,'[','0',';','3','9','m',0};
    char cred[8]     = {0x1b,'[','1',';','3','1','m',0};
    char cgreen[8]   = {0x1b,'[','1',';','3','2','m',0};
    char cyellow[8]  = {0x1b,'[','1',';','3','3','m',0};
    char cblue[8]    = {0x1b,'[','1',';','3','4','m',0};
    char cmagenta[8] = {0x1b,'[','1',';','3','5','m',0};
    char cwhite[8]   = {0x1b,'[','1',';','3','9','m',0};

    normal   = cnormal;
    red      = cred;
    green    = cgreen;
    yellow   = cyellow;
    blue     = cblue;
    magenta  = cmagenta;
    white    = cwhite;

    // ------------
    // Parameters
    //

    // collection
    std::vector<std::string> hcalCollections;
    hcalCollections.push_back(std::string("DHCALRawHits"));

    // Option of output file with clean events
    _outFileName=pParameters->getStringVal("LCIOOutputFile");

    // Option of output file with clean events
    _rootFileName=pParameters->getStringVal("ROOTOutputFile");

    // Option of output file with noise
    _noiseFileName=pParameters->getStringVal("NOISEOutputFile");

    // Energy
    _beamEnergy = 10;
    _beamEnergy=pParameters->getIntVal("beamEnergy");

    // layer cut
    _LayerCut = 10;
    _LayerCut=pParameters->getIntVal("LayerCut");

    // noise cut
    _noiseCut = 10;
    _noiseCut=pParameters->getIntVal("noiseCut");

    // time windows
    _timeWin = 2;
    _timeWin=pParameters->getIntVal("timeWin");

    //maping on XML file
    _geomXML = "setup_geometry.xml";
    _geomXML=pParameters->getStringVal("setup_geometry");

    //maping on txt file
    _mappingfile = "mapping_ps.txt";
    _mappingfile=pParameters->getStringVal("DIFMapping");

    //inter layer
    _layer_gap = 9.;
    _layer_gap=pParameters->getFloatVal("layerGap");

    // electronic noise cut
    _elec_noise_cut = 5000;
    _elec_noise_cut=pParameters->getIntVal("electronic_noise_cut");

    // electronic noise cut
    _time2prev_event_cut = 0;
    _time2prev_event_cut=pParameters->getIntVal("_time2prev_event_cut");

    // Gain Correction
    GAIN_CORRECTION_MODE = false;
    GAIN_CORRECTION_MODE=pParameters->getIntVal("GAIN_CORRECTION_MODE");

    // Cerenkov Window ( Offset between CTag and Physic event, in TimeBin )
    _cerenkWindow = 7;
    _cerenkWindow = pParameters->getIntVal("CerenkovWindow");


    // Cerenkov Signal Length ( In time Bin )
    _cerenkLength = 1;
    _cerenkLength = pParameters->getIntVal("CerenkovLenght");

    // ------------


    _lcWriter = IOIMPL::LCFactory::getInstance()->createLCWriter() ;
    _lcWriter->setCompressionLevel( 0 ) ;
    _lcWriter->open(_outFileName.c_str(),EVENT::LCIO::WRITE_NEW) ;


    XMLReader(_geomXML.c_str());
    printDifGeom();
    evtnum=0; // event number
    previousEvtNbr =0;
}

//==================================================================================
void TriventApp::ClearVector()
{
    x.clear();
    y.clear();
    z.clear();
    fThr.clear();
    fDifId.clear();
    fNevt.clear();
    fEvtReconstructed.clear();
    fPrevTime.clear();
    fTime.clear();
    fTimePeak.clear();
    fPrevTimePeak.clear();
    fDeltaTimePeak.clear();
    fNhit.clear();
    fTriggerNr.clear();
    fCerenkovTag1.clear();
    fCerenkovTag2.clear();
    fCerenkovCount1.clear();
    fCerenkovCount2.clear();
    fCerenkovCountDecember.clear();

    fPion.clear();
    fElectron.clear();
    fMuon.clear();
}

//==================================================================================
// void TriventApp::ClearCerenkovVector(){


// }
//==================================================================================
void TriventApp::processRunHeader( EVENT::LCRunHeader * runHd ) {
    // _nRun++ ; From LcioToRoot
    // _nEvt = 0;
}

//==================================================================================
void TriventApp::processEvent( EVENT::LCEvent * evtP ){

    ClearVector();
    std::vector<IMPL::RawCalorimeterHitImpl*> _cerenkovHits;


    if (evtP != NULL)
    {
        try
        {
            for(unsigned int i=0; i< _hcalCollections.size(); i++){//!loop over collection
                try
                {
                    EVENT::LCCollection * col = NULL;
                    col = evtP ->getCollection(_hcalCollections[i].c_str());
                    //                    int numElements = col->getNumberOfElements();// hit number
                    trig_count++;
                    if ( (trig_count%1000) == 0)
                    {
                        streamlog_out( DEBUG5 ) << "Selected Event so far: " << _selectedNum << std::endl;
                        streamlog_out( MESSAGE ) << yellow << "Trigger number == " << trig_count << normal << std::endl;
                    }

                    if(col == NULL )  {
                        streamlog_out( WARNING )<< red << "TRIGGER SKIPED ...Col == NULL"<< normal <<std::endl;
                        break;
                    }

                    //                if(numElements > _elec_noise_cut)  {
                    //                  streamlog_out( MESSAGE ) << red << "TRIGGER SKIPED ...NoiseCut"<< normal <<std::endl;
                    //                  break;
                    //                }

                    // set raw hits
                    _trigger_raw_hit.clear();
                    std::vector<int> vTrigger;

                    for (int ihit(0); ihit < col->getNumberOfElements(); ++ihit) {// loop over the hits
                        EVENT::RawCalorimeterHit *raw_hit =
                                dynamic_cast<EVENT::RawCalorimeterHit*>( col->getElementAt(ihit)) ;
                        if (NULL != raw_hit)
                        {
                            unsigned int difid=0;
                            difid = raw_hit->getCellID0()&0xFF;
                            if (difid==3) {
                                IMPL::RawCalorimeterHitImpl *pCerenkRawHit = new IMPL::RawCalorimeterHitImpl;

                                pCerenkRawHit->setAmplitude(raw_hit->getAmplitude());
                                pCerenkRawHit->setCellID0(raw_hit->getCellID0());
                                pCerenkRawHit->setCellID1(raw_hit->getCellID1());
                                pCerenkRawHit->setTimeStamp(raw_hit->getTimeStamp());
                                _cerenkovHits.push_back(pCerenkRawHit);
                                _trigger_raw_hit.push_back(pCerenkRawHit);
                                
                            }else   _trigger_raw_hit.push_back(raw_hit);

                            //extract abolute bcid information:
                            if(ihit==0)
                            {
                                if (difid==0) continue;

                                std::stringstream pname("");
                                pname <<"DIF"<<difid<<"_Triggers";

                                col->getParameters().getIntVals(pname.str(),vTrigger);
                                if (vTrigger.size()!=0)
                                {
                                    _bcid1=vTrigger[4] ;
                                    _bcid2=vTrigger[3] ;
                                    unsigned long long Shift=16777216ULL;//to shift the value from the 24 first bits
                                    unsigned long long theBCID_=_bcid1*Shift+_bcid2;
                                    streamlog_out( DEBUG1 ) << "trigger time : " << theBCID_ << std::endl;
                                }
                            }
                        }
                    }
                    getMaxTime();
                    std::vector<int> time_spectrum = getTimeSpectrum();

                    //---------------------------------------------------------------
                    //! Find the condidate event
                    int ibin=0;
                    int bin_c_prev = -2 * _timeWin; //  the previous bin center

                    int time_prev = 0;
                    while(ibin < (_maxtime+1))
                    {
                        if(time_spectrum[ibin] >= _noiseCut
                                && time_spectrum[ibin] >  time_spectrum[ibin+1]
                                && time_spectrum[ibin] >= time_spectrum[ibin+1]
                                && time_spectrum[ibin] >= time_spectrum[ibin+1]
                                && time_spectrum[ibin] >= time_spectrum[ibin+2] )
                        {
                            IMPL::LCEventImpl*  evt = new IMPL::LCEventImpl() ;     // create the event

                            //---------- set event paramters ------
                            const std::string parname_trigger = "trigger";
                            const std::string parname_energy  = "beamEnergy";
                            const std::string parname_bcid1 = "bcid1";
                            const std::string parname_bcid2 = "bcid2";
                            const std::string parname_cerenkov1 = "cerenkov1";
                            const std::string parname_cerenkov2 = "cerenkov2";

                            evt->parameters().setValue(parname_trigger,evtP->getEventNumber());
                            evt->parameters().setValue(parname_energy , _beamEnergy);
                            evt->parameters().setValue(parname_bcid1 , _bcid1);
                            evt->parameters().setValue(parname_bcid2 , _bcid2);
                            evt->setRunNumber( evtP->getRunNumber()) ;

                            //-------------------------------------
                            IMPL::LCCollectionVec* outcol = new IMPL::LCCollectionVec(EVENT::LCIO::CALORIMETERHIT);
                            evtnum++;
                            _cerenkCount1 = 0;
                            _cerenkCount2 = 0;
                            TriventApp::eventBuilder(outcol,ibin,bin_c_prev);
                            evt->parameters().setValue(parname_cerenkov1, _cerenkCount1); // Value determined in the eventBuilder
                            evt->parameters().setValue(parname_cerenkov2, _cerenkCount2); // Value determined in the eventBuilder
                            // ->Need to be after that function
                            // Print if there is more than 1 CTag in the current event
                            if (_cerenkCount1>1) streamlog_out ( MESSAGE )<< red << "CerenCount1: " << _cerenkCount1 << normal << std::endl;
                            if (_cerenkCount1>2) streamlog_out ( MESSAGE )<< red << "CerenCount2: " << _cerenkCount2 << normal << std::endl;

                            streamlog_out( DEBUG1 ) << "zcut.size() = " << zcut.size() << "\t _LayerCut = " << _LayerCut << std::endl;
                            if( (int)zcut.size() >_LayerCut &&                                  // the min layer numb cut
                                    abs(int(ibin)-time_prev) > _time2prev_event_cut){ // time2prev event  cut

                                streamlog_out( DEBUG5 ) << "Evt number : " <<  evtnum << std::endl;
                                streamlog_out( DEBUG5 ) <<green<<" Trivent find event at :==> "<< red << ibin
                                                       <<green<<"\t :Nhit: ==> "<< magenta
                                                      <<outcol->getNumberOfElements() << normal <<std::endl;
                                evt->setEventNumber( evtnum ) ;
                                evt->addCollection(outcol, "SDHCAL_HIT");
                                _lcWriter->writeEvent( evt ) ;
                                evts.push_back(evt);
                                _selectedNum++;

                            }else{
                                _rejectedNum++;

                                // ------ The following part is to delete event in the tree that were rejected by trivent
                                // ------  Probably not very elegant/efficient for now ...
                                // ------

                                // --- To account for multiple hit/cell
                                int sum_of_elems(0);
                                std::vector<int>::iterator i=fNhit.end()-outcol->getNumberOfElements();
                                while(int(*i)>0) {
                                    //  If there is more than one hit in a cell: i would be > 0
                                    *i--;   streamlog_out ( MESSAGE ) << red << "Multiple hit: " << *i << normal << std::endl;
                                }
                                for(*i; i!=fNhit.end();++i)
                                    sum_of_elems += 1;
                                // ---

                                for (int j=0; j < sum_of_elems; j++ ){
                                    x.pop_back();
                                    z.pop_back();
                                    fThr.pop_back();
                                    fDifId.pop_back();
                                    fPrevTime.pop_back();
                                    fTime.pop_back();
                                    fTimePeak.pop_back();
                                    fPrevTimePeak.pop_back();
                                    fDeltaTimePeak.pop_back();
                                    fNevt.pop_back();
                                    fEvtReconstructed.pop_back();
                                    fTriggerNr.pop_back();
                                    fNhit.pop_back();
                                    fCerenkovTag1.pop_back();
                                    fCerenkovTag2.pop_back();
                                    fCerenkovCount1.pop_back();
                                    fCerenkovCount2.pop_back();
                                    fCerenkovCountDecember.pop_back();
                                }
                                // ------

                                 delete outcol;
                            }
                            time_prev = ibin;
                            delete evt; evt=NULL;

                            bin_c_prev = ibin;
                            ibin = ibin+_timeWin;
                        }else{ibin++;}
                    }
                    //                    _tree->Fill();
                }catch (lcio::DataNotAvailableException zero) {}
            }
        }catch (lcio::DataNotAvailableException err) {}
    }


}
//==============================================================
void TriventApp::end()
{
    streamlog_out( MESSAGE )<< "Trivent Rejected "<<_rejectedNum <<" Candidate event"<<std::endl;
    streamlog_out( MESSAGE )<< "Trivent Selected "<<_selectedNum <<" Candidate event"<<std::endl;
    streamlog_out( MESSAGE )<< "Trivent end"<<std::endl;
    streamlog_out( MESSAGE )<< "Total Cerenkov1 Tags: " << _cerenkCountTot1 << std::endl;
    streamlog_out( MESSAGE )<< "Total Cerenkov2 Tags: " << _cerenkCountTot2 << std::endl;

    _lcWriter->close();
}
