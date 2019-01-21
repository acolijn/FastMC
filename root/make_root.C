#include "TFile.h"
#include "TTree.h"

void make_root(string fname) {

   const Int_t kMaxHit = 1000;
   Int_t ievent;
   Float_t weight;
   Float_t edep;
   Float_t x0;
   Float_t y0;
   Float_t z0;
   Int_t nhit;
   Float_t xh[kMaxHit];
   Float_t yh[kMaxHit];
   Float_t zh[kMaxHit];
   Float_t de[kMaxHit];

   string rootfile = fname + ".root";
   TFile f(rootfile.c_str(),"recreate");
   TTree *t3 = new TTree("t","MC");
   t3->Branch("ieve",&ievent,"ieve/I");
   t3->Branch("w",&weight,"w/F");
   t3->Branch("e",&edep,"e/F");
   t3->Branch("x0",&x0,"x0/F");
   t3->Branch("y0",&y0,"y0/F");
   t3->Branch("z0",&z0,"z0/F");
   t3->Branch("nhit",&nhit,"nhit/I");
   t3->Branch("xh",xh,"xh[nhit]/F");
   t3->Branch("yh",yh,"yh[nhit]/F");
   t3->Branch("zh",zh,"zh[nhit]/F");
   t3->Branch("de",de,"de[nhit]/F");

   ifstream in;
   //in.open("../mcdata/testdata.txt");
   in.open(fname.c_str());
   ievent = 0;
   while(ievent >= 0) {    
      in >> ievent >> nhit >>weight >> edep >> x0 >> y0 >> z0;

      for (Int_t i=0; i<nhit; i++){
         in >> xh[i] >> yh[i] >> zh[i] >> de[i];
      }
      t3->Fill();

      if (ievent%100000 == 0) cout <<"processed "<<ievent<< " events"<<endl;
   }
   t3->Print();
   f.cd();
   t3->Write();

}
