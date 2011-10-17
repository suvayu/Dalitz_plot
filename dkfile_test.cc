// dkfile_test.cc -- Dalitz plot from MC info

/*
  This macro generates a Dalitz plot for the 3-body decay of the Ds
  meson into KKpi using the MC information dumped by Gauss when run in
  Generator only mode. The intention is to cross check the Dalitz
  model used by EvtGen for such decays against observations from BaBar
  and other experiments (e.g. arXiv:1011.4190 [hep-ex]).
*/

#include <list>
#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TH2D.h"
#include "TPad.h"
#include "TStyle.h"

int dkfile_test(std::string fileglob="GaussMonitor.root", int DEBUG=0)
{
  gROOT->Reset();

  // To be modified later for TChains
  TFile f( fileglob.c_str() ) ;
  TTree * T = (TTree*) f.Get( "GeneratorFullMonitor/1" ) ;

  Int_t nEvents = T -> GetEntries() ;
  printf( "Nentries = %d\n" , nEvents ) ;

  Int_t nPart, nInter;
  Float_t pdgId[2000] ;
  Float_t nDau[2000] ;
  Float_t pdgIdDau1[2000] ;
  Float_t pdgIdDau2[2000] ;
  Float_t pdgIdDau3[2000] ;
  Float_t pdgIdDau4[2000] ;
  Float_t pdgIdDau5[2000] ;
  Float_t pdgIdDau6[2000] ;
  Float_t pdgIdMother[2000] ;
  Float_t IdMother[2000] ;
  Float_t e[2000] , px[2000] , py[2000], pz[2000] ;

  T -> SetBranchAddress( "NInter" , &nInter ) ;
  T -> SetBranchAddress( "NPart" , &nPart ) ;
  T -> SetBranchAddress( "pdgId" , &pdgId ) ;
  T -> SetBranchAddress( "nDau" , &nDau ) ;
  T -> SetBranchAddress( "pdgIdDau1" , &pdgIdDau1 ) ;
  T -> SetBranchAddress( "pdgIdDau2" , &pdgIdDau2 ) ;
  T -> SetBranchAddress( "pdgIdDau3" , &pdgIdDau3 ) ;
  T -> SetBranchAddress( "pdgIdDau4" , &pdgIdDau4 ) ;
  T -> SetBranchAddress( "pdgIdDau5" , &pdgIdDau5 ) ;
  T -> SetBranchAddress( "pdgIdDau6" , &pdgIdDau6 ) ;
  T -> SetBranchAddress( "pdgIdMother" , &pdgIdMother ) ;
  T -> SetBranchAddress( "indexMother" , &IdMother ) ;
  T -> SetBranchAddress( "e"  , &e  ) ;
  T -> SetBranchAddress( "px" , &px ) ;
  T -> SetBranchAddress( "py" , &py ) ;
  T -> SetBranchAddress( "pz" , &pz ) ;

  // for the Dalitz plot
  TLorentzVector lv_Kp(0,0,0,0), lv_Km(0,0,0,0), lv_pip(0,0,0,0),
    lv_pim(0,0,0,0), lv_12(0,0,0,0), lv_23(0,0,0,0);

  TH2D *hDalitz = new TH2D( "hDalitz", "Dalitz plot", 20, 0, 4, 20, 0, 4);

  // temporary variables used later
  Float_t PDG(0), PDGm(0);
  Int_t sizep(0), sizem(0);
  std::vector<Int_t> Dsp, Dsm;
  bool isDsp(false), isDsm(false), hasKp(false), hasKm(false),
    haspip(false), haspim(false);

  // test with some const (e.g. 2000) instead of nEvents
  for ( Int_t i(0); i < nEvents ; ++i )
    {
      T -> GetEntry( i ) ;

      if (DEBUG >= 4) std::cout << "DEBUG1: NInter : " << nInter
                                << ", NPart : " << nPart << std::endl;

      for ( Int_t j(0); j < nPart; ++j )
        {
          /*
            | particle |   Lb | Ds+ | Ds*+ |  K+ | pi+ |
            |----------+------+-----+------+-----+-----|
            | PDG id   | 5122 | 431 |  433 | 321 | 211 |

            Decay: Ds- -> K+ K- pi-
                   Ds+ -> K- K+ pi+

            1. Find K+/K-/pi+/pi- coming from Ds+/Ds-
            2. Ensure pi and Ds has same charge to eliminate other channels
            3. Make separate list of index of K/pi based on Ds+/Ds-
            5. Check if it is one of the desired decay mode
            6. Reconstruct Ds+ and Ds-
            7. Calculate m12 and m23 for Dalitz plot
          */

          PDG = pdgId[j];
          PDGm = pdgIdMother[j];

          if (( fabs( PDG ) == 321 && fabs( PDGm ) == 431 ) // K from Ds
              || ( fabs( PDG ) == 211 && fabs( PDGm ) == 431
                   && (PDG * PDGm) > 0 )) // pi from Ds with same sign
            {
              if (DEBUG >= 3) printf("  DEBUG3: PDG:%+4.0f, idx:% 4d, PDG mom:%+4.0f\n",
                                     pdgId[j], j, pdgIdMother[j]);

              if ( PDGm ==  431 ) Dsp.push_back(j);
              if ( PDGm == -431 ) Dsm.push_back(j);
              // pdgIdMother[] is giving random results, maybe problem
              // with type casting. abandon for now
            }
        }

      sizep = (Int_t) Dsp.size();
      sizem = (Int_t) Dsm.size();

      // reconstruct Ds+
      if ( sizep == 3 )
        {
          for ( Int_t k(0); k < sizep; ++k)
            {
              if ( pdgId[ Dsp[k] ] ==  321 )
                {
                  hasKp  = true;
                  lv_Kp .SetPxPyPzE( px[ Dsp[k] ] / 1E3, py[ Dsp[k] ] / 1E3,
                                     pz[ Dsp[k] ] / 1E3,  e[ Dsp[k] ] / 1E3 );
                }

              if ( pdgId[ Dsp[k] ] == -321 )
                {
                  hasKm  = true;
                  lv_Km .SetPxPyPzE( px[ Dsp[k] ] / 1E3, py[ Dsp[k] ] / 1E3,
                                     pz[ Dsp[k] ] / 1E3,  e[ Dsp[k] ] / 1E3 );
                }

              if ( pdgId[ Dsp[k] ] ==  211 )
                {
                  haspip = true;
                  lv_pip.SetPxPyPzE( px[ Dsp[k] ] / 1E3, py[ Dsp[k] ] / 1E3,
                                     pz[ Dsp[k] ] / 1E3,  e[ Dsp[k] ] / 1E3 );
                }

              if (DEBUG >= 2) printf("  DEBUG2: PDG:%+4.0f, idx:%4d\n",
                                     pdgId[ Dsp[k] ], Dsp[k]);
            }

          if ( hasKp && hasKm && haspip ) isDsp = true;
          if (DEBUG >= 2) std::cout << "  DEBUG2: Ds+:" << isDsp << std::endl;

          // cleanup
          hasKp = hasKm = haspip = false;
        }

      // reconstruct Ds-
      if ( sizem == 3 )
        {
          for ( Int_t k(0); k < sizem; ++k)
            {
              if ( pdgId[ Dsm[k] ] ==  321 )
                {
                  hasKp  = true;
                  lv_Kp.SetPxPyPzE( px[ Dsm[k] ] / 1E3, py[ Dsm[k] ] / 1E3,
                                    pz[ Dsm[k] ] / 1E3,  e[ Dsm[k] ] / 1E3 );
                }

              if ( pdgId[ Dsm[k] ] == -321 )
                {
                  hasKm  = true;
                  lv_Km.SetPxPyPzE( px[ Dsm[k] ] / 1E3, py[ Dsm[k] ] / 1E3,
                                    pz[ Dsm[k] ] / 1E3,  e[ Dsm[k] ] / 1E3 );
                }

              if ( pdgId[ Dsm[k] ] == -211 )
                {
                  haspim = true;
                  lv_pim.SetPxPyPzE( px[ Dsm[k] ] / 1E3, py[ Dsm[k] ] / 1E3,
                                     pz[ Dsm[k] ] / 1E3,  e[ Dsm[k] ] / 1E3 );
                }

              if (DEBUG >= 2) printf("  DEBUG2: PDG:%+4.0f, idx:%4d\n",
                                     pdgId[ Dsm[k] ], Dsm[k]);
            }

          if ( hasKp && hasKm && haspim ) isDsm = true;
          if (DEBUG >= 2) std::cout << "  DEBUG2: Ds-:" << isDsm << std::endl;

          // cleanup
          hasKp = hasKm = haspim = false;
        }

      if (isDsp)                // Ds+ -> K- K+ pi+
        {
          lv_12 = lv_Km + lv_Kp;
          lv_23 = lv_Kp + lv_pip; // ??
          if (DEBUG >= 1) printf("  DEBUG1: m12: %1.3f GeV, m23: %1.3f GeV\n",
                                 lv_12.M2(), lv_23.M2() );
          hDalitz->Fill( lv_12.M2(), lv_23.M2() );
        }

      if (isDsm)                // Ds- -> K+ K- pi-
        {
          lv_12 = lv_Km + lv_Kp;
          lv_23 = lv_Km + lv_pim; // ??
          if (DEBUG >= 1) printf("  DEBUG1: m12: %1.3f GeV, m23: %1.3f GeV\n",
                                 lv_12.M2(), lv_23.M2() );
          hDalitz->Fill( lv_12.M2(), lv_23.M2() );
        }

      // cleanup
      Dsp.clear();
      Dsm.clear();
      isDsp = isDsm = false;
      lv_12.SetXYZT(0,0,0,0);
      lv_23.SetXYZT(0,0,0,0);
    }

  // gStyle->SetTitleOffset(1, "XY");
  hDalitz->SetXTitle("m^{2}(K^{+}K^{-}) GeV^{2}/c^{4}");
  hDalitz->SetYTitle("m^{2}(K^{+}#pi^{+} / K^{-}#pi^{-}) GeV^{2}/c^{4}");
  hDalitz->Draw("SCAT");
  gPad->Print("Dalitz_plot_Ds.png");
  if (DEBUG >= 1) hDalitz->Print("all");

  return 0 ;
}


/*

  List of TTree branches:

  |---------+-------------+----------------------|
  | Tree    | 1           | MCTruth              |
  |---------+-------------+----------------------|
  | Br    0 | runN        | runN/I               |
  |---------+-------------+----------------------|
  | Br    1 | evtN        | evtN/I               |
  |---------+-------------+----------------------|
  | Br    2 | NInter      | NInter/I             |
  |---------+-------------+----------------------|
  | Br    3 | procId      | procId[NInter]/F     |
  |---------+-------------+----------------------|
  | Br    4 | s_hat       | s_hat[NInter]/F      |
  |---------+-------------+----------------------|
  | Br    5 | t_hat       | t_hat[NInter]/F      |
  |---------+-------------+----------------------|
  | Br    6 | u_hat       | u_hat[NInter]/F      |
  |---------+-------------+----------------------|
  | Br    7 | pt_hat      | pt_hat[NInter]/F     |
  |---------+-------------+----------------------|
  | Br    8 | x1_Bjork    | x1_Bjork[NInter]/F   |
  |---------+-------------+----------------------|
  | Br    9 | x2_Bjork    | x2_Bjork[NInter]/F   |
  |---------+-------------+----------------------|
  | Br   10 | NPart       | NPart/I              |
  |---------+-------------+----------------------|
  | Br   11 | e           | e[NPart]/F           |
  |---------+-------------+----------------------|
  | Br   12 | px          | px[NPart]/F          |
  |---------+-------------+----------------------|
  | Br   13 | py          | py[NPart]/F          |
  |---------+-------------+----------------------|
  | Br   14 | pz          | pz[NPart]/F          |
  |---------+-------------+----------------------|
  | Br   15 | vxProd      | vxProd[NPart]/F      |
  |---------+-------------+----------------------|
  | Br   16 | vyProd      | vyProd[NPart]/F      |
  |---------+-------------+----------------------|
  | Br   17 | vzProd      | vzProd[NPart]/F      |
  |---------+-------------+----------------------|
  | Br   18 | vtProd      | vtProd[NPart]/F      |
  |---------+-------------+----------------------|
  | Br   19 | vxDecay     | vxDecay[NPart]/F     |
  |---------+-------------+----------------------|
  | Br   20 | vyDecay     | vyDecay[NPart]/F     |
  |---------+-------------+----------------------|
  | Br   21 | vzDecay     | vzDecay[NPart]/F     |
  |---------+-------------+----------------------|
  | Br   22 | vtDecay     | vtDecay[NPart]/F     |
  |---------+-------------+----------------------|
  | Br   23 | pdgId       | pdgId[NPart]/F       |
  |---------+-------------+----------------------|
  | Br   24 | nDau        | nDau[NPart]/F        |
  |---------+-------------+----------------------|
  | Br   25 | pdgIdMother | pdgIdMother[NPart]/F |
  |---------+-------------+----------------------|
  | Br   26 | pdgIdDau1   | pdgIdDau1[NPart]/F   |
  |---------+-------------+----------------------|
  | Br   27 | pdgIdDau2   | pdgIdDau2[NPart]/F   |
  |---------+-------------+----------------------|
  | Br   28 | pdgIdDau3   | pdgIdDau3[NPart]/F   |
  |---------+-------------+----------------------|
  | Br   29 | pdgIdDau4   | pdgIdDau4[NPart]/F   |
  |---------+-------------+----------------------|
  | Br   30 | pdgIdDau5   | pdgIdDau5[NPart]/F   |
  |---------+-------------+----------------------|
  | Br   31 | pdgIdDau6   | pdgIdDau6[NPart]/F   |
  |---------+-------------+----------------------|
  | Br   32 | indexMother | indexMother[NPart]/F |
  |---------+-------------+----------------------|
  | Br   33 | indexInter  | indexInter[NPart]/F  |
  |---------+-------------+----------------------|

*/
