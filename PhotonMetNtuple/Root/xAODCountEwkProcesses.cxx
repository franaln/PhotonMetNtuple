#include <TSystem.h>

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <PhotonMetNtuple/xAODCountEwkProcesses.h>
#include "EventLoop/OutputStream.h"
#include <TTreeFormula.h>

// this is needed to distribute the algorithm to the workers
ClassImp(xAODCountEwkProcesses)


xAODCountEwkProcesses::xAODCountEwkProcesses()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}

EL::StatusCode xAODCountEwkProcesses::setupJob(EL::Job &job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  job.useXAOD();
  
  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init("xAODCountEwkProcesses").ignore(); // call before opening first file
  
  // tell EventLoop about our output:
  EL::OutputStream out("output");
  job.outputAdd(out);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses::histInitialize()
{
  h_hp = new TH1D("hp", "hp", 30, 0., 30.);
  h_hp->GetXaxis()->SetBinLabel(1,  "111 (chi^{0}_1 - chi^{0}_1)");
  h_hp->GetXaxis()->SetBinLabel(2,  "112 (chi^{0}_1 - chi^{0}_2)");
  h_hp->GetXaxis()->SetBinLabel(3,  "113 (chi^{0}_1 - chi^{0}_3)");
  h_hp->GetXaxis()->SetBinLabel(4,  "114 (chi^{0}_1 - chi^{0}_4)");
  h_hp->GetXaxis()->SetBinLabel(5,  "115 (chi^{0}_1 - chi^{+}_1)");
  h_hp->GetXaxis()->SetBinLabel(6,  "116 (chi^{0}_1 - chi^{+}_2)");
  h_hp->GetXaxis()->SetBinLabel(7,  "117 (chi^{0}_1 - chi^{-}_1)");
  h_hp->GetXaxis()->SetBinLabel(8,  "118 (chi^{0}_1 - chi^{-}_2)");
  h_hp->GetXaxis()->SetBinLabel(9,  "122 (chi^{0}_2 - chi^{0}_2)");
  h_hp->GetXaxis()->SetBinLabel(10, "123 (chi^{0}_2 - chi^{0}_3)");
  h_hp->GetXaxis()->SetBinLabel(11, "124 (chi^{0}_2 - chi^{0}_4)");
  h_hp->GetXaxis()->SetBinLabel(12, "125 (chi^{0}_2 - chi^{+}_1)");
  h_hp->GetXaxis()->SetBinLabel(13, "126 (chi^{0}_2 - chi^{+}_2)");
  h_hp->GetXaxis()->SetBinLabel(14, "127 (chi^{0}_2 - chi^{-}_1)");
  h_hp->GetXaxis()->SetBinLabel(15, "128 (chi^{0}_2 - chi^{-}_2)");
  h_hp->GetXaxis()->SetBinLabel(16, "133 (chi^{0}_3 - chi^{0}_3)");
  h_hp->GetXaxis()->SetBinLabel(17, "134 (chi^{0}_3 - chi^{0}_4)");
  h_hp->GetXaxis()->SetBinLabel(18, "135 (chi^{0}_3 - chi^{+}_1)");
  h_hp->GetXaxis()->SetBinLabel(19, "136 (chi^{0}_3 - chi^{+}_2)");
  h_hp->GetXaxis()->SetBinLabel(20, "137 (chi^{0}_3 - chi^{-}_1)");
  h_hp->GetXaxis()->SetBinLabel(21, "138 (chi^{0}_3 - chi^{-}_2)");
  h_hp->GetXaxis()->SetBinLabel(22, "144 (chi^{0}_4 - chi^{0}_4)");
  h_hp->GetXaxis()->SetBinLabel(23, "145 (chi^{0}_4 - chi^{+}_1)");
  h_hp->GetXaxis()->SetBinLabel(24, "146 (chi^{0}_4 - chi^{+}_2)");
  h_hp->GetXaxis()->SetBinLabel(25, "147 (chi^{0}_4 - chi^{-}_1)");
  h_hp->GetXaxis()->SetBinLabel(26, "148 (chi^{0}_4 - chi^{-}_2)");
  h_hp->GetXaxis()->SetBinLabel(27, "157 (chi^{+}_1 - chi^{-}_1)");
  h_hp->GetXaxis()->SetBinLabel(28, "158 (chi^{+}_1 - chi^{-}_2)");
  h_hp->GetXaxis()->SetBinLabel(29, "167 (chi^{+}_2 - chi^{-}_1)");
  h_hp->GetXaxis()->SetBinLabel(30, "168 (chi^{+}_2 - chi^{-}_2)");

  TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
  h_hp->SetDirectory(out_dir);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses::fileExecute()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses::changeInput(bool firstFile)
{
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses::initialize()
{
  m_event = wk()->xaodEvent();


  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses :: execute ()
{
  
  const xAOD::TruthParticleContainer* truthP = 0;
  if(!m_event->retrieve(truthP, "TruthParticles").isSuccess()) {
    Error("xAODCountEwkProcesses", "Failed to retrieve truth particles collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // FindSusyHardProc(*truth_particles, pdg1, pdg2);
  int pdg1 = 0;
  int pdg2 = 0;

  const xAOD::TruthParticle* firstsp(0);
  const xAOD::TruthParticle* secondsp(0);

  if (!truthP || truthP->empty()) {
    return false;
  }
  for (const auto& tp : *truthP) {

    //check ifSUSY particle
    if ((abs(tp->pdgId()) > 1000000 && abs(tp->pdgId()) < 1000007) || // squarkL
        (abs(tp->pdgId()) > 1000010 && abs(tp->pdgId()) < 1000017) || // sleptonL
        (abs(tp->pdgId()) > 2000000 && abs(tp->pdgId()) < 2000007) || // squarkR
        (abs(tp->pdgId()) > 2000010 && abs(tp->pdgId()) < 2000017) || // sleptonR
        (abs(tp->pdgId()) > 1000020 && abs(tp->pdgId()) < 1000040)) { // gauginos

      if (tp->nParents() != 0) {
        if ( tp->parent(0)->absPdgId()  < 1000000) {
          if (!firstsp) {
            firstsp = tp;
          } else if (!secondsp) {
            secondsp = tp;
          } else {
            if (firstsp->nChildren() != 0 && tp->barcode() == firstsp->child(0)->barcode()) {
              firstsp = tp;
            }
            else if (secondsp->nChildren() != 0 && tp->barcode() == secondsp->child(0)->barcode()) {
              secondsp = tp;
            }
            else if (firstsp->nChildren() != 0 && firstsp->child(0)->barcode() == secondsp->barcode()) {
              firstsp = secondsp;
              secondsp = tp;
            }
            else if (secondsp->nChildren() != 0 && secondsp->child(0)->barcode() == firstsp->barcode()) {
              secondsp = firstsp;
              firstsp = tp;
            }
          }
        }
      }
    }
  }

  // quit if no sparticles found
  if (!firstsp && !secondsp) return true; // should find none or two

  if (firstsp->nChildren() == 1) {
    for (const auto& tp : *truthP) {
      if (tp->barcode() == firstsp->child(0)->barcode() && tp->pdgId() != firstsp->pdgId()) {
        firstsp = tp;
        break;
      }
    }
  }

  if (secondsp->nChildren() == 1) {
    for (const auto& tp : *truthP) {
      if (tp->barcode() == secondsp->child(0)->barcode() && tp->pdgId() != secondsp->pdgId()) {
        secondsp = tp;
        break;
      }
    }
  }

  if (abs(firstsp->pdgId()) > 1000000) pdg1 = firstsp->pdgId();
  if (abs(secondsp->pdgId()) > 1000000) pdg2 = secondsp->pdgId();


  // Get Final state and fill histogram
  int fs = GetFinalState(pdg1, pdg2);

  if (fs == 111)
    h_hp->Fill(0.5);
  else if (fs == 112)
    h_hp->Fill(1.5);
  else if (fs == 113)
    h_hp->Fill(2.5);
  else if (fs == 114)
    h_hp->Fill(3.5);
  else if (fs == 115)
    h_hp->Fill(4.5);
  else if (fs == 116)
    h_hp->Fill(5.5);
  else if (fs == 117)
    h_hp->Fill(6.5);
  else if (fs == 118)
    h_hp->Fill(7.5);
  else if (fs == 122)
    h_hp->Fill(8.5);
  else if (fs == 123)
    h_hp->Fill(9.5);
  else if (fs == 124)
    h_hp->Fill(10.5);
  else if (fs == 125)
    h_hp->Fill(11.5);
  else if (fs == 126)
    h_hp->Fill(12.5);
  else if (fs == 127)
    h_hp->Fill(13.5);
  else if (fs == 128)
    h_hp->Fill(14.5);
  else if (fs == 133)
    h_hp->Fill(15.5);
  else if (fs == 134)
    h_hp->Fill(16.5);
  else if (fs == 135)
    h_hp->Fill(17.5);
  else if (fs == 136)
    h_hp->Fill(18.5);
  else if (fs == 137)
    h_hp->Fill(19.5);
  else if (fs == 138)
    h_hp->Fill(20.5);
  else if (fs == 144)
    h_hp->Fill(21.5);
  else if (fs == 145)
    h_hp->Fill(22.5);
  else if (fs == 146)
    h_hp->Fill(23.5);
  else if (fs == 147)
    h_hp->Fill(24.5);
  else if (fs == 148)
    h_hp->Fill(25.5);
  else if (fs == 157)
    h_hp->Fill(26.5);
  else if (fs == 158)
    h_hp->Fill(27.5);
  else if (fs == 167)
    h_hp->Fill(28.5);
  else if (fs == 168)
    h_hp->Fill(29.5);
  else
    std::cout << "weird: fs =" << fs << std::endl;


  // std::cout << pdg1 << pdg2 << std::endl;

  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses::postExecute()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses::finalize()
{
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODCountEwkProcesses::histFinalize()
{
  return EL::StatusCode::SUCCESS;
}



unsigned int xAODCountEwkProcesses::GetFinalState(const int SUSY_Spart1_pdgId, const int SUSY_Spart2_pdgId)
{
  int ngluino = 0;
  int nsquark = 0; // (up and down type without bottom/top)
  int nantisquark = 0; // (up and down type without bottom/top)

  int nsbottom = 0;
  int nstop = 0;
  int nsbottom2 = 0;
  int nstop2 = 0;
  int nantisbottom = 0;
  int nantistop = 0;
  int nantisbottom2 = 0;
  int nantistop2 = 0;

  int nchi01 = 0;
  int nchi02 = 0;
  int nchi03 = 0;
  int nchi04 = 0;
  int nch1plus = 0;
  int nch2plus = 0;
  int nch1minus = 0;
  int nch2minus = 0;

  // sleptons
  int nsmuonRplus = 0;
  int nsmuonRminus = 0;
  int nselecRplus = 0;
  int nselecRminus = 0;

  int nsmuonLplus = 0;
  int nsmuonLminus = 0;
  int nselecLplus = 0;
  int nselecLminus = 0;

  int nstau1plus = 0;
  int nstau1minus = 0;
  int nstau2plus = 0;
  int nstau2minus = 0;

  // snutrinos
  int nselnuL = 0;
  int nsmunuL = 0;
  int nstaunuL = 0;

  int nother = 0;

  //Classification of the event follows (gg, sq...):

  if      (abs(SUSY_Spart1_pdgId) == 1000022) nchi01++;
  else if (abs(SUSY_Spart1_pdgId) == 1000023) nchi02++;
  else if (abs(SUSY_Spart1_pdgId) == 1000025) nchi03++;
  else if (abs(SUSY_Spart1_pdgId) == 1000035) nchi04++;
  else if (    SUSY_Spart1_pdgId == 1000024) nch1plus++;
  else if (    SUSY_Spart1_pdgId == -1000024) nch1minus++;
  else if (    SUSY_Spart1_pdgId == 1000037) nch2plus++;
  else if (    SUSY_Spart1_pdgId == -1000037) nch2minus++;
  else if (    SUSY_Spart1_pdgId == 1000021) ngluino++;
  else if ((abs(SUSY_Spart1_pdgId) > 1000000 && abs(SUSY_Spart1_pdgId) <= 1000004) || (abs(SUSY_Spart1_pdgId) > 2000000 && abs(SUSY_Spart1_pdgId) <= 2000004)) {
    if (SUSY_Spart1_pdgId > 0) nsquark++;
    else nantisquark++;
  }
  else if (SUSY_Spart1_pdgId == 1000005) nsbottom++;
  else if (SUSY_Spart1_pdgId == 1000006) nstop++;
  else if (SUSY_Spart1_pdgId == 2000005) nsbottom2++;
  else if (SUSY_Spart1_pdgId == 2000006) nstop2++;
  else if (SUSY_Spart1_pdgId == -1000005) nantisbottom++;
  else if (SUSY_Spart1_pdgId == -1000006) nantistop++;
  else if (SUSY_Spart1_pdgId == -2000005) nantisbottom2++;
  else if (SUSY_Spart1_pdgId == -2000006) nantistop2++;
  else if (SUSY_Spart1_pdgId == 2000011) nselecRminus++;
  else if (SUSY_Spart1_pdgId == -2000011) nselecRplus++;
  else if (SUSY_Spart1_pdgId == 1000011) nselecLminus++;
  else if (SUSY_Spart1_pdgId == -1000011) nselecLplus++;
  else if (abs(SUSY_Spart1_pdgId) == 1000012) nselnuL++;
  else if (SUSY_Spart1_pdgId == 2000013) nsmuonRminus++;
  else if (SUSY_Spart1_pdgId == -2000013) nsmuonRplus++;
  else if (SUSY_Spart1_pdgId == 1000013) nsmuonLminus++;
  else if (SUSY_Spart1_pdgId == -1000013) nsmuonLplus++;
  else if (abs(SUSY_Spart1_pdgId) == 1000014) nsmunuL++;
  else if (SUSY_Spart1_pdgId == 1000015) nstau1minus++;
  else if (SUSY_Spart1_pdgId == -1000015) nstau1plus++;
  else if (SUSY_Spart1_pdgId == 2000015) nstau2minus++;
  else if (SUSY_Spart1_pdgId == -2000015) nstau2plus++;
  else if (abs(SUSY_Spart1_pdgId) == 1000016) nstaunuL++;
  else nother++;




  if (abs(SUSY_Spart2_pdgId) == 1000022) nchi01++;
  else if (abs(SUSY_Spart2_pdgId) == 1000023) nchi02++;
  else if (abs(SUSY_Spart2_pdgId) == 1000025) nchi03++;
  else if (abs(SUSY_Spart2_pdgId) == 1000035) nchi04++;
  else if (SUSY_Spart2_pdgId == 1000024) nch1plus++;
  else if (SUSY_Spart2_pdgId == -1000024) nch1minus++;
  else if (SUSY_Spart2_pdgId == 1000037) nch2plus++;
  else if (SUSY_Spart2_pdgId == -1000037) nch2minus++;

  else if (SUSY_Spart2_pdgId == 1000021) ngluino++;
  else if ((abs(SUSY_Spart2_pdgId) > 1000000 && abs(SUSY_Spart2_pdgId) <= 1000004) || (abs(SUSY_Spart2_pdgId) > 2000000 && abs(SUSY_Spart2_pdgId) <= 2000004)) {
    if (SUSY_Spart2_pdgId > 0) nsquark++;
    else nantisquark++;
  }
  else if (SUSY_Spart2_pdgId == 1000005) nsbottom++;
  else if (SUSY_Spart2_pdgId == 1000006) nstop++;
  else if (SUSY_Spart2_pdgId == 2000005) nsbottom2++;
  else if (SUSY_Spart2_pdgId == 2000006) nstop2++;
  else if (SUSY_Spart2_pdgId == -1000005) nantisbottom++;
  else if (SUSY_Spart2_pdgId == -1000006) nantistop++;
  else if (SUSY_Spart2_pdgId == -2000005) nantisbottom2++;
  else if (SUSY_Spart2_pdgId == -2000006) nantistop2++;

  else if (SUSY_Spart2_pdgId == 2000011) nselecRminus++;
  else if (SUSY_Spart2_pdgId == -2000011) nselecRplus++;
  else if (SUSY_Spart2_pdgId == 1000011) nselecLminus++;
  else if (SUSY_Spart2_pdgId == -1000011) nselecLplus++;
  else if (abs(SUSY_Spart2_pdgId) == 1000012) nselnuL++;
  else if (SUSY_Spart2_pdgId == 2000013) nsmuonRminus++;
  else if (SUSY_Spart2_pdgId == -2000013) nsmuonRplus++;
  else if (SUSY_Spart2_pdgId == 1000013) nsmuonLminus++;
  else if (SUSY_Spart2_pdgId == -1000013) nsmuonLplus++;
  else if (abs(SUSY_Spart2_pdgId) == 1000014) nsmunuL++;
  else if (SUSY_Spart2_pdgId == 1000015) nstau1minus++;
  else if (SUSY_Spart2_pdgId == -1000015) nstau1plus++;
  else if (SUSY_Spart2_pdgId == 2000015) nstau2minus++;
  else if (SUSY_Spart2_pdgId == -2000015) nstau2plus++;
  else if (abs(SUSY_Spart2_pdgId) == 1000016) nstaunuL++;
  else nother++;


  ///Final classification
  // gluino/squark + X
  if (ngluino == 1 && (nsquark == 1 || nantisquark == 1)) return 1;
  else if (ngluino == 2) return 2;
  else if (nsquark == 2 || nantisquark == 2) return 3;
  else if (nsquark == 1 && nantisquark == 1) return 4;

  else if (nsbottom == 1 && nantisbottom == 1) return 51;
  else if (nsbottom2 == 1 && nantisbottom2 == 1) return 52;
  else if (nstop == 1 && nantistop == 1) return 61;
  else if (nstop2 == 1 && nantistop2 == 1) return 62;

  else if (ngluino == 1 && nchi01 == 1) return 71;
  else if (ngluino == 1 && nchi02 == 1) return 72;
  else if (ngluino == 1 && nchi03 == 1) return 73;
  else if (ngluino == 1 && nchi04 == 1) return 74;

  else if (ngluino == 1 && nch1plus == 1) return 75;
  else if (ngluino == 1 && nch2plus == 1) return 76;
  else if (ngluino == 1 && nch1minus == 1) return 77;
  else if (ngluino == 1 && nch2minus == 1) return 78;

  else if ((nsquark == 1 || nantisquark == 1) && nchi01 == 1) return 81;
  else if ((nsquark == 1 || nantisquark == 1) && nchi02 == 1) return 82;
  else if ((nsquark == 1 || nantisquark == 1) && nchi03 == 1) return 83;
  else if ((nsquark == 1 || nantisquark == 1) && nchi04 == 1) return 84;

  else if ((nsquark == 1 || nantisquark == 1) && nch1plus == 1) return 85;
  else if ((nsquark == 1 || nantisquark == 1) && nch2plus == 1) return 86;
  else if ((nsquark == 1 || nantisquark == 1) && nch1minus == 1) return 87;
  else if ((nsquark == 1 || nantisquark == 1) && nch2minus == 1) return 88;


  // Gaugino pair-production
  // chi^{0}_1 + X
  else if (nchi01 == 2) return 111;
  else if (nchi01 == 1 && nchi02 == 1) return 112;
  else if (nchi01 == 1 && nchi03 == 1) return 113;
  else if (nchi01 == 1 && nchi04 == 1) return 114;
  else if (nchi01 == 1 && nch1plus == 1) return 115;
  else if (nchi01 == 1 && nch2plus == 1) return 116;
  else if (nchi01 == 1 && nch1minus == 1) return 117;
  else if (nchi01 == 1 && nch2minus == 1) return 118;

  // chi^{0}_2 + X
  else if (nchi02 == 2) return 122;
  else if (nchi02 == 1 && nchi03 == 1) return 123;
  else if (nchi02 == 1 && nchi04 == 1) return 124;
  else if (nchi02 == 1 && nch1plus == 1) return 125;
  else if (nchi02 == 1 && nch2plus == 1) return 126;
  else if (nchi02 == 1 && nch1minus == 1) return 127;
  else if (nchi02 == 1 && nch2minus == 1) return 128;

  // chi^{0}_3 + X
  else if (nchi03 == 2) return 133;
  else if (nchi03 == 1 && nchi04 == 1) return 134;
  else if (nchi03 == 1 && nch1plus == 1) return 135;
  else if (nchi03 == 1 && nch2plus == 1) return 136;
  else if (nchi03 == 1 && nch1minus == 1) return 137;
  else if (nchi03 == 1 && nch2minus == 1) return 138;

  // chi^{0}_4 + X
  else if (nchi04 == 2) return 144;
  else if (nchi04 == 1 && nch1plus == 1) return 145;
  else if (nchi04 == 1 && nch2plus == 1) return 146;
  else if (nchi04 == 1 && nch1minus == 1) return 147;
  else if (nchi04 == 1 && nch2minus == 1) return 148;

  // chi^{+}_1/2 + chi^{-}_1/2
  else if (nch1plus == 1 && nch1minus == 1) return 157;
  else if (nch1plus == 1 && nch2minus == 1) return 158;

  else if (nch2plus == 1 && nch1minus == 1) return 167;
  else if (nch2plus == 1 && nch2minus == 1) return 168;

  // slepton
  else if (nselecLplus == 1 && nselecLminus == 1) return 201; // sElectronLPair
  else if (nselecRplus == 1 && nselecRminus == 1) return 202; // sElectronRPair
  else if (nselnuL == 2) return 203; // sElectron neutrino pair
  else if (nselecLplus == 1 && nselnuL == 1) return 204; // sElectron+ sNutrino
  else if (nselecLminus == 1 && nselnuL == 1) return 205; // sElectron- sNutrino
  else if (nstau1plus == 1 && nstau1minus == 1) return 206;
  else if (nstau2plus == 1 && nstau2minus == 1) return 207;
  else if ((nstau1plus == 1 || nstau1minus == 1) && (nstau2plus == 1 || nstau2minus == 1)) return 208;
  else if (nstaunuL == 2) return 209; // sTau neutrino pair
  else if (nstau1plus == 1 && nstaunuL == 1) return 210;
  else if (nstau1minus == 1 && nstaunuL == 1) return 211;
  else if (nstau2plus == 1 && nstaunuL == 1) return 212;
  else if (nstau2minus == 1 && nstaunuL == 1) return 213;

  else if (nsmuonLplus == 1 && nsmuonLminus == 1) return 216; // sMuonPair
  else if (nsmuonRplus == 1 && nsmuonRminus == 1) return 217; // sMuonPair
  else if (nsmunuL == 2) return 218; // sMuon neutrino pair
  else if (nsmuonLplus == 1 && nsmunuL == 1) return 219; // sMuon+ sNutrino
  else if (nsmuonLminus == 1 && nsmunuL == 1) return 220; // sMuon- sNutrino

  std::cerr << "ERROR. could not determine finalState for:" << std::endl;
  std::cerr << "  SUSY_Spart1_pdgId: " << SUSY_Spart1_pdgId << std::endl;
  std::cerr << "  SUSY_Spart2_pdgId: " << SUSY_Spart2_pdgId << std::endl;
  std::cerr << "Returning 0" << std::endl;

  return 0;
}
