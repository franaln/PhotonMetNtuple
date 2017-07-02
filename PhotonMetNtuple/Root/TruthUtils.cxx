#include <PhotonMetNtuple/TruthUtils.h>

TruthParticle::TruthParticle(const xAOD::TruthParticle *ptcl) :
  //TLorentzVector(ptcl->px()*0.001, ptcl->py()*0.001, ptcl->pz()*0.001, ptcl->e()*0.001),
  TLorentzVector(),
  index(-1),
  good(true)
{
  SetPtEtaPhiE(ptcl->pt()*0.001, ptcl->eta(), ptcl->phi(), ptcl->e()*0.001);
}

TruthParticle::TruthParticle(const xAOD::Jet *ptcl, Bool_t isbjet) :
  TLorentzVector(), //ptcl->px()*0.001, ptcl->py()*0.001, ptcl->pz()*0.001, ptcl->e()*0.001),
  index(-1),
  good(true),
  isbjet(isbjet)
{
  SetPtEtaPhiE(ptcl->pt()*0.001, ptcl->eta(), ptcl->phi(), ptcl->e()*0.001);
}


Double_t TruthUtils::GetDeltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t  phi = fabs(phi1 - phi2);
  if(phi <= TMath::Pi())  return phi;
  else                    return (2 * TMath::Pi() - phi);
}

Double_t TruthUtils::GetDeltaR(const TruthParticle &p, const TruthParticle &q)
{
  float dphi = GetDeltaPhi(p.Phi(), q.Phi());
  float deta = p.Eta() - q.Eta();
  return TMath::Hypot(deta, dphi);
}

Bool_t TruthUtils::OverlapsOthers(TruthParticle &p, std::vector<TruthParticle> &others, Double_t deltaR_cut) 
{
  std::vector<TruthParticle>::const_iterator p_it;
  for (p_it=others.begin(); p_it!=others.end(); ++p_it) {
    if ((*p_it).good && GetDeltaR(p, (*p_it)) < deltaR_cut) 
      return true;
  }
  
  return false;
}

void TruthUtils::CleanBads(std::vector<TruthParticle> &vector)
{
  std::vector<TruthParticle>::iterator p_it;
  
  for (p_it=vector.begin(); p_it!=vector.end();) {
    if (!(*p_it).good)
      p_it = vector.erase(p_it);
    else
      ++p_it;
  }
}




/// print details about the truth particle to the screen
// void print_truth_ptcl(const xAOD::TruthParticle *ptcl, TString comment, 
//                       int childDepth, int parentDepth, int currentDepth) 
// {
//   // indentation in case we print decay chain. Three spaces per level
//   TString indent(Form(Form("%%%ds",3*currentDepth),""));
//   if (ptcl == NULL) { 
//     printf("%sNULL\n",indent.Data()); 
//     return; 
//   }
//   printf("%sTruth part. ID:%5d, status: %2d, %s  %s\n",
//          indent.Data(), ptcl->pdgId(), ptcl->status(), fourVecAsText(ptcl).Data(), comment.Data());

//   if (childDepth > 0 || parentDepth > 0) {
//     int npar = ptcl->nParents();
//     int nchild = ptcl->nChildren();
//     printf("%s-> %d parent and %d children\n", indent.Data(), npar, nchild);

//     if (parentDepth > 0) 
//       for (int ip=0; ip<npar; ++ip) print_truth_ptcl(ptcl->parent(ip), Form("parent %d of ",ip+1)+comment,
//                                                      childDepth-1, parentDepth-1, currentDepth+1);
//     if (childDepth > 0)
//       for (int ic=0; ic<nchild; ++ic) print_truth_ptcl(ptcl->child(ic), Form("child %d of ",ic+1)+comment,
//                                                        childDepth-1, parentDepth-1, currentDepth+1);
//   }
// }
  
// Bool_t is_stable(const xAOD::TruthParticle *ptcl) 
// {
//   return (ptcl->status() % 1000) == 1 && ptcl->barcode() < 200000;
// }

Bool_t is_stable(const xAOD::TruthParticle* p){ 

  if (p->barcode() > 200000) 
    return false;

  if (p->pdgId() == 21 && p->e() == 0) 
    return false;

  if (p->status() % 1000 == 1) 
    return true;

  if (p->hasDecayVtx()) {
    const xAOD::TruthVertex* dvtx = p->decayVtx();
    if (dvtx){
      if (p->status() == 2 && dvtx->barcode() < -200000) 
        return true;
    }
  }

  if (p->status()% 1000 == 2 && p->status() > 1000 ) 
    return true;

  return false;
}
  
// // Return true if from hadron
// Bool_t is_from_hadron(const xAOD::TruthParticle *ptcl) 
// {

//   int id = ptcl->pdgId();
    
//   // if the particle is a hadron, return true
//   if (ptcl->isHadron(id)) return true;
    
//   // if there are no parents, not from hadron
//   if (ptcl->nParents()==0) return false;
  
//   const xAOD::TruthParticle *parent = ptcl->parent(0);
//   int parent_id = parent->pdgId();
//   if (parent->isHadron(parent_id)) return true; // from hadron!
//   if (parent_id == 15 || parent_id == id) return is_from_hadron(parent);
  
//   // if we get here, all is good
//   return false;
// }
  
Int_t get_photon_type(const xAOD::TruthParticle *ptcl)
{
  int pdg_id = ptcl->absPdgId(); //reader->ph_truth_type->at(d3pd_index));

  int parent = -1;
  if (ptcl->nParents() > 0)
    parent = ptcl->parent(0)->pdgId();

  //std::cout << "pdg_id: " << pdg_id << " parent: " << parent << std::endl;
  if (pdg_id == 22) {  // real photons
    if (parent < 0 || parent == 22) {
      // if (reader->ph_truth_isPhotonFromHardProc->at(d3pd_index)) return 1;     //from hard proc
      // else if(reader->ph_truth_isBrem->at(d3pd_index)) return 2;               //from brems
      // else 
      return 0;                                                   //signal photons (at least for GGM samples)
    }
    else if (parent == 24 || parent == 23) return 3;  //from W/Z decay
    else if (parent > 10 && parent < 17) return 4;  //from lepton/neutrino (i.e. from W/Z but...)
    else if (parent < 7)   return 5;  //radiated off from quarks
    else                   return 6;  //from hadron decay
  }
  else if (pdg_id == 11) { // electrons faking photons
    if (parent == 24)      return 7; //from W decay
    else if (parent == 23) return 8; //from Z decay
    else                   return 9; //from other sources
  }
  else if (pdg_id != 0) { // other particles faking photons
     return 10;
  }
  else { // no truth matching found. Impossible to tell
    return 11;
  }

  return -1;
}

// // adds up 4-vectors of all stable particles:
// //  should always give E=m=sqrt(s) \vec{p}=\vec{0} !
// // unless partciels are missing from the file
// TLorentzVector getStableParticle4VectorSum(const xAOD::TruthParticleContainer *truthPtcls) {
//   TLorentzVector sum;
//   for ( const xAOD::TruthParticle *ptcl : *truthPtcls)
//     if (isStable(ptcl)) sum += ptcl->p4();
//   return sum;
// }
  
// bool isFromHiggs(const xAOD::TruthParticle *ptcl) {
//   if (MC::PID::isHiggs(ptcl->pdgId())) return true;
//   if (ptcl->parent()==nullptr) return false;
//   return isFromHiggs(ptcl->parent());
// }

// bool isFromBhadron(const xAOD::TruthParticle *ptcl) {
//   if (MC::PID::isBottomHadron(ptcl->pdgId())) return true;
//   if (ptcl->parent()==nullptr) return false;
//   return isFromBhadron(ptcl->parent());
// }

// bool isGoodTruthPhoton(const xAOD::TruthParticle *ptcl) {
//   return isStable(ptcl) && MC::PID::isPhoton(ptcl->pdgId()) && notFromHadron(ptcl);
// }

// bool isGoodTruthElectron(const xAOD::TruthParticle *ptcl) {
//   return isStable(ptcl) && MC::PID::isElectron(ptcl->pdgId()) && notFromHadron(ptcl);
// }

// bool isGoodTruthMuon(const xAOD::TruthParticle *ptcl) {
//   return isStable(ptcl) && MC::PID::isMuon(ptcl->pdgId()) && notFromHadron(ptcl);
// }

// std::vector<const xAOD::TruthParticle*> getGoodTruthPhotonsOld(const xAOD::TruthParticleContainer *truthPtcls) {
//   ::std::vector<const xAOD::TruthParticle*> truthPhotons;
//   for (const xAOD::TruthParticle *ptcl : *truthPtcls)
//     if (isGoodTruthPhoton(ptcl)) truthPhotons.push_back(ptcl);
//   return truthPhotons;
// }

// //! /brief returns all stable electrons that do not originate from hadrons
// TruthPtcls getGoodTruthPhotons( const xAOD::TruthParticleContainer *truthPtcls ) {
//   TruthPtcls ys(SG::VIEW_ELEMENTS);
//   for (auto ptcl : *truthPtcls) if ( isGoodTruthPhoton(ptcl) ) ys.push_back(ptcl);
//   return ys;
// }

// //! /brief returns all stable electrons that do not originate from hadrons
// TruthPtcls getGoodTruthElectrons( const xAOD::TruthParticleContainer * truthPtcls ) {
//   TruthPtcls es(SG::VIEW_ELEMENTS);
//   for (auto ptcl : *truthPtcls) if (isGoodTruthElectron(ptcl)) es.push_back(ptcl);
//   return es;
// }

// //! /brief returns all stable electrons that do not originate from hadrons
// TruthPtcls getGoodTruthMuons( const xAOD::TruthParticleContainer * truthPtcls ) {
//   TruthPtcls mus(SG::VIEW_ELEMENTS);
//   for (auto ptcl : *truthPtcls) if (isGoodTruthMuon(ptcl)) mus.push_back(ptcl);
//   return mus;
// }

// TruthPtcls getHadronsAndTheirDecay( const xAOD::TruthParticleContainer * truthPtcls ) {
//   TruthPtcls hadrons(SG::VIEW_ELEMENTS);
//   for (auto ptcl : *truthPtcls) {
//     if (isGoodTruthPhoton(ptcl)) continue;
//     if (isGoodTruthElectron(ptcl)) continue;
//     if (isGoodTruthMuon(ptcl)) continue;
//     hadrons.push_back(ptcl);
//   }
//   return hadrons;
// }


// TruthParticleStruct identifyTruthParticles(const xAOD::TruthParticleContainer *truthPtcls,
//                                            const xAOD::JetContainer *truthJets,
//                                            double jet_pTcut) {
//   TruthParticleStruct tp;
//   TruthPtcls ys  = getGoodTruthPhotons(truthPtcls);
//   TruthPtcls es  = getGoodTruthElectrons(truthPtcls);
//   TruthPtcls mus = getGoodTruthMuons(truthPtcls);
//   TruthPtcls hads =  getHadronsAndTheirDecay(truthPtcls);
  
//   // TO-DO
//   // Dressing should happen here !
//   tp.electrons = es; tp.muons = mus;
//   // some ys should probably go to hads here
//   tp.photons = ys; tp.hadrons = hads;
  
//   tp.photonsFromHiggs = getPhotonsFromHiggs(truthPtcls);
//   // this one might be slow ... ?
//   tp.HiggsDecay = getHiggsDecayProducts(truthPtcls);
  
//   tp.Bhadrons    = getBHadrons(truthPtcls);
//   tp.Dhadrons    = getDHadrons(truthPtcls);
//   tp.muonsFromBs = getMuonsFromBs(truthPtcls);
  
//   TruthJets jets(SG::VIEW_ELEMENTS);
//   TruthJets bjets(SG::VIEW_ELEMENTS);
//   TruthJets cjets(SG::VIEW_ELEMENTS);
//   TruthJets lightJets(SG::VIEW_ELEMENTS);
  
//   // Here applying a 5 GeV cut for the jet labelling
//   TruthPtcls Bs    = getBHadrons(truthPtcls,5.0*HG::GeV);
//   TruthPtcls Ds    = getDHadrons(truthPtcls,5.0*HG::GeV);
//   for (const xAOD::Jet *tjet : *truthJets) {
//     // apply a pT cut, if requested
//     if ( jet_pTcut>0 && tjet->pt()<jet_pTcut ) continue;
    
//     // ignore jets overlapping with good photons, electrons or muons
//     if (HG::minDRrap(tjet,ys)<0.4) continue;
//     if (HG::minDRrap(tjet,es)<0.4) continue; // <<== WZ jets should not do this
//     // if (HG::minDRrap(tjet,mus)<0.4) continue; ??
//     jets.push_back(tjet);
//     // classify all jets into b, c or light
//     if      (HG::minDRrap(tjet,Bs)<0.4) bjets.push_back(tjet);
//     else if (HG::minDRrap(tjet,Ds)<0.4) cjets.push_back(tjet);
//     else lightJets.push_back(tjet);
//   }
  
//   // later: further split light jets into: LQ, gluon, unmatched
//   tp.jets  = jets;
//   tp.bJets = bjets;
//   tp.cJets = cjets;
//   tp.lightJets = lightJets;
//   return tp;
// }
