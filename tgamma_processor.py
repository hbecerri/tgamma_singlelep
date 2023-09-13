# This file illustrates how to implement a processor, realizing the selection
# steps and outputting histograms and a cutflow with efficiencies.
# Here we create a very simplified version of the ttbar-to-dilep processor.
# One can run this processor using
# 'python3 -m pepper.runproc --debug example_processor.py example_config.json'
# Above command probably will need a little bit of time before all cuts are
# applied once. This is because a chunk of events are processed simultaneously.
# You change adjust the number of events in a chunk and thereby the memory
# usage by using the --chunksize parameter (the default value is 500000).

import pepper
import awkward as ak
from functools import partial
import math
import numpy as np
from coffea.nanoevents.methods import vector
from pepper import utils
from pepper.utils import pxpypz_from_ptetaphi
import pepper.top_reco as top_reco
import sys
np.set_printoptions(threshold=sys.maxsize)

def xyzE(pt,eta,phi,mass):

        pt = abs(pt)
        px = pt*np.cos(phi)
        py = pt*np.sin(phi)
        pz = pt*np.sinh(eta)
        E = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
        return px, py, pz, E
  
# All processors should inherit from pepper.ProcessorBasicPhysics
class Processor(pepper.ProcessorBasicPhysics):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigTTbarLL

    def __init__(self, config, eventdir):
        # Initialize the class, maybe overwrite some config variables and
        # load additional files if needed
        # Can set and modify configuration here as well
#        config["histogram_format"] = "root"
        # Need to call parent init to make histograms and such ready
        super().__init__(config, eventdir)
        config["histogram_format"] = "root"

        # It is not recommended to put anything as member variable into a
        # a Processor because the Processor instance is sent as raw bytes
        # between nodes when running on HTCondor.

    def process_selection(self, selector, dsname, is_mc, filler):
        # Implement the selection steps: add cuts, define objects and/or
        # compute event weights
        # Add a cut only allowing events according to the golden JSON
        # The good_lumimask method is specified in pepper.ProcessorBasicPhysics
        # It also requires a lumimask to be specified in config
        era = self.get_era(selector.data, is_mc)
        if not is_mc:
            selector.add_cut("Lumi", partial(
                self.good_lumimask, is_mc, dsname))

        #        if is_mc:
        #            selector.set_column("GenNeutrino", partial(self.build_genneutrino_column, is_mc))


        # Only allow events that pass triggers specified in config
        # This also takes into account a trigger order to avoid triggering
        # the same event if it's in two different data datasets.
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname, is_mc, self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))

        if is_mc and self.config["year"] in ("2016", "2017", "ul2016pre",
                                             "ul2016post", "ul2017"):
            selector.add_cut("L1 prefiring", self.add_l1_prefiring_weights)

        selector.add_cut("MET filters", partial(self.met_filters, is_mc))

        # Pick electrons satisfying our criterias
        selector.set_multiple_columns(self.pick_electrons)
        # Pick muons satisfying our criterias
        selector.set_multiple_columns(self.pick_muons)

        selector.set_cat("channel",{"ele", "muon"})
        selector.set_multiple_columns(self.lepton_categories)

        # Combine electron and muon to lepton
        selector.set_column("Lepton", partial(self.build_lepton_column, is_mc, selector.rng))

        # Only accept events that have one lepton and not have extra leptons
        selector.add_cut("oneLepton", self.one_good_lepton)

        # Pick photons satisfying our criterias
        selector.set_column("Photon", self.pick_medium_photons)

        # Only accept events that have at least one lepton
        selector.add_cut("atLeastOnePhoton",self.one_good_photon)
        selector.set_column("mlg",self.mass_lg) 

        # Pick Jets satisfying our criterias
        selector.set_column("Jet", partial(self.build_jet_column, is_mc))        
        selector.set_column("bJet", self.build_bjet_column)        
        selector.set_column("nbtag", self.num_btags)

        # Only accept events that have at least two jets and one bjet
        selector.add_cut("atLeast2jets",self.has_jets)

        # HEM issue cut
        if (self.config["hem_cut_if_ele"] or self.config["hem_cut_if_muon"]
                or self.config["hem_cut_if_jet"]):
            selector.add_cut("HEM cut", self.hem_cut)

        # Build MET column
        smear_met = "smear_met" in self.config and self.config["smear_met"]
        variation = self.get_jetmet_nominal_arg()
        selector.set_column("OrigJet", selector.data["Jet"])
        selector.set_column(
            "MET", partial(self.build_met_column, is_mc, variation.junc,
                           variation.jer if smear_met else None, selector.rng,
                           era, variation=variation.met))
        # Only accept events with MET pt more than 20 GeV
        selector.add_cut("Req MET", self.met_requirement)
        selector.add_cut("atLeast1bjet",partial(self.btag_cut, is_mc))

        # Build different categories according to the number of jets
#        selector.set_cat("jet_btag", {"j2+_b0", "j2_b1", "j3+_b1", "j2_b2+", "j3+_b2+"})
#        selector.set_multiple_columns(self.btag_categories)

        #ttbar semileptnic reconstruction
        selector.set_column("Neutrino1",self.neutrino_reco)  
        #selector.set_column("Top", self.ttbar_reco)
#        selector.set_multiple_columns(self.ttbar_reco)
        topVars = top_reco.topreco(selector.data)

        selector.set_column("tophad_m", topVars["mtophad"])           
        selector.set_column("tophad_pt", topVars["pttophad"])
        selector.set_column("tophad_eta", topVars["etatophad"])
        selector.set_column("tophad_phi", topVars["phitophad"])

        selector.set_column("toplep_m",   topVars["mtoplep"])
        selector.set_column("toplep_pt",  topVars["pttoplep"])
        selector.set_column("toplep_eta", topVars["etatoplep"])
        selector.set_column("toplep_phi", topVars["phitoplep"])

        chi2_flat =  [y for x in topVars["chisquare"] for y in x]
        chi2 = [i for i in chi2_flat if i is not None] 
        selector.set_column("chi2",chi2)       

        selector.set_column("GenTopPos", self.gentoppos)
        selector.set_column("GenTopNeg", self.gentopneg)

        charge = ak.concatenate(
            [selector.data["Electron"].charge, selector.data["Muon"].charge], axis=1)
        selector.set_column("Lepton_charge", charge)
 
  
    def pick_electrons(self, data):
        ele = data["Electron"]

        # We do not want electrons that are between the barrel and the end cap
        # For this, we need the eta of the electron with respect to its
        # supercluster
        sc_eta_abs = abs(ele.eta + ele.deltaEtaSC)
        is_in_transreg = (1.444 < sc_eta_abs) & (sc_eta_abs < 1.566)

        # Electron ID, as an example we use the MVA one here
#        has_id = ele.mvaFall17V2Iso_WP90
        has_id = ele.cutBased >= 3

        # Finally combine all the requirements
        is_good = (
            has_id
            & (~is_in_transreg)
            & (self.config["ele_eta_min"] < ele.eta)
            & (ele.eta < self.config["ele_eta_max"])
            & (self.config["good_ele_pt_min"] < ele.pt))
 
        veto_id = ele.cutBased >=1

        is_veto = (
                veto_id
              & (~is_in_transreg)
              & (self.config["ele_eta_min"] < ele.eta)
              & (ele.eta < self.config["ele_eta_max"])
              & (self.config["veto_ele_pt_min"] < ele.pt))

        # Return all electrons with are deemed to be good
        return {"Electron": ele[is_good], "VetoEle": ele[is_veto]}

    def pick_muons(self, data):
        muon = data["Muon"]
        etacuts = (self.config["muon_eta_min"] < muon.eta) & (muon.eta < self.config["muon_eta_max"])

        good_id = muon.tightId
        good_iso = muon.pfIsoId > 3
        is_good = (
            good_id
            & good_iso
            & etacuts
            & (self.config["good_muon_pt_min"] < muon.pt))
        
        veto_id = muon.looseId
        veto_iso = muon.pfIsoId > 1
        is_veto = (
            veto_id
            & veto_iso
            & etacuts
            & (self.config["veto_muon_pt_min"] < muon.pt))

        return {"Muon": muon[is_good], "VetoMuon": muon[is_veto]}

    def pick_medium_photons(self, data):
        photon = data["Photon"]
        leptons = data["Lepton"]
        has_id = photon.cutBased>1
        pass_psv = (photon.pixelSeed==False)

        eta_abs = abs(photon.eta)
        is_in_transreg = (1.4442 < eta_abs) & (eta_abs < 1.566)
    
        etacuts = (abs(photon["eta"])<2.5)
        ptcuts = (photon.pt>15)

        has_lepton_close = ak.any(
            photon.metric_table(leptons) < 0.4, axis=2)

        is_good = (
                has_id
	      & (~has_lepton_close)
              & pass_psv
              & etacuts
              & (~is_in_transreg)
              & ptcuts)

        return photon[is_good]

    def good_jet(self, data):
        """Apply some basic jet quality cuts."""
        jets = data["Jet"]
        leptons = data["Lepton"]
        photons = data["Photon"]

        j_id, j_puId, lep_dist, pho_dist, eta_min, eta_max, pt_min = self.config[[
            "good_jet_id", "good_jet_puId", "good_jet_lepton_distance", "good_jet_photon_distance",
            "good_jet_eta_min", "good_jet_eta_max", "good_jet_pt_min"]]

        if j_id == "skip":
            has_id = True
        elif j_id == "cut:loose":
            has_id = jets.isLoose
            # Always False in 2017 and 2018
        elif j_id == "cut:tight":
            has_id = jets.isTight
        elif j_id == "cut:tightlepveto":
            has_id = jets.isTightLeptonVeto
        else:
            raise pepper.config.ConfigError(
                    "Invalid good_jet_id: {}".format(j_id))

        if j_puId == "skip":
            has_puId = True
        elif j_puId == "cut:loose":
            has_puId = ak.values_astype(jets["puId"] & 0b100, bool)
        elif j_puId == "cut:medium":
            has_puId = ak.values_astype(jets["puId"] & 0b10, bool)
        elif j_puId == "cut:tight":
            has_puId = ak.values_astype(jets["puId"] & 0b1, bool)
        else:
            raise pepper.config.ConfigError(
                    "Invalid good_jet_id: {}".format(j_puId))

        # Only apply PUID if pT < 50 GeV
        has_puId = has_puId | (jets.pt >= 50)

        j_pt = jets.pt
        if "jetfac" in ak.fields(data):
            jets["pt"] = jets["pt"] * data["jetfac"][is_good_jet]
            jets["mass"] = jets["mass"] * data["jetfac"][is_good_jet]
            jets = jets[ak.argsort(jets["pt"], ascending=False)]

        has_lepton_close = ak.any(
            jets.metric_table(leptons) < lep_dist, axis=2)
        has_photon_close = ak.any(
            jets.metric_table(photons) < pho_dist, axis=2)     

        return (has_id & has_puId
                & (~has_lepton_close)
                & (~has_photon_close)
                & (eta_min < jets.eta)
                & (jets.eta < eta_max)
                & (pt_min < j_pt))
 
    def build_bjet_column(self,data):

        jets = data["Jet"]
        bjets = jets[data["Jet"].btagged]

        return bjets

    def one_good_muon(self, data):
        return (ak.num(data["Muon"]) == 1) & (ak.num(data["VetoMuon"]) == 1)
    
    def one_good_ele(self, data):
        return (ak.num(data["Electron"]) == 1) & (ak.num(data["VetoEle"]) == 1)
    
    def one_good_lepton(self, data):
        return ((ak.num(data["Muon"]) == 1) & (ak.num(data["VetoMuon"]) == 1)) | ((ak.num(data["Electron"]) == 1) & (ak.num(data["VetoEle"]) == 1))

    def one_good_photon(self,data):
        return ak.num(data["Photon"])>0
   
    def num_btags(self, data):
        return ak.num(data['bJet'])

    def lepton_categories(self,data):
        
        cat = {}
        nele = ak.num(data['VetoEle'])
        nmuon = ak.num(data['VetoMuon'])

        cat['ele'] = (nele==1) & (nmuon==0)
        cat['muon'] = (nele==0) & (nmuon==1)

        return cat

    def btag_categories(self,data):
        
        cats = {}
        
        num_btagged = data["nbtag"]
        njet = ak.num(data["Jet"])

        cats["j2+_b0"] = (num_btagged == 0) & (njet == 2)
        cats["j2_b1"] = (num_btagged == 1) & (njet == 2)
        cats["j3+_b1"] = (num_btagged == 1) & (njet > 2)
        cats["j2_b2+"] = (num_btagged >= 2) & (njet == 2)
        cats["j3+_b2+"] = (num_btagged >= 2) & (njet > 2)


        return cats

    def met_requirement(self, data):
        met = data["MET"].pt
        return met > self.config["met_min_met"]

    def mass_lg(self, data):
        """Return invariant mass of lepton plus photon"""
        return (data["Lepton"][:, 0] + data["Photon"][:, 0]).mass

    def opposite_sign_lepton_pair(self, data):
        # At this point we only have events with exactly two leptons, but now
        # we want only events where they have opposite charge

        # First concatenate the charge of our electron(s) and our muon(s)
        # into one array
        charge = ak.concatenate(
            [data["Electron"].charge, data["Muon"].charge], axis=1)

        # Now in this array we can simply compare the first and the second
        # element. Note that this is done on axis 1, axis 0 is always used for
        # event indexing, e.g. you would compare charges from event 0 and 1 if
        # you do charge[0] != charge[1]
        return charge[:, 0] != charge[:, 1]

    def lepton_pair(self, data):
        # We only want events with excatly two leptons, thus look at our
        # electron and muon counts and pick events accordingly
        return ak.num(data["Electron"]) + ak.num(data["Muon"]) == 2
   
    def neutrino_reco(self,data):

        met = data["MET"]  
        lepton = data["Lepton"][:,0] 
        pxl, pyl, pzl = pxpypz_from_ptetaphi(lepton.pt, lepton.eta, lepton.phi)
        pxnu, pynu, pznu = pxpypz_from_ptetaphi(met.pt, lepton.eta, met.phi)        
     
        Enu = pxnu**2 + pynu**2
 
        mWT = np.sqrt((lepton.pt + met.pt)**2 - (pxl + pxnu)**2 -
                  (pyl + pynu)**2)        
 
        mW = ak.Array(np.full(len(mWT), 80.4, dtype=float))
        mask_mWT_GT_mW = mWT > mW
   
        dummy_mask = ak.full_like(mask_mWT_GT_mW, True)
        mask_mWT_GT_mW = ak.singletons(ak.mask(mask_mWT_GT_mW, dummy_mask))
        pxnu = ak.singletons(ak.mask(pxnu, dummy_mask))
        pynu = ak.singletons(ak.mask(pynu, dummy_mask))
        pznu = ak.singletons(ak.mask(pznu, dummy_mask))
        ptW = ak.singletons(ak.mask(lepton.pt, dummy_mask))
        pxl = ak.singletons(ak.mask(pxl, dummy_mask))
        pyl = ak.singletons(ak.mask(pyl, dummy_mask))
        pzl = ak.singletons(ak.mask(pzl, dummy_mask))   
 
        k = met.pt * lepton.pt - pxnu * pxl - pynu * pyl

        k = ak.fill_none(ak.pad_none(k[k > 0.0001], 1), 0.0001)
        scf = 0.5 * (mW * mW) / k
        pxnu = ak.concatenate([pxnu[mask_mWT_GT_mW] * scf[mask_mWT_GT_mW],
                           pxnu[~mask_mWT_GT_mW]], axis=1)
        pynu = ak.concatenate([pynu[mask_mWT_GT_mW] * scf[mask_mWT_GT_mW],
                           pynu[~mask_mWT_GT_mW]], axis=1)
        Etnu = np.sqrt(pxnu**2 + pynu**2)

        Lambda = (mW**2)/2. + pxl * pxnu + pyl * pynu
        discr = ((Lambda * pzl)**2)/(ptW**4) - (((lepton.energy * Etnu)**2) -
                                            (Lambda**2))/(ptW**2) 

        mask_posDiscr = discr > 0
  
        sol = Lambda * pzl/(lepton.pt * lepton.pt)
        pxnu_neg = pxnu[~mask_posDiscr]
        pynu_neg = pynu[~mask_posDiscr]
        pznu_neg = sol[~mask_posDiscr]
        Enu_neg = np.sqrt(pxnu_neg**2 + pynu_neg**2 + pznu_neg**2)
 
        sol1 = (
            Lambda[mask_posDiscr] * pzl[mask_posDiscr]/(ptW[mask_posDiscr]**2) +
            np.sqrt(discr[mask_posDiscr])
        )
        sol2 = (
            Lambda[mask_posDiscr] * pzl[mask_posDiscr]/(ptW[mask_posDiscr]**2) -
            np.sqrt(discr[mask_posDiscr])
        )

        pxnu_pos = ak.concatenate([pxnu[mask_posDiscr],
                              pxnu[mask_posDiscr]], axis=1)
        pynu_pos = ak.concatenate([pynu[mask_posDiscr],
                               pynu[mask_posDiscr]], axis=1)
        pznu_pos = ak.concatenate([sol1, sol2], axis=1)
        Enu_pos = np.sqrt(pxnu_pos**2 + pynu_pos**2 + pznu_pos**2)

        pxnu = ak.concatenate([pxnu_pos, pxnu_neg], axis=1)
        pynu = ak.concatenate([pynu_pos, pynu_neg], axis=1)
        pznu = ak.concatenate([pznu_pos, pznu_neg], axis=1)
        Enu = ak.concatenate([Enu_pos, Enu_neg], axis=1)
  
        mnu2 = np.minimum(10e20, np.maximum((Enu*Enu - (pxnu*pxnu + pynu*pynu + pznu*pznu)), 0))
        mnu = np.sqrt(mnu2)
        ptnu = np.hypot(pxnu, pynu)
        phinu = np.arctan2(pynu, pxnu)
        etanu = np.arcsinh(pznu/ptnu)
        neutrinos = ak.zip({"pt": ptnu, "eta": etanu, "phi": phinu, "mass": mnu},
                       with_name="PtEtaPhiMLorentzVector")
        return neutrinos

    def gentoppos(self, data):
        part = data["GenPart"]
        part = part[~ak.is_none(part.parent, axis=1)]
        part = part[part.hasFlags("isLastCopy")]
        part = part[(part.pdgId) == 6]
        return part

    def gentopneg(self, data):
        part = data["GenPart"]
        part = part[~ak.is_none(part.parent, axis=1)]
        part = part[part.hasFlags("isLastCopy")]
        part = part[(part.pdgId) == -6]
        return part
