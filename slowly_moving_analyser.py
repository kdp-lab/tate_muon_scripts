import pyLCIO
import ROOT
import glob
import json
from math import *
import numpy as np


# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = 1000 # Set to -1 to run over all events

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
# fnames = glob.glob("/local/d1/lrozanov/mucoll-tutorial-2023/reco_Hbb_bib/134_0.1_reco_bib.slcio") 
fnames = glob.glob("/local/d1/mu+mu-/reco_v3/100_150_0_timingchange_32ns64ns/4500_10_reco.slcio") 
print("Found %i files."%len(fnames))

#Make lists

#MCPs
mcp_pt = []
mcp_eta = []
mcp_phi = []
mcp_theta = []
pdgid = []
status = []
prod_vertex = []
prod_time = []
prod_traveldist = []
id = []

#DAUGHTERS
daughter_pdgid = []

# TRACKS
d0_res = []
z0_res = []
pt_res = []

track_pt = []
track_eta = []
track_theta = []
track_phi = []
nhits = []
ndf = []
chi2 = []
           
#LC RELATION                   
#confirmed staus
    #mcp values
mcp_stau_pt = []
mcp_stau_eta = []
mcp_stau_phi = []
mcp_stau_theta = []
mcp_stau_traveldist = []
st_pt_match = []
st_eta_match = []
st_phi_match = []
st_theta_match = []
    #track values
st_track_pt = []
st_track_theta = []
st_track_phi = []
st_ndf = []
st_chi2 = []
st_d0 = []
st_z0 = []
    #hit values
st_nhits = []
st_pt_res = []
st_pix_nhits = []
st_inn_nhits = []
st_out_nhits = []
st_x = []
st_y = []
st_z = []
#st_hit_pdg = []
st_time = []
st_corrected_time = []
st_hit_detector = []
st_hit_layer = []
st_hit_side = []

# HITS
x = []
y = []
z = []
hit_pdgid = []
time = []
corrected_time = []
hit_layer = []
hit_detector = []
hit_side = []
pixel_nhits = []
inner_nhits = []
outer_nhits = []

sim_VB_x, sim_VB_y, sim_VB_z, sim_VB_time, sim_VB_pdg, sim_VB_mcpid, sim_VB_layer = [], [], [], [], [], [], []
sim_VE_x, sim_VE_y, sim_VE_z, sim_VE_time, sim_VE_pdg, sim_VE_mcpid, sim_VE_layer = [], [], [], [], [], [], []
sim_IB_x, sim_IB_y, sim_IB_z, sim_IB_time, sim_IB_pdg, sim_IB_mcpid, sim_IB_layer = [], [], [], [], [], [], []
sim_IE_x, sim_IE_y, sim_IE_z, sim_IE_time, sim_IE_pdg, sim_IE_mcpid, sim_IE_layer = [], [], [], [], [], [], []
sim_OB_x, sim_OB_y, sim_OB_z, sim_OB_time, sim_OB_pdg, sim_OB_mcpid, sim_OB_layer = [], [], [], [], [], [], []
sim_OE_x, sim_OE_y, sim_OE_z, sim_OE_time, sim_OE_pdg, sim_OE_mcpid, sim_OE_layer = [], [], [], [], [], [], []

reco_VB_x, reco_VB_y, reco_VB_z, reco_VB_time, reco_VB_layer = [], [], [], [], []
reco_VE_x, reco_VE_y, reco_VE_z, reco_VE_time, reco_VE_layer = [], [], [], [], []
reco_IB_x, reco_IB_y, reco_IB_z, reco_IB_time, reco_IB_layer = [], [], [], [], []
reco_IE_x, reco_IE_y, reco_IE_z, reco_IE_time, reco_IE_layer = [], [], [], [], []
reco_OB_x, reco_OB_y, reco_OB_z, reco_OB_time, reco_OB_layer = [], [], [], [], []
reco_OE_x, reco_OE_y, reco_OE_z, reco_OE_time, reco_OE_layer = [], [], [], [], []

speedoflight = 299792458/1000000  # mm/ns

# Make counter variables
n_mcp_stau = 0
no_inner_hits = 0
i = 0
num_matched_tracks = 0
num_dupes = 0
num_fake_tracks = 0
# reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
# reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks", "SiTracks_Refitted", "MCParticle_SiTracks", "MCParticle_SiTracks_Refitted", "IBTrackerHits", "IETrackerHits", "OBTrackerHits", "OETrackerHits", "VBTrackerHits", "VETrackerHits"])
# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
for f in fnames:
    if max_events > 0 and i >= max_events: break
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)
    for ievt,event in enumerate(reader): 
        if max_events > 0 and i >= max_events: break
        if i%1 == 0: print("Processing event %i."%i)
        print(" ")

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        trackCollection = event.getCollection("SiTracks")
        relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
        relation = pyLCIO.UTIL.LCRelationNavigator(relationCollection)

        hit_collections = []
        IBTrackerHits = event.getCollection('ITBarrelHits')
        hit_collections.append(IBTrackerHits)
        IETrackerHits = event.getCollection('ITEndcapHits')
        hit_collections.append(IETrackerHits)
        OBTrackerHits = event.getCollection('OTBarrelHits')
        hit_collections.append(OBTrackerHits)
        OETrackerHits = event.getCollection('OTEndcapHits')
        hit_collections.append(OETrackerHits)
        VBTrackerHits = event.getCollection('VXDBarrelHits')
        hit_collections.append(VBTrackerHits)
        VETrackerHits = event.getCollection('VXDEndcapHits')
        hit_collections.append(VETrackerHits)
        
        # Relations
        VBrelationCollection = event.getCollection('VXDBarrelHitsRelations')
        VBrelation = pyLCIO.UTIL.LCRelationNavigator(VBrelationCollection)
        
        VErelationCollection = event.getCollection('VXDEndcapHitsRelations')
        VErelation = pyLCIO.UTIL.LCRelationNavigator(VErelationCollection)
        
        IBrelationCollection = event.getCollection('ITBarrelHitsRelations')
        IBrelation = pyLCIO.UTIL.LCRelationNavigator(IBrelationCollection)
        
        IErelationCollection = event.getCollection('ITEndcapHitsRelations')
        IErelation = pyLCIO.UTIL.LCRelationNavigator(IErelationCollection)
        
        OBrelationCollection = event.getCollection('OTBarrelHitsRelations')
        OBrelation = pyLCIO.UTIL.LCRelationNavigator(OBrelationCollection)
        
        OErelationCollection = event.getCollection('OTEndcapHitsRelations')
        OErelation = pyLCIO.UTIL.LCRelationNavigator(OErelationCollection)

        # MCPs
        imcp_pt = []
        imcp_eta = []
        imcp_phi = []
        imcp_theta = []
        ipdgid = []
        istatus = []
        iprod_vertex = []
        iprod_time = []
        iprod_traveldist = []
        iid = []

        #DAUGHTERS
        idaughter_pdgid = []

        # TRACKS
        id0_res = []
        iz0_res = []
        ipt_res = []

        itrack_pt = []
        itrack_eta = []
        itrack_theta = []
        itrack_phi = []
        inhits = []
        indf = []
        ichi2 = []
           
        #LC RELATION            
        #confirmed staus (hopefully)
            #mcp values
        imcp_stau_pt = []
        imcp_stau_eta = []
        imcp_stau_phi = []
        imcp_stau_theta = []
        imcp_stau_traveldist = []

        ist_pt_match = []    
        ist_eta_match = []
        ist_phi_match = []
        ist_theta_match = []
            #track values
        ist_track_pt = []
        ist_track_theta = []
        ist_track_phi = []
        ist_ndf = []
        ist_chi2 = []
        ist_d0 = []
        ist_z0 = []
            #hit values
        ist_nhits = []
        ist_pt_res = []
        ist_x = []
        ist_y = []
        ist_z = []
        #ist_hit_pdg = []
        ist_time = []
        ist_corrected_time = []
        ist_hit_detector = []
        ist_hit_layer = []
        ist_hit_side = []

        # HITS
        ix = []
        iy = []
        iz = []
        ihit_pdg = []
        itime = []
        icorrected_time = []
        ihit_layer = []
        ihit_detector = []
        ihit_side = []
       

        
        
        # Sim Hits
        isim_VB_x, isim_VB_y, isim_VB_z, isim_VB_time, isim_VB_pdg, isim_VB_mcpid, isim_VB_layer = [], [], [], [], [], [], []
        isim_VE_x, isim_VE_y, isim_VE_z, isim_VE_time, isim_VE_pdg, isim_VE_mcpid, isim_VE_layer = [], [], [], [], [], [], []
        isim_IB_x, isim_IB_y, isim_IB_z, isim_IB_time, isim_IB_pdg, isim_IB_mcpid, isim_IB_layer = [], [], [], [], [], [], []
        isim_IE_x, isim_IE_y, isim_IE_z, isim_IE_time, isim_IE_pdg, isim_IE_mcpid, isim_IE_layer = [], [], [], [], [], [], []
        isim_OB_x, isim_OB_y, isim_OB_z, isim_OB_time, isim_OB_pdg, isim_OB_mcpid, isim_OB_layer = [], [], [], [], [], [], []
        isim_OE_x, isim_OE_y, isim_OE_z, isim_OE_time, isim_OE_pdg, isim_OE_mcpid, isim_OE_layer = [], [], [], [], [], [], []

        ireco_VB_x, ireco_VB_y, ireco_VB_z, ireco_VB_time, ireco_VB_layer = [], [], [], [], []
        ireco_VE_x, ireco_VE_y, ireco_VE_z, ireco_VE_time, ireco_VE_layer = [], [], [], [], []
        ireco_IB_x, ireco_IB_y, ireco_IB_z, ireco_IB_time, ireco_IB_layer = [], [], [], [], []
        ireco_IE_x, ireco_IE_y, ireco_IE_z, ireco_IE_time, ireco_IE_layer = [], [], [], [], []
        ireco_OB_x, ireco_OB_y, ireco_OB_z, ireco_OB_time, ireco_OB_layer = [], [], [], [], []
        ireco_OE_x, ireco_OE_y, ireco_OE_z, ireco_OE_time, ireco_OE_layer = [], [], [], [], []
        
        # Matching using LCRelation:
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
            pdg = mcp.getPDG()

            #print(mcp_tlv)
            imcp_pt.append(mcp_tlv.Perp())
            imcp_eta.append(mcp_tlv.Eta())
            imcp_phi.append(mcp_tlv.Phi())
            imcp_theta.append(mcp_tlv.Theta())
            ipdgid.append(pdg)
            istatus.append(mcp.getGeneratorStatus())
            iprod_vertex.append([mcp.getVertex()[i] for i in range(3)])
            iprod_time.append(mcp.getTime())
            iid.append(mcp.id())

            travel_dist = sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2) - sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2 + mcp.getVertex()[2]**2)
            iprod_traveldist.append(travel_dist)
                
            #print("PID, Status, pT, eta, phi: ", pdg, mcp.getGeneratorStatus(), mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())

            daughters = mcp.getDaughters()
            
            for daughter in daughters:
                idaughter_pdgid.append(daughter.getPDG())

            #print(idaughter_pdgid[0:10])
            stau_ids = [1000015, 2000015]
            if abs(mcp.getPDG()) in stau_ids: 
                stau_tracks = relation.getRelatedToObjects(mcp)

                imcp_stau_pt.append(mcp_tlv.Perp())
                imcp_stau_eta.append(mcp_tlv.Eta())
                imcp_stau_phi.append(mcp_tlv.Phi())
                imcp_stau_theta.append(mcp_tlv.Theta())
                imcp_stau_traveldist.append(sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2) - sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2 + mcp.getVertex()[2]**2))
                n_mcp_stau += 1

                
                #print("Stau Truth pt, eta, phi:", mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())
                #print("Stau truth pt, eta, phi:", imcp_stau_pt, imcp_stau_eta, imcp_stau_phi)

                for st in stau_tracks:
                    Bfield = 5 #T, 3.57 for legacy
                    theta = np.pi/2- np.arctan(st.getTanLambda())
                    phi = st.getPhi()
                    eta = -np.log(np.tan(theta/2))
                    pt  = 0.3 * Bfield / fabs(st.getOmega() * 1000.)
                    ptres = (mcp_tlv.Perp() - pt) / mcp_tlv.Perp() 
                    track_tlv = ROOT.TLorentzVector()
                    track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
                    st_nhitz = st.getTrackerHits().size()
                    d0 = st.getD0()
                    z0 = st.getZ0()
                    
                    
                    #hits
                    iist_x = []
                    iist_y = []
                    iist_z = []
                    #iist_hit_pdg = []
                    iist_time = []
                    iist_corrected_time = []
                    iist_hit_layer = []
                    iist_hit_detector = []
                    iist_hit_side = []
                    st_pix_nhit = 0
                    st_inn_nhit = 0
                    st_out_nhit = 0
                    lastLayer = -1
                    nLayersCrossed = 0
                    for hit in st.getTrackerHits():
                        encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                        decoder = pyLCIO.UTIL.BitField64(encoding)
                        cellID = int(hit.getCellID0())
                        decoder.setValue(cellID)
                        detector = decoder["system"].value()
                        layer = decoder['layer'].value()
                        side = decoder["side"].value()
                        if (lastLayer != layer):
                            if (detector <= 2): # if vdx
                                nLayersCrossed += 0.5
                                #if (detector == 1 or detector == 2):

                            #print("increment nlayerscrossed")
                            nLayersCrossed += 1 ### NOTE counting each of vertex doublet layers as individual layer
                           
                        if detector == 1 or detector == 2:
                            st_pix_nhit += 1
                        if detector == 3 or detector == 4:
                            st_inn_nhit += 1
                        if detector == 5 or detector == 6:
                            st_out_nhit += 1
                        lastLayer = layer

                        st_position = hit.getPosition()
                        st_pos_x = st_position[0]
                        st_pos_y = st_position[1]
                        st_pos_z = st_position[2]

                        nTotalHits = (st_pix_nhit)/2.0 + st_inn_nhit + st_out_nhit
    
                        if (nTotalHits < 3.5): # account for if 1 part of vertex doublet misses hit
                            #print("doesn't pass nhits cut, skip stau track")
                            continue 


                        d = sqrt(st_position[0]*st_position[0] + st_position[1]*st_position[1] + st_position[2]*st_position[2])
                        tof = d/speedoflight

                        resolution = 0.03
                        if detector > 2:
                            resolution = 0.06

                        st_corrected_t = hit.getTime()*(1.+ROOT.TRandom3(ievt).Gaus(0., resolution)) - tof

                        if (nLayersCrossed >= 3.5):
                    
                            iist_x.append(st_pos_x) 
                            iist_y.append(st_pos_y)
                            iist_z.append(st_pos_z)
                            #iist_hit_pdg.append(hit_pdg)
                            iist_time.append(hit.getTime())
                            iist_corrected_time.append(st_corrected_t)
                            iist_hit_detector.append(detector)
                            iist_hit_layer.append(layer)
                            iist_hit_side.append(side)  

                            #print(iist_time)

                    if (nLayersCrossed >= 3.5):

                        #appending matched mcp values
                        ist_pt_match.append(imcp_stau_pt)
                        ist_eta_match.append(imcp_stau_eta)
                        ist_phi_match.append(imcp_stau_phi)
                        ist_theta_match.append(imcp_stau_theta)

                        #appending matched track values
                        ist_track_pt.append([pt])
                        ist_track_theta.append([theta])
                        ist_track_phi.append([phi])
                        ist_ndf.append([st.getNdf()])
                        ist_chi2.append([st.getChi2()])
                        ist_d0.append([d0])
                        ist_z0.append([z0])
                        ist_nhits.append([st_nhitz])
                        ist_pt_res.append([pt_res])

                        st_pix_nhits.append([st_pix_nhit])
                        st_inn_nhits.append([st_inn_nhit])
                        st_out_nhits.append([st_out_nhit])
                        ist_x.append(iist_x)
                        ist_y.append(iist_y)
                        ist_z.append(iist_z)
                        #ist_hit_pdg.append(iihit_pdg)
                        ist_time.append(iist_time)
                        ist_corrected_time.append(iist_corrected_time)
                        ist_hit_detector.append(iist_hit_detector)
                        ist_hit_layer.append(iist_hit_layer)
                        ist_hit_side.append(iist_hit_side)

                        #print(ist_time)


        # Loop over the track objects
        for track in trackCollection:
            Bfield = 5 #T, 3.57 for legacy
            theta = np.pi/2- np.arctan(track.getTanLambda())
            phi = track.getPhi()
            eta = -np.log(np.tan(theta/2))
            pt  = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
            track_tlv = ROOT.TLorentzVector()
            track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
            # dr = mcp_tlv.DeltaR(track_tlv) # I don't think this works
            nhitz = track.getTrackerHits().size()
            d0 = track.getD0()
            z0 = track.getZ0()
            inhits.append(nhitz)
            itrack_pt.append(pt)
            itrack_eta.append(eta)
            itrack_theta.append(theta)
            itrack_phi.append(track.getPhi())

            indf.append(track.getNdf())
            ichi2.append(track.getChi2())
            # print("Reco pt, eta, phi, nhits, dr:", pt, eta, phi, nhitz, dr)

            # for j, particle_pt in enumerate(imcp_pt):
            #     particle_eta = imcp_eta[j]
            #     particle_phi = imcp_phi[j]
            #     mcp_tlv = ROOT.TLorentzVector()
            #     mcp_tlv.SetPtEtaPhiE(particle_pt, particle_eta, particle_phi, 0)
            #     dr = mcp_tlv.DeltaR(track_tlv)
            #     if dr < 0.005:
            #         #print(particle_pt, pt)
            #         ptres = (particle_pt - pt) / particle_pt
            #         #print(j, dr)
            #         num_matched_tracks += 1

            #         id0_res_vs_pt.append([particle_pt, d0])
            #         id0_res_vs_eta.append([particle_eta, d0])
            #         iz0_res_vs_pt.append([particle_pt, z0])
            #         iz0_res_vs_eta.append([particle_eta, z0])
            #         ipt_res_vs_eta.append([particle_eta, ptres])
            #         ipt_res_vs_pt.append([particle_pt, ptres])
            #         ipt_match.append(particle_pt)
            #         ieta_match.append(particle_eta)
            #         id0_res_match.append(d0)
            #         iz0_res_match.append(z0)
            #         ipt_res.append(ptres)
                    #itheta_match.append(mcp_stau_theta[j])

            # Hit per track

            iix = []
            iiy = []
            iiz = []
            iihit_pdg = []
            iitime = []
            iicorrected_time = []
            iihit_layer = []
            iihit_detector = []
            iihit_side = []
            
            pixel_nhit = 0
            inner_nhit = 0
            outer_nhit = 0
            for hit in track.getTrackerHits():
                    try:
                        mcp = hit.getMCParticle()
                        hit_pdg = mcp.getPDG()
                    except:
                        hit_pdg = 0
                # now decode hits
                    encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                    decoder = pyLCIO.UTIL.BitField64(encoding)
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    detector = decoder["system"].value()
                    layer = decoder['layer'].value()
                    side = decoder["side"].value()
                    if detector == 1 or detector == 2:
                        pixel_nhit += 1
                    if detector == 3 or detector == 4:
                        inner_nhit += 1
                    if detector == 5 or detector == 6:
                        outer_nhit += 1
                    position = hit.getPosition()
                    pos_x = position[0]
                    pos_y = position[1]
                    pos_z = position[2]

                    d = sqrt(position[0]*position[0] + position[1]
                     * position[1] + position[2]*position[2])
                    tof = d/speedoflight

                    resolution = 0.03
                    if detector > 2:
                        resolution = 0.06

                    corrected_t = hit.getTime()*(1.+ROOT.TRandom3(ievt).Gaus(0., resolution)) - tof
                    
                    iix.append(pos_x)
                    iiy.append(pos_y)
                    iiz.append(pos_z)
                    iihit_pdg.append(hit_pdg)
                    iitime.append(hit.getTime())
                    iicorrected_time.append(corrected_t)
                    iihit_detector.append(detector)
                    iihit_layer.append(layer)
                    iihit_side.append(side)


            pixel_nhits.append([pixel_nhit])
            inner_nhits.append([inner_nhit])
            outer_nhits.append([outer_nhit])
            ix.append(iix)
            iy.append(iiy)
            iz.append(iiz)
            ihit_pdg.append(iihit_pdg)
            itime.append(iitime)
            icorrected_time.append(iicorrected_time)
            ihit_detector.append(iihit_detector)
            ihit_layer.append(iihit_layer)
            ihit_side.append(iihit_side)
                        
        # print("End of event \n")
        # This is here to check that we never reconstruct multiple muons
        # If we did, we'd have to match the correct muon to the MCP object to do eff/res plots
        # But since we don't, we can skip that step

        #MCP append
        mcp_pt.append(imcp_pt)
        mcp_eta.append(imcp_eta)
        mcp_phi.append(imcp_phi)
        mcp_theta.append(imcp_theta)
        pdgid.append(ipdgid)
        status.append(istatus)
        prod_vertex.append(iprod_vertex)
        prod_time.append(iprod_time)
        prod_traveldist.append(iprod_traveldist)
        id.append(iid)

        sim_VB_x.append(isim_VB_x); sim_VB_y.append(isim_VB_y); sim_VB_z.append(isim_VB_z); sim_VB_time.append(isim_VB_time); sim_VB_pdg.append(isim_VB_pdg); sim_VB_mcpid.append(isim_VB_mcpid); sim_VB_layer.append(isim_VB_layer)
        sim_VE_x.append(isim_VE_x); sim_VE_y.append(isim_VE_y); sim_VE_z.append(isim_VE_z); sim_VE_time.append(isim_VE_time); sim_VE_pdg.append(isim_VE_pdg); sim_VE_mcpid.append(isim_VE_mcpid); sim_VE_layer.append(isim_VE_layer)
        sim_IB_x.append(isim_IB_x); sim_IB_y.append(isim_IB_y); sim_IB_z.append(isim_IB_z); sim_IB_time.append(isim_IB_time); sim_IB_pdg.append(isim_IB_pdg); sim_IB_mcpid.append(isim_IB_mcpid); sim_IB_layer.append(isim_IB_layer)
        sim_IE_x.append(isim_IE_x); sim_IE_y.append(isim_IE_y); sim_IE_z.append(isim_IE_z); sim_IE_time.append(isim_IE_time); sim_IE_pdg.append(isim_IE_pdg); sim_IE_mcpid.append(isim_IE_mcpid); sim_IE_layer.append(isim_IE_layer)
        sim_OB_x.append(isim_OB_x); sim_OB_y.append(isim_OB_y); sim_OB_z.append(isim_OB_z); sim_OB_time.append(isim_OB_time); sim_OB_pdg.append(isim_OB_pdg); sim_OB_mcpid.append(isim_OB_mcpid); sim_OB_layer.append(isim_OB_layer)
        sim_OE_x.append(isim_OE_x); sim_OE_y.append(isim_OE_y); sim_OE_z.append(isim_OE_z); sim_OE_time.append(isim_OE_time); sim_OE_pdg.append(isim_OE_pdg); sim_OE_mcpid.append(isim_OE_mcpid); sim_OE_layer.append(isim_OE_layer)
        reco_VB_x.append(ireco_VB_x); reco_VB_y.append(ireco_VB_y); reco_VB_z.append(ireco_VB_z); reco_VB_time.append(ireco_VB_time); reco_VB_layer.append(ireco_VB_layer)
        reco_VE_x.append(ireco_VE_x); reco_VE_y.append(ireco_VE_y); reco_VE_z.append(ireco_VE_z); reco_VE_time.append(ireco_VE_time); reco_VE_layer.append(ireco_VE_layer)
        reco_IB_x.append(ireco_IB_x); reco_IB_y.append(ireco_IB_y); reco_IB_z.append(ireco_IB_z); reco_IB_time.append(ireco_IB_time); reco_IB_layer.append(ireco_IB_layer)
        reco_IE_x.append(ireco_IE_x); reco_IE_y.append(ireco_IE_y); reco_IE_z.append(ireco_IE_z); reco_IE_time.append(ireco_IE_time); reco_IE_layer.append(ireco_IE_layer)
        reco_OB_x.append(ireco_OB_x); reco_OB_y.append(ireco_OB_y); reco_OB_z.append(ireco_OB_z); reco_OB_time.append(ireco_OB_time); reco_OB_layer.append(ireco_OB_layer)
        reco_OE_x.append(ireco_OE_x); reco_OE_y.append(ireco_OE_y); reco_OE_z.append(ireco_OE_z); reco_OE_time.append(ireco_OE_time); reco_OE_layer.append(ireco_OE_layer)
       
        #DAUGHTERS append
        daughter_pdgid.append(idaughter_pdgid)

        #TRACK append
        track_pt.append(itrack_pt)
        track_eta.append(itrack_eta)
        track_theta.append(itrack_theta)
        track_phi.append(itrack_phi)
        ndf.append(indf)
        chi2.append(ichi2)

        
        #LCRELATION append
        #confirmed staus
            #mcp values
        if n_mcp_stau > 0:
            mcp_stau_pt.append(imcp_stau_pt)
            mcp_stau_eta.append(imcp_stau_eta)
            mcp_stau_phi.append(imcp_stau_phi)
            mcp_stau_theta.append(imcp_stau_theta)
            mcp_stau_traveldist.append(imcp_stau_traveldist)

            #print("Outer list stau pt, eta, phi:", mcp_stau_pt, mcp_stau_eta, mcp_stau_phi)


        st_pt_match.append(ist_pt_match)
        st_eta_match.append(ist_eta_match)
        st_phi_match.append(ist_phi_match)
        st_theta_match.append(ist_theta_match)
            #track values
        st_track_pt.append(ist_track_pt)
        st_track_theta.append(ist_track_theta)
        st_track_phi.append(ist_track_phi)
        st_ndf.append(ist_ndf)
        st_chi2.append(ist_chi2)
        st_d0.append(ist_d0)
        st_z0.append(ist_z0)
            #hit values
        st_nhits.append(ist_nhits)
        st_pt_res.append(ist_pt_res)
        st_x.append(ist_x)
        st_y.append(ist_y)
        st_z.append(ist_z)
        #st_hit_pdg = []
        st_time.append(ist_time)
        st_corrected_time.append(ist_corrected_time)
        st_hit_detector.append(ist_hit_detector)
        st_hit_layer.append(ist_hit_layer)
        st_hit_side.append(ist_hit_side)

        #print(st_time)
        
        #HIT append
        i+=1
        nhits.append(inhits)
        x.append(ix)
        y.append(iy)
        z.append(iz)
        hit_pdgid.append(ihit_pdg)
        time.append(itime)
        corrected_time.append(icorrected_time)
        hit_layer.append(ihit_layer)
        hit_detector.append(ihit_detector)
        hit_side.append(ihit_side)
        
        
        # if len(id0_res_vs_pt) > 0:
        #     #pt_res_hits.append(ipt_res_hits)
        #     d0_res_vs_pt.append(id0_res_vs_pt)
        #     d0_res_vs_eta.append(id0_res_vs_eta)
        #     z0_res_vs_pt.append(iz0_res_vs_pt)
        #     z0_res_vs_eta.append(iz0_res_vs_eta)
        #     pt_res_vs_eta.append(ipt_res_vs_eta)
        #     pt_res_vs_pt.append(ipt_res_vs_pt)
        #     pt_res.append(ipt_res)
        #     pt_match.append(ipt_match)
        #     eta_match.append(ieta_match)
        #     theta_match.append(itheta_match)
        #     d0_res_match.append(id0_res_match)
        #     z0_res_match.append(iz0_res_match)
        
    reader.close()

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
print("\nSummary statistics:")
print("Ran over %i events."%i)
print("Found:")
# print("\t%i MCPs"%len(np.ravel(mcp_pt)))
# print("\t%i Stau MCPs"%n_mcp_stau)
print("\tSanity check mcpCollection length:", len(imcp_pt))
print("\tSanity check trackCollection length:", len(itrack_pt))
# print("\tSanity check # hits:", len(np.ravel(x)))
# print("\t%i PFOs"%hists["pfo_pt"].GetEntries())
# print("\t%i mu PFOs"%hists["pfo_mu_pt"].GetEntries())
# print('\t%i matched muon tracks'%(num_matched_tracks))
# print('\t%i duplicates eliminated'%num_dupes)
# print('\t%i hard radiations discarded'%hard_rad_discard)
# print('\t%i fake tracks'%num_fake_tracks)
# print('\t%i GeV'%np.max(mcp_stau_pt))


#Make a list of all the data you want to save
data_list = {
    "sim_VB_x": sim_VB_x, "sim_VB_y": sim_VB_y, "sim_VB_z": sim_VB_z, "sim_VB_time": sim_VB_time, "sim_VB_pdg": sim_VB_pdg, "sim_VB_mcpid": sim_VB_mcpid, "sim_VB_layer": sim_VB_layer,
    "sim_VE_x": sim_VE_x, "sim_VE_y": sim_VE_y, "sim_VE_z": sim_VE_z, "sim_VE_time": sim_VE_time, "sim_VE_pdg": sim_VE_pdg, "sim_VE_mcpid": sim_VE_mcpid, "sim_VE_layer": sim_VE_layer,
    "sim_IB_x": sim_IB_x, "sim_IB_y": sim_IB_y, "sim_IB_z": sim_IB_z, "sim_IB_time": sim_IB_time, "sim_IB_pdg": sim_IB_pdg, "sim_IB_mcpid": sim_IB_mcpid, "sim_IB_layer": sim_IB_layer,
    "sim_IE_x": sim_IE_x, "sim_IE_y": sim_IE_y, "sim_IE_z": sim_IE_z, "sim_IE_time": sim_IE_time, "sim_IE_pdg": sim_IE_pdg, "sim_IE_mcpid": sim_IE_mcpid, "sim_IE_layer": sim_IE_layer,
    "sim_OB_x": sim_OB_x, "sim_OB_y": sim_OB_y, "sim_OB_z": sim_OB_z, "sim_OB_time": sim_OB_time, "sim_OB_pdg": sim_OB_pdg, "sim_OB_mcpid": sim_OB_mcpid, "sim_OB_layer": sim_OB_layer,
    "sim_OE_x": sim_OE_x, "sim_OE_y": sim_OE_y, "sim_OE_z": sim_OE_z, "sim_OE_time": sim_OE_time, "sim_OE_pdg": sim_OE_pdg, "sim_OE_mcpid": sim_OE_mcpid, "sim_OE_layer": sim_OE_layer,
    "reco_VB_x": reco_VB_x, "reco_VB_y": reco_VB_y, "reco_VB_z": reco_VB_z, "reco_VB_time": reco_VB_time, "reco_VB_layer": reco_VB_layer,
    "reco_VE_x": reco_VE_x, "reco_VE_y": reco_VE_y, "reco_VE_z": reco_VE_z, "reco_VE_time": reco_VE_time, "reco_VE_layer": reco_VE_layer,
    "reco_IB_x": reco_IB_x, "reco_IB_y": reco_IB_y, "reco_IB_z": reco_IB_z, "reco_IB_time": reco_IB_time, "reco_IB_layer": reco_IB_layer,
    "reco_IE_x": reco_IE_x, "reco_IE_y": reco_IE_y, "reco_IE_z": reco_IE_z, "reco_IE_time": reco_IE_time, "reco_IE_layer": reco_IE_layer,
    "reco_OB_x": reco_OB_x, "reco_OB_y": reco_OB_y, "reco_OB_z": reco_OB_z, "reco_OB_time": reco_OB_time, "reco_OB_layer": reco_OB_layer,
    "reco_OE_x": reco_OE_x, "reco_OE_y": reco_OE_y, "reco_OE_z": reco_OE_z, "reco_OE_time": reco_OE_time, "reco_OE_layer": reco_OE_layer
}

#data_list = {}

#MCP
data_list["mcp_pt"] = mcp_pt
data_list["mcp_eta"] = mcp_eta
data_list["mcp_phi"] = mcp_phi
data_list["mcp_theta"] = mcp_theta
data_list["pdgid"] = pdgid
data_list["status"] = status
data_list["prod_vertex"] = prod_vertex
data_list["prod_time"] = prod_time
data_list["prod_traveldist"] = prod_traveldist
data_list["id"] = id

#DAUGHTERS
data_list["daughter_pdgid"] = daughter_pdgid

#TRACK
data_list["d0_res"] = d0_res
data_list["z0_res"] = z0_res
data_list["pt_res"] = pt_res
data_list["track_pt"] = track_pt
data_list["track_eta"] = track_eta
data_list["track_theta"] = track_theta
data_list["track_phi"] = track_phi
data_list["ndf"] = ndf
data_list["chi2"] = chi2
data_list["nhits"] = nhits

#LCRELATION
    #stau pdgid matching
data_list["mcp_stau_pt"] = mcp_stau_pt
data_list["mcp_stau_eta"] = mcp_stau_eta
data_list["mcp_stau_phi"] = mcp_stau_phi
data_list["mcp_stau_theta"] = mcp_stau_theta
data_list["mcp_stau_traveldist"] = mcp_stau_traveldist
data_list["st_pt_match"] = st_pt_match
data_list["st_eta_match"] = st_eta_match
data_list["st_phi_match"] = st_phi_match
data_list["st_theta_match"] = st_theta_match
data_list["st_track_pt"] = st_track_pt
data_list["st_track_theta"] = st_track_theta
data_list["st_track_phi"] = st_track_phi
data_list["st_ndf"] = st_ndf
data_list["st_d0"] = st_d0
data_list["st_z0"] = st_z0
data_list["st_nhits"] = st_nhits
data_list["st_pt_res"] = st_pt_res
data_list["st_pix_nhits"] = st_pix_nhits
data_list["st_inn_nhits"] = st_inn_nhits
data_list["st_out_nhits"] = st_out_nhits
data_list["st_x"] = st_x
data_list["st_y"] = st_y
data_list["st_z"] = st_z
data_list["st_time"] = st_time
data_list["st_corrected_time"] = st_corrected_time
data_list["st_hit_detector"] = st_hit_detector
data_list["st_hit_layer"] = st_hit_layer
data_list["st_hit_side"] = st_hit_side

#HIT
data_list["pixel_nhits"] = pixel_nhits
data_list["inner_nhits"] = inner_nhits
data_list["outer_nhits"] = outer_nhits
# data_list["pt_res_hits"] = pt_res_hits
data_list["x"] = x
data_list["y"] = y
data_list["z"] = z
data_list["hit_pdgid"] = hit_pdgid
data_list["time"] = time
data_list["corrected_time"] = corrected_time
data_list["hit_layer"] = hit_layer
data_list["hit_detector"] = hit_detector
data_list["hit_side"] = hit_side

# After the loop is finished, save the data_list to a .json file
output_json = "/home/tateflicker/newest_files/4500_10_reco_32ns64ns.json"
with open(output_json, 'w') as fp:
    json.dump(data_list, fp)


#print("Sanity check st_time event 0 length:", len(st_time[0][0]))
