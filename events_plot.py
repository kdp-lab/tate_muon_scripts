import uproot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import awkward as ak

#pythia file
events = uproot.open("~/MG5_aMC_v3_5_1/muon_to_stau_decay/Events/run_07/tag_1_pythia8_events.root:Events")
events.show()
events.keys()
events.arrays()

branches = events.arrays()

#print(branches["Particle_pid"])
#print(branches["Particle_status"])

staus = abs(branches["Particle_pid"]) == 1000015
print(staus)

#status_2 = branches["Particle_status"] == 2
#print(status_2)

stau_px = branches["Particle_px"][staus]
stau_py = branches["Particle_py"][staus]
stau_pz = branches["Particle_pz"][staus]
stau_x = branches["Particle_x"][staus]
stau_y = branches["Particle_y"][staus]
stau_z = branches["Particle_z"][staus]
stau_M = branches["Particle_mass"][staus]
stau_t = (branches["Particle_ctau"][staus])/(300)       #proper lifetime
#stau_tof       time of flight in lab frame

#print(stau_px)
#print(stau_t)
#print(stau_M)

#print(branches["Particle_ctau"][staus])

#plt.hist(ak.flatten(stau_t), bins=100, range=(0, 5))
#plt.title("Stau Proper Lifetime (ns)")
#plt.show()

#eta
stau_p_squared = np.square(stau_px) + np.square(stau_py) + np.square(stau_pz)
stau_p = np.sqrt(stau_p_squared)
stau_eta = np.arctanh(stau_pz/stau_p)
print(stau_eta)

plt.hist(ak.flatten(stau_eta), bins=100, range=(-1,1))
plt.title("Stau Pseudorapidity")
plt.show()

#theta
stau_theta = 2*np.arctan((np.e)**-(stau_eta))
print(stau_theta)

plt.hist(ak.flatten(stau_theta), bins=100, range=(0,3))
plt.title("Stau Theta")
plt.show()

#pT
stau_pt = stau_pz*np.arcsinh(stau_eta)
print(stau_pt)

plt.hist(ak.flatten(stau_pt), bins=100, range=(0,3e+03))
plt.title("Stau pT")
plt.show()

#phi
stau_phi = np.arccos(stau_px/stau_pt)
print(stau_phi)

#px = events["Particle_px"].array(library="np")
#print(px)
#if px.type != np.dtype("float"):
#       px = px.astype("float")
#py = events["Particle_py"].array(library="np")
#print(py)
#if py.type != np.dtype("float"):
#       py = py.astype("float")
#pz = events["Particle_pz"].array(library="np")
#print(pz)
#if pz.type != np.dtype("float"):
#       pz = pz.astype("float")

#p_squared = np.square(px) + np.square(py) + np.square(pz)
#p = np.sqrt(p_squared)
#eta = np.arctanh(pz/p)
#print(eta)
#theta = 2*np.arctan((np.e)**-eta)

#phi =
#print(phi)

#pT = py*sin^-1(phi)

#calc_pT = py*np.arcsin(phi)

#madgraph file
#unweighted_events = uproot.open("~/MG5_aMC_v3_5_1/muon_to_stau_decay/Events/run_07/unweighted_events.root:LHEF")
#unweighted_events.show()

#pT_mg = unweighted_events["Particle.PT"].array(library="np")
#print(pT_mg)

#computing velocity:

#computing time of flight in lab frame:

#computing lifetime:

#decay position:

#number of charged status 1 particles:

                                                                                                                                                              105,1         78%



                                                                                                                                                              25,1          Top








