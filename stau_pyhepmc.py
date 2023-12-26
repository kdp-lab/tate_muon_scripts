import pyhepmc
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import os
import json
import argparse
from particle import literals as lp
from IPython.display import display

parser = argparse.ArgumentParser(description="Input configuration for running")
parser.add_argument("--mass", type=str, default="400", help="slepton mass in GeV")
parser.add_argument("--lifetime", type=str, default="1ns", help="lifetime in ns")
parser.add_argument("--nevents", type=int, default=-1, help="nevents to run")
parser.add_argument("--doTest", type=bool, default=False, help="run a test")

args, unknown = parser.parse_known_args()
mass = args.mass
lifetime = args.lifetime
doTest = args.doTest
nevents = args.nevents

pid = 1000015
status = 2
toGeV = 5000.

staus = []
events = []

with pyhepmc.open("/home/tate/MG5_aMC_v3_5_1/muon_to_stau_decay/Events/run_11/tag_1_pythia8_events.hepmc") as f:
    while True:
        event = f.read()
        if not event:
            break
        if doTest and event.event_number > 10:
            break
        if nevents > 0 and event.event_number > nevents:
            break
        if event.event_number % (nevents/10) == 0:
            print("Event", event.event_number)
        if doTest:
            print("In event", event.event_number)

        event = {}

        def sum_stau_energy(event):
            esum = 0.0
            for p in event.particles:
                if abs(p.pid) != lp.stau.pdgid:
                    continue
                if p.status != 2:
                    continue
                esum += p.momentum.e
            return esum

        sum_stau_energy(event)

