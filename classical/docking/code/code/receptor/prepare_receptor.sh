#!/bin/bash

# Prepare the Chain A PfDHFR-TS receptor for molecular docking
# Input: 6A2L_chainA.pdbqt
# Tool: Meeko mk_prepare_receptor.py

mk_prepare_receptor.py -i 6A2L_chainA.pdbqt
