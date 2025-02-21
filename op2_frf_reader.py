# -*- coding: utf-8 -*-
"""This script reads a selected OP2 result file from an FRF (SEMFREQ) NASTRAN simulation, and
computes and plots the CoG transfer function obtained from the CBUSH forces (assumed to ) and
the specified model mass.
"""
import yaml
import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askopenfilename

import pyNastran
from pyNastran.op2.op2 import OP2

# Load general configuration YAML file
with open('frf_config.yaml', 'r') as stream:
    frf_config = yaml.safe_load(stream)

DIRS = frf_config['dirs']
EXCIT_DIR = frf_config['excit_dir']
COLORS = frf_config['colors']
SUBCASE_LABEL = frf_config['subcase_label']
SUBSYSTEM_MASS = frf_config['subsystem_mass']
BASE_ACCEL = frf_config['base_accel']

Tk().withdraw()
accepted_filetypes = [('OP2 files', '*.OP2')]
OP2_FILENAME = askopenfilename(
    title='Select an OP2 file...',
    filetypes=accepted_filetypes,
)
model = OP2()
model.read_op2(OP2_FILENAME)
# TODO: read base acceleration from BDF file

CBUSH_force_data = model.op2_results.force.cbush_force[SUBCASE_LABEL]
acceleration_data = model.accelerations[SUBCASE_LABEL]
elements = CBUSH_force_data.element
nodes = acceleration_data.node_gridtype[:,0]

# node_idx = np.where(nodes==node)[0]
# node_TF = abs((accel[1].data[:, node_idx, DIRS[dir]]))
freqs = CBUSH_force_data.freqs
FRF_accel= dict()
for dir in DIRS.keys():
    FRF_force = np.abs(CBUSH_force_data.data[:,:,DIRS[dir]].sum(axis=1))
    FRF_accel[dir] = FRF_force / SUBSYSTEM_MASS  / BASE_ACCEL
    peak_accel = np.max(FRF_accel[dir])
    peak_freq = np.where(FRF_accel[dir] == peak_accel)[0][0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    label_str = f'Peak acceleration of {peak_accel:.1f}g at {peak_freq:.0f} Hz'
    ax.loglog(freqs, FRF_accel[dir], color=COLORS[dir], label=label_str)
    ax.legend(loc='lower right')
    plt.title(f'{EXCIT_DIR} axis excitation - Antenna Set CoG {dir} Acceleration Response')
    plt.ylim([np.min(FRF_accel[dir]), peak_accel])
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Acceleration Transfer Function [g/g]')
    plt.grid()
plt.show()