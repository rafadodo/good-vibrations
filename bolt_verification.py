# -*- coding: utf-8 -*-
"""This script performs the pre- and post-processing for the ECSS-style structural verification of
a set of bolted joints. It uses inputs (material and thread properties, and bolt loads) from sheet
files whose paths are read from a YAML configuration file, and prints bolt margins in an output
sheet whose path is also read from the configuration file. The rest of the parameters needed are
defined as constants inside the script.
"""
import pandas as pd
import numpy as np
import yaml

import bolt_utils


# Load general configuration YAML file
with open('bolt_general_config.yaml', 'r') as stream:
    general_config = yaml.safe_load(stream)

# Load safety factors configuration YAML file
with open('bolt_safety_config.yaml', 'r') as stream:
    safety_config = yaml.safe_load(stream)

# Load tightening configuration YAML file
with open('bolt_tightening_config.yaml', 'r') as stream:
    tightening_config = yaml.safe_load(stream)

# Extract bolt parameters from YAML
NUMBER_OF_BOLTS = general_config['number_of_bolts']
CLAMPED_MAT_NAME = general_config['clamped_mat_name']
BOLT_MAT_NAME = general_config['bolt_mat_name']
THREAD_DIAM = general_config['thread_diam']
HOLE_DIAM = general_config['hole_diam']
HALF_THREAD_ANGLE = np.deg2rad(general_config['half_thread_angle'])
BOLT_BENDING_LOGICAL = general_config['bolt_bending_logical']

# Conversion factors
FORCE_UNITS_FACTOR = general_config['force_units_factor']
MOMENT_UNITS_FACTOR = general_config['moment_units_factor']

# Load case params
LOAD_CASES = general_config['load_cases']
LOADCASE_IDX = general_config['loadcase_idx']

# Thread params
THREAD_SHEET_FILENAME = general_config['thread_sheet_filename']
THREAD_SHEET_COL_RANGE = general_config['thread_sheet_col_range']
THREAD_SHEET_HEADERS_ROW = general_config['thread_sheet_headers_row']

# Load params
LOADS_SHEET_FILENAME = general_config['loads_sheet_filename']
LOAD_SHEET_TABNAME = LOAD_CASES[LOADCASE_IDX]
LOADS_SHEET_COL_RANGE = general_config['loads_sheet_col_range']
LOADS_SHEET_HEADERS_ROW = general_config['loads_sheet_headers_row']

# Material params
MATERIALS_SHEET_FILENAME = general_config['materials_sheet_filename']
MATERIALS_SHEET_COL_RANGE = general_config['materials_sheet_col_range']
MATERIALS_SHEET_HEADERS_ROW = general_config['materials_sheet_headers_row']

# Output
OUTPUT_SHEET_FILENAME = general_config['output_sheet_filename']
OUTPUT_DF_HEADERS = general_config['output_df_headers']


# Tightening parameters
TOOL_SCATTER = tightening_config['tool_scatter']
PRELOAD_LOSS = tightening_config['preload_loss']
M_PREV_MIN = tightening_config['m_prev_min']
M_PREV_MAX = tightening_config['m_prev_max']
NUT_FACTOR = tightening_config['nut_factor']
NUT_F_UNCERT = tightening_config['nut_f_uncert']
HEAD_MIN_FRIC = tightening_config['head_min_fric']
CLAMP_MIN_FRIC = tightening_config['clamp_min_fric']
FORCE_RATIO_LOW = tightening_config['force_ratio_low']
FORCE_RATIO_UP = tightening_config['force_ratio_up']

# Safety Factors
K_P = safety_config['k_p']
K_M = safety_config['k_m']
FOS_Y = safety_config['fos_y']
FOS_U = safety_config['fos_u']
FOS_SEPARATION = safety_config['fos_separation']
FOS_SLIP = safety_config['fos_slip']
K_F = safety_config['k_f']
FOS_INSTALL_Y = safety_config['fos_install_y']
FOS_INSTALL_U = safety_config['fos_install_u']
ENG_LEN_FACTOR = safety_config['eng_len_factor']

mat_properties = bolt_utils.get_mat_properties(sheet_filename=MATERIALS_SHEET_FILENAME,
                                               sheet_col_range=MATERIALS_SHEET_COL_RANGE,
                                               headers_row=MATERIALS_SHEET_HEADERS_ROW
                                               )
clamped_mat_properties = mat_properties[CLAMPED_MAT_NAME]
clamp_E = clamped_mat_properties['E']
clamp_G = clamped_mat_properties['G']
clamp_Sy = clamped_mat_properties['Tensile yield strength']
clamp_Su = clamped_mat_properties['Tensile ultimate strength']
clamp_Tu = clamped_mat_properties['Shear ultimate strength']
bolt_mat_properties = mat_properties[BOLT_MAT_NAME]
bolt_E = bolt_mat_properties['E']
bolt_G = bolt_mat_properties['G']
bolt_Sy = bolt_mat_properties['Tensile yield strength']
bolt_Su = bolt_mat_properties['Tensile ultimate strength']
bolt_Ty = bolt_mat_properties['Shear yield strength']
bolt_Tu = bolt_mat_properties['Shear ultimate strength']

metric_thread_size = f"M{THREAD_DIAM}"
thread_props = bolt_utils.get_thread_properties(sheet_filename=THREAD_SHEET_FILENAME,
                                                sheet_col_range=THREAD_SHEET_COL_RANGE,
                                                headers_row=THREAD_SHEET_HEADERS_ROW,
                                                thread_size=metric_thread_size)
pitch = thread_props['P']
d_head = thread_props['dk_min']
Torque_nom = 1700
d1 = THREAD_DIAM - 1.0825 * pitch # Thread basic minor diameter [mm] %Verified
d2 = THREAD_DIAM - 0.64952 * pitch # Pitch Diameter [mm] %Verified
d3 = THREAD_DIAM - 1.22687 * pitch # Minor Diamter [mm] %Verified
d_stress = 0.5 * (d2+d3) # Stress Diameter [mm] %Verified
underhead_diam = 0.5*(d_head + HOLE_DIAM) # Effective diameter for under head friction [mm]
A_nom = 0.25 * np.pi * THREAD_DIAM**2 # Nominal Area [mm2]  %Verified
A_stress = 0.25 * np.pi * d_stress**2 # Stress Area [mm2] %Verified
A_stiff = 0.25 * np.pi * d3**2 # Stiffness Area [mm2] %Verified
engaged_length = THREAD_DIAM * ENG_LEN_FACTOR # Engaged length [mm]

#Preload calculation
preload_max =  ((1+TOOL_SCATTER)*Torque_nom - M_PREV_MIN) / (NUT_FACTOR/(1+NUT_F_UNCERT)) / THREAD_DIAM
preload_min = (1-PRELOAD_LOSS) * ((1-TOOL_SCATTER)*Torque_nom - M_PREV_MAX) / (NUT_FACTOR/(1-NUT_F_UNCERT)) / THREAD_DIAM


beam_forces = bolt_utils.get_beam_forces(sheet_filename=LOADS_SHEET_FILENAME,
                                         sheet_tabname=LOAD_SHEET_TABNAME,
                                         sheet_col_range=LOADS_SHEET_COL_RANGE,
                                         headers_row=LOADS_SHEET_HEADERS_ROW)

output_array = np.zeros(shape=(NUMBER_OF_BOLTS, len(OUTPUT_DF_HEADERS)))
for bolt_index in range(NUMBER_OF_BOLTS):
    F_axial = abs(beam_forces['Axial_F'][bolt_index])/FORCE_UNITS_FACTOR*K_P*K_M # Axial Load Modulus [N]
    F_shear = beam_forces['Shear_F'][bolt_index]/FORCE_UNITS_FACTOR*K_P*K_M
    if BOLT_BENDING_LOGICAL:
        M_bending = beam_forces['Bending_M'][bolt_index]/MOMENT_UNITS_FACTOR*K_P*K_M # Bending moment if aplicable [N.mm]
    else:
        M_bending = 0
    FaxCorr = F_axial + 8*M_bending/d_stress # Corrected force with bending moment

    # Tightening MoS
    TauMax = (Torque_nom*(1+TOOL_SCATTER) - underhead_diam/2*preload_max*HEAD_MIN_FRIC) / (np.pi*d_stress**3/16) # Max Shear Stress
    SigmaMax = preload_max/A_stress # Max. Normal Stress (neglecting tamperature change)
    Svm = np.sqrt(SigmaMax**2 + 3*TauMax**2) # Von Mises Stress
    Tight_MoSy = bolt_Sy/(Svm*FOS_INSTALL_Y) -1 # Tightening MoS Yield
    Tight_MoSu = bolt_Su/(Svm*FOS_INSTALL_U) -1 # Tightening MoS Ultimate %ESA recomineda no aplicar SF y K_M por que ya estan incluidas las incertezas

    # Slip MoS Stand Alone
    Slip_MoS = ((preload_min-(1-FORCE_RATIO_LOW)*F_axial)*CLAMP_MIN_FRIC)/(F_shear*FOS_SLIP*K_F) - 1
    Slip_LSF = preload_min/((1-FORCE_RATIO_LOW)*F_axial + (F_shear*FOS_SLIP*K_F)/CLAMP_MIN_FRIC)

    # Separation
    Sep_MoS = preload_min/((1-FORCE_RATIO_LOW)*F_axial*FOS_SEPARATION*K_F) - 1
    Sep_LSF = preload_min/((1-FORCE_RATIO_LOW)*F_axial*FOS_SEPARATION*K_F)

    # Shear: Bolt margins of safety due to the pure shear
    Shear_MoSy = (bolt_Ty*A_stress) / (F_shear*FOS_Y*K_F) - 1 # Bolt Shear MoS Yield
    Shear_MoSu = (bolt_Tu*A_stress) / (F_shear*FOS_U*K_F) - 1 # Bolt Shear MoS Ultimate
    Shear_LSFy = (bolt_Ty*A_stress) / (F_shear*FOS_Y*K_F) # Bolt Shear MoS Yield
    Shear_LSFu = (bolt_Tu*A_stress) / (F_shear*FOS_U*K_F) # Bolt Shear MoS Ultimate

    # Axial
    Ax_MoSy = A_stress*bolt_Sy / (preload_max + FORCE_RATIO_UP*FaxCorr*FOS_Y*K_F) - 1
    Ax_MoSu = A_stress*bolt_Su / (preload_max + FORCE_RATIO_UP*FaxCorr*FOS_U*K_F) - 1
    Ax_LSFy = (A_stress*bolt_Sy - preload_max) /(FORCE_RATIO_UP*FaxCorr*FOS_Y*K_F)
    Ax_LSFu = (A_stress*bolt_Su - preload_max) /(FORCE_RATIO_UP*FaxCorr*FOS_U*K_F)

    # Combined Stress
    Ra_Y = (preload_max + FORCE_RATIO_UP*FaxCorr*FOS_Y*K_F) / (bolt_Sy*A_stress)
    Ra_U = (preload_max + FORCE_RATIO_UP*FaxCorr*FOS_U*K_F) / (bolt_Su*A_stress)
    Rq_Y = (F_shear*FOS_Y*K_F) / (bolt_Ty*A_stress)
    Rq_U = (F_shear*FOS_U*K_F) / (bolt_Tu*A_stress)
    Comb_MoSy = (Ra_Y**2 + Rq_Y**2)**(-1/2)-1 # Combined MoS to Yield
    Comb_MoSu = (Ra_U**2 + Rq_U**2)**(-1/2)-1 # Combined MoS to Ultimate

    # Combined LSF
    aFy =  (FORCE_RATIO_UP*FaxCorr*FOS_Y*K_F / (bolt_Sy*A_stress))**2 + Rq_Y**2
    bFy = 2*preload_max*FORCE_RATIO_UP*FaxCorr*FOS_Y*K_F / (bolt_Sy*A_stress)**2
    cFy = preload_max**2 / (bolt_Sy*A_stress)**2 - 1

    aFu = (FORCE_RATIO_UP*FaxCorr*FOS_U*K_F / (bolt_Su*A_stress))**2 + Rq_U**2
    bFu = 2*preload_max*FORCE_RATIO_UP*FaxCorr*FOS_U*K_F / (bolt_Su*A_stress)**2
    cFu = preload_max**2 / (bolt_Su*A_stress)**2 - 1

    Comb_LSFy = (-bFy + np.sqrt(bFy**2-4*aFy*cFy)) / (2*aFy)
    Comb_LSFu = (-bFu + np.sqrt(bFu**2-4*aFu*cFu)) / (2*aFu)

    # Thread Shear Pull Out
    effect_length = engaged_length - 0.8*pitch # Effective engaged thread [mm]

    A_thread_Fem = np.pi*THREAD_DIAM*effect_length / pitch * (0.5*pitch+(THREAD_DIAM-d2)*np.tan(HALF_THREAD_ANGLE)) # Female Thread Area [mm2]
    A_thread_Male = np.pi*d1*effect_length / pitch * (0.5*pitch+(d2-d1)*np.tan(HALF_THREAD_ANGLE)) # Male Thread Area [mm2]

    FtotMax = preload_max + FORCE_RATIO_UP*F_axial*FOS_U*K_F # %Maximum Load for ultimate calculation [N]

    Rs = (clamp_Tu*A_thread_Fem)/ (bolt_Tu*A_thread_Male) # Shear Strength Ratio [~]

    c1 = 1 # For threaded hole or...
    # c1 = 3.8*Sw/D - Sw**2/D**2 - 2.61 # For Nut Sw = Wrench size [mm]

    # c2
    c2b = 0.897
    c2n = 0.728 + 1.769*Rs - 2.896*Rs**2 + 1.296*Rs**2

    ThrBMoSu = (bolt_Tu*A_thread_Male*c1*c2b) / FtotMax - 1 # Thread Shear Pull Out MALE  (Bolt) MoS
    ThrNMoSu = (clamp_Tu*A_thread_Fem*c1*c2n) / FtotMax - 1 # Thread Shear Pull Out FEMALE  (Nut/Thread) MoS
    ThrMoSu = min(ThrBMoSu,ThrNMoSu)

    ThrB_LSFu = (bolt_Tu*A_thread_Male*c1*c2b - preload_max) / (FORCE_RATIO_UP*F_axial*FOS_U*K_F) # Thread Shear Pull Out MALE  (Bolt) LS
    ThrN_LSFu = (clamp_Tu*A_thread_Fem*c1*c2n - preload_max) / (FORCE_RATIO_UP*F_axial*FOS_U*K_F) # Thread Shear Pull Out FEMALE  (Nut/Thread) LS
    Thr_LSFu = min(ThrB_LSFu,ThrN_LSFu)

    ThrBminLu = ((bolt_Tu*A_thread_Male/effect_length*c1*c2b)/FtotMax)**-1 + 0.8*pitch # Thread Shear Pull Out MALE  (Bolt) MoS=0
    ThrNminLu = ((clamp_Tu*A_thread_Fem/effect_length*c1*c2n)/FtotMax)**-1 + 0.8*pitch # Thread Shear Pull Out FEMALE  (Nut/Thread) MoS=0

    FtotMaxBL = A_stress*bolt_Su # Maximum Load for ultimate calculation [N] - Bolt Axial Limited AxLSFu

    # Thread Shear Pull Out MALE  (Bolt) - Bolt Filure Before Thread Failure
    ThrBminLuREQ = ((bolt_Tu*A_thread_Male/effect_length*c1*c2b)/FtotMaxBL)**-1 + 0.8*pitch
    # Thread Shear Pull Out FEMALE  (Nut/Thread) - Bolt Filure Before Thread Failure
    ThrNminLuREQ = ((clamp_Tu*A_thread_Fem/effect_length*c1*c2n)/FtotMaxBL)**-1 + 0.8*pitch
    ThrminLuREQ = max(ThrBminLuREQ,ThrNminLuREQ)

    output_array[bolt_index,:] = np.array([Tight_MoSy, Slip_LSF, Sep_LSF, Shear_LSFu, Ax_LSFu, Comb_LSFu, Thr_LSFu])


output_df_indexes = [f"Bolt N°{bolt_num+1}" for bolt_num in range(NUMBER_OF_BOLTS)]
output_df = pd.DataFrame(data=np.round(output_array,2),
                         index=output_df_indexes,
                         columns=OUTPUT_DF_HEADERS)
output_df.to_excel(OUTPUT_SHEET_FILENAME, sheet_name=LOAD_SHEET_TABNAME)