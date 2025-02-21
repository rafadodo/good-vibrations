# -*- coding: utf-8 -*-
"""This module contains functions auxiliary to the pre- and post-processing of bolted joint analysis
"""
import pandas as pd

def get_thread_properties(sheet_filename, sheet_col_range, headers_row, thread_size):
    df_thread_props = pd.read_excel(sheet_filename,
                                    usecols=sheet_col_range,
                                    header=headers_row
                                    )
    thread_props = dict(zip(df_thread_props['Thread'], df_thread_props[thread_size]))
    for key in thread_props:
        thread_props[key] = float(thread_props[key])
    return thread_props


def get_beam_forces(sheet_filename, sheet_tabname, sheet_col_range, headers_row):
    "Returns a dict containing "
    df_beam_forces = pd.read_excel(sheet_filename,
                                   sheet_name=sheet_tabname,
                                   usecols=sheet_col_range,
                                   header=headers_row
                                   )
    df_beam_forces_stripped = df_beam_forces.dropna(how='all')
    beam_forces_dict = df_beam_forces_stripped.to_dict()

    return beam_forces_dict


def get_mat_properties(sheet_filename, sheet_col_range, headers_row):
    "Returns a dict containing "
    df_mat_properties = pd.read_excel(sheet_filename,
                                      usecols=sheet_col_range,
                                      header=headers_row
                                      )
    df_mat_properties = df_mat_properties.T
    df_mat_properties.columns = df_mat_properties.iloc[0]
    mat_properties = df_mat_properties.to_dict()

    return mat_properties
