#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    Script to perform water bridges analysis between two atom
#    selections in a molecular dynamics simulation as implemented 
#    in the MDAnalysis package.
#    More information about water bridges analysis can be found at:
#
#    https://www.mdanalysis.org/docs/documentation_pages/analysis/
#    wbridge_analysis.html
#
#    Copyright (C) 2020 Valentina Sora <sora.valentina1@gmail.com>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  
#    If not, see <http://www.gnu.org/licenses/>.

import argparse
import sys

import MDAnalysis as mda
from MDAnalysis.analysis.hbonds.wbridge_analysis import WaterBridgeAnalysis
import pandas as pd

class WaterBridgeAnalysis_Custom(WaterBridgeAnalysis):
    
    """Class for water bridges analysis including atom
    names for additional forcefields (for CHARMM27 and GLYCAM06
    they are already present in the WaterBridgeAnalysis base
    class, but the attribute cannot be simply updated in
    subclasses).
    """
    
    # forcefield-specific donors
    DEFAULT_DONORS = \
        {"CHARMM27" : \
            tuple(set(["N", "OH2", "OW", "NE", "NH1", "NH2", \
                       "ND2", "SG", "NE2", "ND1", "NZ", "OG", \
                       "OG1", "NE1", "OH"])), \
        
        "GLYCAM06" : \
            tuple(set(["N", "NT", "N3", "OH", "OW"])), \
        
        "CHARMM22ST_PHOSPHO" : \
            tuple(set(["NH1", "NE2", "NE", "OW", "OH2", "OG1", \
                       "NE1", "N", "OH", "SG", "ND1", "NH2", \
                       "OG", "NZ", "ND2"]))}
    
    # forcefield-specific acceptors
    DEFAULT_ACCEPTORS = \
        {"CHARMM27" : \
            tuple(set(["O", "OC1", "OC2", "OH2", "OW", "OD1", \
                       "OD2", "SG", "OE1", "OE1", "OE2", "ND1", \
                       "NE2", "SD", "OG", "OG1", "OH"])), \
        "GLYCAM06" : \
            tuple(set(["N", "NT", "O", "O2", "OH", \
                       "OS", "OW", "OY", "SM"])), \

        "CHARMM22ST_PHOSPHO" : \
            tuple(set(["OE1", "SD", "OE2", "OC2", "OD1", "NE2", \
                       "OD2", "OH2", "OW", "OG1", "OH", "SG", \
                       "ND1", "O", "OC1", "OG", "O1P", "O2P", "OJ"]))}


class AnalysisFunctions:
    """Class of custom functions for analyzing water bridges.

    Notes
    -----
    All functions in the public API return `None` and 
    have three parameters:
    
    current : `list`
        List of hydrogen bonds from selection 1 to selection 2.
    output : `dict`
        A dictionary which is modified in-place where the key
        is the type of water bridge and the value is the weight 
        of this type of water bridge.
    u : `MDAnalysis.universe`
        The current Universe for looking up atoms.
    """

    @staticmethod
    def _get_wb_attributes(u, current):
        """Get the atom-, residue- and segment-level attributes
        of the two selections for the current water bridge.
        """

        # decompose the first hydrogen bond of the water bridge
        sele1_index, sele1_heavy_index, atom2, \
            heavy_atom2, dist, angle = current[0]
        # decompose the last hydrogen bond of the water bridge
        atom1, heavy_atom1, sele2_index, \
            sele2_heavy_index, dist, angle = current[-1]
        # expand the atom index to the resname, resid, atom names
        sele1 = u.atoms[sele1_index]
        sele2 = u.atoms[sele2_index]
        # get the order of the current water bridge
        order_of_wb = len(current) - 1
        # return water bridge attributes
        return ((sele1.segid, sele1.resname, sele1.resid, \
                 sele1.name, sele1_index), \
                (sele2.segid, sele2.resname, sele2.resid, \
                 sele2.name, sele2_index), \
                order_of_wb)

    @staticmethod
    def wb_per_atom(current, output, u):
        """Count the water bridges per atom, differentiating 
        them according to their order.
        """

        # retrieve attributes for the current water bridge
        s1, s2, order_of_wb = \
            AnalysisFunctions._get_wb_attributes(u, current)
        s1_segid, s1_resname, s1_resid, s1_name, s1_index = s1
        s2_segid, s2_resname, s2_resid, s2_name, s2_index = s2
        # atom attributes are included in the key
        key = (s1_segid, s1_resname, s1_resid, s1_name, s1_index, \
               s2_segid, s2_resname, s2_resid, s2_name, s2_index, \
               order_of_wb)
        # Update the count
        output[key] += 1


    @staticmethod
    def wb_per_residue(current, output, u):
        """Count the water bridges per residue, differentiating 
        them according to their order (multiple water bridges 
        of the same order between two residues are counted as 
        one per frame).
        """

        # retrieve attributes for the current water bridge
        s1, s2, order_of_wb = \
            AnalysisFunctions._get_wb_attributes(u, current)
        s1_segid, s1_resname, s1_resid, s1_name, s1_index = s1
        s2_segid, s2_resname, s2_resid, s2_name, s2_index = s2 
        # atom-level attributes are not included in the key, since
        # we are only interested in counting per-residue; the order
        # of the water bridge is included in the key
        key = (s1_segid, s1_resname, s1_resid, \
               s2_segid, s2_resname, s2_resid, \
               order_of_wb)
        # Update the count
        output[key] += 1


if __name__ == "__main__":

    # create an argument parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    f_helpstr = "Name of/path to the input trajectory."
    parser.add_argument("-f", "--traj", \
                        dest = "traj", \
                        type = str, \
                        required = True, \
                        help = f_helpstr)

    s_helpstr = "Name of/path to the input topology."
    parser.add_argument("-s", "--top", \
                        dest = "top", \
                        type = str, \
                        required = True, \
                        help = s_helpstr)

    o_helpstr = "Name of/path to the output file."
    parser.add_argument("-o", "--output-file", \
                        dest = "output_file", \
                        type = str, \
                        required = True, \
                        help = o_helpstr)

    t_helpstr = \
        "How to analyze the water bridges found. " \
        "Choices are: {:s} (default: {:s}). "
    t_choices = ["atom", "residue"]
    t_default = "atom"
    parser.add_argument("-t", "--analysis-type", \
                        dest = "analysis_type", \
                        type = str, \
                        choices = t_choices, \
                        required = False, \
                        default = t_default, \
                        help = t_helpstr.format(", ".join(t_choices), \
                                                t_default))

    s1_helpstr = \
        "First selection (MDAnalysis syntax with underscores " \
        "instead of spaces)."
    parser.add_argument("-s1", "--selection1", \
                        dest = "selection1", \
                        type = str, \
                        required = True, \
                        help = s1_helpstr)

    s2_helpstr = \
        "Second selection (MDAnalysis syntax with underscores " \
        "instead of spaces)."
    parser.add_argument("-s2", "--selection2", \
                        dest = "selection2", \
                        type = str, \
                        required = True, \
                        help = s2_helpstr)

    sw_helpstr = \
        "Water selection (MDAnalysis syntax with underscores " \
        "instead of spaces) (default: {:s})."
    sw_default = "resname_SOL"
    parser.add_argument("-sw", "--water-selection", \
                        dest = "water_selection", \
                        type = str, \
                        required = False, \
                        default = sw_default, \
                        help = sw_helpstr.format(sw_default))

    forcefield_helpstr = \
        "Force field from which to take atom names. Choices " \
        "are {:s} (default: {:s})."
    forcefield_choices = ["CHARMM27", "GLYCAM06", "CHARMM22ST_PHOSPHO"]
    forcefield_default = "CHARMM27"
    parser.add_argument("--forcefield", \
                        dest = "forcefield", \
                        type = str, \
                        default = forcefield_default, \
                        required = False, \
                        help = forcefield_helpstr.format(\
                                ", ".join(forcefield_choices), \
                                forcefield_default))

    order_default = 1
    order_helpstr = \
        "Maximum water bridge order (default: {:d})."
    parser.add_argument("--order", \
                        dest = "order", \
                        type = int, \
                        default = order_default, \
                        required = False, \
                        help = order_helpstr.format(order_default))

    # parse the arguments
    args = parser.parse_args()
    top = args.top
    traj = args.traj
    output_file = args.output_file
    analysis_type = args.analysis_type
    selection1 = args.selection1.replace("_", " ")
    selection2 = args.selection2.replace("_", " ")
    water_selection = args.water_selection.replace("_", " ")
    forcefield = args.forcefield
    order = args.order
    
    # create the Universe
    u = mda.Universe(top, traj)
    
    # set up the water bridges analyis
    analysis = \
        WaterBridgeAnalysis_Custom(universe = u, \
                                   selection1 = selection1, \
                                   selection2 = selection2, \
                                   water_selection = water_selection, \
                                   forcefield = forcefield, \
                                   order = order)
    
    # run the analysis
    analysis.run()
    
    # analyze the water bridges found
    if analysis_type == "atom":
        # default counting, each water bridge treated
        # separately.
        analysis_func = AnalysisFunctions.wb_per_atom
        cols = \
            ["s1_segid", "s1_resname", "s1_resid", \
             "s1_name", "s1_index", \
             "s2_segid", "s2_resname", "s2_resid", \
             "s2_name", "s2_index", \
             "order_of_wb", "persistence"]
    
    elif analysis_type == "residue":
        # water bridges per each pair of residues, 
        # differentiated by order
        analysis_func = AnalysisFunctions.wb_per_residue_order
        cols = \
            ["s1_segid", "s1_resname", "s1_resid", \
             "s2_segid", "s2_resname", "s2_resid", \
             "order_of_wb", "persistence"]

    # count the water bridges by type according to the selected
    # criteria
    wb_count = analysis.count_by_type(analysis_func = analysis_func)  
    
    # convert each item of the list into a flattened tuple (just
    # leave the tuple as it is if the first item is not a tuple, as
    # it happens when the default counting is used).
    wb_count_flat = [(*item[0], *item[1:]) for item in wb_count]
    
    # convert the output to a dataframe
    df = pd.DataFrame(data = wb_count_flat, \
                      columns = cols)
    
    # sort the water bridges by increasing order and decreasing
    # persistence (more persistent water bridges will come first)
    sort_keys = ["order_of_wb", "persistence"]
    sort_ascending = [True, False]
    df.sort_values(by = sort_keys, \
                   ascending = sort_ascending, \
                   inplace = True)
    
    # save the output to a CSV file
    df.to_csv(output_file, \
              sep = ",", \
              float_format = "%.5f", \
              index = False)
