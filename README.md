# analyze_water_bridges.py

This script provides a convenient command-line tool to perform water bridges analysis between two atom selections in a molecular dynamics simulations, taking advantage of the water bridges analysis functionalities implemented in the MDAnalysis package.

For a detailed description of how water bridges analysis is performed in MDAnalysis, please visit [this page](https://www.mdanalysis.org/docs/documentation_pages/analysis/wbridge_analysis.html).

## Requirements

`python3.7`

`MDAnalysis`

`pandas`

The script has been tested with python3.7.4, MDAnalysis version 0.20.1 and pandas version 1.0.1.

## Installation

No installation is required, just put the script into your working directory or into another location of your choice and run it from your working directory.

## Usage

The script should be run as:

`python3.7 analyze_water_bridges.py [-h] -f TRAJ -s TOP -o OUTPUT_FILE [-t {atom,residue}] -s1 SELECTION1 -s2 SELECTION2 [-sw WATER_SELECTION] [--forcefield FORCEFIELD] [--order ORDER]`

Available arguments and options are:

* `-h` or `--help` to show the help message and exit.
* `-f` or `--traj` followed by the name of/path to the input trajectory.
* `-s` or `--top` followed by name of/path to the input topology.
* `-o` or `--output-file` followed by the name of/path to the CSV file where the results of the analysis will be stored as a data frame.
* `-t` or `--analysis-type` followed by either `atom`or `residue` according to the level at which you want the water bridges to be analyzed and reported. If `atom`, water bridges between each donor/acceptor pair found to form a water bridge at least once along the trajectory will be reported separately. If `residue`, water bridges will be counted between pairs of residues, and multiple water bridges between two residues will be counted as one when counting along the simulation.
* `-s1` or `--selection1` followed by the first atom selection, defined according to the MDAnalysis syntax, but with whitespaces replaced by underscores (e.g. `resname ALA` becomes `resname_ALA`).
* `-s2` or `--selection2` followed by the second atom selection, defined according to the MDAnalysis syntax, but with whitespaces replaced by underscores (e.g. `resname ALA` becomes `resname_ALA`).
* `-sw` or `--water-selection` followed by the water selection, defined according to the MDAnalysis syntax, but with whitespaces replaced by underscores (e.g. `resname SOL` becomes `resname_SOL`). The default is `resname_SOL`.
* `--forcefield` followed by the  name of the force field from which atom names should be taken. `CHARMM27`, `GLYCAM06` and `CHARMM22ST_PHOSPHO` are currently supported. The default is `CHARMM27`.
* `--order` followed by the maximum order for a water bridge to be reported. The order of a water bridge is the number of water molecules bridging two atoms. Water bridges of order 0 are simply hydrogen bonds between two atoms. The default is `1`. 
