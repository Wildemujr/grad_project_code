#!/usr/bin/env python3
# coding=utf8


import numpy as np
import pandas as pd
import argparse
from pprint import pprint as pp



if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="rmsf_pymol_colorizer.py")
    parser.add_argument(
        "-c",
        "--color",
        dest="color",
        help="Isolating residues from the RMSF .dat file for coloring.",
    )
    parser.add_argument(
        "-s",
        "--save",
        dest="save_file",
        help="Select the name that you would like to save pymol session as."
    )

    options = parser.parse_args()

    if options.color:
        b_data = np.genfromtxt(f"{options.color}")
        df = pd.DataFrame(b_data)
        x = pd.DataFrame([i for i in range(0, len(b_data))], columns=["Residue Index"])
        y = df.iloc[:, 1:2]
        y.columns = ["Atomic Fluctuations"]

        df = pd.merge(x, y, left_index=True, right_index=True)

        """
        I've decided that, after looking at the RMSF plot, that this should be
        the cutoff values: 

        0.8 -> inf == red
        0.6 - 0.8 -> color == green
        0.0 -> 0.6 -> color == blue
        """

        s = pd.Series(np.linspace(0, 2, 20, endpoint=True))
        y_filtered_blue = y[(y["Atomic Fluctuations"] < 0.6)]
        y_filtered_green = y[
            (y["Atomic Fluctuations"] > 0.6) & (y["Atomic Fluctuations"] < 0.75)
        ]
        y_filtered_red = y[(y["Atomic Fluctuations"] > 0.75)]

        ### For the blue color ###
        y_filtered_blue_string = pd.DataFrame.to_string(
            y_filtered_blue, columns=["Atomic Fluctuations"]
        )
        y_filtered_blue_list = y_filtered_blue.index.tolist()
        y_filtered_blue_list = list(map(lambda x: str(x+1), y_filtered_blue_list))
        blue_residue_nums = "+".join(y_filtered_blue_list)
        f = open("y_filtered_blue.txt", "w+")
        f.write(y_filtered_blue_string)
        f.close()

        ### For the green color ###
        y_filtered_green_string = pd.DataFrame.to_string(
            y_filtered_green, columns=["Atomic Fluctuations"]
        )
        y_filtered_green_list = y_filtered_green.index.tolist()
        y_filtered_green_list = list(map(lambda x: str(x+1), y_filtered_green_list))
        green_residue_nums = "+".join(y_filtered_green_list)
        f = open("y_filtered_green.txt", "w+")
        f.write(y_filtered_green_string)
        f.close()

        ###  For the red color ###
        y_filtered_red_string = pd.DataFrame.to_string(
            y_filtered_red, columns=["Atomic Fluctuations"]
        )
        y_filtered_red_list = y_filtered_red.index.tolist()
        y_filtered_red_list = list(map(lambda x: str(x+1), y_filtered_red_list))
        red_residue_nums = "+".join(y_filtered_red_list)
        f = open("y_filtered_red.txt", "w+")
        f.write(y_filtered_red_string)
        f.close()

        # This generates the PyMol script. Load in PyMol via: @script_name.pml
        f = open("pyMol_rmsf_coloring.pml", "w+")
        f.write("fetch 1fdn \n")
        f.write("remove solvent \n")
        f.write("set seq_view, 1 \n")
        f.write(f"select nterm, resi {blue_residue_nums} \n")
        f.write("color cyan, nterm \n")
        f.write(f"select nterm, resi {green_residue_nums} \n")
        f.write("color green, nterm \n")
        f.write(f"select nterm, resi {red_residue_nums} \n")
        f.write("color red, nterm \n")
        f.write(f"save {options.save}.pse")
        f.close()

        ### Here if anyone wants output to the command line

        # print("\n")
        # print("Copy the following into pymol to color these residues cyan:", "\n")
        # print(f"select nterm, resi {blue_residue_nums}")
        # print("color cyan, nterm", "\n")

        # print("\n")
        # print("Copy the following into pymol to color these residues green:", "\n")
        # print(f"select nterm, resi {green_residue_nums}")
        # print("color green, nterm", "\n")

        # print("\n")
        # print("Copy the following into pymol to color these residues red:", "\n")
        # print(f"select nterm, resi {red_residue_nums}")
        # print("color red, nterm", "\n")
