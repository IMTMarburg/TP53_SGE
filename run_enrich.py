#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
run_enrich.py: Contains a script to run an enrichment analysis using Enrich2,
written by  Alan F Rubin ORCID_icon http://orcid.org/0000-0003-1474-605X.

The software is freely available from https://github.com/FowlerLab/Enrich2/
under the GPLv3 license.
"""
from pathlib import Path
from typing import Optional
import pandas as pd

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"

import pandas as pd
import subprocess
import argparse
import json
from pathlib import Path
from typing import Optional


def create_count_files_for_enrich2(out: Path, count_file: Path):
    """
    create_count_files_for_enrich2 writes the raw count files for the enrich2 analyis
    """
    # create output directory if it does not exist
    out.mkdir(exist_ok=True, parents=True)
    # the replicates used in the analysis
    replicates = {
        "DMSO_1": "read_count_dmso_1",
        "DMSO_2": "read_count_dmso_2",
        "DMSO_3": "read_count_dmso_3",
        "N3a_1": "read_count_n3a_1",
        "N3a_2": "read_count_n3a_2",
        "N3a_3": "read_count_n3a_3",
    }
    # load the raw count
    df_rfs = pd.read_csv(count_file, dtype={"mut_ID": "str"})
    # write out the count files for enrich2, change the mut_ID to be non-numeric
    for replicate in replicates:
        out_file = out / f"identifier_count_{replicate}.tsv"
        if not out_file.exists():
            df = df_rfs[["mut_ID", replicates[replicate]]]
            df.loc[:, "mut_ID"] = df["mut_ID"].apply(lambda x: f"M{x}")
            df = df.rename(
                columns={"mut_ID": "identifier", replicates[replicate]: "count"}
            )
            df.to_csv(out / f"identifier_count_{replicate}.tsv", sep="\t", index=False)


def run_enrich(
    enrich_conda_python: str,
    config_file: str,
    path_to_enrich2: Path,
):
    """runs the enrich2 analysis"""
    subprocess.check_call(
        [
            str(enrich_conda_python),
            str(path_to_enrich2 / "enrich2" / "main.py"),
            config_file,
            "ratios",
            "full",
        ]
    )


def get_main_scores(path: Path):
    """Read and process the main score tsv file from enrich2"""
    df_ex5_in = pd.read_csv(path, sep="\t", header=[1]).loc[1:]
    df_ex5_in["value"] = df_ex5_in["value"].apply(lambda x: x[1:]).astype(int)
    df_ex5_in = df_ex5_in.rename(columns={"value": "mut_ID"}).set_index("mut_ID")
    return df_ex5_in


def get_main_shared_full(path: Path):
    """Read and process the main shared full tsv file from enrich2"""
    df_ex5_in_shared = pd.read_csv(path, sep="\t", header=[1, 2]).loc[1:]
    df_ex5_in_shared.columns = [
        f"{'_'.join(x.split('_')[1:])}_{y}"
        for (x, y) in df_ex5_in_shared.columns.values
    ]
    df_ex5_in_shared["mut_ID"] = (
        df_ex5_in_shared["_value"].apply(lambda x: x[1:]).astype(int)
    )
    del df_ex5_in_shared["_value"]
    df_ex5_in_shared = df_ex5_in_shared.set_index("mut_ID")
    return df_ex5_in_shared


def get_output_directory_from_config(config_file: str) -> Path:
    with open(config_file, "r") as file:
        data = json.load(file)
    return Path("../" + data["output directory"])


def consolidate_results(
    count_file: Path, config_file: str, outfolder: Optional[Path] = None
):
    enrich_output_folder = get_output_directory_from_config(config_file)

    # read in the input file
    df_rfs = pd.read_csv(count_file)
    # read the main scores from enrich2 results
    main_scores = (
        enrich_output_folder
        / "tsv"
        / f"{enrich_output_folder.name}_exp"
        / "main_identifiers_scores.tsv"
    )
    df_main_scores = get_main_scores(main_scores)
    # read the main shared full scores from enrich2 results
    main_scores_shared_full = (
        enrich_output_folder
        / "tsv"
        / f"{enrich_output_folder.name}_exp"
        / "main_identifiers_scores_shared_full.tsv"
    )
    df_shared_full = get_main_shared_full(main_scores_shared_full)
    # consolidate the results
    df_all_in_one = df_rfs.copy().set_index("mut_ID")
    df_all_in_one = df_all_in_one.join(df_main_scores)
    df_all_in_one = df_all_in_one.join(df_shared_full)
    df_all_in_one = df_all_in_one.reset_index()
    outpath = outfolder
    if outfolder is None:
        outpath = enrich_output_folder
    outpath.mkdir(exist_ok=True, parents=True)
    df_all_in_one.to_csv(
        outpath / f"{enrich_output_folder.name}_scores.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--enrich_conda_python",
        type=Path,
        default="/talizorah/opt/environments/anaconda/envs/enrich2/bin/python",
        help="Path to the python binary in the conda environment",
    )
    parser.add_argument(
        "--config_file",
        type=str,
        default="enrich2_files/Exon5678_identifier",
        help="The configuration file for Enrich2",
    )
    parser.add_argument(
        "--path_to_enrich2",
        type=Path,
        default=Path("../code/Enrich2"),
        help="The path to the Enrich2 github code",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("enrich2_files"),
        help="The path to the output directory for Enrich2 count files",
    )
    parser.add_argument(
        "--count_file",
        type=Path,
        default=Path("out/Exon5-8_RFS.csv"),
        help="The path to the raw count file for all replicates",
    )
    args = parser.parse_args()
    create_count_files_for_enrich2(args.out, args.count_file)
    run_enrich(args.enrich_conda_python, args.config_file, args.path_to_enrich2)
    consolidate_results(args.count_file, args.config_file, args.out)
