import pypipegraph2 as ppg

ppg.replace_ppg1()
import pandas as pd
import re
import mbf
import subprocess
import sys
import katha_tools
import mmdemultiplex
import counting_sequences
from typing import List, Dict, Callable
from pypipegraph import Job
from pathlib import Path
from counting_sequences import ngmerge

########################
# instantiate pipegraph
########################

ppg.new()


#######################
# initialize variables
#######################

incoming_dir = "incoming"
result_dir = "results"
info_data_file = f"{incoming_dir}/NovaSeq_incoming_sequences_Exon5678.xlsx"
barcode_sheet = "Übersicht"
barcode_df = pd.read_excel(info_data_file, sheet_name=barcode_sheet, skiprows=2)
barcode_df = barcode_df.fillna(method="ffill", axis=0)
pool_col_name = "Library Name"
sample_col_name = "Probenname"
fw_col_name = "Barcode fw"  # the
rv_col_name = "Barcode rev"
trim_start = "trim after start"  # we want to trim the constant part at the start
trim_end = "trim before end"  # we want to trim the constant part at the end
summary_output_file = results / "summary_counts.tsv"
# The exons to be analyzed. Keys are file prefixes of the fastq files, values
# are the names of the sequencing libraries.
exons = {
    "Exon5_DMSO": "Exon 5",
    "Exon5_N3a": "Exon 5",
    "Exon6_DMSO": "Exon 6",
    "Exon6_N3a": "Exon 6",
    "Exon8_DMSO": "Exon 8",
    "Exon8_N3a": "Exon 8",
    "cDNA_Exon5": "Exon 5 cDNA",
    "cDNA_Exon6": "Exon 6 cDNA",
    "cDNA_Exon8": "Exon 8 cDNA",
}

################################
# declare the jobs for pipegraph
################################

dependencies_for_summary = []
per_sample_count_files = {}

# analyze each exon separately
for exon in exons:
    # demultiplex each exon library into replicates
    raw_pool = mmdemultiplex.samples.DemultiplexInputSample(
        exon,
        input_strategy=mbf.align.strategies.FASTQsFromPrefix(
            f"{incoming_dir}/fastqs/{exon}_"
        ),
        reverse_reads=False,
        pairing="paired",
    )
    demultiplexer = mmdemultiplex.demultiplex.Demultiplexer(
        raw_pool,
        katha_tools.tables.get_barcode_df(
            exon,
            barcode_df,
            pool_col_name,
            sample_col_name,
            fw_col_name,
            rv_col_name,
            trim_start,
            trim_end,
        ),
        output_folder=Path(f"{result_dir}/demultiplexing"),
    )
    demultiplexer.do_demultiplex()

    # get a dict of demultiplexed samples/replicates from the demultiplexer
    samples = demultiplexer.get_samples()
    sequence_counter = counting_sequences.SequenceCounter(
        sequence_file_path=Path(info_data_file),  # file with predefined sequences
        name=f"count_{exon}",
        seqs_to_trim_reads=None,
        seqs_to_trim_predefined=None,
        trimmed_length=None,
        result_folder=f"{result_dir}/seq_counts/{exon}",
        sequence_df_filter=None,
        sheet_name=exons[exon],
    )  # Initialize the sequence counter for the exon

    # for each sample/replicate
    for sample_name in samples:
        #  merge overlapping R1 and R2 reads into a single read and create a new merged single-end fastq file
        sample = samples[sample_name]
        sample.save_input()
        r1, r2 = sample.get_aligner_input_filenames()
        merged_path = f"{result_dir}/merged"
        merged_fastq = Path(f"{merged_path}/{sample_name}_merged.fastq.gz")
        merge_job = ngmerge(
            merged_fastq,
            f"{merged_path}/{sample_name}_merged_failed",
            r1,
            r2,
            dependencies=[sample.prepare_input()],
            options={
                "-l": f"{merged_path}/{sample_name}_merged.log",
                "-z": "",
                "-d": "",
            },
        )
        merged_sample = mbf.align.raw.Sample(
            f"{sample_name}_merged",
            input_strategy=mbf.align.strategies.FASTQsFromJob(merge_job),
            reverse_reads=False,
            pairing="single",
        )

        # count the number of reads remaining in the merged fastq file
        counting_sequences.get_reads_for_lanes_df(
            Path(result_dir)
            / "read_counts"
            / f"{sample_name}_read_counts_after_merge.tsv",
            {f"{sample_name}_merged": merged_sample},
            dependencies=[merged_sample.prepare_input()],
        )
        r1, *_ = merged_sample.get_aligner_input_filenames()

        # count the most common sequences
        counting_sequences.CutadaptMatch.count_most_common_sequences(
            r1=r1,
            r2=None,
            output_file=Path(result_dir)
            / "most_common_sequences"
            / "merged"
            / f"{sample_name}_merged_most_common.tsv",
            max=100000,
            dependencies=[merged_sample.prepare_input()],
        )

        # count the number of occurences of each sequence in the merged reads
        # by exact matching
        job = sequence_counter.write_count_table(merged_sample)
        per_sample_count_files[sample_name] = job.filenames[0]
        dependencies_for_summary.append(job)
# combine all the count tables into one
sequence_counter.combine(
    per_sample_count_files, summary_output_file, dependencies_for_summary
)


###################
# run the pipegraph
###################

ppg.run()
