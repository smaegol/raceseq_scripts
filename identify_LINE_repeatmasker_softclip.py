#!/usr/bin/env python

#######################################################################################
###                                                                                 ###
###     Copyright (C) 2017  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)                  ###
###                                                                                 ###
###     This program is free software: you can redistribute it and/or modify        ###
###     it under the terms of the GNU General Public License as published by        ###
###     the Free Software Foundation, either version 3 of the License, or           ###
###     (at your option) any later version.                                         ###
###                                                                                 ###
###     This program is distributed in the hope that it will be useful,             ###
###     but WITHOUT ANY WARRANTY; without even the implied warranty of              ###
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               ###
###     GNU General Public License for more details.                                ###
###                                                                                 ###
###     You should have received a copy of the GNU General Public License           ###
###     along with this program. If not, see <http://www.gnu.org/licenses/>.        ###
###                                                                                 ###
#######################################################################################

import argparse

# parse command line arguments
parser = argparse.ArgumentParser(
    description='Analyze repatmasker output file - identify soft clipped fragments')

parser.add_argument('--input', dest='inputfile', action='store',
                    help='Input file from Parsing-RepeatMasker-Outputs/parseRM.pl (required)', required=True)
parser.add_argument('--fasta', dest='fastafile', action='store',
                    help='Fasta file with sequences provided as an input for repeatmasker (required)', required=True)
parser.add_argument('--output', dest='output', action='store',
                    help='Output tsv file (required)', required=True)
args = parser.parse_args()


from Bio import SeqIO
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# read parseRM output file:
data = pd.DataFrame.from_csv(args.inputfile, sep='\t')

# index sequences file:
sequences = SeqIO.index(args.fastafile, "fasta")


regex_for_R5_CTGAC = r"^(?P<possible_tail>.*?CTGAC)(?P<adapter_15N_etc>.*)"
regex_for_LINE_CTGAC = r"^(?P<LINE_seq>.*?)(?P<clip_seq>CTGAC.*)"


seq_records = []
seq_names = []

for row in data.itertuples():
    LINE_name = row.Rname  # get reference name
    overlap = row.IfOverlap  # get overlap
    seq_name = row.Gname  # get name of the sequence
    # check if sequence was already processed
    # to choose only the best (first) match of repeatmasker for further analyses
    if (seq_name in seq_names):
        print("sequence already processed: " + seq_name)
    else:
        seq_names.append(seq_name)  # append seq name to processed seq_names
        # get LINE positions in analyzed sequence:
        LINE_start_in_seq = row.Gstart
        LINE_end_in_seq = row.Gend
        LINE_left_in_seq = row.Gleft
        LINE_start_in_ref = row.Rstart
        LINE_end_in_ref = row.Rend
        LINE_left_in_ref = row.Rleft
        # get the full sequence
        record2 = sequences[str(seq_name)]
        whole_seq = record2.seq
        whole_seq_id = record2.id
        # get the LINE seq based on coordinates from repeatmasker
        LINE_seq = whole_seq[LINE_start_in_seq - 1:LINE_end_in_seq]
        # get clipped fragments:
        clip_seq = whole_seq[LINE_end_in_seq:]
        clip5_seq = whole_seq[0:LINE_start_in_seq]

        # look for broken CTGAC tailseq delimiter sequence:
        merged_seq = LINE_seq[-4:] + clip_seq[0:4]
        match_CTGAC_R5_merged_seq = re.search(
            regex_for_R5_CTGAC, str(merged_seq))

        # if CTGAC was broken:
        if(match_CTGAC_R5_merged_seq):
            clipped_R5 = match_CTGAC_R5_merged_seq.group("possible_tail")
            length_matched = len(clipped_R5) - 5
            matched_CTGAC = clipped_R5[length_matched:4]
            LINE_seq = LINE_seq[:-len(matched_CTGAC)]
            LINE_end_in_seq = LINE_end_in_seq - len(matched_CTGAC)
            LINE_end_in_ref = LINE_end_in_ref - len(matched_CTGAC)
            clip_seq = matched_CTGAC + clip_seq

        match_CTGAC_R5_LINE = re.search(regex_for_LINE_CTGAC, str(LINE_seq))

        # look for CTGAC inside LINE sequence (possible for highly heteroegenous LINE ends, without tailing)
        # if present - correct identified LINE
        if(match_CTGAC_R5_LINE):
            true_LINE_seq = match_CTGAC_R5_LINE.group("LINE_seq")
            additional_clip_seq = match_CTGAC_R5_LINE.group("clip_seq")
            additional_clip_length = len(LINE_seq) - len(true_LINE_seq)
            LINE_end_in_seq = LINE_end_in_seq - additional_clip_length
            LINE_end_in_ref = LINE_end_in_ref - additional_clip_length
            LINE_seq = true_LINE_seq
            clip_seq = additional_clip_seq + clip_seq

        # produce output sequence:
        seq_descr = "\tclip5: " + str(clip5_seq) + "\tclip3: " + clip_seq + \
            "\tpos: " + str(LINE_end_in_ref) + "\tref: " + str(LINE_name)
        processed_seq = SeqRecord(
            Seq(str(LINE_seq)), id=whole_seq_id, description=str(seq_descr))
        seq_records.append(processed_seq)


# include sequences for which LINE sequences were not found
for record in SeqIO.parse(args.fastafile, "fasta"):
    seq_id = record.id  # get id of read
    if (seq_id in seq_names):
        next
    else:
        whole_seq = record.seq
        #clip_seq = str(whole_seq.reverse_complement())
        seq_descr = "\tclip5: \tclip3: \tpos: -1\tref: -1"
        processed_seq = SeqRecord(
            Seq(str(whole_seq)), id=seq_id, description=str(seq_descr))
        seq_records.append(processed_seq)

# write fasta output file
SeqIO.write(seq_records, args.output, "fasta")


### END ###
