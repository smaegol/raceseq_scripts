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
    description='Analyze repatmasker output file - identify soft clipped fragments (R3 reads)')

parser.add_argument('--input', dest='inputfile', action='store',
                    help='Input file from Parsing-RepeatMasker-Outputs/parseRM.pl (required)', required=True)
parser.add_argument('--fasta', dest='fastafile', action='store',
                    help='Fasta file with sequences provided as an input for repeatmasker (required)', required=True)
parser.add_argument('--output', dest='output', action='store',
                    help='Output tsv file (required)', required=True)
args = parser.parse_args()

from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

data = pd.DataFrame.from_csv(args.inputfile, sep='\t')

sequences = SeqIO.index(args.fastafile, "fasta")

seq_records = []
seq_names = []


for row in data.itertuples():
    LINE_name = row.Rname
    overlap = row.IfOverlap
    seq_name = row.Gname
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
        LINE_seq_obj = Seq(str(LINE_seq))
        # reverse complement the LINE sequence
        LINE_seq = str(LINE_seq_obj.reverse_complement())
        clip_seq = whole_seq[:LINE_start_in_seq - 1]
        clip_seq_obj = Seq(str(clip_seq))
        # reverse complement clipped fragment
        clip_seq = str(clip_seq_obj.reverse_complement())

        # produce output sequence:
        seq_descr = "\tclip5: \tclip3: " + clip_seq + "\tpos: " + \
            str(LINE_end_in_ref) + "\tref: " + str(LINE_name)
        processed_seq = SeqRecord(
            Seq(str(LINE_seq)), id=whole_seq_id, description=str(seq_descr))
        seq_records.append(processed_seq)


#include sequences for which LINE sequences were not found
for record in SeqIO.parse(args.fastafile, "fasta"):
    seq_id = record.id  # get id of read
    if (seq_id in seq_names):
        next
    else:
        whole_seq = record.seq
        clip_seq = str(whole_seq.reverse_complement())
        seq_descr = "\tclip5: \tclip3: \tpos: -1\tref: -1"
        processed_seq = SeqRecord(
            Seq(str(clip_seq)), id=seq_id, description=str(seq_descr))
        seq_records.append(processed_seq)


# write fasta output file
SeqIO.write(seq_records, args.output, "fasta")


### END ###
