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



import os
import sys
# srcipt path is required to find the location of files required for analysis (indexes and other scripts)
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

###        !!!!!!!!!!!! Set paths and options before you start !!!!!!!!!!!!!!!!!! 	###

bowtie2_path = "bowtie2"

bowtie_threads = '10'
get_sofclipped_script_path = script_path + "/get_softclipped_reads_from_sam.pl"
identify_LINEs_script_path = script_path + "/identify_LINE_repeatmasker.py"
transcript_genomes = {'GAPDH': script_path + '/indexes/GAPDH_noA',
                      'REPORTERL1KD': script_path + '/indexes/reporter_L1_sirna',
                      'REPORTERL1_overexp': script_path + '/indexes/reporter_L1_overexp',
                      'REPORTERL1': script_path + '/indexes/reporterL1',
                      'LEAP': script_path + '/indexes/LEAP',
                      'LEAP_AU': script_path + '/indexes/LEAP_AU',
                      'ACTB': '/home/smaegol/storage/analyses/tail_seq_3/genome/ActB',
                      'ETV4': '/home/smaegol/storage/analyses/tail_seq_3/genome/ETV4_noA',
                      'GAPDH': '/home/smaegol/storage/analyses/tail_seq_3/genome/GAPDH_noA',
                      'SOGA2': '/home/smaegol/storage/analyses/tail_seq_3/genome/SOGA2',
                      'PABPC4': '/home/smaegol/storage/analyses/tail_seq_3/genome/PABPC4_noA',
                      }


samplesheet_location = script_path + \
    '/flowcell2/flowcell2_analysis_samplesheet.csv'


###        !!!!!!!!!!!!                 					   !!!!!!!!!!!!!!!!!! 	###

import argparse

# parse command line arguments
parser = argparse.ArgumentParser(
    description='Analyze paired output of tailseeker to identify tails')

parser.add_argument('--inputdir', dest='inputdir',
                    action='store', help='Input dir(required)', required=True)
parser.add_argument('--output', dest='output', action='store',
                    help='Output tsv file (required)', required=True)
parser.add_argument('--glob', dest='glob', action='store',
                    help='Custom specification of files to analyze (optional)', required=False)
parser.add_argument('--samplesheet', dest='samplesheet', action='store',
                    help='Alternative samplesheet (optional)', required=False)
parser.add_argument('--window', dest='window', action='store',
                    help='Window size [nucleotides] for 3prime end terminal nucleotides analysis (optional, default = 7)', required=False, default=7)
parser.add_argument('--min_R5_length', dest='min_R5_length', action='store',
                    help='Minimum size of R5 read to include R5 clipping in tails prediction', required=False, default=100)

args = parser.parse_args()

from Bio import SeqIO
import re
import glob
import pandas as pd
from Bio.Seq import Seq
import subprocess
from Bio.SeqRecord import SeqRecord

# get samplesheet from command-line
if(args.samplesheet):
    samplesheet_location = args.samplesheet
# read samplesheet into pandas dataframe
data = pd.DataFrame.from_csv(samplesheet_location, sep='\t')


def get_3end_nucleotides(sequence, window_size):
    """Analyze terminal nucleotides."""
    terminal_nucleotides = sequence[-int(window_size):]
    return(terminal_nucleotides)


def analyze_tails(R1, R2, transcript, sample_name, localization, replicate, condition, cell_line, primer_name, person, project_name, exp_type, R5_sequenced_length):
    """Main processing function.

    Takes localization of softclipped fasta and description of different features of samples (provided in the samplesheet)
    Returns dict with data concerning all sequences analyzed.
    Output contains tails sequences, lengths and classification
    """
    # for L1 genomic sequences - use repeatmasker files
    if((transcript == 'ENDOL1')):
        R1 = R1 + '.rmasker.fasta'
        R2 = R2 + '.rmasker.fasta'

    # index R2 (R3 - tailseeker output) reads
    R2_reads = SeqIO.index(R2, "fasta")
    # dict storing results which will be saved in tsv file
    tails_results = {}
    final_results = {}

    # regex used in the further analyses
    regex_for_heuristic_tail_identification_R3 = r"(?P<tail>^A{2}A*[TGC]?A*[TGC]?A{2}A*T*|^A{2}A*[TGC]?A*T*|^A+T*|^T+|^T+[ACG]?T+|^T+A*)(?P<rest_of_clip>.*)"
    regex_for_R5_CTGAC = r"^(?P<possible_tail>.*CTGAC)(?P<adapter_15N_etc>.*)"
    regex_for_R5_CTGAC_mism = r"^(?P<possible_tail>.*(CTGAC){e<=1})(?P<adapter_15N_etc>.*)"
    regex_for_tail_R5_CTGAC_mism = r"^(?P<possible_tail>.*?)(^TGAC|.TGAC|C.GAC|CT.AC|CTG.C|CTGA.)(?P<other>.*)"
    regex_for_tail_R5_CTGAC = r"^(?P<possible_tail>.*?)(CTGAC)(?P<other>.*)"
    regex_for_genome_encoded_Atail = r"(?P<genome_encoded_tail>A+)$"
    regex_for_genome_encoded_Ttail = r"(?P<genome_encoded_tail>T+)$"
    regex_for_plasmid_seq = r"^(?P<plasmid_seq>GGGGTGGGCG.*)"
    regex_for_heterogenous_end_tail = r"^(?P<genomic_fragment>.*?)(?P<tail>AA+|TT+|AA+T+)"
    regex_for_A_only_in_seq = r"^A*$"

    #read all fasta records from R1(R5 - tailseeker) file
    for record in SeqIO.parse(R1, "fasta"):

        seq_id = record.id  # get id of read
        # create dict for storing temp results for pair
        tails_results[seq_id] = {}
        tails_results[seq_id]['CTGAC_R5'] = 0  # set the initial value to 0


        tails_results[seq_id]['heterogenous_end'] = ''
        tails_results[seq_id]['heterogenous_end_R3'] = ''
        # create dict for storing final results for pair
        final_results[seq_id] = {}
        R5_seq = record.seq  # get seq of R5 read (after clipping)

        # check if mate is present in the R2 reads file (can be absent in case of rmasker output)
        if(str(seq_id) in R2_reads):
            record2 = R2_reads[str(seq_id)]  # get R5 read from rmasker output
            R3_seq = record2.seq  # get seq of R3 read (after clipping)
        else:
            # if mate is absent - create the empty read for compatibility with further processing steps
            record2 = SeqRecord(Seq(
                ''), description=">a0000:00000000:0000:0:0:\tclip5: \tclip3: \tpos: -1\tref: -1")
            R3_seq = ''

        #get terminal nucleotides information:
        terminal_nucleotides=get_3end_nucleotides(R3_seq, args.window)
        tails_results[seq_id]['terminal_nucleotides'] = terminal_nucleotides

        # get all the information from the fasta header of both files (including tailseeker tail, softclipping, mapping position, template)
        regex_match_for_R1 = re.search(
            '(?P<tile>.{5})\:(?P<position>\d{8})\:(?P<tailseq_score>\d{4})\:(?P<PCRduplicates>\d+)\:(?P<Atail_length>.*)\:(?P<additional_bases>.*)\tclip5: (?P<clip5>[ACGTN]*)\tclip3: (?P<clip3>[ACGTNacgtn]*)\tpos: (?P<pos>.*)\tref: (?P<ref>.*)', record.description)
        regex_match_for_R2 = re.search(
            '(?P<tile>.{5})\:(?P<position>\d{8})\:(?P<tailseq_score>\d{4})\:(?P<PCRduplicates>\d+)\:(?P<Atail_length>.*)\:(?P<additional_bases>.*)\tclip5: (?P<clip5>[ACGTN]*)\tclip3: (?P<clip3>[ACGTNacgtn]*)\tpos: (?P<pos>.*)\tref: (?P<ref>.*)', record2.description)

        #get sequences of 3' clipping (potential tails)
        if (regex_match_for_R1):
            clipped_R5 = regex_match_for_R1.group('clip3')
        else:
            clipped_R5 = ""
        if (regex_match_for_R2):
            clipped_R3 = regex_match_for_R2.group('clip3')
        else:
            clipped_R3 = ""

        #get lengths of clipping:
        clip3_R5_length = len(clipped_R5)
        clip3_R3_length = len(clipped_R3)

        tails_results[seq_id]["original_clipped_R5"]=clipped_R5
        tails_results[seq_id]["original_clipped_R3"]=clipped_R3

        # get tailseq predictions from read id:
        A_tail_length = regex_match_for_R1.group('Atail_length')
        T_tail_length = len(regex_match_for_R1.group('additional_bases')) - 1
        additional_bases = regex_match_for_R1.group('additional_bases')
        ref_name_R5 = regex_match_for_R1.group('ref')
        ref_name_R3 = regex_match_for_R2.group('ref')
        additional_bases = re.sub(r'\s', '', additional_bases)
        tailseq_tail_length = int(A_tail_length) + T_tail_length
        PCRduplicates = regex_match_for_R1.group('PCRduplicates')
        tailseq_score = regex_match_for_R1.group('tailseq_score')

        tailseq_score_bitwise = bin(int(tailseq_score))

        tailseq_polyA_detected = 0
        tailseq_delimiter_mismatch = 0
        tailseq_delimiter_shifted = 0
        tailseq_delimiter_not_found = 0

        if (len(tailseq_score_bitwise)-2>=1):
            tailseq_polyA_detected = int(tailseq_score_bitwise[-1])
        if (len(tailseq_score_bitwise)-2>=2):
            tailseq_delimiter_mismatch = int(tailseq_score_bitwise[-2])
        if (len(tailseq_score_bitwise)-2>=6):
            tailseq_delimiter_shifted = int(tailseq_score_bitwise[-6])
        if (len(tailseq_score_bitwise)-2>=8):
            tailseq_delimiter_not_found = int(tailseq_score_bitwise[-8])


        R5_mapping_pos = regex_match_for_R1.group('pos')
        R3_mapping_pos = regex_match_for_R2.group('pos')

        # if R5 read sequence is composed only of A - discard
        if (re.search(regex_for_A_only_in_seq, str(R5_seq))):
            #print("found: " + seq_id)
            R5_seq = ''
            clipped_R5 = ''
            R5_mapping_pos = "-1"

        # identify incorrectly predicted ends - if CTGAC spans mapping site
        # can occur if the end of transcript is similar to the tailseq delimiter sequence(CTGAC)
        tails_results[seq_id]["mapping_spanning_delimiter"]=0
        merged_seq = R5_seq[-4:] + clipped_R5[0:4]
        match_CTGAC_R5_merged_seq = re.search(regex_for_R5_CTGAC, str(merged_seq))
        if(match_CTGAC_R5_merged_seq):
            clipped_from_R5_seq = match_CTGAC_R5_merged_seq.group("possible_tail")
            length_matched = len(clipped_from_R5_seq) - 5
            matched_CTGAC = clipped_from_R5_seq[length_matched:4]
            R5_seq = R5_seq[:-len(matched_CTGAC)]
            R5_mapping_pos = int(R5_mapping_pos) - len(matched_CTGAC)
            clipped_R5 = matched_CTGAC + clipped_R5
            tails_results[seq_id]["mapping_spanning_delimiter"]=1


        # check for the presence of CTGAC in clipped R5 fragments
        # if yes - remove CTGAC and following bases, leave clipped fragment and correct its length
        tails_results[seq_id]['CTGAC_R5_mismatched']=0
        tails_results[seq_id]['heterogenous_end_length']=0
        #match_CTGAC_R5 = ''
        #print("pre",clipped_R5)
        match_CTGAC_R5 = re.search(regex_for_tail_R5_CTGAC_mism, clipped_R5)
        if ((tailseq_delimiter_not_found==0) & (tailseq_delimiter_mismatch==0)):
            match_CTGAC_R5 = re.search(regex_for_tail_R5_CTGAC_mism, clipped_R5)
        elif ((tailseq_delimiter_not_found==0) & (tailseq_delimiter_mismatch>0)):
            match_CTGAC_R5 = re.search(regex_for_tail_R5_CTGAC_mism, clipped_R5)
        #    print("mismatch")
        else:
            match_CTGAC_R5 = re.search(regex_for_tail_R5_CTGAC, clipped_R5)
        #print(match_CTGAC_R5)
        if(match_CTGAC_R5):
            if (match_CTGAC_R5.group("possible_tail")==''):
                clipped_R5 = ''
            else:
                clipped_R5 = match_CTGAC_R5.group("possible_tail")
            # mark the presence of CTGAC in clipped fragment (prerequisite of analysis)
            #print("mm",clipped_R5)
            tails_results[seq_id]['CTGAC_R5'] = 1
            if(tailseq_delimiter_mismatch>0):
                tails_results[seq_id]['CTGAC_R5_mismatched']=1
            clip3_R5_length = len(clipped_R5)
#################### TO BE REWRITTEN #######################33
            # check if clipped fragment contains other stuff than A/U/AU tails (possible heterogenity of 3'end of LINE1)
            match_heterogenous_end_tail = re.search(
                regex_for_heterogenous_end_tail, clipped_R5)
            if(match_heterogenous_end_tail):
                # if such heterogenity was found - append possible genomic fragment to the mapped sequence and levae only identified tail
                heterogenous_end = match_heterogenous_end_tail.group(
                    "genomic_fragment")
                number_heterogenous_nucleotides = len(heterogenous_end)
                tails_results[seq_id]['heterogenous_end'] = heterogenous_end
                tails_results[seq_id]['heterogenous_end_length'] = number_heterogenous_nucleotides
                #clipped_R5 = match_heterogenous_end_tail.group("tail")
                #print("het",clipped_R5)
                #clip3_R5_length = len(clipped_R5)
                #R5_seq = R5_seq + heterogenous_end
                #R5_mapping_pos = int(R5_mapping_pos) + \
                ##    number_heterogenous_nucleotides
        ##else:
        #    print("nn",clipped_R5)
        #print("aa",clipped_R5)

        # if R3 read got soft-clipped
        tails_results[seq_id]['heterogenous_end_R3_length']=0
        if(clip3_R3_length > 0):
            # check if clipped fragment contains other stuff than A/U/AU tails (possible heterogenity of 3'end of LINE1)
            match_heterogenous_end_tail_R3 = re.search(
                regex_for_heterogenous_end_tail, clipped_R3)
            if(match_heterogenous_end_tail_R3):
                heterogenous_end = match_heterogenous_end_tail_R3.group(
                    "genomic_fragment")
                number_heterogenous_nucleotides = len(heterogenous_end)
                tails_results[seq_id]['heterogenous_end_R3'] = heterogenous_end
                tails_results[seq_id]['heterogenous_end_R3_length'] = number_heterogenous_nucleotides
                #clipped_R3 = match_heterogenous_end_tail_R3.group("tail")
                #R3_seq = R3_seq + heterogenous_end
                #R3_mapping_pos = int(R3_mapping_pos) + \
                #    number_heterogenous_nucleotides


##########################################################################

        # get last mapped nucleotides - will be further used to identify A or U nucleotides which should be included in the tail sequence but were also present in the reference
        if (R5_mapping_pos != "-1"):
            # if R5 read was mapped
            R5_last_mapped_nucleotide = R5_seq[-1]
            # if possible - get first clipped nucleotide
            if(clip3_R5_length > 0):
                R5_first_clipped_nucleotide = clipped_R5[0]
            else:
                R5_first_clipped_nucleotide = 'NA'
        else:
            # if read was not mapped - return NA
            R5_last_mapped_nucleotide = 'NA'
            R5_first_clipped_nucleotide = 'NA'

        # process R3 reads analogously to above:
        if (R3_mapping_pos != "-1"):
            R3_last_mapped_nucleotide = R3_seq[-1]
            if(clip3_R3_length > 0):
                R3_first_clipped_nucleotide = clipped_R3[0]
            else:
                R3_first_clipped_nucleotide = 'NA'
        else:
            R3_last_mapped_nucleotide = 'NA'

        # perform correction of obtained clipping sequences - in case if the first clipped nucletoide was A or U(T)
        # then if last mapped nucleotide was A or U(T) - include them in the soft-clipped fragment
        tails_results[seq_id]['corrected_genomic_A'] = 0
        tails_results[seq_id]['corrected_genomic_T'] = 0
        if (R5_mapping_pos != -1):
            # check if the end of mapping was A (it is hard to determine if those are coded in genome or added later) (if first clipped is also A)
            if ((R5_last_mapped_nucleotide == 'A') & ((R5_first_clipped_nucleotide == 'A') | (R5_first_clipped_nucleotide == 'NA'))):
                # if last mapped nucleotide is 'A' - check if all A mapped in the end of read could be tail
                # taking into account R5 read as it is more probable to have proper nucleotides than in the R3 read which has a large number of homopolymers
                match_read_nucleotides = re.search(
                    regex_for_genome_encoded_Atail, str(R5_seq))
                if(match_read_nucleotides):
                    # if there is a match - correct clipping
                    mapped_tail_nucleotides = match_read_nucleotides.group(
                        "genome_encoded_tail")
                    number_mapped_tail_nucleotides = len(mapped_tail_nucleotides)
                    clipped_R5 = mapped_tail_nucleotides + clipped_R5
                    clip3_R5_length = number_mapped_tail_nucleotides + clip3_R5_length
                    # modify mapping position
                    R5_mapping_pos = int(R5_mapping_pos) - \
                        number_mapped_tail_nucleotides
                    R5_seq = R5_seq[:-number_mapped_tail_nucleotides]
                    if(len(str(R5_seq)) > 0):
                        R5_last_mapped_nucleotide = R5_seq[-1]
                    else:
                        R5_last_mapped_nucleotide = 'NA'
                        clipped_R5 = ''
                        R5_seq = ''
                        R5_mapping_pos = '-1'
                    tails_results[seq_id]['corrected_genomic_A'] = 1
                else:
                    print("Error - should have a match")
                    sys.exit()

        # check if the end of mapping was T (it is hard to determine if those are coded in genome or added later) (if first clipped is also T)
        if ((R5_last_mapped_nucleotide == 'T') & ((R5_first_clipped_nucleotide == 'T') | (R5_first_clipped_nucleotide == 'NA'))):
            # if last mapped nucleotide is 'T' - check if all T mapped in the end of read could be tail
            # taking into account R5 read as it is more probable to have proper nucleotides than in the R3 read which has a large number of homopolymers
            match_read_nucleotides = re.search(
                regex_for_genome_encoded_Ttail, str(R5_seq))
            if(match_read_nucleotides):
                # if there is a match - correct clipping
                mapped_tail_nucleotides = match_read_nucleotides.group(
                    "genome_encoded_tail")
                number_mapped_tail_nucleotides = len(mapped_tail_nucleotides)
                clipped_R5 = mapped_tail_nucleotides + clipped_R5
                clip3_R5_length = number_mapped_tail_nucleotides + clip3_R5_length
                # modify mapping position
                R5_mapping_pos = int(R5_mapping_pos) - \
                    number_mapped_tail_nucleotides
                R5_seq = R5_seq[:-number_mapped_tail_nucleotides]
                if(len(str(R5_seq)) > 0):
                    R5_last_mapped_nucleotide = R5_seq[-1]
                else:
                    R5_last_mapped_nucleotide = 'NA'
                    clipped_R5 = ''
                    R5_seq = ''
                    R5_mapping_pos = '-1'
                tails_results[seq_id]['corrected_genomic_T'] = 1
            else:
                print("Error - should have a match")
                sys.exit()

        # store initial results in temp dict:
        tails_results[seq_id]['tailseq_A_tail_length'] = A_tail_length
        tails_results[seq_id]['tailseq_additional_bases'] = additional_bases
        tails_results[seq_id]['additional_bases_length'] = T_tail_length
        tails_results[seq_id]['tailseq_tail_length'] = tailseq_tail_length
        tails_results[seq_id]['R5_read_3prime_clip_sequence'] = clipped_R5
        tails_results[seq_id]['R3_read_3prime_clip_sequence'] = clipped_R3
        tails_results[seq_id]['R5_mapping_pos'] = R5_mapping_pos
        tails_results[seq_id]['R3_mapping_pos'] = R3_mapping_pos
        tails_results[seq_id]['R5_read_3prime_clip_length'] = clip3_R5_length
        tails_results[seq_id]['R3_read_3prime_clip_length'] = clip3_R3_length
        tails_results[seq_id]['processed_R5_clip'] = clipped_R5
        tails_results[seq_id]['processed_R3_clip'] = clipped_R3
        tails_results[seq_id]['ref_name_R5'] = ref_name_R5
        tails_results[seq_id]['ref_name_R3'] = ref_name_R3
        tails_results[seq_id]['tailseq_delimiter_not_found'] = tailseq_delimiter_not_found
        tails_results[seq_id]['tailseq_delimiter_shifted'] = tailseq_delimiter_shifted
        tails_results[seq_id]['tailseq_delimiter_mismatch'] = tailseq_delimiter_mismatch
        tails_results[seq_id]['tailseq_predicted_tail'] = ''

        if (R5_sequenced_length<args.min_R5_length):
            tails_results[seq_id]['CTGAC_R5'] = 1


        # check for the presence of overexpression plasmid in the clipped fragment (in case of reporter LINE1 analyses)
        match_plasmid_R5 = re.search(regex_for_plasmid_seq, clipped_R5)
        if(match_plasmid_R5):
            tails_results[seq_id]['tail_source'] = 'plasmid_match_no_tail_plasmid'
            tails_results[seq_id]['tail_sequence'] = ''
            tails_results[seq_id]['mapping_position'] = -1

        #if no plasmid was identified - perform processing
        else:
            # create representation of tailseq-identified tail
            tailseq_tail = ''
            if (int(A_tail_length) > 0):
                for i in range(0, int(A_tail_length)):
                    tailseq_tail = tailseq_tail + "A"
                tailseq_tail = tailseq_tail + additional_bases

            tails_results[seq_id]['tailseq_predicted_tail'] = tailseq_tail

            if (tailseq_tail_length > 0):
                tails_results[seq_id]['CTGAC_R5'] = 1 #set to treat all tailseq_tails as true tails
                # if tailseq identified tail is > 0 bp
                # treat this as a true tail but try to find this tails in the softclipping
                if (regex_match_for_R1.group('pos') != "-1"):  # mapping of R5 read
                    # -1 means that read was unmapped
                    if (regex_match_for_R2.group('pos') != "-1"):  # mapping of R3 read
                        # paired mapping - best situation
                        tails_results[seq_id]['mapping'] = 'both'
                        if (clip3_R5_length == clip3_R3_length):
                            # if clipped sequences from both R5 and R3 read has the same length
                            #print("match!! " + clipped_R5 + " (" + str(clip3_R3_length) + ")")
                            #matched_pairs = matched_pairs + 1
                            if (tailseq_tail_length == clip3_R5_length):
                                # if tailseq-identified tails has the same length as sequence clipped from the 3' end of read R5:
                                #print("R5 matching tailseq length")
                                tails_results[seq_id]['tailseq_tail_match_clipped_tail'] = 1
                                tails_results[seq_id]['tail_source'] = 'tailseq_clip_match_R5_clip_match_R3_clip'
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                if (R5_mapping_pos == R3_mapping_pos):
                                    # if sequences were mapped in the same position
                                    tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                else:
                                    # else - treat R3 mapping as the proper one
                                    tails_results[seq_id]['mapping_position'] = R3_mapping_pos

                            else:
                                # clipped sequences lengths don't match tailseq_length
                                # treat tailseeker identified tail as the proper one
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                if (clip3_R5_length > tailseq_tail_length):
                                    tails_results[seq_id]['tail_source'] = 'tailseq_clip_clip_R3_R5_longer_than_tailseq'
                                else:
                                    tails_results[seq_id]['tail_source'] = 'tailseq_clip_clip_R3_R5_shorter_than_tailseq'
                                if (R5_mapping_pos == R3_mapping_pos):
                                    # if sequences were mapped in the same position
                                    tails_results[seq_id]['mapping_position'] = R3_mapping_pos
                                else:
                                    # else - treat R5 mapping as the proper one
                                    tails_results[seq_id]['mapping_position'] = R5_mapping_pos

                        else:
                            # clipped sequences dont have the same length
                            # treat the tailseeker identified tail as the proper one
                            if (tailseq_tail_length == clip3_R3_length):
                                # sequence clipped from R3 has proper length
                                tails_results[seq_id]['tail_source'] = 'tailseq_clip_R3_match_length'
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                tails_results[seq_id]['mapping_position'] = R3_mapping_pos
                            elif (tailseq_tail_length == clip3_R5_length):
                                # sequence clipped from R5 has proper length
                                tails_results[seq_id]['tail_source'] = 'tailseq_clip_R5_match_length'
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                            else:
                                # both R3 and R5 have different length than tailseq identified tail
                                # treat tailseeker identified tail as the proper one
                                # but try to find tail in the softclipped output (for additional analyses)
                                # find R3 tail in R5 tail
                                my_regex = r"^(?P<tail>" + re.escape(clipped_R3) + \
                                    r")(?P<rest_of_seq>.*)"
                                tail_match = re.match(my_regex, clipped_R5)
                                # find R5 tail in R5 tail
                                my_regex2 = r"^(?P<tail>" + re.escape(clipped_R5) + \
                                    r")(?P<rest_of_seq>.*)"
                                tail_match2 = re.match(my_regex2, clipped_R3)
                                if (tail_match):
                                    # if match found (R3 soft clipped fragment in R5 soft-clipped):
                                    tails_results[seq_id]['matched_R3_tail_in_R5'] = 1
                                    refined_R5_clip = tail_match.group('tail')
                                    tails_results[seq_id]['processed_R5_clip'] = refined_R5_clip
                                    if (len(refined_R5_clip) > tailseq_tail_length):
                                        tails_results[seq_id]['tail_source'] = 'tailseq_clip_clip_R3_longer_than_tailseq'
                                        tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                        tails_results[seq_id]['mapping_position'] = int(
                                            R5_mapping_pos)
                                    else:
                                        tails_results[seq_id]['tail_source'] = 'tailseq_clip_clip_R3_shorter_than_tailseq'
                                        tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                        tails_results[seq_id]['mapping_position'] = int(
                                            R5_mapping_pos)
                                elif (tail_match2):
                                    # if match found:
                                    tails_results[seq_id]['matched_R5_tail_in_R3'] = 1
                                    refined_R3_clip = tail_match2.group('tail')
                                    tails_results[seq_id]['processed_R3_clip'] = refined_R3_clip
                                    tails_results[seq_id]['mapping_position'] = int(R5_mapping_pos)
                                    if (len(refined_R3_clip) > tailseq_tail_length):
                                        tails_results[seq_id]['tail_source'] = 'tailseq_clip_clip_R5_longer_than_tailseq'
                                        tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                    else:
                                        tails_results[seq_id]['tail_source'] = 'tailseq_clip_clip_R5_shorter_than_tailseq'
                                        tails_results[seq_id]['tail_sequence'] = tailseq_tail

                                else:
                                    tails_results[seq_id]['tail_source'] = 'tailseq_clip_clipping_different_lengths'
                                    tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                    # cannot determine correct mapping position
                                    tails_results[seq_id]['mapping_position'] = -1

                    else:
                        # R3 read was unmapped:
                        if (tailseq_tail_length == clip3_R5_length):
                            tails_results[seq_id]['tailseq_tail_match_clipped_tail'] = 1
                            tails_results[seq_id]['tail_source'] = 'tailseq_clip_no_R3'
                            tails_results[seq_id]['tail_sequence'] = tailseq_tail
                            tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                        else:
                            # tailseq tail length have different length than clipping
                            # try to identify tail
                            tail_match = re.match(
                                regex_for_heuristic_tail_identification_R3, str(clipped_R5))
                            # changed R3 to R5 -in regex cause no CTGAC should occur in R5 clip
                            if (tail_match):
                                tails_results[seq_id]['processed_R5_clip'] = tail_match.group(
                                    'tail')
                                tails_results[seq_id]['tail_source'] = 'tailseq_clip_heuristic_R5_clip'
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                            else:
                                tails_results[seq_id]['tail_source'] = 'tailseq_only_no_R3'
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                tails_results[seq_id]['mapping_position'] = -1

                else:
                    # R5 read was unmapped:
                    if (regex_match_for_R2.group('pos') != "-1"):  # mapping of R3 read
                        if (tailseq_tail_length == clip3_R3_length):
                            # if tailseq-identified tails has the same length as sequence clipped from the 3' end of read R5:
                            tails_results[seq_id]['tailseq_tail_match_clipped_tail'] = 1
                            tails_results[seq_id]['tail_source'] = 'tailseq_clip_noR5'
                            tails_results[seq_id]['tail_sequence'] = tailseq_tail
                            tails_results[seq_id]['mapping_position'] = R3_mapping_pos
                        else:
                            # try to identify tail based on heuristics (regex matching)
                            # tailseq tail length have different length than clipping
                            tail_match = re.match(
                                regex_for_heuristic_tail_identification_R3, str(clipped_R3))
                            if (tail_match):
                                tails_results[seq_id]['processed_R3_clip'] = tail_match.group(
                                    'tail')
                                tails_results[seq_id]['tail_source'] = 'tailseq_clip_heuristic_R3_clip'
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                tails_results[seq_id]['mapping_position'] = R3_mapping_pos
                            else:
                                tails_results[seq_id]['tail_source'] = 'tailseq_only_no_R5'
                                tails_results[seq_id]['tail_sequence'] = tailseq_tail
                                tails_results[seq_id]['mapping_position'] = -1
                    else:
                        tails_results[seq_id]['tail_source'] = 'tailseq_clip_no_mapping_R5_R3'
                        tails_results[seq_id]['mapping_position'] = -1
                        tails_results[seq_id]['tail_sequence'] = tailseq_tail

            else:
                # tailseq tail not found
                # have to identify tails based on clipping only
                # process mapped reads:
                if (regex_match_for_R1.group('pos') != "-1"):  # mapping of R5 read
                    # -1 means that read was unmapped
                    if (regex_match_for_R2.group('pos') != "-1"):  # mapping of R3 read
                        # paired mapping - best situation
                        tails_results[seq_id]['mapping'] = 'both'
                        if (clip3_R5_length == clip3_R3_length):
                            # if clipped sequences from both R5 and R3 read has the same length
                            if (R5_sequenced_length>args.min_R5_length): #check if R5 had proper length, enough for looking for tails
                                if (clip3_R5_length > 0):  # if there is any clipping
                                    if(tails_results[seq_id]['CTGAC_R5'] > 0):
                                        # check if clipping of R5 ends with tailseq delimiter CTGAC (identified at the beginning of analysis)
                                        tails_results[seq_id]['tail_source'] = 'no_tailseq_clip_R5_R3_CTGAC'
                                        tails_results[seq_id]['tail_sequence'] = clipped_R5
                                        tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                    else:
                                        # if no CTGAC was identified - store the possible tail sequence (but it will not be used in further analysis)
                                        tails_results[seq_id]['tail_source'] = 'no_tailseq_clip_R5_R3'
                                        #for reporter (short) reads - take R3 sequence (as even those which dont have CTGAC in R5 will be treated as possible tails)
                                        #if ((transcript=="REPORTERL1") or (transcript=="REPORTERL1_overexp")):
                                    #        tails_results[seq_id]['tail_sequence'] = clipped_R3
                                        #else:
                                        tails_results[seq_id]['tail_sequence'] = clipped_R5
                                        tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                else:
                                    # softclipping fragment got 0 length - treat as no_tail
                                    if (tails_results[seq_id]['CTGAC_R5'] > 0):
                                        # check if clipping of R5 ends with tailseq delimiter CTGAC (identified at the beginning of analysis)
                                        tails_results[seq_id]['tail_source'] = 'no_tailseq_no_tail_CTGAC'
                                        tails_results[seq_id]['tail_sequence'] = ''
                                        if (R5_mapping_pos == R3_mapping_pos):
                                            # if sequences were mapped in the same position
                                            tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                        else:
                                            # else - treat R5 mapping as the proper one
                                            tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                    else:
                                        tails_results[seq_id]['tail_source'] = 'no_tail'
                                        tails_results[seq_id]['tail_sequence'] = ''
                                        tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                            else:
                                print ("R5 too short")
                                #if R5 length is too low, look for possible tails in R3 read
                                tail_match = re.match(
                                    regex_for_heuristic_tail_identification_R3, str(clipped_R3))
                                if (tail_match):
                                    tails_results[seq_id]['processed_R3_clip'] = tail_match.group(
                                        'tail')
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R5_heuristic_R3_clip'
                                    tails_results[seq_id]['tail_sequence'] = tail_match.group(
                                        'tail')
                                    tails_results[seq_id]['mapping_position'] = R3_mapping_pos
                                else:
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R5_no_tail_pattern'
                                    tails_results[seq_id]['tail_sequence'] = ''
                                    tails_results[seq_id]['mapping_position'] = -1

                        else:
                            # clipped sequences dont have the same length
                            # both R3 and R5 have different length than tailseq identified tail
                            # treat the R5 clipping as more reliable (especially if CTGAC delimiter sequence was identified in clipping)
                            if (R5_sequenced_length>args.min_R5_length): #check if R5 had proper length, enough for looking for tails
                                if (tails_results[seq_id]['CTGAC_R5'] > 0):
                                    # check if clipping of R5 ends with tailseq delimiter CTGAC (identified at the beginning of analysis)
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_clip_clipping_different_lengths_R5_CTGAC'
                                    tails_results[seq_id]['tail_sequence'] = clipped_R5
                                    tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                else:
                                    # if not CTGAC was identified - store the possible tail sequence (but it will not be used in further analysis)
                                    tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_clip_clipping_different_lengths_R5'
                                    #for reporter (short) reads - take R3 sequence (as even those which dont have CTGAC in R5 will be treated as possible tails)
                                    #if ((transcript=="REPORTERL1") or (transcript=="REPORTERL1_overexp")):
                                    #    tails_results[seq_id]['tail_sequence'] = clipped_R3
                                    #else:
                                    tails_results[seq_id]['tail_sequence'] = clipped_R5
                            else: #look for tail in R3
                                tail_match = re.match(
                                    regex_for_heuristic_tail_identification_R3, str(clipped_R3))
                                if (tail_match):
                                    tails_results[seq_id]['processed_R3_clip'] = tail_match.group(
                                        'tail')
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R5_heuristic_R3_clip'
                                    tails_results[seq_id]['tail_sequence'] = tail_match.group(
                                        'tail')
                                    tails_results[seq_id]['mapping_position'] = R3_mapping_pos
                                else:
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R5_no_tail_pattern'
                                    tails_results[seq_id]['tail_sequence'] = clipped_R3
                                    tails_results[seq_id]['mapping_position'] = -1

                    else:
                        # no R3 read was mapped, try to find tail in R5 read
                        if (R5_sequenced_length>args.min_R5_length): #check if R5 had proper length, enough for looking for tails
                            if (tails_results[seq_id]['CTGAC_R5'] > 0):
                                # check if clipping of R5 ends with tailseq delimiter CTGAC (identified at the beginning of analysis)
                                if(clip3_R5_length == 0):
                                    # softclipping fragment got 0 length - treat as no_tail
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R3_R5_no_tail_CTGAC'
                                    tails_results[seq_id]['tail_sequence'] = ''
                                    tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                                else:
                                    # else - store clipped fragment as a possible tail
                                    tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R3_R5_CTGAC'
                                    tails_results[seq_id]['tail_sequence'] = clipped_R5
                                    tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                            else:
                                # if not CTGAC was identified - store the possible tail sequence (but it will not be used in further analysis)
                                tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R3_R5_clip'
                                tails_results[seq_id]['tail_sequence'] = clipped_R5
                                tails_results[seq_id]['mapping_position'] = R5_mapping_pos
                        else:
                            tails_results[seq_id]['tail_source'] = 'no_tailseq_no_mapping_too_short_R5'
                            tails_results[seq_id]['tail_sequence'] = ''
                            tails_results[seq_id]['mapping_position'] = -1

                else:
                    # if R5 was unmapped, check if R3 was mapped:
                    if (regex_match_for_R2.group('pos') != "-1"):
                        # use heuristics (regex matching) to identify possible tail
                        tail_match = re.match(
                            regex_for_heuristic_tail_identification_R3, str(clipped_R3))
                        if (tail_match):
                            tails_results[seq_id]['processed_R3_clip'] = tail_match.group(
                                'tail')
                            tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R5_heuristic_R3_clip'
                            tails_results[seq_id]['tail_sequence'] = tail_match.group(
                                'tail')
                            tails_results[seq_id]['mapping_position'] = R3_mapping_pos
                        else:
                            tails_results[seq_id]['tail_source'] = 'no_tailseq_no_R5_no_tail_pattern'
                            tails_results[seq_id]['tail_sequence'] = ''
                            tails_results[seq_id]['mapping_position'] = -1
                    else:
                        # if both reads were unmapped - no tail:
                        tails_results[seq_id]['tail_source'] = 'no_tailseq_no_mapping'
                        tails_results[seq_id]['tail_sequence'] = ''
                        tails_results[seq_id]['mapping_position'] = -1

        #for reporter sequences (which got short reads for R5 and CTGAC presence is not expected)
        #for compatibility with genomic sequences treat all as containing CTGAC
        #if ((transcript == "REPORTERL1") or (transcript == "REPORTERL1_overexp")):
        #    tails_results[seq_id]['CTGAC_R5']=1

        # perform final processing of tail data
        tail_sequence = tails_results[seq_id]['tail_sequence']
        # store final results for sequence
        final_results[seq_id]['tail_sequence'] = tail_sequence
        # initialize values:
        Atail = ''
        Atail_length = 0
        Utail = ''
        Utail_length = 0
        Gtail = ''
        Gtail_length = 0
        tail_type = ''

        # define regex for ientification of tail types:
        A_only_tail_match = re.match("^(?P<Atail>^A+$)", tail_sequence)
        AU_tail_match = re.match("^(?P<Atail>^A+)(?P<Utail>T+)$", tail_sequence)
        AG_tail_match = re.match("^(?P<Atail>^A+)(?P<Utail>G+)$", tail_sequence)
        U_only_tail_match = re.match("^(?P<Utail>^T+$)", tail_sequence)
        UA_tail_match = re.match("^(?P<Utail>^T+)(?P<Atail>A+)$", tail_sequence)
        UG_tail_match = re.match("^(?P<Utail>^T+)(?P<Atail>G+)$", tail_sequence)
        Amixed_match = re.match("^(?P<Atail>A+[TGCA]{0,5}?A+$)", tail_sequence)
        Umixed_match = re.match("^(?P<Utail>T+[TGCA]{0,5}?T+$)", tail_sequence)
        AmixedU_match = re.match(
            "^(?P<Atail>A+[TGCA]{0,5}?A+.)(?P<Utail>T+[TGCA]{0,2}T+.)$", tail_sequence)
        AmixedG_match = re.match(
            "^(?P<Atail>A+[TGCA]{0,5}?A+.)(?P<Utail>G+.)$", tail_sequence)
        UmixedA_match = re.match(
            "^(?P<Utail>T+[TGCA]{0,5}?T+.)(?P<Atail>A+[TGCA]{0,2}A+.)$", tail_sequence)
        UmixedG_match = re.match(
            "^(?P<Utail>T+[TGCA]{0,5}?T+.)(?P<Atail>G+.)$", tail_sequence)
        AmixedUmixed_match = re.match(
            "^(?P<Atail>A+[TGCA]{0,5}?A+.)(?P<Utail>T+[TGCA]{0,2}T+.$)", tail_sequence)

        T_count = tail_sequence.count("T")
        A_count = tail_sequence.count("A")
        C_count = tail_sequence.count("C")
        G_count = tail_sequence.count("G")
        A_threshold_heteroA = 0.75
        T_threshold_heteroT = 0.75

        uridylated=0

        # analyze tail types
        if (tail_sequence == ''):
            Atail = ''
            Utail = ''
            Atail_length = 0
            Utail_length = 0
            Gtail = ''
            tail_type_for_anal = 'no_tail'
            tail_type = 'no_tail'
            #below - switched off, as it can miss some tailseq notails, if R5 was not properly mapped
            # require a CTGAC sequence in clipped fragment or R5 read to consider sequence without tail for further analyses
#            if(re.match("(.*)CTGAC$", tails_results[seq_id]['tail_source'])):
#                tail_type_for_anal = 'no_tail'
#                tail_type = 'no_tail'
#            else:
#                tail_type_for_anal = 'false_no_tail_no_CTGAC'
#                tail_type = 'false_no_tail_no_CTGAC'
        elif (A_only_tail_match):
            Atail = A_only_tail_match.group("Atail")
            Atail_length = len(A_only_tail_match.group("Atail"))
            Utail = ''
            Utail_length = 0
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'A_only'
            tail_type = 'A_only'
        elif (AU_tail_match):
            Atail = AU_tail_match.group("Atail")
            Atail_length = len(AU_tail_match.group("Atail"))
            Utail = AU_tail_match.group("Utail")
            Utail_length = len(AU_tail_match.group("Utail"))
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'AU'
            tail_type = 'AU'
            uridylated=1
        elif (AG_tail_match):
            Atail = AG_tail_match.group("Atail")
            Atail_length = len(AG_tail_match.group("Atail"))
            Gtail = AG_tail_match.group("Utail")
            Gtail_length = len(AG_tail_match.group("Utail"))
            Utail = ''
            Utail_length = 0
            tail_type_for_anal = 'AG'
            tail_type = 'AG'
        elif (U_only_tail_match):
            Atail = ""
            Atail_length = 0
            Utail = U_only_tail_match.group("Utail")
            Utail_length = len(U_only_tail_match.group("Utail"))
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'U_only'
            tail_type = 'U_only'
            uridylated=1
        elif (UA_tail_match):
            Atail = UA_tail_match.group("Atail")
            Atail_length = len(UA_tail_match.group("Atail"))
            Utail = UA_tail_match.group("Utail")
            Utail_length = len(UA_tail_match.group("Utail"))
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'UA'
            tail_type = 'UA'
            uridylated=1
        elif (UG_tail_match):
            Gtail = UG_tail_match.group("Atail")
            Gtail_length = len(UG_tail_match.group("Atail"))
            Utail = UG_tail_match.group("Utail")
            Utail_length = len(UG_tail_match.group("Utail"))
            Atail = ''
            tail_type_for_anal = 'UG'
            tail_type = 'UG'
            uridylated=1
        elif (Umixed_match):
            Atail = ''
            Atail_length = 0
            Utail = Umixed_match.group("Utail")
            Utail_length = len(Umixed_match.group("Utail"))
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'U_mixed'
            tail_type = 'U_heterogenous'
            uridylated=1
        elif (AmixedU_match):
            Atail = AmixedU_match.group("Atail")
            Atail_length = len(AmixedU_match.group("Atail"))
            Utail = AmixedU_match.group("Utail")
            Utail_length = len(AmixedU_match.group("Utail"))
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'AU_mixed'
            tail_type = 'A_heterogenous'
            uridylated=1
        elif (AmixedG_match):
            Atail = AmixedG_match.group("Atail")
            Atail_length = len(AmixedG_match.group("Atail"))
            Gtail = AmixedG_match.group("Utail")
            Gtail_length = len(AmixedG_match.group("Utail"))
            Utail = ''
            Utail_length = 0
            tail_type_for_anal = 'AG_mixed'
            tail_type = 'A_heterogenous'
        elif (UmixedA_match):
            Atail = UmixedA_match.group("Atail")
            Atail_length = len(UmixedA_match.group("Atail"))
            Utail = UmixedA_match.group("Utail")
            Utail_length = len(UmixedA_match.group("Utail"))
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'UA_mixed'
            tail_type = 'U_heterogenous'
            uridylated=1
        elif (UmixedG_match):
            Gtail = UmixedG_match.group("Atail")
            Gtail_length = len(UmixedG_match.group("Atail"))
            Utail = UmixedG_match.group("Utail")
            Utail_length = len(UmixedG_match.group("Utail"))
            Atail = ''
            Atail_length = 0
            tail_type_for_anal = 'UG_mixed'
            tail_type = 'U_heterogenous'
            uridylated=1
        elif (AmixedUmixed_match):
            Atail = AmixedUmixed_match.group("Atail")
            Atail_length = len(AmixedUmixed_match.group("Atail"))
            Utail = AmixedUmixed_match.group("Utail")
            Utail_length = len(AmixedUmixed_match.group("Utail"))
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'A_mixed_U_mixed'
            tail_type = 'heterogenous'
            uridylated=1
        elif (Amixed_match):
            Atail = Amixed_match.group("Atail")
            Atail_length = len(Amixed_match.group("Atail"))
            Utail = ''
            Utail_length = 0
            Gtail = ''
            Gtail_length = 0
            tail_type_for_anal = 'A_mixed'
            tail_type = 'A_heterogenous'
        else:
            # if no regex was matched - assign to class "other"
            Atail = ''
            Atail_length = 0
            Utail = ''
            Utail_length = 0
            Gtail = ''
            Gtail_length = 0
            tail_type = 'other'
            tail_type_for_anal = 'other'

        if(tails_results[seq_id]['tail_source'] == 'plasmid_match_no_tail'):
            tail_type = 'plasmid_match_no_tail'

        if (tail_type=='other'):
            tail_length = len(tails_results[seq_id]['tail_sequence'])
            count_As = tails_results[seq_id]['tail_sequence'].count("A")
            count_Ts = tails_results[seq_id]['tail_sequence'].count("T")
            ratio_As = count_As/tail_length
            ratio_Ts = count_Ts/tail_length
            sumAT=ratio_As+ratio_Ts
            #print(seq_id,count_As,ratio_As,count_Ts,ratio_Ts,sumAT,tails_results[seq_id]['tail_sequence'])
            if((tail_length>=4) & (sumAT>=0.75)):
                match_AUtail = re.match("^(?P<Atail>^A+).*(T.|T)$", tail_sequence)
                match_Utail = re.match("^TT.*",tail_sequence)
                match_Atail = re.match("^AA.*A$",tail_sequence)
                if (ratio_Ts>0.6):
                    if(match_AUtail):
                        tail_type_for_anal = 'A_mixed_U_mixed'
                        tail_type = 'heterogenous'
                    else:
                        tail_type_for_anal = 'U_mixed'
                        tail_type = 'U_heterogenous'
                    uridylated=1
                elif(ratio_Ts<0.05):
                    if(match_AUtail):
                        tail_type_for_anal = 'A_mixed_U_mixed'
                        tail_type = 'heterogenous'
                        uridylated=1
                    #elif(match_Utail):
                    #    tail_type_for_anal = 'U_mixed'
                    #    tail_type = 'U_heterogenous'
                    #    uridylated=1
                    else:
                        tail_type_for_anal = 'A_mixed'
                        tail_type = 'A_heterogenous'
                else:
                    if(match_AUtail):
                        tail_type_for_anal = 'A_mixed_U_mixed'
                        tail_type = 'heterogenous'
                        uridylated=1
                    elif(match_Atail):
                        tail_type_for_anal = 'A_mixed'
                        tail_type = 'A_heterogenous'



                #print(tail_type)


        # store final results
        final_results[seq_id]['tail_type'] = tail_type
        final_results[seq_id]['tail_type_mixed'] = tail_type_for_anal
        final_results[seq_id]['Atail'] = Atail
        final_results[seq_id]['Atail_length'] = Atail_length
        final_results[seq_id]['Utail'] = Utail
        final_results[seq_id]['Utail_length'] = Utail_length
        final_results[seq_id]['Gtail'] = Gtail
        final_results[seq_id]['Gtail_length'] = Gtail_length
        final_results[seq_id]['tail_length'] = len(
            final_results[seq_id]['tail_sequence'])
        final_results[seq_id]['tail_source'] = tails_results[seq_id]['tail_source']
        final_results[seq_id]['transcript'] = transcript
        final_results[seq_id]['cell_line'] = cell_line
        final_results[seq_id]['person'] = person
        final_results[seq_id]['localization'] = localization
        final_results[seq_id]['condition'] = condition
        final_results[seq_id]['replicate'] = replicate
        final_results[seq_id]['sample_name'] = sample_name
        final_results[seq_id]['primer_name'] = primer_name
        final_results[seq_id]['project_name'] = project_name
        final_results[seq_id]['mapping_position'] = tails_results[seq_id]['mapping_position']
        final_results[seq_id]['exp_type'] = exp_type
        final_results[seq_id]['R5_sequenced_length'] = R5_sequenced_length
        final_results[seq_id]['R3_mapping_position'] = tails_results[seq_id]['R3_mapping_pos']
        final_results[seq_id]['R5_mapping_position'] = tails_results[seq_id]['R5_mapping_pos']
        final_results[seq_id]['R5_seq'] = record.seq
        final_results[seq_id]['R3_seq'] = record2.seq
        final_results[seq_id]['R5_clip'] = clipped_R5
        final_results[seq_id]['R3_clip'] = clipped_R3
        final_results[seq_id]['R5_last_mapped_nucleotide'] = R5_last_mapped_nucleotide
        final_results[seq_id]['R3_last_mapped_nucleotide'] = R3_last_mapped_nucleotide
        final_results[seq_id]['PCR_duplicates'] = PCRduplicates
        final_results[seq_id]['ref_name_R5'] = ref_name_R5
        final_results[seq_id]['ref_name_R3'] = ref_name_R3
        #final_results[seq_id]['CTGAC_R5'] = 1
        final_results[seq_id]['CTGAC_R5'] = tails_results[seq_id]['CTGAC_R5']
        final_results[seq_id]['CTGAC_R5_mismatched'] = tails_results[seq_id]['CTGAC_R5_mismatched']
        final_results[seq_id]['corrected_genomic_T'] = tails_results[seq_id]['corrected_genomic_T']
        final_results[seq_id]['corrected_genomic_A'] = tails_results[seq_id]['corrected_genomic_A']
        final_results[seq_id]['heterogenous_end'] = tails_results[seq_id]['heterogenous_end']
        final_results[seq_id]['heterogenous_end_R3'] = tails_results[seq_id]['heterogenous_end_R3']
        final_results[seq_id]['terminal_nucleotides'] = tails_results[seq_id]['terminal_nucleotides']
        final_results[seq_id]["tailseq_predicted_tail"]=tails_results[seq_id]['tailseq_predicted_tail']
        final_results[seq_id]["heterogenous_end_length"]=tails_results[seq_id]['heterogenous_end_length']
        final_results[seq_id]["heterogenous_end_R3_length"]=tails_results[seq_id]['heterogenous_end_R3_length']
        final_results[seq_id]["mapping_spanning_delimiter"]=tails_results[seq_id]['mapping_spanning_delimiter']
        final_results[seq_id]['tailseq_delimiter_mismatch']=tails_results[seq_id]['tailseq_delimiter_mismatch']
        final_results[seq_id]['tailseq_delimiter_not_found']=tails_results[seq_id]['tailseq_delimiter_not_found']
        final_results[seq_id]['uridylated']=uridylated

    return final_results

# end


analyzed = 0
os.chdir(args.inputdir)


# define default files to search
files_to_search = "*R5.fastq"

# get files to search from command-line (if present)
if (args.glob):
    files_to_search = args.glob

for R5_file in glob.glob(files_to_search):
    # iterate through R5 files
    analyzed = analyzed + 1  # increment number of analyzed files
    print(R5_file)
    file_parts = re.search("(?P<path>.*//|)(?P<basename>.*fastq)", R5_file)
    file_basename = file_parts.group("basename")  # get file basename
    file_path = file_parts.group("path")
    file_parts = re.search("(?P<prefix>.*)_R5(?P<suffix>.*)", file_basename)
    file_prefix = file_parts.group("prefix")
    file_suffix = file_parts.group("suffix")
    R5_fasta_file = file_prefix + "_.fasta"
    R3_file = file_prefix + "_R3" + file_suffix  # create R3 file name
    SAM_file_R5 = R5_file + ".sam"  # bowtie output files (SAM)
    SAM_file_R3 = R3_file + ".sam"  # bowtie output files (SAM)
    bowtie_output_R5 = R5_file + ".bowtie.out"  # bowtie output files (log)
    bowtie_output_R3 = R3_file + ".bowtie.out"  # bowtie output files (log)
    softclipped_fasta_R5 = R5_file + ".sam.clipped.fasta"
    softclipped_fasta_R3 = R3_file + ".sam.clipped.fasta"
    # get information about sample pair from samplesheet:
    temp = data[data["R5_file"] == R5_file]
    transcript = temp['transcript'][0]  # transcript name
    cell_line = temp['cell_line'][0]  # cell line
    condition = temp['condition'][0]  # condition
    person = temp['person'][0]  # person - who prepared the library
    # subcellular localization (fractionation experiment)
    localization = temp['localization'][0]
    replicate = temp['replicate'][0]  # replicate number
    sample_name = temp['Sample_Name'][0]  # sample name
    primer_name = temp['primer_name'][0]  # name of the primer used in RACE
    project_name = temp.index.values[0] #get Project name
    exp_type = temp['type'][0] #xperiment yype - OVR, KD, NT, LEAP ...
    R5_sequenced_length = temp['R5_length'][0]
    # if bowtie was not run before - run bowtie2 on R5 and R3 files:
    # check if there is any sequence in R5 file, if not - skip analysis for this pair
    R5_records = list(SeqIO.parse(R5_file, "fastq"))

    if(len(R5_records) == 0):
        print("skipping analysis because no sequences found in the R5 file")
    else:
        if (transcript != 'ENDOL1'):
            # get genome - index to be used for mapping
            genome = transcript_genomes[transcript]
            print ("running bowtie on file " + R5_file +
                   " with transcript " + transcript)
            if os.path.isfile(SAM_file_R5):
                if os.stat(SAM_file_R5).st_size > 0:
                    print("bowtie output " + SAM_file_R5 +
                          " exists. Will reuse previous results.")
                else:
                    print("bowtie output file " + SAM_file_R5 +
                          " exists but is zero-size. Will attempt to rerun the mapping.")
                    subprocess.call(bowtie2_path + " -x " + genome + " -U " + R5_file + " -S " + SAM_file_R5 +
                                    " --very-sensitive-local -L 25 -D 10 -p " + bowtie_threads + " 2> " + bowtie_output_R5, shell=True)
            else:
                subprocess.call(bowtie2_path + " -x " + genome + " -U " + R5_file + " -S " + SAM_file_R5 +
                                " --very-sensitive-local -L 25 -D 10 -p " + bowtie_threads + " 2> " + bowtie_output_R5, shell=True)
            if os.path.isfile(SAM_file_R3):
                if os.stat(SAM_file_R3).st_size > 0:
                    print("bowtie output " + SAM_file_R3 +
                          " exists. Will reuse previous results.")
                else:

                    print("bowtie output file " + SAM_file_R3 +
                          " exists but is zero-size. Will attempt to rerun the mapping.")
                    subprocess.call(bowtie2_path + " -x " + genome + " -U " + R3_file + " -S " + SAM_file_R3 +
                                    " --very-sensitive-local -L 25 -D 10 -p " + bowtie_threads + " 2> " + bowtie_output_R3, shell=True)
            else:
                subprocess.call(bowtie2_path + " -x " + genome + " -U " + R3_file + " -S " + SAM_file_R3 +
                                " --very-sensitive-local -L 25 -D 10 -p " + bowtie_threads + " 2> " + bowtie_output_R3, shell=True)

            # if softclipping output does not exist or R5 and R3 reads - perform identification of soft clipped fragments:
            if os.path.isfile(softclipped_fasta_R5):
                if os.stat(softclipped_fasta_R5).st_size > 0:
                    print("softclipping output " + softclipped_fasta_R5 +
                          " exists. Will reuse previous results.")
                else:
                    print("softclipping output file " + softclipped_fasta_R5 +
                          " exists but is zero-size. Will attempt to rerun the analysis.")
                    subprocess.call(get_sofclipped_script_path + " -input " +
                                    SAM_file_R5 + " -output " + softclipped_fasta_R5, shell=True)
            else:
                subprocess.call(get_sofclipped_script_path + " -input " +
                                SAM_file_R5 + " -output " + softclipped_fasta_R5, shell=True)
            if os.path.isfile(softclipped_fasta_R3):
                if os.stat(softclipped_fasta_R3).st_size > 0:
                    print("softclipping output " + softclipped_fasta_R3 +
                          " exists. Will reuse previous results.")
                else:
                    print("softclipping output file " + softclipped_fasta_R3 +
                          " exists but is zero-size. Will attempt to rerun the analysis.")
                    subprocess.call(get_sofclipped_script_path + " -input " +
                                    SAM_file_R3 + " -output " + softclipped_fasta_R3 + " -R3 1", shell=True)
            else:
                subprocess.call(get_sofclipped_script_path + " -input " +
                                SAM_file_R3 + " -output " + softclipped_fasta_R3 + " -R3 1", shell=True)
        else:
            print(
                "Skipping bowtie part because genomic LINE sequences are processed using RepeatMasker")

        # Run the analysis of tails:
        paired_results = analyze_tails(softclipped_fasta_R5, softclipped_fasta_R3, transcript,
                                       sample_name, localization, replicate, condition, cell_line, primer_name, person, project_name, exp_type, R5_sequenced_length)
        # Create pandas data frame with result
        tails_df = pd.DataFrame.from_dict(paired_results, orient='index')
#        print(tails_df)
        # make sure data are saved after each library processed:
        print("Saving data for sample " + sample_name)
        if (analyzed > 1):
            tails_df.to_csv(args.output, mode='a', sep='\t', header=False)
        else:
            tails_df.to_csv(args.output, mode='w', sep='\t', header=True)


print("all " + str(analyzed) + " samples analyzed succesfully\n")


### END ###
