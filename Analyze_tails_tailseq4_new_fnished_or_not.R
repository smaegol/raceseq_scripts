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



# load libraries
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggseqlogo)
library(RColorBrewer)

melt_data_localization <- "/home/smaegol/storage/analyses/tail_seq_5/ALL/processing_out/all_samples_2_3_4_5_anal.tsv"
# read tailing information to data.frame using data.table
tails_data_melt <- fread(melt_data_localization, sep = "\t", header = T, stringsAsFactors = T,
                         data.table = F, showProgress = TRUE)  # read data





#tails_data_melt[tails_data_melt$localization=='JADRO',]$localization<-"TOTAL"
tails_data_melt$tailed <- tails_data_melt$tail_length > 0  #mark tailed reads
tails_data_melt$mapped <- tails_data_melt$mapping_position != -1  #mark mapped reads
tails_data_mapped <- tails_data_melt[tails_data_melt$mapped != 0, ]  #discard unmapped reads

# fix of flase_no_tail - should be assigned no tail cause CTGAC was found in the
# clipped fragment - to be fixed in python script
#tails_data_mapped$tail_type <- as.character(tails_data_mapped$tail_type)
#tails_data_mapped[tails_data_mapped$tail_type == "false_no_tail_no_CTGAC", ]$tail_type <- "no_tail"

tails_data_mapped$ref_name_R5 <- as.character(tails_data_mapped$ref_name_R5)
tails_data_mapped$ref_name_R3 <- as.character(tails_data_mapped$ref_name_R3)
# tails_data_mapped_same_ref<-tails_data_mapped[tails_data_mapped$ref_name_R5==tails_data_mapped$ref_name_R3,]

# mark uridylated reads
tails_data_mapped$uridylated2 <- FALSE
tails_data_mapped[tails_data_mapped$Utail_length > 0, ]$uridylated2 = TRUE

#convert T to U in terminal nucleotides (for seqlogo)
tails_data_mapped$terminal_nucleotides<-as.character(tails_data_mapped$terminal_nucleotides)
tails_data_mapped$terminal_nucleotides<-gsub("T","U",tails_data_mapped$terminal_nucleotides)

summary(as.factor(tails_data_mapped$CTGAC_R5))

colnames(tails_data_mapped)
head(tails_data_mapped[,c(1,2,3,12,24,25,26,27,33,34,37,38)])

# in further analyses use only those read which got CTGAC delimiter identified in
# the clipped fragment
# for reporter analyses all mapped reads got CTGAC_R5 variable = 1 (because of short reads) so they will be included in the analysis
tails_data_mapped_true <- tails_data_mapped[tails_data_mapped$CTGAC_R5 > 0, ]
tails_data_mapped_true$ref_name = tails_data_mapped_true$ref_name_R5  #use ref_name_R5 as ref_name

head(tails_data_mapped_true[,c(1,2,3,12,24,25,26,27,33,34,37,38)])

head(tails_data_mapped[tails_data_mapped$tail_source=='tailseq_clip_clip_R3_shorter_than_tailseq',c(1,2,26,27,40)],30)

tails_data_mapped_true$tail_type = as.character(tails_data_mapped_true$tail_type)

head(tails_data_mapped_true[grep("A_mix",tails_data_mapped_true$tail_type_mixed),])

#treat all mixed (heterogenous) as true tails
#tails_data_mapped_true[grep("AU_mixed",tails_data_mapped_true$tail_type_mixed),]$tail_type = 'AU'
#tails_data_mapped_true[grep("^A_mixed",tails_data_mapped_true$tail_type_mixed),]$tail_type = 'A_only'
#tails_data_mapped_true[grep("^U_mixed",tails_data_mapped_true$tail_type_mixed),]$tail_type = 'U_only'
tails_data_mapped_true[tails_data_mapped_true$tail_type=='AU',]$uridylated2 = TRUE
tails_data_mapped_true[tails_data_mapped_true$tail_type=='U_only',]$uridylated2 = TRUE
#if clipping was shorter than tailseq - treat as other, as it is probably stretch of As in the reference seq
#tails_data_mapped_true[tails_data_mapped_true$tail_source=='tailseq_clip_clip_R3_shorter_than_tailseq',]$tail_type = 'no_tail'
#tails_data_mapped_true[tails_data_mapped_true$tail_source=='tailseq_clip_clip_R3_shorter_than_tailseq',]$uridylated = FALSE


# remove heterogenous tails from analysis
tails_data_mapped_true_no_hetero = tails_data_mapped_true[-grep("hetero", tails_data_mapped_true$tail_type),
                                                          ]

head(tails_data_mapped_true_no_hetero,20)
# remove other type tails (for which we can suspect they are not tails but rather origin from improper mapping/repeatmasker) from the analysis
tails_data_mapped_true_no_hetero_no_other = tails_data_mapped_true_no_hetero[-grep("other",
                                                                                   tails_data_mapped_true_no_hetero$tail_type), ]


# treat all AG,UG or UA tails as other_no_tail

tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type ==
                                            "AG", ]$tail_type <- "other_no_tail"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type ==
                                            "UG", ]$tail_type <- "other_no_tail"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$tail_type ==
                                            "UA", ]$tail_type <- "other_no_tail"
tails_data_mapped_true_no_hetero_no_other <- tails_data_mapped_true_no_hetero_no_other[-grep("other",
                                                                                             tails_data_mapped_true_no_hetero_no_other$tail_type), ] #remove all other from analysis


# create classes for A-tail lengths (0,1,2-5,6-10,11-20,21-30,30+)
tails_data_mapped_true_no_hetero_no_other$A_length = ""
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length ==
                                            0, ]$A_length = "0"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length ==
                                            1, ]$A_length = "1"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
                                            seq(2, 5, 1), ]$A_length = "2-5"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
                                            seq(6, 10, 1), ]$A_length = "6-10"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
                                            seq(11, 20, 1), ]$A_length = "11-20"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length %in%
                                            seq(21, 30, 1), ]$A_length = "21-30"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Atail_length >
                                            30, ]$A_length = "30+"
tails_data_mapped_true_no_hetero_no_other$A_length <- factor(tails_data_mapped_true_no_hetero_no_other$A_length,
                                                             levels = c("0", "1", "2-5", "6-10", "11-20", "21-30", "30+"))

# create classes for U-tail lengths (0,1,2,3-5,6-10,10+)
tails_data_mapped_true_no_hetero_no_other$U_length = ""
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length ==
                                            0, ]$U_length = "0"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length ==
                                            1, ]$U_length = "1"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length ==
                                            2, ]$U_length = "2"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length %in%
                                            seq(3, 5, 1), ]$U_length = "3-5"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length %in%
                                            seq(6, 10, 1), ]$U_length = "6-10"
tails_data_mapped_true_no_hetero_no_other[tails_data_mapped_true_no_hetero_no_other$Utail_length >
                                            10, ]$U_length = "10+"
tails_data_mapped_true_no_hetero_no_other$U_length <- factor(tails_data_mapped_true_no_hetero_no_other$U_length,
                                                             levels = c("10+", "6-10", "3-5", "2", "1", "0"))


# modify levels of tail_types to have U_only,A-only,AU or no_tail
tails_data_mapped_true_no_hetero_no_other_tails <- tails_data_mapped_true_no_hetero_no_other
tails_data_mapped_true_no_hetero_no_other_tails$tail_type <- as.character(tails_data_mapped_true_no_hetero_no_other_tails$tail_type)
tails_data_mapped_true_no_hetero_no_other_tails$tail_type <- factor(tails_data_mapped_true_no_hetero_no_other_tails$tail_type,
                                                                    levels = c("U_only", "AU", "no_tail", "A_only"))




# create dataframe with PA1 cell_line data
tails_data_mapped_true_no_hetero_no_other_PA1 <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$cell_line ==
                                                                                                   "PA1", ]
# create dataframe with PA1 cell_line data for knockdown conditions
tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD <- tails_data_mapped_true_no_hetero_no_other_PA1[tails_data_mapped_true_no_hetero_no_other_PA1$condition !=
                                                                                                          "NT", ]
# filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD <- tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$tail_length <=
                                                                                                                   64, ]

# create dataframe with PA1 cell_line data for untreated
tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT <- tails_data_mapped_true_no_hetero_no_other_PA1[tails_data_mapped_true_no_hetero_no_other_PA1$condition ==
                                                                                                          "NT" & tails_data_mapped_true_no_hetero_no_other_PA1$localization %in% c("CYTO","NUC") & tails_data_mapped_true_no_hetero_no_other_PA1$primer_name %in% c("L1NGS0","GAPDH"), ]
# filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT <- tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT[tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT$tail_length <=
                                                                                                                   64, ]




tails_data_mapped_true_no_hetero_no_other_tails_NT1 <- tails_data_mapped_true_no_hetero_no_other_tails[(tails_data_mapped_true_no_hetero_no_other_tails$condition ==
                                                                                                          "NT" & tails_data_mapped_true_no_hetero_no_other_tails$localization %in% c("TOTAL","JADRO")),
                                                                                                       ]
# filter out tails longer than 64nt - due to a limit of spike-ins used
tails_data_mapped_true_no_hetero_no_other_tails_NT1 <- tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$tail_length <=
                                                                                                             64, ]



tails_data_for_analysis_mysz <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$cell_line=="MYSZ",]


tails_data_for_analysis_H9 <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$cell_line=="H9",]

tails_data_for_analysis_LEAP147 <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$cell_line=="LEAP147",]

tails_data_for_analysis_LEAP317 <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$cell_line=="LEAP317",]


tails_data_for_analysis_reporter_overexp <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript=="REPORTERL1_overexp",]
tails_data_for_analysis_reporter_overexp <- tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$tail_length <=
                                                                                       64, ]

tails_data_for_analysis_reporter_UTR_all_mapped <- tails_data_mapped_true[tails_data_mapped_true$mapping_position>8950 & tails_data_mapped_true$mapping_position<9010,]



tails_data_for_analysis_reporter_UTR <- tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$mapping_position>8950 & tails_data_for_analysis_reporter_overexp$mapping_position<9010,]
tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$uridylated==1,]$uridylated=TRUE

tails_data_for_analysis_reporter_SV40body <- tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$mapping_position<9140 & tails_data_for_analysis_reporter_overexp$mapping_position>=9010,]

tails_data_for_analysis_reporter_SV40end <- tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$mapping_position>9140,]


tails_data_for_analysis_reporter_KD <- tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript=="REPORTERL1",]
tails_data_for_analysis_reporter_KD <- tails_data_for_analysis_reporter_KD[tails_data_for_analysis_reporter_KD$tail_length <=
                                                                             64, ]


### functions for summarizing data ###

summarize_tails_by_tail_type <- function(tail_data2) {
  # function summarizing tails by tail type
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
                                             localization, primer_name, tail_type), summarise, N = length(primer_name))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, localization, primer_name), transform, freq = N/sum(N))
  return(tail_data_summarized)
}


summarize_tails_by_tail_type_replicates <- function(tail_data2) {
  #function summarizing tail length data for given dataset
  tail_data<-tail_data2
  #tail_data[tail_data$Utail_length>=10,]$Utail_length<-'10+'
  #tail_data$Utail_length<-as.factor(tail_data$Utail_length)
  #tail_data$Utail_length<-factor(tail_data$Utail_length,levels=c('0','1','2','3','4','5','6','7','8','9','10+'))
  tail_data_summarized<-ddply(tail_data,.(transcript,condition,cell_line,localization,primer_name,tail_type,replicate),summarise,N=length(primer_name))
  tail_data_summarized<-ddply(tail_data_summarized,.(transcript,condition,cell_line,localization,primer_name,replicate),transform, freq=N/sum(N))
  return(tail_data_summarized)
}


summarize_tails_by_tail_type_replicates <- function(tail_data2) {
  # function summarizing tails by tail type
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
                                             localization, primer_name, tail_type), summarise, N = length(primer_name))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, localization, primer_name), transform, freq = N/sum(N))
  return(tail_data_summarized)
}



summarize_tails_by_tail_type_loc <- function(tail_data2) {
  # function summarizing tails by tail type
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
                                             localization, primer_name, tail_type), summarise, N = length(primer_name))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, primer_name), transform, freq = N/sum(N))
  return(tail_data_summarized)
}

summarize_Utails_calculate_means <- function(tail_data2) {
  # function summarizing Utail length data for given dataset
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
                                             primer_name), summarise, N = length(Utail_length), mean_Utail = mean(Utail_length,
                                                                                                                  na.rm = T), sd_utail = sd(Utail_length))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, primer_name), transform, se_utail = sd_utail/sqrt(N))
  return(tail_data_summarized)
}

summarize_Utails_calculate_means_repl <- function(tail_data2) {
  # function summarizing Utail length data for given dataset
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,replicate,
                                             primer_name), summarise, N = length(Utail_length), mean_Utail_rep = mean(Utail_length,
                                                                                                                      na.rm = T), sd_utail_rep = sd(Utail_length))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, primer_name,replicate), transform, se_utail_rep = sd_utail_rep/sqrt(N))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, primer_name), transform, mean_Utail = mean(mean_Utail_rep),sd_utail = sd(mean_Utail_rep))
  return(tail_data_summarized)
}

summarize_Atails_calculate_means <- function(tail_data2) {
  # function summarizing Atail length data for given dataset
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
                                             primer_name), summarise, N = length(Atail_length), mean_Atail = mean(Atail_length,
                                                                                                                  na.rm = T), sd_Atail = sd(Atail_length))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, primer_name), transform, se_Atail = sd_Atail/sqrt(N))
  return(tail_data_summarized)
}

summarize_Atails_calculate_medians <- function(tail_data2) {
  # function summarizing Atail length data for given dataset
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
                                             primer_name), summarise, N = length(Atail_length), mean_Atail = median(Atail_length,
                                                                                                                    na.rm = T), sd_Atail = sd(Atail_length))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, primer_name), transform, se_Atail = sd_Atail/sqrt(N))
  return(tail_data_summarized)
}


summarize_Utails_lengths <- function(tail_data2) {
  # function summarizing tail length data for given dataset
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(transcript, condition, cell_line,
                                             localization, primer_name, U_length), summarise, N = length(U_length))
  tail_data_summarized <- ddply(tail_data_summarized, .(transcript, condition,
                                                        cell_line, localization, primer_name), transform, freq = N/sum(N))
  return(tail_data_summarized)
}

summarize_tails_by_tail_type_collapse_short_relative <- function(tail_data2) {
  #function summarizing tails based on their length and type
  tail_data<-tail_data2
  tail_data$tail_length<-as.numeric(tail_data$tail_length)
  #calculate maximum possible tail length in given dataset
  max_tail_length<-summary(tail_data2$tail_length)[6]
  #bin tails based on their length
  if (max_tail_length>=60) {
    tail_data[tail_data$tail_length>=60,]$tail_length<-'60+'
  }
  if (max_tail_length>=50) {
    for (temp_len in seq(50,59,1)) {
      if(any(tail_data$tail_length==temp_len)) {
        tail_data[tail_data$tail_length==temp_len,]$tail_length<-'50-59'

      }
    }
  }
  if (max_tail_length>=40) {
    for (temp_len in seq(40,49,1)) {
      if(any(tail_data$tail_length==temp_len)) {
        tail_data[tail_data$tail_length==temp_len,]$tail_length<-'40-49'
      }
    }
  }
  if (max_tail_length>=30) {
    for (temp_len in seq(30,39,1)) {
      if(any(tail_data$tail_length==temp_len)) {
        tail_data[tail_data$tail_length==temp_len,]$tail_length<-'30-39'
      }
    }
  }
  if (max_tail_length>=20) {
    for (temp_len in seq(20,29,1)) {
      if(any(tail_data$tail_length==temp_len)) {
        tail_data[tail_data$tail_length==temp_len,]$tail_length<-'20-29'
      }
    }
  }
  for (temp_len in seq(10,19,1)) {
    if(any(tail_data$tail_length==temp_len)) {
      tail_data[tail_data$tail_length==temp_len,]$tail_length<-'10-19'
    }
  }

  for (temp_len in seq(1,9,1)) {
    if(any(tail_data$tail_length==temp_len)) {
      tail_data[tail_data$tail_length==temp_len,]$tail_length<-'1-9'
    }
  }
  #reorder tail length factor
  tail_data$tail_length<-as.factor(tail_data$tail_length)
  tail_data$tail_length<-factor(tail_data$tail_length,levels=c('0','1-9','10-19','20-29','30-39','40-49','50-59','60+'))
  #summarize:
  tail_data_summarized<-ddply(tail_data,.(transcript,condition,tail_length,cell_line,localization,primer_name,tail_type),summarise,N=length(tail_length))
  tail_data_summarized<-ddply(tail_data_summarized,.(transcript,condition,cell_line,localization,primer_name,tail_length),transform, freq=N/sum(N))
  return(tail_data_summarized)
}


summarize_uridylation_replicates2 <- function(tail_data2) {
  # function summarizing tails by tail type
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(condition, primer_name, uridylated,replicate), summarise, N = length(primer_name))
  tail_data_summarized <- ddply(tail_data_summarized, .(condition, primer_name,replicate), transform, freq = N/sum(N))
  tail_data_summarized <- ddply(tail_data_summarized, .(condition, primer_name, uridylated), transform, mean_freq = mean(freq), sd = sd(freq))
  return(tail_data_summarized)
}

summarize_uridylation_replicates3 <- function(tail_data2) {
  # function summarizing tails by tail type
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(condition,  uridylated,replicate), summarise, N = length(primer_name))
  tail_data_summarized <- ddply(tail_data_summarized, .(condition,replicate), transform, freq = N/sum(N))
  tail_data_summarized <- ddply(tail_data_summarized, .(condition,uridylated), transform, mean_freq = mean(freq), sd = sd(freq))
  return(tail_data_summarized)
}


summarize_uridylation_replicates4 <- function(tail_data2) {
  # function summarizing tails by tail type
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(condition,  uridylated2,replicate), summarise, N = length(primer_name))
  tail_data_summarized <- ddply(tail_data_summarized, .(condition,replicate), transform, freq = N/sum(N))
  tail_data_summarized <- ddply(tail_data_summarized, .(condition,uridylated2), transform, mean_freq = mean(freq), sd = sd(freq))
  return(tail_data_summarized)
}


summarize_uridylation_replicates5 <- function(tail_data2) {
  # function summarizing tails by tail type
  tail_data <- tail_data2
  tail_data_summarized <- ddply(tail_data, .(project_name,condition,  uridylated2,replicate), summarise, N = length(primer_name))
  tail_data_summarized <- ddply(tail_data_summarized, .(project_name,condition,replicate), transform, freq = N/sum(N))
  tail_data_summarized <- ddply(tail_data_summarized, .(condition,uridylated2), transform, mean_freq = mean(freq), sd = sd(freq))
  return(tail_data_summarized)
}


## FIGURES ##


tails_data_mapped_reporter_overexp_for_terminal <- tails_data_mapped[tails_data_mapped$transcript ==
                                                                       "REPORTERL1_overexp", ]
tails_data_mapped_reporter_overexp_for_terminal$condition <- as.character(tails_data_mapped_reporter_overexp_for_terminal$condition)
tails_data_mapped_reporter_overexp_for_terminal$condition <- as.factor(tails_data_mapped_reporter_overexp_for_terminal$condition)


tails_data_mapped_gapdh_overexp_for_terminal <- tails_data_mapped[tails_data_mapped$transcript ==
                                                                    "GAPDH" & tails_data_mapped$condition %in% c("CNTRL", "TUT7WT", "MOV10"), ]
tails_data_mapped_gapdh_overexp_for_terminal$condition <- as.character(tails_data_mapped_gapdh_overexp_for_terminal$condition)
tails_data_mapped_gapdh_overexp_for_terminal$condition <- as.factor(tails_data_mapped_gapdh_overexp_for_terminal$condition)


for (cond in levels(tails_data_mapped_reporter_overexp_for_terminal$condition)) {
  print(cond)
  terminal_nucleotides_cond_name = paste("terminal_nucleotides_reporter_overexp",
                                         cond, sep = "_")
  # temp2=eval(as.symbol(summary_tail_lengths_table_name))
  assign(terminal_nucleotides_cond_name, tails_data_mapped_reporter_overexp_for_terminal[tails_data_mapped_reporter_overexp_for_terminal$condition ==
                                                                                           cond, ]$terminal_nucleotides)
}

terminal_nucleotides_reporter_overexp = list(CTRL = terminal_nucleotides_reporter_overexp_MBP,
                                             TUT7WT = terminal_nucleotides_reporter_overexp_TUT7WT, MOV10 = terminal_nucleotides_reporter_overexp_MOV10,
                                             TUT4WT = terminal_nucleotides_reporter_overexp_TUT4WT, TUT7MUT = terminal_nucleotides_reporter_overexp_TUT7MT,
                                             TUT4MUT = terminal_nucleotides_reporter_overexp_TUT4MT, FAM46CMT = terminal_nucleotides_reporter_overexp_FAM46CMT, FAM46CWT = terminal_nucleotides_reporter_overexp_FAM46CWT,
                                             FTUT7WT = terminal_nucleotides_reporter_overexp_FTUT7WT,
                                             GLD2 = terminal_nucleotides_reporter_overexp_GLD2,
                                             PAPD5 = terminal_nucleotides_reporter_overexp_PAPD5,
                                             U6TUT = terminal_nucleotides_reporter_overexp_U6TUT)


terminal_nucleotides_reporter_overexp = list(CTRL = terminal_nucleotides_reporter_overexp_CNTRL,
                                             TUT7WT = terminal_nucleotides_reporter_overexp_TUT7WT, MOV10 = terminal_nucleotides_reporter_overexp_MOV10,
                                             TUT4WT = terminal_nucleotides_reporter_overexp_TUT4WT, TUT7MUT = terminal_nucleotides_reporter_overexp_TUT7MUT,
                                             TUT4MUT = terminal_nucleotides_reporter_overexp_TUT4MUT)

for (cond in levels(tails_data_mapped_gapdh_overexp_for_terminal$condition)) {
  print(cond)
  terminal_nucleotides_cond_name = paste("terminal_nucleotides_gapdh_overexp",
                                         cond, sep = "_")
  # temp2=eval(as.symbol(summary_tail_lengths_table_name))
  assign(terminal_nucleotides_cond_name, tails_data_mapped_gapdh_overexp_for_terminal[tails_data_mapped_gapdh_overexp_for_terminal$condition ==
                                                                                        cond, ]$terminal_nucleotides)
}

terminal_nucleotides_gapdh_overexp = list(CTRL = terminal_nucleotides_gapdh_overexp_CNTRL,
                                          TUT7WT = terminal_nucleotides_gapdh_overexp_TUT7WT, MOV10 = terminal_nucleotides_gapdh_overexp_MOV10)


pdf("fig_4a_reporter.pdf")
print(ggseqlogo(terminal_nucleotides_reporter_overexp, ncol = 2, method = "bits"))
dev.off()

pdf("fig_4a_gapdh.pdf")
print(ggseqlogo(terminal_nucleotides_gapdh_overexp, ncol = 2, method = "prob"))
dev.off()


pdf("fig_4B_old_style_all.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                              c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "FTUT7WT","TUT7WT", "TUT7MT", "TUT4WT", "TUT4MT", "MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION"))
print(plot_tail_lengths)
dev.off()

pdf("fig_4B_old_style_TUTMOV.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                              c("MBP", "TUT7WT", "TUT4WT", "MOV10","TUT7MT","TUT4MT","FTUT7WT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "FTUT7WT","TUT7WT", "TUT7MT", "TUT4WT", "TUT4MT", "MOV10"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION"))
print(plot_tail_lengths)
dev.off()

pdf("fig_4B_reporter_all.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                   c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_reporter_UTR_less_than_9010.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                               c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),y=mean_freq,
                                                                                                                                              fill = uridylated, colours = uridylated)) + geom_bar(position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()
ggsave(filename = "fig_4B_reporter_UTR_less_than_9010_new.pdf", plot = plot_tail_lengths)


pdf("fig_4B_reporter_UTR_less_than_9010_no_fam_mut_no_GLD.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                               c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CWT","FTUT7WT","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CWT","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()
ggsave(filename = "fig_4B_reporter_UTR_less_than_9010_no_fam_mut_no_GLD.pdf", plot = plot_tail_lengths)


#Create fig for UTR only LINE reads for all conditions except FTUT7, GLD2, FAM46MT

summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                               c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CWT","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CWT","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==1 & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()
ggsave(filename = "fig_4B_reporter_UTR_less_than_9010_no_fam_mut_no_GLD_no_FTUT7WT.pdf", plot = plot_tail_lengths)



pdf("fig_4B_old_style_reporter_UTR_less_than_9010_new.pdf",width=2048,height=2048)
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name,summarize_tails_by_tail_type(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                         c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))

plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION, UTR only LINE"))
print(plot_tail_lengths)
dev.off()
ggsave(filename = "fig_4B_old_style_reporter_UTR_less_than_9010_new.pdf", plot = plot_tail_lengths)

pdf("fig_4B_reporter_SV40bod_9010_9140.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_SV40body[tails_data_for_analysis_reporter_SV40body$condition %in%
                                                                                                                    c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_reporter_SV40end_more_than9140.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_SV40end[tails_data_for_analysis_reporter_SV40end$condition %in%
                                                                                                                   c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()



pdf("fig_4B_reporter_all.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                   c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_reporter_all_test.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_test2))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_reporter_all_hetero.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_test_reporter_all_tails[tails_test_reporter_all_tails$condition %in%
                                                                                                        c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_MYSZ.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_mysz))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
#summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
#summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
#                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(primer_name),
                                                                                                                                              fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_4B_H9.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_H9))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
#summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
#summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
#                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(primer_name),
                                                                                                                                              fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()



pdf("fig_4B_reporter_all_all_mapped.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true[tails_data_mapped$condition %in%
                                                                                                 c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true$transcript=="REPORTERL1_overexp", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$primer=='L1NGS1' & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                                                            fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                                                                 stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_reporter_all_all_mapped_ACTB.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="ACTB", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_4B_reporter_all_all_mapped_GAPDH.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="GAPDH", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


#fig_4B_reporter_all_all_mapped_GAPDH.pdf - no GLD2 no FTUT7WT no FAM46CMT
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CWT","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="GAPDH", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("GAPDH OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
ggsave(filename = "fig_4B_GAPDH_noGLD2_noFTUT7WT_noFAM46CMT.pdf", plot = plot_tail_lengths)


#fig_4B_reporter_all_all_mapped_ACTB.pdf - no GLD2 no FTUT7WT no FAM46CMT
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CWT","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="ACTB", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("GAPDH OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
ggsave(filename = "fig_4B_ACTB_noGLD2_noFTUT7WT_noFAM46CMT.pdf", plot = plot_tail_lengths)



#fig_4B_reporter_all_all_mapped_PABPC4.pdf - no GLD2 no FTUT7WT no FAM46CMT
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CWT","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="PABPC4", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("GAPDH OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
ggsave(filename = "fig_4B_PABPC4_noGLD2_noFTUT7WT_noFAM46CMT.pdf", plot = plot_tail_lengths)





pdf("fig_4B_reporter_all_all_mapped_ETV4.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="ETV4", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_reporter_all_all_mapped_SOGA2.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="SOGA2", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_4B_reporter_all_all_mapped_PABPC4.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_mapped_true_no_hetero_no_other_tails$transcript=="PABPC4", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_4B_reporter_TUTMOV.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                   c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_4B_rebut2.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_mapped_true_no_hetero_no_other_tails[(tails_data_mapped_true_no_hetero_no_other_tails$transcript=='REPORTERL1_overexp') & tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                          c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT", "FTUT7WT", "TUT7MT","TUT4WT","TUT4MT", "MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate=='1',], aes(x = as.factor(condition),
                                                                                                                                                fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                     stat = "identity") + scale_y_continuous() + xlab("condition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_4C.pdf")
summary_tail_types_table_name <- paste("repo_KD_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_KD[tails_data_for_analysis_reporter_KD$condition %in%
                                                                                                         c("CNTRLKD", "TUT7KD", "TUT4KD", "TUT4TUT7KD"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("CNTRLKD", "TUT7KD", "TUT4KD", "TUT4TUT7KD"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("DEPLETION"))
print(plot_tail_lengths)
dev.off()


pdf("fig_4D.pdf")
summary_tail_types_table_name <- paste("genomic_L1_summarized_tails_types_by_cell_lines")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$transcript ==
                                                                                                                         "ENDOL1" & tails_data_mapped_true_no_hetero_no_other_tails_NT1$cell_line %in%
                                                                                                                         c("293T", "H9", "HELAHA", "NPC", "PA1", "MYSZ"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$cell_line <- as.character(summary_tail_types_table$cell_line)
summary_tail_types_table[summary_tail_types_table$cell_line == "MYSZ", ]$cell_line = "MOUSE_TESTIS"
summary_tail_types_table$cell_line <- factor(summary_tail_types_table$cell_line,
                                             levels = c("293T", "HELAHA", "NPC", "PA1", "H9", "MOUSE_TESTIS"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(cell_line),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("cell_line") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("DEPLETION"))
print(plot_tail_lengths)
dev.off()

pdf("fig_4F.pdf")
summary_tail_types_table_name <- paste("genomic_L1_summarized_tails_types_by_condition_KD")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$condition %in%
                                                                                                                            c("CNTRLKD", "MOV10KD", "TUT4TUT7KD") & tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$transcript ==
                                                                                                                            "ENDOL1", ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("CNTRLKD", "TUT4TUT7KD", "MOV10KD"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("condition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("DEPLETION (PA-1)"))
print(plot_tail_lengths)
dev.off()


pdf("fig_4G.pdf")
summary_tail_types_table_name <- paste("genomic_L1_and_GAPDH_summarized_tails_types_by_subcellular_localization")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type_loc(tails_data_mapped_true_no_hetero_no_other_PA1_tails_NT))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(localization),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("localization") + ylab("fraction of transcripts") +
  scale_fill_grey() + facet_grid(. ~ transcript)
print(plot_tail_lengths)
dev.off()


pdf("fig_S5A.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                  c("MBP", "FTUT7WT","TUT7WT", "TUT7MT", "MOV10","TUT4WT","TUT4MT") & tails_data_for_analysis_reporter_overexp$uridylated ==
                                                                                                                  TRUE, ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "FTUT7WT", "TUT7WT", "TUT7MT", "MOV10","TUT4WT","TUT4MT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          se_utail, ymax = mean_Utail + se_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_S5A.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means_repl(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                       c("MBP", "FTUT7WT","TUT7WT", "TUT7MT", "MOV10") & tails_data_for_analysis_reporter_overexp$uridylated ==
                                                                                                                       TRUE, ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "FTUT7WT", "TUT7WT", "TUT7MT", "MOV10"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          se_utail, ymax = mean_Utail + se_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_S5A_2.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                  c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_for_analysis_reporter_overexp$uridylated ==
                                                                                                                  TRUE, ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT", "TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          se_utail, ymax = mean_Utail + se_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                              c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_for_analysis_reporter_UTR$uridylated ==
                                                                                                              TRUE, ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT", "TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          se_utail, ymax = mean_Utail + se_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)


pdf("fig_S5A_ACTB.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means_repl(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                              c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_mapped_true_no_hetero_no_other_tails$uridylated == TRUE & tails_data_mapped_true_no_hetero_no_other_tails$transcript=='ACTB',]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
#summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- as.factor(summary_mean_Utail_table$condition)
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          sd_utail, ymax = mean_Utail + sd_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_S5A_GAPDH.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means_repl(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                              c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_mapped_true_no_hetero_no_other_tails$uridylated == TRUE & tails_data_mapped_true_no_hetero_no_other_tails$transcript=='GAPDH',]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
#summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- as.factor(summary_mean_Utail_table$condition)
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          sd_utail, ymax = mean_Utail + sd_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_S5A_PABPC4.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means_repl(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                              c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_mapped_true_no_hetero_no_other_tails$uridylated == TRUE & tails_data_mapped_true_no_hetero_no_other_tails$transcript=='PABPC4',]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
#summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- as.factor(summary_mean_Utail_table$condition)
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          sd_utail, ymax = mean_Utail + sd_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_S5A_SOGA2.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means_repl(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                              c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_mapped_true_no_hetero_no_other_tails$uridylated == TRUE & tails_data_mapped_true_no_hetero_no_other_tails$transcript=='SOGA2',]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
#summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- as.factor(summary_mean_Utail_table$condition)
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          sd_utail, ymax = mean_Utail + sd_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_S5A_ETV4.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means_repl(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$condition %in%
                                                                                                                              c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_mapped_true_no_hetero_no_other_tails$uridylated == TRUE & tails_data_mapped_true_no_hetero_no_other_tails$transcript=='ETV4',]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
#summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- as.factor(summary_mean_Utail_table$condition)
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          sd_utail, ymax = mean_Utail + sd_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_S5A_repor_UTR.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Utails_calculate_means_repl(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                                   c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_for_analysis_reporter_UTR$uridylated == TRUE,]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
#summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- as.factor(summary_mean_Utail_table$condition)
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Utail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Utail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Utail -
                                                                                          sd_utail, ymax = mean_Utail + sd_utail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_SJ.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Atails_calculate_means(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                  c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") , ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT", "TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Atail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Atail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Atail -
                                                                                          se_Atail, ymax = mean_Atail + se_Atail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_SJ_Aonly.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Atails_calculate_means(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                  c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_for_analysis_reporter_overexp$Atail_length>0 , ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT", "TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Atail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Atail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Atail -
                                                                                          se_Atail, ymax = mean_Atail + se_Atail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_SJ_Aonly_median.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Atails_calculate_medians(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                    c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_for_analysis_reporter_overexp$Atail_length>0 , ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT", "TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Atail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Atail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Atail -
                                                                                          se_Atail, ymax = mean_Atail + se_Atail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()


pdf("fig_SJ_Aonly_median_UTR_only.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Atails_calculate_medians(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                                c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_for_analysis_reporter_UTR$Atail_length>0 , ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT", "TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Atail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Atail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Atail -
                                                                                          se_Atail, ymax = mean_Atail + se_Atail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()



pdf("fig_SJ_uridylated_only.pdf")
summary_mean_Utail_table_name <- paste("reporter_mean_Utails")
assign(summary_mean_Utail_table_name, summarize_Atails_calculate_means(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                  c("MBP", "TUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT") & tails_data_for_analysis_reporter_overexp$uridylated == TRUE , ]))
summary_mean_Utail_table <- eval(as.symbol(summary_mean_Utail_table_name))
summary_mean_Utail_table$condition <- as.character(summary_mean_Utail_table$condition)
summary_mean_Utail_table$condition <- factor(summary_mean_Utail_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT", "TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))
plot_tail_lengths <- ggplot(summary_mean_Utail_table, aes(x = as.factor(condition))) +
  geom_bar(aes(y = mean_Atail), position = position_dodge(), stat = "identity") +
  xlab("condition") + ylab("Atail mean length") + scale_fill_grey() + geom_errorbar(aes(ymin = mean_Atail -
                                                                                          se_Atail, ymax = mean_Atail + se_Atail), colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()

pdf("fig_S5B.pdf")
summary_Utail_lengths_table_name <- paste("reporter_overexp_summarized_Utails_lengths_by_condition")
assign(summary_Utail_lengths_table_name, summarize_Utails_lengths(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in% c("MBP", "MOV10", "TUT7WT", "TUT7MT", "TUT4WT", "TUT4MT","FTUT7WT"),]))
summary_Utail_lengths_table <- eval(as.symbol(summary_Utail_lengths_table_name))
summary_Utail_lengths_table$condition <- as.character(summary_Utail_lengths_table$condition)
summary_Utail_lengths_table$condition <- factor(summary_Utail_lengths_table$condition,
                                                levels = c("MBP", "FTUT7WT", "MOV10", "TUT7WT", "TUT7MT", "TUT4WT", "TUT4MT"))
plot_tail_lengths <- ggplot(summary_Utail_lengths_table, aes(x = as.factor(condition),
                                                             fill = U_length, colours = U_length)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                              stat = "identity") + xlab("cell line") + ylab("fraction of transcripts") + scale_fill_grey()
print(plot_tail_lengths)
dev.off()


pdf("fig_S5C-G.pdf")
for (transcript in c("ENDOL1")) {
  for (cell_line in c("H9", "HELAHA", "293T", "NPC", "PA1")) {
    summary_tail_lengths_table_name <- paste("_summarized_tails_lengths_by_type")
    assign(summary_tail_lengths_table_name, summarize_tails_by_tail_type_collapse_short_relative(tails_data_mapped_true_no_hetero_no_other_tails_NT1[tails_data_mapped_true_no_hetero_no_other_tails_NT1$transcript ==
                                                                                                                                                       transcript & tails_data_mapped_true_no_hetero_no_other_tails_NT1$cell_line ==
                                                                                                                                                       cell_line & tails_data_mapped_true_no_hetero_no_other_tails_NT1$primer_name ==
                                                                                                                                                       "L1NGS0" & tails_data_mapped_true_no_hetero_no_other_tails_NT1$tail_length >
                                                                                                                                                       0, ]))
    summary_tail_lengths_table <- eval(as.symbol(summary_tail_lengths_table_name))
    summary_tail_lengths_table$tail_type<-as.character(summary_tail_lengths_table$tail_type)
    summary_tail_lengths_table$tail_type<-factor(summary_tail_lengths_table$tail_type,levels=c("U_only","AU","A_only"))
    plot_tail_lengths <- ggplot(summary_tail_lengths_table, aes(x = as.factor(tail_length),
                                                                fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("Atail_length") + ylab("% of 3' ends") +
      facet_grid(. ~ cell_line) + scale_fill_manual(values = brewer.pal(11, "Spectral")[c(2,4,10)])
    print(plot_tail_lengths)
  }
}
dev.off()



pdf("fig_S5I.pdf")
for (condition in c("CNTRLKD", "TUT4TUT7KD", "MOV10KD")) {
  summary_tail_lengths_table_name <- paste("_summarized_tails_lengths_by_type")
  assign(summary_tail_lengths_table_name, summarize_tails_by_tail_type_collapse_short_relative(tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$condition ==
                                                                                                                                                        condition & tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$primer_name ==
                                                                                                                                                        "L1NGS0" & tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$tail_length >
                                                                                                                                                        0, ]))
  summary_tail_lengths_table <- eval(as.symbol(summary_tail_lengths_table_name))
  summary_tail_lengths_table$tail_type <- as.character(summary_tail_lengths_table$tail_type)
  summary_tail_lengths_table$tail_type <- factor(summary_tail_lengths_table$tail_type,
                                                 levels = c("U_only", "AU", "A_only"))
  plot_tail_lengths <- ggplot(summary_tail_lengths_table, aes(x = as.factor(tail_length),
                                                              fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                                 stat = "identity") + scale_y_continuous() + xlab("Atail_length") + ylab("% of 3' ends") +
    facet_grid(. ~ condition) + scale_fill_manual(values = brewer.pal(11, "Spectral")[c(2,
                                                                                        4, 10)])
  print(plot_tail_lengths)
}

dev.off()



uridylated_TUT7WT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT7WT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_CNTRL<-summary_tail_types_table[summary_tail_types_table$condition=='MBP' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_TUT4WT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT4WT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_MOV10<-summary_tail_types_table[summary_tail_types_table$condition=='MOV10' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_TUT7MT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT7MT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_TUT4MT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT4MT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_FAM46CWT<-summary_tail_types_table[summary_tail_types_table$condition=='FAM46CWT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_FAM46CMT<-summary_tail_types_table[summary_tail_types_table$condition=='FAM46CMT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_PAPD5<-summary_tail_types_table[summary_tail_types_table$condition=='PAPD5' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_U6TUT<-summary_tail_types_table[summary_tail_types_table$condition=='U6TUT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_GLD2<-summary_tail_types_table[summary_tail_types_table$condition=='GLD2' & summary_tail_types_table$uridylated==TRUE,6]
#uridylated_TUT4TUT7KD<-summary_tail_types_table[summary_tail_types_table$condition=='TUT4TUT7KD' & summary_tail_types_table$uridylated==TRUE,5]
#uridylated_OE_df<-list(CTRL=uridylated_CNTRL,TUT4WT=uridylated_TUT4WT,TUT7WT=uridylated_TUT7WT,MOV10=uridylated_MOV10,TUT7MT=uridylated_TUT7MT,TUT4MT=uridylated_TUT4MT,GLD2=uridylated_GLD2,PAPD5=uridylated_PAPD5,FAM46CWT=uridylated_FAM46CWT,FAM46CMT=uridylated_FAM46CMT,U6TUT=uridylated_U6TUT,check.rows=F)
uridylated_OE_df<-list(CTRL=uridylated_CNTRL,TUT4WT=uridylated_TUT4WT,TUT7WT=uridylated_TUT7WT,MOV10=uridylated_MOV10,TUT7MT=uridylated_TUT7MT,TUT4MT=uridylated_TUT4MT,GLD2=uridylated_GLD2,PAPD5=uridylated_PAPD5,FAM46CWT=uridylated_FAM46CWT,U6TUT=uridylated_U6TUT)
attributes(uridylated_OE_df) = list(names = names(uridylated_OE_df),
                                    row.names=1:max(length(uridylated_TUT4WT), length(uridylated_CNTRL)), class='data.frame')
uridylated_OE_df_melt<-melt(uridylated_OE_df)
aov_test<-aov(value ~ variable,uridylated_OE_df_melt)
TukeyHSD(aov_test)
tukey_4C<-TukeyHSD(aov_test)


#reporter parts analysis
tails_data_for_analysis_reporter_UTR <- tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$mapping_position>8950 & tails_data_for_analysis_reporter_overexp$mapping_position<9010,]


reporter_start=8800
reporter_end=9200
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_overexp[tails_data_for_analysis_reporter_overexp$condition %in%
                                                                                                                   c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CMT","FAM46CWT","FTUT7WT","GLD2","PAPD5","U6TUT","TUT7MT","TUT4MT") & tails_data_for_analysis_reporter_overexp$mapping_position>reporter_start & tails_data_for_analysis_reporter_overexp$mapping_position<reporter_end & tails_data_for_analysis_reporter_overexp$PCR_duplicates<2, ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CMT","FAM46CWT","GLD2","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
uridylated_TUT7WT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT7WT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_CNTRL<-summary_tail_types_table[summary_tail_types_table$condition=='MBP' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_TUT4WT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT4WT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_MOV10<-summary_tail_types_table[summary_tail_types_table$condition=='MOV10' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_TUT7MT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT7MT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_TUT4MT<-summary_tail_types_table[summary_tail_types_table$condition=='TUT4MT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_FAM46CWT<-summary_tail_types_table[summary_tail_types_table$condition=='FAM46CWT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_FAM46CMT<-summary_tail_types_table[summary_tail_types_table$condition=='FAM46CMT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_PAPD5<-summary_tail_types_table[summary_tail_types_table$condition=='PAPD5' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_U6TUT<-summary_tail_types_table[summary_tail_types_table$condition=='U6TUT' & summary_tail_types_table$uridylated==TRUE,6]
uridylated_GLD2<-summary_tail_types_table[summary_tail_types_table$condition=='GLD2' & summary_tail_types_table$uridylated==TRUE,6]
#uridylated_TUT4TUT7KD<-summary_tail_types_table[summary_tail_types_table$condition=='TUT4TUT7KD' & summary_tail_types_table$uridylated==TRUE,5]
#uridylated_OE_df<-list(CTRL=uridylated_CNTRL,TUT4WT=uridylated_TUT4WT,TUT7WT=uridylated_TUT7WT,MOV10=uridylated_MOV10,TUT7MT=uridylated_TUT7MT,TUT4MT=uridylated_TUT4MT,GLD2=uridylated_GLD2,PAPD5=uridylated_PAPD5,FAM46CWT=uridylated_FAM46CWT,FAM46CMT=uridylated_FAM46CMT,U6TUT=uridylated_U6TUT,check.rows=F)
uridylated_OE_df<-list(CTRL=uridylated_CNTRL,TUT4WT=uridylated_TUT4WT,TUT7WT=uridylated_TUT7WT,MOV10=uridylated_MOV10,TUT7MT=uridylated_TUT7MT,TUT4MT=uridylated_TUT4MT,GLD2=uridylated_GLD2,PAPD5=uridylated_PAPD5,FAM46CWT=uridylated_FAM46CWT,FAM46CMT=uridylated_FAM46CMT,U6TUT=uridylated_U6TUT)
attributes(uridylated_OE_df) = list(names = names(uridylated_OE_df),
                                    row.names=1:max(length(uridylated_TUT4WT), length(uridylated_CNTRL)), class='data.frame')
uridylated_OE_df_melt<-melt(uridylated_OE_df)
aov_test<-aov(value ~ variable,uridylated_OE_df_melt)
TukeyHSD(aov_test)
tukey_4C<-TukeyHSD(aov_test)




tails_dplyr_test <- tails_data_for_analysis_reporter_overexp %>% group_by(condition,tail_type)
tails_dplyr_aaa <- tails_dplyr_test %>% summarize(ile = n()) %>% mutate(freq = ile/sum(ile))
ggplot(tails_dplyr_aaa %>% filter(tail_type=='U_only'),aes(x=condition,y=freq,fill=tail_type)) + geom_bar(stat="identity")

tails_dplyr_test <- tails_data_for_analysis_reporter_overexp %>% group_by(condition,tail_sequence)
tails_dplyr_aaa <- tails_dplyr_test %>% summarize(ile = n()) %>% mutate(freq = ile/sum(ile))
ggplot(tails_dplyr_aaa %>% filter(tail_sequence=='AT'),aes(x=condition,y=freq)) + geom_bar(stat="identity")

tails_dplyr_test <- tails_data_for_analysis_reporter_UTR %>% group_by(condition,tail_type)
tails_dplyr_aaa <- tails_dplyr_test %>% summarize(ile = n()) %>% mutate(freq = ile/sum(ile))
ggplot(tails_dplyr_aaa %>% filter(tail_type=='no_tail'),aes(x=condition,y=freq,fill=tail_type)) + geom_bar(stat="identity")


for (cond in levels(tails_data_for_analysis_reporter_UTR$condition)) {
  print(cond)
  terminal_nucleotides_cond_name = paste("tails_data_for_analysis_reporter_UTR",
                                         cond, sep = "_")
  # temp2=eval(as.symbol(summary_tail_lengths_table_name))
  assign(terminal_nucleotides_cond_name, tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition ==
                                                                                           cond, ])
}



pdf("fig_4B_reporter_UTR_less_than_9010_no_fam_mut_no_GLD.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates2(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$condition %in%
                                                                                                               c("MBP", "TUT7WT", "TUT4WT", "MOV10","FAM46CWT","FTUT7WT","PAPD5","U6TUT","TUT7MT","TUT4MT"), ]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CWT","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = uridylated, colours = uridylated)) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()
ggsave(filename = "fig_4B_reporter_UTR_less_than_9010_no_fam_mut_no_GLD.pdf", plot = plot_tail_lengths)

pdf("fig_4B_reporter_UTR_less_than_9010_no_fam_mut_no_GLD.pdf")
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates3(tails_data_for_analysis_reporter_UTR))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
#summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
#summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
#                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CWT","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = as.factor(uridylated), colours = as.factor(uridylated))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                   stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)
dev.off()
ggsave(filename = "fig_4B_reporter_UTR_less_than_9010_no_fam_mut_no_GLD.pdf", plot = plot_tail_lengths)

#old style
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_tails_by_tail_type(tails_data_for_analysis_reporter_UTR))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
#summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
#summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
#                                             levels = c("MBP", "FTUT7WT","TUT7WT", "TUT7MT", "TUT4WT", "TUT4MT", "MOV10"))
plot_tail_lengths <- ggplot(summary_tail_types_table, aes(x = as.factor(condition),
                                                          fill = tail_type, colours = tail_type)) + geom_bar(aes(y = freq), position = position_stack(),
                                                                                                             stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION"))
print(plot_tail_lengths)


tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='U_heterogenous',]$uridylated=TRUE
tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='AU',]$uridylated=TRUE
tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='U_only',]$uridylated=TRUE
tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='UA',]$uridylated=TRUE
tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='A_only',]$uridylated=FALSE
tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='A_heterogenous',]$uridylated=FALSE
tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='other',]$uridylated=FALSE
tails_data_for_analysis_reporter_UTR_all_mapped[tails_data_for_analysis_reporter_UTR_all_mapped$tail_type=='no_tail',]$uridylated=FALSE

#tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='U_heterogenous',]$uridylated=TRUE
tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='AU',]$uridylated=TRUE
tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='U_only',]$uridylated=TRUE
#tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='UA',]$uridylated=TRUE
tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='A_only',]$uridylated=FALSE
#tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='A_heterogenous',]$uridylated=FALSE
#tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='other',]$uridylated=FALSE
tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$tail_type=='no_tail',]$uridylated=FALSE



summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates3(tails_data_for_analysis_reporter_overexp))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
#summary_tail_types_table$condition <- as.character(summary_tail_types_table$condition)
#summary_tail_types_table$condition <- factor(summary_tail_types_table$condition,
#                                             levels = c("MBP", "TUT7WT","FTUT7WT","TUT7MT","TUT4WT", "TUT4MT","MOV10","FAM46CWT","PAPD5","U6TUT"))



plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = as.factor(uridylated), colours = as.factor(uridylated))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                                         stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("OVEREXPRESSION, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)


tails_data_for_analysis_reporter_KD2<-tails_data_for_analysis_reporter_KD
tails_data_for_analysis_reporter_KD2[tails_data_for_analysis_reporter_KD2$tail_length==1,]$tail_length<-0
tails_data_for_analysis_reporter_KD2[tails_data_for_analysis_reporter_KD2$Utail_length<2 & tails_data_for_analysis_reporter_KD2$tail_type=='U_only',]$uridylated<-FALSE
tails_data_for_analysis_reporter_KD2[tails_data_for_analysis_reporter_KD2$tail_length==0,]$uridylated<-FALSE


summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates3(tails_data_for_analysis_reporter_KD2))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))

plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = as.factor(uridylated), colours = as.factor(uridylated))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                                         stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)

summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates4(tails_data_mapped_true[tails_data_mapped_true$transcript=='REPORTERL1KD',]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))

plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                                         stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)

summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates4(tails_data_for_analysis_reporter_KD))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))

plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                                         stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)


tails_data_for_analysis_reporter_KD3<-tails_data_for_analysis_reporter_KD
tails_data_for_analysis_reporter_KD3[tails_data_for_analysis_reporter_KD3$condition=='CNTRLKD_HGC',]$condition<-'CNTRLKD'


summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates4(tails_data_for_analysis_reporter_KD[tails_data_for_analysis_reporter_KD$project_name=='Tailsq_2',]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))

plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated==TRUE & summary_tail_types_table$replicate==1,], aes(x = as.factor(condition),
                                                                                                                                              fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                                         stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
print(plot_tail_lengths)

test <- tails_data_for_analysis_reporter_KD  %>% select(V1,tail_sequence,tail_type,Atail,Atail_length,Utail_length,tail_length,tail_source,transcript,cell_line,localization,condition,replicate,primer_name,project_name,mapping_position,exp_type,CTGAC_R5,terminal_nucleotides,uridylated,uridylated2,tailed,mapped) %>% group_by(project_name,condition,replicate,transcript,primer_name)

test <- test %>% filter(mapping_position>8500 & mapping_position<8700)



test2 <- ggplot(test,aes(x=mapping_position,group=condition,colour=condition)) + stat_bin(binwidth=1,aes(x=mapping_position,y=..ncount..),geom="line")
test2



test_ovr <- tails_data_for_analysis_reporter_UTR  %>% select(V1,tail_sequence,tail_type,Atail,Atail_length,Utail_length,tail_length,tail_source,transcript,cell_line,localization,condition,replicate,primer_name,project_name,mapping_position,exp_type,CTGAC_R5,terminal_nucleotides,uridylated,uridylated2,tailed,mapped) %>% group_by(project_name,condition,replicate,transcript,primer_name)

test_ovr <- test_ovr %>% filter(mapping_position>8900 & mapping_position<9100)



test_ovr2 <- tails_data_for_analysis_reporter_overexp  %>% select(tail_type,Utail_length,tail_length,transcript,replicate,condition,replicate,primer_name,project_name,mapping_position,uridylated,uridylated2) %>% group_by(project_name,condition,transcript,primer_name,uridylated,replicate)
test_ovr2 <- test_ovr2 %>% filter(mapping_position>8900 & mapping_position<9200)
test_ovr2 <- test_ovr2 %>% filter(condition == 'CNTRL')
#%>% filter(project_name=='Tailseq_2')

#test_ovr2 <- test_ovr2 %>% summarize(uridyl = n()) %>% ungroup() %>% group_by(condition,replicate) %>% mutate(freq_repl = uridyl/sum(uridyl))
#test_ovr2 <- test_ovr2 %>%  group_by(condition,uridylated) %>% summarize(mean_freq=mean(freq_repl),sd=sd(freq_repl))
test_ovr2

test_ovr_plot <- ggplot(test_ovr2[test_ovr2$uridylated==1,],aes(x=condition,y=mean_freq)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
test_ovr_plot



test2ovr <- ggplot(test_ovr2,aes(x=mapping_position,group=condition,colour=condition)) + stat_bin(binwidth=1,aes(x=mapping_position,y=..count..),geom="line")
test2ovr

test2ovr + facet_grid(project_name ~ .)


test <- tails_data_for_analysis_reporter_KD  %>% select(V1,tail_sequence,tail_type,Atail,Atail_length,Utail_length,tail_length,tail_source,transcript,cell_line,localization,condition,replicate,primer_name,project_name,mapping_position,exp_type,CTGAC_R5,terminal_nucleotides,uridylated,uridylated2,tailed,mapped) %>% group_by(project_name,condition,replicate,transcript,primer_name)

test <- tails_data_for_analysis_reporter_KD_filtr  %>% select(V1,tail_sequence,tail_type,Atail,Atail_length,Utail_length,tail_length,tail_source,transcript,cell_line,localization,condition,replicate,primer_name,project_name,R3_mapping_position,mapping_position,exp_type,CTGAC_R5,terminal_nucleotides,uridylated,uridylated2,tailed,mapped) %>% group_by(project_name,condition,replicate,transcript,primer_name)



test <- test %>% filter(R3_mapping_position>8500 & R3_mapping_position<8700)
test <- test %>% filter(condition=='TUT4KD')

test2 <- ggplot(test,aes(group=condition,colour=condition)) + stat_bin(binwidth=1,aes(x=R3_mapping_position,y=..ncount..),geom="line")
test2

test2 + facet_grid(project_name ~ .)



ovr_plot <- ggplot(test_ovr,aes(x=))

summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates4(tails_data_for_analysis_reporter_UTR[tails_data_for_analysis_reporter_UTR$project_name=='Tailseq_5',]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                                                              fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                                                           stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
#print(plot_tail_lengths)
plot_tail_lengths


summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_data_for_analysis_reporter_UTR))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
#summary_tail_types_table <- summary_tail_types_table %>% group_by(project_name,condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                       fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
#print(plot_tail_lengths)
plot_tail_lengths


tails_data_for_analysis_reporter_KD_filtr <- tails_data_for_analysis_reporter_KD[tails_data_for_analysis_reporter_KD$project_name=='Tailseq_5' & (tails_data_for_analysis_reporter_KD$R3_mapping_position>8660 | is.na(tails_data_for_analysis_reporter_KD$R3_mapping_position)),]

tails_data_for_analysis_reporter_KD_filtr[tails_data_for_analysis_reporter_KD_filtr$condition %in% c("CNTRLKD_MGC") & tails_data_for_analysis_reporter_KD_filtr$replicate==1,]$replicate<-'4'
tails_data_for_analysis_reporter_KD_filtr[tails_data_for_analysis_reporter_KD_filtr$condition %in% c("CNTRLKD_MGC") & tails_data_for_analysis_reporter_KD_filtr$replicate==2,]$replicate<-'5'
tails_data_for_analysis_reporter_KD_filtr[tails_data_for_analysis_reporter_KD_filtr$condition %in% c("CNTRLKD_MGC") & tails_data_for_analysis_reporter_KD_filtr$replicate==3,]$replicate<-'6'

tails_data_for_analysis_reporter_KD_filtr[tails_data_for_analysis_reporter_KD_filtr$condition %in% c("CNTRLKD_HGC") & tails_data_for_analysis_reporter_KD_filtr$replicate==1,]$replicate<-'7'
tails_data_for_analysis_reporter_KD_filtr[tails_data_for_analysis_reporter_KD_filtr$condition %in% c("CNTRLKD_HGC") & tails_data_for_analysis_reporter_KD_filtr$replicate==2,]$replicate<-'8'
tails_data_for_analysis_reporter_KD_filtr[tails_data_for_analysis_reporter_KD_filtr$condition %in% c("CNTRLKD_HGC") & tails_data_for_analysis_reporter_KD_filtr$replicate==3,]$replicate<-'9'
tails_data_for_analysis_reporter_KD_filtr[tails_data_for_analysis_reporter_KD_filtr$condition %in% c("CNTRLKD_HGC","CNTRLKD_MGC"),]$condition<-'CNTRLKD'


summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_data_for_analysis_reporter_KD_filtr))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                       fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
#print(plot_tail_lengths)
plot_tail_lengths


tails_test <- tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD[tails_data_mapped_true_no_hetero_no_other_PA1_tails_KD$transcript=='ENDOL1',]
summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_test))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
#summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                       fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN PA1 ENDO")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
#print(plot_tail_lengths)
plot_tail_lengths


summary_tail_types_table


summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript=='ACTB' & tails_data_mapped_true_no_hetero_no_other_tails$exp_type=='OVR',]))
summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                       fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                    stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
  scale_fill_grey() + ggtitle(paste("KNOCKDOWN mod, LINE1 UTR only reporter")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
#print(plot_tail_lengths)
plot_tail_lengths


pdf("transcripts_overexp.pdf",width=12)
for (transcript in levels(as.factor(tails_data_mapped_true_no_hetero_no_other_tails$transcript))) {
  print(transcript)
  summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
  assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript==transcript & tails_data_mapped_true_no_hetero_no_other_tails$exp_type=='OVR',]))
  summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
  if (nrow(summary_tail_types_table)>1) {
    print(summary_tail_types_table)
    summary_test <- summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,] %>% select(condition,freq,replicate,project_name) %>% reshape(timevar="condition",direction="wide",idvar=c("replicate","project_name"))
    write.table(summary_test,file=paste(transcript,"_ovr_uridylation_summarized.tsv",sep=""),row.names=F,sep="\t")
    summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


    plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                           fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                        stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
      scale_fill_grey() + ggtitle(paste("OVR",transcript,sep=",")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
    print(plot_tail_lengths)
    #plot_tail_lengths
  }
  for (project in levels(as.factor(tails_data_mapped_true_no_hetero_no_other_tails$project_name))) {
    summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
    assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript==transcript & tails_data_mapped_true_no_hetero_no_other_tails$exp_type=='OVR' & tails_data_mapped_true_no_hetero_no_other_tails$project_name==project,]))
    summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
    if (nrow(summary_tail_types_table)>1) {
      print(summary_tail_types_table)
      summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


      plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                             fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                          stat = "identity") + scale_y_continuous() + xlab("ondition") + ylab("fraction of transcripts") +
        scale_fill_grey() + ggtitle(paste("OVR",transcript,project,sep=",")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
      print(plot_tail_lengths)
      #plot_tail_lengths
    }

  }

}

dev.off()

pdf("transcripts_kd.pdf",width=12)
for (transcript in levels(as.factor(tails_data_mapped_true_no_hetero_no_other_tails$transcript))) {
  print(transcript)
  summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
  assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript==transcript & tails_data_mapped_true_no_hetero_no_other_tails$exp_type=='KD',]))
  summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
  if (nrow(summary_tail_types_table)>1) {
    print(summary_tail_types_table)
    summary_test <- summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,] %>% select(condition,freq,replicate,project_name) %>% reshape(timevar="condition",direction="wide",idvar=c("replicate","project_name"))
    write.table(summary_test,file=paste(transcript,"_kd_uridylation_summarized.tsv",sep=""),row.names=F,sep="\t")
    summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


    plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                           fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                        stat = "identity") + scale_y_continuous() + xlab("condition") + ylab("fraction of transcripts") +
      scale_fill_grey() + ggtitle(paste("KD",transcript,"all seqs",sep=",")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
    print(plot_tail_lengths)
    #plot_tail_lengths
  }
  for (project in levels(as.factor(tails_data_mapped_true_no_hetero_no_other_tails$project_name))) {
    summary_tail_types_table_name <- paste("repo_over_summarized_tails_types_by_condition")
    assign(summary_tail_types_table_name, summarize_uridylation_replicates5(tails_data_mapped_true_no_hetero_no_other_tails[tails_data_mapped_true_no_hetero_no_other_tails$transcript==transcript & tails_data_mapped_true_no_hetero_no_other_tails$exp_type=='KD' & tails_data_mapped_true_no_hetero_no_other_tails$project_name==project,]))
    summary_tail_types_table <- eval(as.symbol(summary_tail_types_table_name))
    if (nrow(summary_tail_types_table)>1) {
      print(summary_tail_types_table)

      summary_tail_types_table <- summary_tail_types_table %>% group_by(condition,uridylated2) %>% summarise(mean_freq=mean(mean_freq),sd=mean(sd))


      plot_tail_lengths <- ggplot(summary_tail_types_table[summary_tail_types_table$uridylated2==TRUE,], aes(x = as.factor(condition),
                                                                                                             fill = as.factor(uridylated2), colours = as.factor(uridylated2))) + geom_bar(aes(y = mean_freq), position = position_stack(),
                                                                                                                                                                                          stat = "identity") + scale_y_continuous() + xlab("condition") + ylab("fraction of transcripts") +
        scale_fill_grey() + ggtitle(paste("KD",transcript,project,sep=",")) + geom_errorbar(aes(ymin = mean_freq -sd, ymax = mean_freq + sd),colour = "black", width = 0.1, position = position_dodge(0.9))
      print(plot_tail_lengths)
      #plot_tail_lengths
    }

  }

}

dev.off()

all_data <- tails_data_mapped_true_no_hetero_no_other_tails %>% select(V1,tail_sequence,tail_type,Atail,Atail_length,Utail_length,tail_length,tail_source,transcript,cell_line,localization,condition,replicate,primer_name,project_name,mapping_position,exp_type,CTGAC_R5,terminal_nucleotides,uridylated,uridylated2,tailed,mapped) %>% group_by(project_name,condition,replicate,transcript,primer_name,cell_line,uridylated)
all_data <- all_data %>% filter(tail_length <= 64)


analyze_uridylation <- function(dataset,transcript2,exp_type2,mapping_position_min=NA,mapping_position_max=NA,project=NA,facet_projects=FALSE,conditions=NA) {
 
  print(mapping_position_min)
  print(mapping_position_max)
  print(transcript2)
  print(exp_type2)
  test_trans <- dataset %>% filter(transcript==transcript2,exp_type==exp_type2)
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% filter(mapping_position>mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% filter(mapping_position<mapping_position_max)
  }
  if(!is.na(project)) {
    test_trans <- test_trans %>% filter(project_name==project)
  }
  if(length(conditions)>0) {
    test_trans <- test_trans %>% filter(condition %in% conditions)
  }
  print(test_trans)
  test_trans <- test_trans %>% group_by(condition,replicate,project_name,uridylated2)
  print(test_trans)
  if (facet_projects == FALSE) {
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }
  else {
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }
  print(test_trans)
  plot <- test_trans %>% filter(uridylated2==TRUE) %>% ggplot(aes(x=condition,y=mean_freq_urid)) + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin =  mean_freq_urid - sd_urid, ymax = mean_freq_urid + sd_urid),colour = "black", width = 0.1, position = position_dodge(0.9)) + geom_jitter(aes(y=freq_urid)) + ggtitle(paste(exp_type2,transcript2,"uridylation"))
  if (facet_projects==TRUE) {
    plot <- plot + facet_grid (project_name ~ .)
  }
  print(plot) 
  
}

all_data_reporter_overexp_ACTB <- all_data %>% filter(transcript=='ACTB',exp_type=="OVR")

all_data_reporter_overexp2_ACTB <- all_data_reporter_overexp_ACTB %>% group_by(condition,replicate,project_name,uridylated2)

all_data_reporter_overexp3_ACTB <- all_data_reporter_overexp2_ACTB %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
all_data_reporter_overexp3_ACTB %>% filter(uridylated2==TRUE) %>% ggplot(aes(x=condition)) + geom_bar(aes(y=mean_freq_urid),stat="identity",position="dodge") + geom_errorbar(aes(ymin =  mean_freq_urid - sd_urid, ymax = mean_freq_urid + sd_urid),colour = "black", width = 0.1, position = position_dodge(0.9)) + geom_jitter(aes(y=freq_urid)) + facet_grid(project_name ~ .)


pdf("mapping_trans.pdf",width=12)
test_trans_select <- tails_data_mapped_true_no_hetero_no_other_tails  %>% select(tail_type,Utail_length,tail_length,transcript,replicate,condition,replicate,primer_name,project_name,mapping_position,R3_mapping_position,uridylated,uridylated2,exp_type) %>% group_by(project_name,condition,transcript,primer_name,uridylated,replicate)
for (trans in levels(as.factor(tails_data_mapped_true_no_hetero_no_other_tails$transcript))) {
  print(trans)
  test_trans2 <- test_trans_select %>% filter(transcript==trans)
  test_trans2 <- test_trans2 %>% filter(exp_type=="OVR")
  mapping_summary <- summary(test_trans2$R3_mapping_position)
  min_map = mapping_summary[2]
  max_map = mapping_summary[5]
  test_trans2 <- test_trans2  %>% filter(mapping_position>min_map & mapping_position<max_map)
  test_trans2$condition <- as.character(test_trans2$condition)
  test_trans2$condition <- as.factor(test_trans2$condition)
  for (cond in levels(as.factor(test_trans2$condition))) {
    print(cond)
    test_trans2_cond <- test_trans2 %>% filter(condition == cond)
    test_plot <- ggplot(test_trans2_cond,aes(x=mapping_position,group=condition,colour=condition)) + stat_bin(binwidth=1,aes(x=mapping_position,y=..count..),geom="line") +  ggtitle(paste("OVR",trans,cond,sep=","))
    test_plot <- test_plot + facet_grid(project_name ~ .)
    print(test_plot)
  }
}
dev.off()

pdf("mapping_trans_KD.pdf",width=12)
test_trans_select <- tails_data_mapped_true_no_hetero_no_other_tails  %>% select(tail_type,Utail_length,tail_length,transcript,replicate,condition,replicate,primer_name,project_name,mapping_position,R3_mapping_position,uridylated,uridylated2,exp_type) %>% group_by(project_name,condition,transcript,primer_name,uridylated,replicate)
for (trans in levels(as.factor(tails_data_mapped_true_no_hetero_no_other_tails$transcript))) {
  print(trans)
  test_trans2 <- test_trans_select %>% filter(transcript==trans)
  test_trans2 <- test_trans2 %>% filter(exp_type=="KD")
  mapping_summary <- summary(test_trans2$R3_mapping_position)
  min_map = mapping_summary[2]
  max_map = mapping_summary[5]
  test_trans2 <- test_trans2  %>% filter(mapping_position>min_map & mapping_position<max_map)
  test_trans2$condition <- as.character(test_trans2$condition)
  test_trans2$condition <- as.factor(test_trans2$condition)
  for (cond in levels(as.factor(test_trans2$condition))) {
    print(cond)
    test_trans2_cond <- test_trans2 %>% filter(condition == cond)
    test_plot <- ggplot(test_trans2_cond,aes(x=mapping_position,group=condition,colour=condition)) + stat_bin(binwidth=1,aes(x=mapping_position,y=..count..),geom="line") +  ggtitle(paste("KD",trans,cond,sep=","))
    test_plot <- test_plot + facet_grid(project_name ~ .)
    print(test_plot)
  }
}
dev.off()
