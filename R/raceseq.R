#' Title
#'
#' @param dataset 
#'
#' @return
#' @export
#'
#' @examples
process_data <- function(dataset) {
  processed_dataset <- dataset 
  processed_dataset$tailed <- processed_dataset$tail_length > 0  #mark tailed reads
  processed_dataset$mapped <- processed_dataset$mapping_position != -1  #mark mapped reads
  processed_dataset_mapped <- processed_dataset %>% filter(mapped != 0)  #discard unmapped reads
  
  # fix of flase_no_tail - should be assigned no tail cause CTGAC was found in the
  # clipped fragment - to be fixed in python script
  #tails_data_mapped$tail_type <- as.character(tails_data_mapped$tail_type)
  #tails_data_mapped[tails_data_mapped$tail_type == "false_no_tail_no_CTGAC", ]$tail_type <- "no_tail"
  
  processed_dataset_mapped$ref_name_R5 <- as.character(processed_dataset_mapped$ref_name_R5)
  processed_dataset_mapped$ref_name_R3 <- as.character(processed_dataset_mapped$ref_name_R3)
  # tails_data_mapped_same_ref<-tails_data_mapped[tails_data_mapped$ref_name_R5==tails_data_mapped$ref_name_R3,]
  
  # mark uridylated reads
  processed_dataset_mapped$uridylated2 <- FALSE
  processed_dataset_mapped[processed_dataset_mapped$Utail_length > 0, ]$uridylated2 = TRUE
  
  #convert T to U in terminal nucleotides (for seqlogo)
  processed_dataset_mapped$terminal_nucleotides<-as.character(processed_dataset_mapped$terminal_nucleotides)
  processed_dataset_mapped$terminal_nucleotides<-gsub("T","U",processed_dataset_mapped$terminal_nucleotides)
  
  #summary(as.factor(processed_dataset_mapped$CTGAC_R5))
  
 # colnames(processed_dataset_mapped)
  #head(processed_dataset_mapped[,c(1,2,3,12,24,25,26,27,33,34,37,38)])
  
  # in further analyses use only those read which got CTGAC delimiter identified in
  # the clipped fragment
  # for reporter analyses all mapped reads got CTGAC_R5 variable = 1 (because of short reads) so they will be included in the analysis
  processed_dataset_mapped_true <- processed_dataset_mapped[processed_dataset_mapped$CTGAC_R5 > 0, ]
  processed_dataset_mapped_true$ref_name = processed_dataset_mapped_true$ref_name_R5  #use ref_name_R5 as ref_name
  
  #head(processed_dataset_mapped_true[,c(1,2,3,12,24,25,26,27,33,34,37,38)])
  
  #head(tails_data_mapped[tails_data_mapped$tail_source=='tailseq_clip_clip_R3_shorter_than_tailseq',c(1,2,26,27,40)],30)
  
  processed_dataset_mapped_true$tail_type = as.character(processed_dataset_mapped_true$tail_type)
  
  #head(processed_dataset_mapped_true[grep("A_mix",processed_dataset_mapped_true$tail_type_mixed),])
  
  #treat all mixed (heterogenous) as true tails
  #tails_data_mapped_true[grep("AU_mixed",tails_data_mapped_true$tail_type_mixed),]$tail_type = 'AU'
  #tails_data_mapped_true[grep("^A_mixed",tails_data_mapped_true$tail_type_mixed),]$tail_type = 'A_only'
  #tails_data_mapped_true[grep("^U_mixed",tails_data_mapped_true$tail_type_mixed),]$tail_type = 'U_only'
  processed_dataset_mapped_true[processed_dataset_mapped_true$tail_type=='AU',]$uridylated2 = TRUE
  processed_dataset_mapped_true[processed_dataset_mapped_true$tail_type=='U_only',]$uridylated2 = TRUE
  #if clipping was shorter than tailseq - treat as other, as it is probably stretch of As in the reference seq
  #tails_data_mapped_true[tails_data_mapped_true$tail_source=='tailseq_clip_clip_R3_shorter_than_tailseq',]$tail_type = 'no_tail'
  #tails_data_mapped_true[tails_data_mapped_true$tail_source=='tailseq_clip_clip_R3_shorter_than_tailseq',]$uridylated = FALSE
  
  
  # remove heterogenous tails from analysis
  processed_dataset_mapped_true_no_hetero = processed_dataset_mapped_true[-grep("hetero", processed_dataset_mapped_true$tail_type),
                                                            ]
  
  #head(processed_dataset_mapped_true_no_hetero,20)
  # remove other type tails (for which we can suspect they are not tails but rather origin from improper mapping/repeatmasker) from the analysis
  processed_dataset_mapped_true_no_hetero_no_other = processed_dataset_mapped_true_no_hetero[-grep("other",
                                                                                                   processed_dataset_mapped_true_no_hetero$tail_type), ]
  
  
  # treat all AG,UG or UA tails as other_no_tail
  
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$tail_type ==
                                              "AG", ]$tail_type <- "other_no_tail"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$tail_type ==
                                              "UG", ]$tail_type <- "other_no_tail"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$tail_type ==
                                              "UA", ]$tail_type <- "other_no_tail"
  processed_dataset_mapped_true_no_hetero_no_other <- processed_dataset_mapped_true_no_hetero_no_other[-grep("other",
                                                                                                             processed_dataset_mapped_true_no_hetero_no_other$tail_type), ] #remove all other from analysis
  
  
  # create classes for A-tail lengths (0,1,2-5,6-10,11-20,21-30,30+)
  processed_dataset_mapped_true_no_hetero_no_other$A_length = ""
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length ==
                                              0, ]$A_length = "0"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length ==
                                              1, ]$A_length = "1"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(2, 5, 1), ]$A_length = "2-5"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(6, 10, 1), ]$A_length = "6-10"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(11, 20, 1), ]$A_length = "11-20"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length %in%
                                              seq(21, 30, 1), ]$A_length = "21-30"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Atail_length >
                                              30, ]$A_length = "30+"
  processed_dataset_mapped_true_no_hetero_no_other$A_length <- factor(processed_dataset_mapped_true_no_hetero_no_other$A_length,
                                                               levels = c("0", "1", "2-5", "6-10", "11-20", "21-30", "30+"))
  
  # create classes for U-tail lengths (0,1,2,3-5,6-10,10+)
  processed_dataset_mapped_true_no_hetero_no_other$U_length = ""
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length ==
                                              0, ]$U_length = "0"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length ==
                                              1, ]$U_length = "1"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length ==
                                              2, ]$U_length = "2"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length %in%
                                              seq(3, 5, 1), ]$U_length = "3-5"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length %in%
                                              seq(6, 10, 1), ]$U_length = "6-10"
  processed_dataset_mapped_true_no_hetero_no_other[processed_dataset_mapped_true_no_hetero_no_other$Utail_length >
                                              10, ]$U_length = "10+"
  processed_dataset_mapped_true_no_hetero_no_other$U_length <- factor(processed_dataset_mapped_true_no_hetero_no_other$U_length,
                                                               levels = c("10+", "6-10", "3-5", "2", "1", "0"))
  
  
  # modify levels of tail_types to have U_only,A-only,AU or no_tail
  processed_dataset_mapped_true_no_hetero_no_other_tails <- processed_dataset_mapped_true_no_hetero_no_other
  processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type <- as.character(processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type)
  processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type <- factor(processed_dataset_mapped_true_no_hetero_no_other_tails$tail_type,
                                                                      levels = c("U_only", "AU", "no_tail", "A_only"))
  
  
  
  
  all_data <- processed_dataset_mapped_true_no_hetero_no_other_tails %>% select(V1,tail_sequence,tail_type,Atail,Atail_length,Utail_length,tail_length,tail_source,transcript,cell_line,localization,condition,replicate,primer_name,project_name,mapping_position,exp_type,CTGAC_R5,terminal_nucleotides,uridylated,uridylated2,tailed,mapped) %>% dplyr::group_by(project_name,condition,replicate,transcript,primer_name,cell_line,uridylated)
  all_data <- all_data %>% filter(tail_length <= 64)
  return(all_data)
}

#' Analyze uridylation
#'
#' @param dataset               - dataset to analyze
#' @param transcript2           - transcript to analyze
#' @param exp_type2             - experiment type (OVR,KD,LEAP,NT)
#' @param mapping_position_min  - min position of mapping in the reference
#' @param mapping_position_max  - max position of mapping in the reference
#' @param project               - project (project_name) to include in the analysis (character of character vector)
#' @param facet_projects        - make facets based on projects (T/F)
#' @param conditions            - filter for selected conditions (character vector)
#' @param include_jitter        - include jitter dots on bar plot (T/F)
#'
#' @return                      - list containing calculated uridylation data (calculated output), Dunn test results and lot 
#' @export
#'
#' @examples
#' analyze_uridylation(all_data,"ACTB","OVR",project=c("Tailseq_4","Tailseq_5"),facet_projects==TRUE,conditions=c("CNTRL","TUT4WT","TUT7WT","MOV10"))
analyze_uridylation <- function(dataset,transcript2,exp_type2,mapping_position_min=NA,mapping_position_max=NA,project=NA,facet_projects=FALSE,conditions=NA,include_jitter=TRUE) {
  
  
  output = list() #list for storing output
  #print(mapping_position_min)
 # print(mapping_position_max)
  #print(transcript2)
  #print(exp_type2)
  #get filtered data from input dataset
  #first, filter by transcript and experiment type (OVR,KD,LEAP)
  test_trans <- dataset %>% dplyr::filter(transcript==transcript2,exp_type==exp_type2) 
  plot_title <- paste(exp_type2,transcript2,"uridylation")
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% dplyr::filter(mapping_position>mapping_position_min)
    plot_title <- paste(plot_title,"min map pos: ",mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% dplyr::filter(mapping_position<mapping_position_max)
    plot_title <- paste(plot_title,"max map pos: ",mapping_position_max)
  }
  #if project is specified - used for filtering
  if(!is.na(project)) {
    if (length(project==1)) {
      test_trans <- test_trans %>% dplyr::filter(project_name==project)
    } else {
      test_trans <- test_trans %>% dplyr::filter(project_name %in% project)
    }
  }
  #if conditions are specified - use for filtering
  if(length(conditions)>0 & !is.na(conditions)) {
    test_trans <- test_trans %>% dplyr::filter(condition %in% conditions)
  }
  #print(test_trans)
  #group data by co
  test_trans <- test_trans %>% dplyr::group_by(condition,replicate,project_name,uridylated2)
  #print(test_trans)
  #summarize uridylation data using dplyr
  if (facet_projects == FALSE) {
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }
  else {
    #if using facets - group by projects also
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }
  #print(test_trans)
  #leave only uridylation values (exclude uridylation == FALSE)
  test_trans <- test_trans %>% dplyr::filter(uridylated2==TRUE) 
  #store calculated values
  output$calculated_values = test_trans
  #select data for plot
  #test_trans %>% ungroup() %>% select(replicate,freq_urid,project_name,condition) 
  #create plot
  plot_out <- test_trans %>% ggplot(aes(x=condition)) + geom_bar(stat="identity",position="dodge",aes(y=mean_freq_urid))
  plot_out <- plot_out  + geom_errorbar(aes(ymin =  mean_freq_urid - sd_urid, ymax = mean_freq_urid + sd_urid),colour = "black", width = 0.1, position = position_dodge(0.9)) 
  plot_out <- plot_out + ggtitle(plot_title)
  plot_out <- plot_out + xlab("condition") + ylab("fraction of transcripts")
  if (include_jitter==TRUE) {
    plot_out <- plot_out + geom_jitter(aes(y=freq_urid))
  }
  
  #create facets 
  if (facet_projects==TRUE) {
    plot_out <- plot_out + facet_grid (project_name ~ .)
  }
  #perform dunn test (one to many, first condition is used as control)
  output$dunn_test <- list()
  #if using facets - perform dunn test on each project independently
  if (facet_projects==TRUE) {
    dunn_test_pval2<-c() #vector for storing dunn test pvalues 
    for (proj in levels(test_trans$project_name)) {
      #iterate tgrough projects
      test_proj <- test_trans %>% dplyr::filter(project_name==proj) #filter data for single project
      if(nrow(test_proj)>0) { 
        #if there are any data after filtering, calculate Dunn's test
        #print(proj)
        #calculate Dunn using kwManyOneDunnTest from PCMCRPlus package
        dunn_test <- kwManyOneDunnTest(test_proj$freq_urid,test_proj$condition,p.adjust.method = "BH")
       # print(dunn_test)
        #store Dunn test results in a list
        dunn_output<-dunn_test$p.value
        colnames(dunn_output)<-c("p.value (Dunn test)")
        output$dunn_test[[proj]] <- dunn_output
        dunn_test_pval <- as.data.frame(dunn_test$p.value) # get pvalues
        #dunn_test_pval2 <- c(0,dunn_test_pval[,1]) 
        dunn_test_pval2<-c(dunn_test_pval2,dunn_test_pval[,1]) #store pvalues
      }
    }
  } else #if no facets, calculate Dunn for all projects together
  {
    dunn_test <- kwManyOneDunnTest(test_trans$freq_urid,test_trans$condition,p.adjust.method = "BH")
  #  print(dunn_test)
    dunn_output<-dunn_test$p.value
    colnames(dunn_output)<-c("p.value (Dunn test)")
    output$dunn_test <- dunn_output

    dunn_test_pval <- as.data.frame(dunn_test$p.value)
    dunn_test_pval2 <- dunn_test_pval[,1] #store pvalues
  }
  
  #get plot data to build info about significance
  pg <- ggplot_build(plot_out)
  #get max position of error bars from plot data (pg$data[2]), group by PANEL (facets) and group (condition)
  plot_data<-as.data.frame(pg$data[2]) %>% dplyr::group_by(group,PANEL) %>% summarise(max_pos = max(ymax)) %>% ungroup() %>% mutate(max_pos2=max(max_pos)) %>% mutate(y=max_pos2+0.05 + max_pos2*(0.25*group)) %>% dplyr::filter(group!=1) %>% arrange(PANEL,group)
#  print(plot_data)
  #create_annotation
  if (facet_projects==TRUE) 
  {
    annotation <- test_trans %>% dplyr::filter(uridylated2==TRUE) %>% dplyr::group_by(project_name,condition) %>% summarise(replicate=mean(replicate)) %>% arrange(project_name,condition) %>% slice(-1)
  } else 
  {
    annotation <- test_trans %>% dplyr::filter(uridylated2==TRUE) %>% dplyr::group_by(condition) %>%summarise(replicate=mean(replicate)) %>% slice(-1)
  }
  #calculate ylim
  ylim_value = max(plot_data$y) +0.15 * max(plot_data$y)
  plot_out = plot_out + expand_limits(y=c(0,ylim_value))
  #get position of significance marks from plot data
  annotation$y=plot_data$y
  #store pvalues from Dunn's test
  annotation$pvalue=dunn_test_pval2
  #store start position for significance bars - always first condition - control
  annotation$start=test_trans$condition[1]
  #create significance asterisks:
  annotation$sig=''
  if(any(annotation$pvalue<=0.05)) {
    annotation[annotation$pvalue<=0.05,]$sig <- "*"
  }
  if(any(annotation$pvalue<=0.01)) {
    annotation[annotation$pvalue<=0.05,]$sig <- "**"
  }
  if(any(annotation$pvalue<=0.001)) {
    annotation[annotation$pvalue<=0.05,]$sig <- "***"
  }
  #exclude from annotation all conditions without signficance
  annotation <- annotation[annotation$sig!='',]
 # print(annotation)
  #if there atre any significant differences
  if (nrow(annotation)>0) {
    plot_out <- plot_out + geom_signif(data=annotation,aes(xmin=start, xmax=condition, annotations=sig, y_position=y),textsize = 10, vjust = 0.6,manual=TRUE)
  }
  
  #print(plot_out) #print plot
  output$plot <- plot_out #store plot in output
  return(output) #return output
  
}




#' Plot mapping positions
#'
#' @param dataset                - input dataset 
#' @param trans                  - transcript to analyze
#' @param exp_type2              - experiment type (OVR,KD,LEAP,NT)
#' @param conditions             - conditions to include in the analysis
#' @param mapping_position_min   - minimal mapping position to include on plot
#' @param mapping_position_max   - maximum mapping position to include on plot
#' @param facet_projects         - show facets  on projects
#' @param facet_replicates       - show facets on replicates 
#' @param project 
#'
#' @return                       - plot with mapping positions
#' @export
#'
#' @examples
plot_mapping_positions <- function(dataset,trans,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,project = NA, facet_projects = FALSE,facet_replicates = FALSE) {
  #filter input dataset to include only required transcript and experiment type
  test_trans <- dataset %>% dplyr::filter(transcript==trans)
  test_trans <- test_trans %>% dplyr::filter(exp_type==exp_type2)
  #create plot title
  plot_title <- paste(exp_type2,trans,sep=",")
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% dplyr::filter(mapping_position>mapping_position_min)
    plot_title <- paste(plot_title,"min map pos: ",mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% dplyr::filter(mapping_position<mapping_position_max)
    plot_title <- paste(plot_title,"max map pos: ",mapping_position_max)
  }
  #if project is specified - used for filtering
  if (length(project)>0) {
    test_trans <- test_trans %>% dplyr::filter(project_name %in% project)
  }
  #if conditions are specified - use for filtering
  if(length(conditions)>0 & !is.na(conditions)) {
    test_trans <- test_trans %>% dplyr::filter(condition %in% conditions)
    plot_title <- paste(plot_title,paste(conditions,collapse=", "),sep="; ")
  }
  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  
  #create plot
  plot_out <- ggplot(test_trans,aes(x=mapping_position,group=condition,colour=condition)) 
  plot_out <- plot_out + stat_bin(binwidth=1,aes(x=mapping_position,y=..ncount..),geom="line") 
  plot_out <- plot_out +  ggtitle(plot_title)
  plot_out <- plot_out + xlab("mapping position")
  plot_out <- plot_out + ylab("count")
  #create facets (if specified)
  if (facet_projects==TRUE) {
    plot_out <- plot_out + facet_grid (project_name ~ .)
  }
  #if not faceting by projects, facet by replicates (if specified)
  else if (facet_replicates==TRUE) {
    plot_out <- plot_out + facet_grid (replicate ~ .)
  }
  #print(plot_out)
  return(plot_out)
}



#' Plot tail length disribution
#'
#' @param dataset               - input dataset
#' @param trans                 - transcript to analyze
#' @param exp_type2             - experiment type (OVR,KD,LEAP,NT)
#' @param conditions            - conditions to include (vector)
#' @param mapping_position_min  - minimal mapping position to include
#' @param mapping_position_max  - maximal mapping position to include
#' @param facet_conditions      - facet by conditions
#' @param facet_replicates      - facet by replicates
#' @param max_tail_length       - max tail length to include
#' @param min_tail_length       - min tail length to include
#' @param facet_tail_types      - facet by tail types
#' @param AUtail                - specify A or U tail part only ('Atail' or 'Utail', other values ignored)
#'
#' @return                      - plot with tail lengths distribution
#' @export
#'
#' @examples
plot_tail_length_distribution <- function(dataset,trans,exp_type2,conditions = NA,mapping_position_min = NA,mapping_position_max = NA,facet_conditions = FALSE,facet_replicates = FALSE,max_tail_length = NA, min_tail_length = NA, facet_tail_types = FALSE,AUtail = NA) {
  test_trans <- dataset %>% filter(transcript==trans)
  test_trans <- test_trans %>% filter(exp_type==exp_type2)
  
  plot_title <- trans
  if (!is.na(AUtail)) {
    if (AUtail == 'Atail') {
      test_trans <- test_trans %>% mutate(tail_length = Atail_length)
      plot_title <- paste(plot_title,"A tail length only",sep=",")
    } else if (AUtail =='Utail') {
      test_trans <- test_trans %>% mutate(tail_length = Utail_length)
      plot_title <- paste(plot_title,"U tail length only",sep=",")
    }
  }
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% filter(mapping_position>mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% filter(mapping_position<mapping_position_max)
  }
  if (!is.na(min_tail_length)) {
    test_trans <- test_trans %>% filter(tail_length>min_tail_length)
  }
  if (!is.na(max_tail_length)) {
    test_trans <- test_trans %>% filter(tail_length<max_tail_length)
  }
  #if conditions are specified - use for filtering
  if(length(conditions)>0 & !is.na(conditions)) {
    test_trans <- test_trans %>% filter(condition %in% conditions)
  }

  test_trans$condition <- as.character(test_trans$condition)
  test_trans$condition <- as.factor(test_trans$condition)
  test_trans$tail_type <- as.character(test_trans$tail_type)
  test_trans$tail_type <- as.factor(test_trans$tail_type)
  

  
  plot_out <- ggplot(test_trans,aes(x=tail_length,group=tail_type,colour=tail_type)) 
  plot_out <- plot_out + stat_bin(binwidth=1,aes(x=tail_length,y=..count..,colour=tail_type),geom="line") 
  plot_out <- plot_out +  ggtitle(plot_title)
  #create facets 
  if (facet_conditions==TRUE) {
    plot_out <- plot_out + facet_grid (condition ~ .)
  }
  else if (facet_replicates==TRUE) {
    plot_out <- plot_out + facet_grid (replicate ~ .)
  }
  else if (facet_tail_types==TRUE) {
    plot_out <- plot_out + facet_grid (tail_type ~ .)
  }
 # print(plot_out)
  return(plot_out)
}


#analyze_uridylation()
plot_tail_length_distribution(all_data,"LEAP_AU","LEAP",conditions=c("19A","19A3U"),max_tail_length = 30, min_tail_length=10,facet_conditions = T,AUtail = 'Ataila')

plot_tail_length_distribution(all_data,"LEAP_AU","LEAP",conditions=c("19A","19A3U"),max_tail_length = 30, min_tail_length=10,facet_conditions = T,AUtail = 'Ataila')

