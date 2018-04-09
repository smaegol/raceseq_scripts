library(dplyr)
library(ggplot2)
library(PMCMRplus)



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
#'
#' @return                      - list containing calculated uridylation data (calculated output), Dunn test results and lot 
#' @export
#'
#' @examples
#' analyze_uridylation(all_data,"ACTB","OVR",project=c("Tailseq_4","Tailseq_5"),facet_projects==TRUE,conditions=c("CNTRL","TUT4WT","TUT7WT","MOV10"))
analyze_uridylation <- function(dataset,transcript2,exp_type2,mapping_position_min=NA,mapping_position_max=NA,project=NA,facet_projects=FALSE,conditions=NA) {
  
  output = list() #list for storing output
  print(mapping_position_min)
  print(mapping_position_max)
  print(transcript2)
  print(exp_type2)
  #get filtered data from input dataset
  #first, filter by transcript and experiment type (OVR,KD,LEAP)
  test_trans <- dataset %>% filter(transcript==transcript2,exp_type==exp_type2) 
  #if min and max mapping positions are specified - used for filtering
  if (!is.na(mapping_position_min)) {
    test_trans <- test_trans %>% filter(mapping_position>mapping_position_min)
  }
  if (!is.na(mapping_position_max)) {
    test_trans <- test_trans %>% filter(mapping_position<mapping_position_max)
  }
  #if project is specified - used for filtering
  if(!is.na(project)) {
    if (length(project==1)) {
      test_trans <- test_trans %>% filter(project_name==project)
    } else {
      test_trans <- test_trans %>% filter(project_name %in% project)
    }
  }
  #if conditions are specified - use for filtering
  if(length(conditions)>0 & !is.na(conditions)) {
    test_trans <- test_trans %>% filter(condition %in% conditions)
  }
  print(test_trans)
  #group data by co
  test_trans <- test_trans %>% group_by(condition,replicate,project_name,uridylated2)
  print(test_trans)
  #summarize uridylation data using dplyr
  if (facet_projects == FALSE) {
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }
  else {
    #if using facets - group by projects also
    test_trans <- test_trans %>% dplyr::summarize(n_urid=n()) %>% ungroup() %>% dplyr::group_by(condition,replicate,project_name) %>% dplyr::mutate(freq_urid = n_urid/sum(n_urid)) %>% dplyr::group_by(condition,uridylated2,project_name) %>% dplyr::mutate(mean_freq_urid = mean(freq_urid), sd_urid = sd(freq_urid))
  }
  print(test_trans)
  #leave only uridylation values (exclude uridylation == FALSE)
  test_trans <- test_trans %>% filter(uridylated2==TRUE) 
  #store calculated values
  output$calculated_values = test_trans
  #select data for plot
  #test_trans %>% ungroup() %>% select(replicate,freq_urid,project_name,condition) 
  #create plot
  plot_out <- test_trans %>% ggplot(aes(x=condition)) + geom_bar(stat="identity",position="dodge",aes(y=mean_freq_urid)) + geom_errorbar(aes(ymin =  mean_freq_urid - sd_urid, ymax = mean_freq_urid + sd_urid),colour = "black", width = 0.1, position = position_dodge(0.9)) + geom_jitter(aes(y=freq_urid)) + ggtitle(paste(exp_type2,transcript2,"uridylation"))
  
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
      test_proj <- test_trans %>% filter(project_name==proj) #filter data for single project
      if(nrow(test_proj)>0) { 
        #if there are any data after filtering, calculate Dunn's test
        #print(proj)
        #calculate Dunn using kwManyOneDunnTest from PCMCRPlus package
        dunn_test <- kwManyOneDunnTest(test_proj$freq_urid,test_proj$condition,p.adjust.method = "BH")
        print(dunn_test)
        #store Dunn test results in a list
        output$dunn_test[[proj]] <- dunn_test
        dunn_test_pval <- as.data.frame(dunn_test$p.value) # get pvalues
        #dunn_test_pval2 <- c(0,dunn_test_pval[,1]) 
        dunn_test_pval2<-c(dunn_test_pval2,dunn_test_pval[,1]) #store pvalues
      }
    }
  } else #if no facets, calculate Dunn for all projects together
  {
    dunn_test <- kwManyOneDunnTest(test_trans$freq_urid,test_trans$condition,p.adjust.method = "BH")
    print(dunn_test)
    output$dunn_test <- dunn_test
    dunn_test_pval <- as.data.frame(dunn_test$p.value)
    dunn_test_pval2 <- dunn_test_pval[,1] #store pvalues
  }
  
  #get plot data to build info about significance
  pg <- ggplot_build(plot_out)
  #get max position of error bars from plot data (pg$data[2]), group by PANEL (facets) and group (condition)
  plot_data<-as.data.frame(pg$data[2]) %>% group_by(group,PANEL) %>% summarise(max_pos = max(ymax)) %>% ungroup() %>% mutate(max_pos2=max(max_pos)) %>% mutate(y=max_pos2+0.05 + max_pos2*(0.05*group)) %>% filter(group!=1) %>% arrange(PANEL,group)
  print(plot_data)
  #create_annotation
  if (facet_projects==TRUE) 
  {
    annotation <- test_trans %>% filter(uridylated2==TRUE) %>% group_by(project_name,condition) %>% summarise(replicate=mean(replicate)) %>% arrange(project_name,condition) %>% slice(-1)
  } else 
  {
    annotation <- test_trans %>% filter(uridylated2==TRUE) %>% group_by(condition) %>%summarise(replicate=mean(replicate)) %>% slice(-1)
  }
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
  print(annotation)
  #if there atre any significant differences
  if (nrow(annotation)>0) {
    plot_out <- plot_out + geom_signif(data=annotation,aes(xmin=start, xmax=condition, annotations=sig, y_position=y),textsize = 10, vjust = 0.3,manual=TRUE)
  }
  
  print(plot_out) #print plot
  output$plot <- plot_out #store plot in output
  return(output) #return output
  
}