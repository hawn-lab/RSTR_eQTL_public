plot_eQTL <- function(to_plot, snp_anno, cis_anno=NULL,
                      type=c("cis","trans"), rename_rsID=FALSE,
                      FDR = FALSE){
  plot_ls <- list()
  
  for(i in 1:nrow(to_plot)){
    gene_id <- to_plot[i,"gene"] %>% unlist()
    snp_id <- to_plot[i,"snp"] %>% unlist()
    
    # Expression data
    data <- GENE %>% 
      rownames_to_column("symbol") %>% 
      filter(symbol==gene_id) %>% 
      pivot_longer(-symbol, values_to = "expression") %>% 
      #SNP data
      inner_join(SNP %>% rownames_to_column("snp") %>% 
                   filter(snp == snp_id) %>% 
                   pivot_longer(-snp, values_to ="geno"),
                 by = join_by(name)) %>% 
      #covariates
      inner_join(cov_n160 %>% t() %>% as.data.frame() %>% rownames_to_column("name"),
                 by = join_by(name)) %>% 
      mutate(M0_KCVSEX = factor(M0_KCVSEX)) %>%
      separate(name, into = c("ID","status"), sep="_", remove=FALSE) 
    
    # Get reference and alternative allele of the SNP
    ref_alt <- unlist(snp_anno[snp_anno$rsID == snp_id, c("ref", "alt")])
    
    #Error if no annotation found
    if(length(ref_alt)==0){ stop(paste(snp_id, "not found in snp_anno data.")) }
    
    # Prepare the genotype labels
    gt_states <- c(paste(ref_alt[1], ref_alt[1], sep="/"), 
                   paste(ref_alt[1], ref_alt[2], sep="/"), 
                   paste(ref_alt[2], ref_alt[2], sep="/"))
    gt_states <- factor(gt_states, levels=gt_states)
    
    # Assign the labels based on SNP genotype
    data$gt <- gt_states[data$geno + 1]
    
    # Filter out NA genotypes
    data <- data %>% filter(!geno %in% NA)
    
    # Perform linear regression for TB status
    ###Fit lines don't look correct, esp fo TB
    lm_out1 <- data %>% 
      filter(status == "MEDIA") %>% 
      lm(data=., expression ~ geno + M0_KCVSEX + M0_KCVAGE)
    new_data1 <- data.frame("intercept" = lm_out1$coefficients[1] - lm_out1$coefficients[2],
                            "slope" = lm_out1$coefficients[2]) %>%
      mutate(status = "MEDIA")
    
    # Perform linear regression for MEDIA status
    lm_out2 <- data %>% 
      filter(status == "TB") %>% 
      lm(data=., expression ~ geno + M0_KCVSEX + M0_KCVAGE)
    new_data2 <- data.frame("intercept" = lm_out2$coefficients[1] - lm_out2$coefficients[2],
                            "slope" = lm_out2$coefficients[2]) %>%
      mutate(status = "TB")
    
    data_lm <- full_join(data, bind_rows(new_data1, new_data2),
                         by = join_by(status)) %>% 
      mutate(status = recode_factor(status, "MEDIA"="Media","TB"="+Mtb"))
    
    # Create the plot
    #Rename with rsID if selected
    if(rename_rsID){
      if(is.null(to_plot$rsID)) { stop("Please provide rsID values for renaming.") }
      
      snp_rsid <- to_plot %>% 
        filter(gene==gene_id & snp==snp_id) %>% 
        pull(rsID)
      
      if(is.na(snp_rsid)){ snp_rsid <- snp_id }
    } else{
      snp_rsid <- snp_id
    }
    
    #Make titles
    if(type == "cis"){
      plot_title <- gene_id
      x_lab <- snp_rsid
    } else if(type == "trans"){
      #find cis gene for SNP
      gene_id2 <- cis_anno %>% filter(snp == snp_rsid) %>% pull(gene) %>% unique()
      plot_title <- gene_id
      x_lab <- paste(snp_rsid, "(", gene_id2, ")")
    }
    
    #Add FDR
    if(FDR & type=="cis") {
      data_lm_fdr <- to_plot %>% 
        select(snp, gene, Media, `+Mtb`) %>% 
        pivot_longer(Media:`+Mtb`, names_to = "status", values_to = "FDR") %>% 
        right_join(data_lm) %>% 
        mutate(status = recode(status, "Media"="Uninfected")) %>% 
        mutate(facet.lab = paste(status, FDR, sep="\n"))
      
      media_lvl <- data_lm_fdr$facet.lab[grepl("Uninfected",data_lm_fdr$facet.lab)] %>% 
        unique()
      data_lm_fdr <- data_lm_fdr %>% 
        mutate(facet.lab=fct_relevel(facet.lab, media_lvl, after=0))
    
      } else {
      data_lm_fdr <- data_lm %>% 
        mutate(facet.lab=factor(status, levels=c("Media","+Mtb"))) %>% 
        mutate(facet.lab=fct_recode(status, "Uninfected"="Media"))
    }
    
    plot_ls[[paste(gene_id,snp_id,sep="_")]] <- ggplot(
      data_lm_fdr, aes(gt, expression, color = as.factor(geno), fill = as.factor(geno)))  +
      #Raw data points
      geom_point(position=position_jitterdodge(), color="grey80", show.legend = FALSE) +
      #Stdev bar
      stat_summary(fun.data=mean_sdl, 
                   fun.args = list(mult=1), 
                   geom="errorbar", width=0.25) +
      #mean bar
      stat_summary(fun=mean, geom="errorbar",
                   aes(ymax=after_stat(y), ymin=after_stat(y)),
                   width=0.5) +
      labs(y=paste(gene_id,"Normalized log2 expression",sep="\n"),
           x=x_lab,
           color="Genotype",
           title=plot_title) + 
      facet_grid(~ facet.lab)+
      theme_classic() +  
      #Fit lines
      geom_abline(aes(slope = slope, intercept = intercept),
                  col = "black")  +
      theme(axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5)) +
      #Recode legend
      scale_color_manual(values = c("#117733","#88CCEE","#AA4499"),
                         labels=c("REF/REF","REF/ALT","ALT/ALT"))
  }
  return(plot_ls)
}
