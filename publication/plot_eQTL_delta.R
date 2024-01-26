plot_eQTL_delta <- function(to_plot, snp_anno, cis_anno=NULL,
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
      pivot_longer(-symbol) %>% 
      #Calc delta
      separate(name, into=c("ptID","name"), sep="_") %>% 
      pivot_wider() %>% 
      mutate(delta=TB-MEDIA) %>% 
      select(-MEDIA,-TB) %>% 
      #SNP data
      inner_join(SNP %>% rownames_to_column("snp") %>% 
                   select(snp, contains("MEDIA")) %>% 
                   rename_all(~gsub("_MEDIA","",.)) %>% 
                   filter(snp == snp_id) %>% 
                   pivot_longer(-snp, values_to ="geno", names_to="ptID"),
                 by = join_by(ptID)) %>% 
      #covariates
      inner_join(cov_n160 %>% select(contains("MEDIA")) %>% rename_all(~gsub("_MEDIA","",.))%>% t() %>% as.data.frame() %>% rownames_to_column("ptID"),
                 by = join_by(ptID)) %>% 
      mutate(M0_KCVSEX = factor(M0_KCVSEX))
    
    
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
    lm_delta <- data %>% 
      lm(data=., delta ~ geno + M0_KCVSEX + M0_KCVAGE)
    
    data_lm <-data %>% 
      mutate("intercept" = lm_delta$coefficients[1],
             "slope" = lm_delta$coefficients[2])
    
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
    
    #Add P
    if(FDR & type=="trans"){
      data_lm_fdr <- to_plot %>% 
        mutate(facet.lab=paste("P =", signif(p.value, digits=2))) %>% 
        right_join(data_lm, by = c("gene"="symbol","snp"="snp"))
    }else {
      data_lm_fdr <- data_lm %>% 
        mutate(facet.lab=status)
    }
    
    plot_ls[[paste(gene_id,snp_id,sep="_")]] <-
      ggplot(
      data_lm_fdr, aes(gt, delta, color = as.factor(geno), fill = as.factor(geno)))  +
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
      labs(y=paste(gene_id,"Change in normalized log2 expression","+Mtb - Uninfected",sep="\n"),
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
