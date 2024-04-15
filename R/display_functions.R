#' @description
#' Summary function
#' @examples
#' tar_load(mcmc_output)
#' jags_sample_to_summary_tibble(mcmc_output)|>View()
jags_sample_to_summary_tibble <- function(jags_sample) {
  jags_sample |>
    summary() |>
    (function(ss) {
      tibble(parname = rownames(ss)) |>
        bind_cols(as_tibble(ss))
    })()
}


#' @description
#' color palette
color_palette <- c(RColorBrewer::brewer.pal(n = 12, name = "Paired"),"purple","orange","blue","green","pink","brown","turquoise","grey")

#' @description
#' Afficher l'abondance relative des variants
#' Ce serait bien de pouvoir extraire les noms d'échantillon directement dans la fonction plutôt que de devoir 
#' le passer en argument mais je ne suis pas du tout a l'aise avec ces objets.
#' Idem ce serait top d'avoir les noms de variant directement dans le tau_vga, et qu'ils puissent être extraits 
#' du pi_summary au lieu de les donner en argument.
#' @examples
#' results <- head(jags_sample_to_summary_tibble(mcmc_output),-2)
#' metadata=metadata_file |> read.csv()
#' metafields <- c(colnames(metadata),"cond")
#' metadata <- cbind(metadata,paste0(metadata$Rearing,metadata$Gen))
#' colnames(metadata) <- metafields
#' resultscond <- add_metadata(results, metadata, samples, length(variants),"cond","Age")

plot_pi_summary <- function(pi_summary,samples,variants,alpha,min_abundance,labelcond2) {
  V <- get_abundant_variants(pi_summary,min_abundance)
  vnames <- names(variants)[V]
  correct_abundances <- filter_abundances(pi_summary,V)
  
  correct_abundances |>
    dplyr::mutate(variant_index = parname |>
             lapply(FUN = function(nm) {
               nm |>
                 str_split_i(pattern = "\\[", i = 2) |>
                 str_split_i(pattern = ",", i = 1)
             }) |>
             unlist() |>
             as.numeric()) |>
    dplyr::mutate(sample = parname |>
             lapply(FUN = function(nm) {
               nm |>
                 str_split_i(pattern = ",", i = 2) |>
                 gsub("\\]", "", .)
             }) |>
             unlist() |>
             as.numeric()) |>
    ggplot(aes(x = sample)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(~cond1,scales="free_x",space="free_x") +
    geom_hline(yintercept = 0.05, linetype="dashed", colour="red", linewidth = 0.5) + 
    geom_segment(aes(xend = sample, y = Lower95, yend = Upper95, colour = as.factor(variant_index))) +
    geom_point(aes(y = Mean, shape = factor(cond2), colour = as.factor(variant_index)), size=4) +
    scale_shape_manual(values=c(16:18,15,0:6), name=labelcond2) +
    scale_color_manual(values=color_palette[1:length(vnames)], labels=vnames, name="variants") +
    scale_x_continuous(breaks = c(1:length(samples)),labels=samples) +
    ggtitle(paste0("Variant abundance estimation (alpha=",alpha,")"))
}

#' @description
#' get variants with relative abundance >= P in at least 1 sample
get_abundant_variants <- function(results,P){
  R=results[which(results$Mean >= P),]
  R$parname|>str_split_i(pattern = "\\[", i = 2) |>str_split_i(pattern = ",", i = 1)|>unlist() |>as.numeric()|>unique()|>sort()
} 
#' @description
#' Extract lines from R corresponding to variants with abundance >= P
filter_abundances <- function(results,V){
  Allvar=results$parname|>
    str_split_i(pattern = "\\[", i = 2) |>
    str_split_i(pattern = ",", i = 1)|>
    unlist()|>as.numeric()
  results[which(Allvar %in% V),]
}
#' @description
#' Add two metadata fields to the result
add_metadata <- function(results, metadata, samples, nbvar,m1,m2){
  allS=samples[results$parname|>str_split_i(pattern = ",", i = 2) |>
                 str_split_i(pattern = "]", i = 1)|>unique()|>as.numeric()]
  cond1=rep(metadata[which(metadata$ID %in% allS),"cond"],times=1,each=length(variants))
  cond2=rep(metadata[which(metadata$ID %in% allS),"Age"],times=1,each=length(variants))
  cbind(results,cond1,cond2)
}



if(FALSE){
  function(variants,samples,alpha_best,mcmc_output,mcmc_output_bestalpha,
           metadata_file)
    
results <- head(jags_sample_to_summary_tibble(mcmc_output),-2)
resultba <- head(jags_sample_to_summary_tibble(mcmc_output_bestalpha),-2)


metadata=metadata_file |> read.csv()

metafields <- c(colnames(metadata),"cond")
metadata <- cbind(metadata,paste0(metadata$Rearing,metadata$Gen))
colnames(metadata) <- metafields

resultscond <- add_metadata(results, metadata, samples, length(variants),"cond","Age")
resultbacond <- add_metadata(resultba, metadata, samples, length(variants),"cond","Age")

list(plot1=plot_pi_summary(resultscond,samples,variants,0.1,0,"age (days)"),
     plot2=plot_pi_summary(resultscond,samples,variants,0.1,0.05,"age (days)"),
     plot3=plot_pi_summary(resultbacond,samples,variants,alpha_best,0.05,"age (days)"))}