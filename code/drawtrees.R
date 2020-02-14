#! /usr/bin/env Rscript

library(magrittr)

HHtre <- wilbur::load_tree("HH.tre", filetype = "newick")
meta <- readr::read_delim("influenza_na.dat", delim="\t", 
                          col_names = c("genbank", "host", "segnum", "subtype",
                                        "country", "date", "len", "strain", "x1", "x2",
                                        "complete"))
cmeta <- meta %>%
  subset(date < 1966) %>%
  subset(!grepl("Influenza B", strain)) %>%
  subset(!grepl("Influenza C", strain)) %>%
  subset(!grepl("unidentified influenza virus", strain)) %>%
  dplyr::mutate(
    x1 = NULL,
    x2 = NULL,
    len = NULL,
    segnum = NULL,
    complete = NULL,
    strain = gsub("Influenza A virus \\(","",strain), 
    strain = gsub("\\(.*","", strain) %>% gsub("\\)", "", .)
  )

row.names(cmeta) = paste(cmeta$genbank, cmeta$strain, cmeta$subtype, cmeta$host, cmeta$country, cmeta$date, sep="|") %>% 
  gsub(" ", "_", .) %>%
  gsub("-", "_", .)

row.names(cmeta) = cmeta$genbank

# tre, str -> node id[int] of mrca
mrcaOfStr <- function(tre, node_str = "1931"){
  nodeids <- grep(node_str, tre@data$name)
  shared_mrca <- tidytree::MRCA(tre, nodeids)
  return(shared_mrca)
}


# tre, str -> tre that is rooted based on string
rootBy <- function(tre, node_str = "1931"){
  old_ids <- grep(node_str, tre@data$name)
  if(is.null(old_ids)) {
    cat("No matches")
    return
  }
  if(length(old_ids)>1){
    new_root <- tidytree::MRCA(tre, old_ids)
  }else{
    new_root=tidytree::parent(tre, old_ids)
  }
  outtree <- tidytree::as.treedata(ape::root.phylo(ape::as.phylo(tre), 
                                        node=new_root))
  all_names <- c(outtree@phylo$tip.label, outtree@phylo$node.label)
  outtree@data <- tibble::tibble(node = 1:length(all_names), 
                                 name = all_names)
  outtree
  }

cHHtre <- HHtre %>% 
  #rootBy(., node_str="1931") %>%
  rootBy(., node_str="1918") %>%
  {
    temp_tbl_df <- .@data
    temp_tbl_df %>%
      dplyr::mutate(
        
      )
    .@data$host = dplyr::case_when(grepl("wine", .@data$name) ~ "swine",
                                 grepl("uman", .@data$name) ~ "human",
                                 grepl("vian", .@data$name) ~ "avian",
                                 grepl("quine", .@data$name) ~ "equine",
                                 grepl("ouse", .@data$name) ~ "mouse")
  .
  }

addDataBy <- function(tre, tre_idcol, mdata){
  for (col in names(mdata)){
    cat("added @data$")
    cat(col)
    cat("\n")
    # tre@data[[col]] = unname(unlist(mdata[tre@data[[tre_idcol]] %>% 
    # gsub("\\|\\|", "|NA|", .), col]))
    tre@data[[col]]=unname(unlist(mdata[tre@data[[tre_idcol]], col]))
  }
  tre
}

# tree function
plot_tree <- function(tree){
  tree %>% ggtree::ggplot(.) +
    ggtree::geom_tree() + 
    ggtree::theme_tree2() 
  #    ggtree::geom_tiplab(size=1.5) 
}

cHHtre@data$id = cHHtre@data$name %>% gsub("\\|.+", "", .)
aa <- addDataBy(cHHtre, "id", cmeta)

# == color internal branches
aa@data$color="swine"
aa@data$color[
  aa@data %>% subset(host=="Human" & subtype=="H1N1") %>%
    {.$node %>% tidytree::MRCA(aa, .) %>% tidytree::offspring(aa, .)}
] = "human"

aa@data$color[
  aa@data %>% subset(host=="Swine") %>%
    {.$node %>% tidytree::MRCA(aa, .) %>% tidytree::offspring(aa, .)}
  ] = "swine"

aa@data$color[
  aa@data %>% subset(host=="Mouse") %>%
    {.$node}
  ] = "mouse"

aa@data$color[
  aa@data %>% subset(host=="Avian") %>%
    {.$node %>% tidytree::MRCA(aa, .) %>% tidytree::offspring(aa, .)}
  ] = "avian"

aa@data$color[
  aa@data %>% subset(host=="Human" & subtype=="H2N2") %>%
    {.$node %>% tidytree::MRCA(aa, .) %>% tidytree::offspring(aa, .)}
  ] = "human"

aa@data$color[
  aa@data %>% subset(host=="Equine" & subtype=="H3N8") %>%
    {.$node %>% tidytree::MRCA(aa, .) %>% tidytree::offspring(aa, .)}
  ] = "equine"

aa@data$color[
  aa@data %>% subset(host=="Equine" & subtype=="H7N7") %>%
    {.$node %>% tidytree::MRCA(aa, .) %>% tidytree::offspring(aa, .)}
  ] = "equine"

# clade labels
h2clade <- aa@data %>% subset(host=="Human" & subtype=="H2N2") %>% 
  {.$node %>% tidytree::MRCA(aa, .)}
h1clade <- aa@data %>% subset(host=="Human" & subtype=="H1N1" & date!="1918") %>% 
  {.$node %>% tidytree::MRCA(aa, .)}

hostcolors = c("avian"="#B2182B", "equine"="#EF8A62",
               "human"="#67A9CF", "swine"="#AF8Dc3",
               "mouse"="#555555")

(pp<- aa %>% plot_tree(.) + 
    ggplot2::aes(color=color)+
    ggtree::geom_tiplab(ggplot2::aes(label=paste(subtype, strain, sep="_")), 
                        size=1.5)+
    ggplot2::expand_limits(x=2.3) +
    ggplot2::theme(legend.position = "left") +
    ggplot2::scale_color_manual(values=hostcolors)+
    ggtree::geom_cladelabel(node=h2clade, label="  Human H2N2 (post-1957)", 
                            offset=0.31, color="#555555", alpha=1)+
    ggtree::geom_cladelabel(node=h1clade, label="  Human H1N1 (pre-1957)",
                            offset=0.37, color="#555555", alpha=1)+
    ggplot2::labs(title="Hemagglutinin (HA)")
) 
ggplot2::ggsave("HH2.pdf", plot=pp, dpi=100, width=8, height=8)
ggplot2::ggsave("HH.png", plot=pp, dpi=300, width=8, height=8)
