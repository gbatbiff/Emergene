rm(list = ls())

packages <- c(
  "ggplot2",
  "phangorn",
  "stringr",
  "data.table",
  "dplyr",
  "R.utils",
  "openxlsx",
  "ape",
  "treeio",
  "ggtree",
  "phytools",
  "TreeTools"
)

invisible(lapply(packages, function(pkg) {
  suppressMessages(suppressWarnings(library(pkg, character.only = TRUE)))
}))
#setwd("/home/ghepard/AMR_MODEL/")

args=commandArgs(trailingOnly = TRUE)

help_message <- "
Usage: Rscript scriptname.R [options]

Options:
  -tree <tree_file>          Path to the phylogenetic tree file
  -amr <amr_table_file>      Path to the AMR table file (for amrfinderplus)
  -pangenome <pangenome_file> Path to the gene presence-absence matrix(Panaroo/Roary like)
  --help                     Show this help message and exit
"
# Show help message if '--help' is passed or arguments are missing
if ("--help" %in% args || length(args) == 0) {
  cat(help_message)
  quit(save = "no", status = 0)
}

# Initialize variables for the arguments
amr_table_file <- NULL
pangenome_file <- NULL
tree_file <- NULL
input_type <- NULL

# select for '-amr', '-pangenome', and '-tree'
for (i in 1:length(args)) {
  if (args[i] == "-amr" && i + 1 <= length(args)) {
    amr_table_file <- args[i + 1]
    input_type <- "amrfinder"
  }
  if (args[i] == "-pangenome" && i + 1 <= length(args)) {
    pangenome_file <- args[i + 1]
    input_type <- "panaroo"
  }
  if (args[i] == "-tree" && i + 1 <= length(args)) {
    tree_file <- args[i + 1]
  }
}

# Check if tree file is provided 
if (is.null(tree_file)) {
  cat("Error: Missing '-tree' argument for the phylogenetic tree file.\n")
  cat(help_message)
  quit(save = "no", status = 1)
}

# Ensure that only one input type is provided
if (is.null(input_type)) {
  cat("Error: Either '-amr' or '-pangenome' must be specified.\n")
  cat(help_message)
  quit(save = "no", status = 1)
}

if (input_type == "amrfinder" && is.null(amr_table_file)) {
  cat("Error: Missing '-amr' argument for AMR table file.\n")
  cat(help_message)
  quit(save = "no", status = 1)
}

if (input_type == "panaroo" && is.null(pangenome_file)) {
  cat("Error: Missing '-pangenome' argument for pangenome file.\n")
  cat(help_message)
  quit(save = "no", status = 1)
}

cat("Loading phylogenetic tree...\n")
tree <- read.tree(tree_file)
cat("Tree Loaded!\n")

if (input_type == "amrfinder") {
  cat("Processing AMR Table (amrfinder)...\n")
  amr <- read.delim(amr_table_file, sep = "\t")
  cat("AMR Table Loaded!\n")
  #print(head(amr))  # Just to check if the file is loaded correctly
}

if (input_type == "panaroo") {
  cat("Processing Pangenome (panaroo)...\n")
  pangenome <- read.delim(pangenome_file, sep = "\t")
  cat("Pangenome File Loaded!\n")
  #print(head(pangenome))  # Just to check if the file is loaded correctly
}


timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")  # Get the current timestamp
output_dir <- file.path("results", timestamp)  # Create folder name based on timestamp

# Create the timestamped folder
cat("Creating results directory: ", output_dir, "\n")
dir.create(output_dir, recursive = TRUE)


#phylotree<-args[1]

tr <- read.tree(tree_file)


#amrtable<-args[2]

amr<-read.delim(amr_table_file, sep = "\t")


amr<-amr[which(amr$X..Coverage.of.reference.sequence==100 & amr$X..Identity.to.reference.sequence==100),]


meta <- amr[amr$Name %in% tree$tip.label,]


rates<-data.frame()


nmost<-as.data.frame(table(meta$Gene.symbol))


### take minimum occurrence of strains harbouring the resistance

nmost<-nmost[nmost$Freq>=10,"Var1"]


for (amr in unique(nmost)){
  
  
  
  amr<-"blaTEM-1"
  
  df<-meta[which(meta$Gene.symbol == amr),]
  
  res_df<-table(df$Name, df$Gene.symbol)
  
  res<-data.frame(row.names = tr$tip.label)
  res$amr<-ifelse(rownames(res) %in% names(res_df[,1]), "R", "S")
  xx<-setNames(as.factor(res$amr), rownames(res))
  
  if (length(levels(xx))==1) {
    
    next
  }  
  
  
  
  if (table(xx)[[1]]<table(xx)[[2]]) {
    
    
   tree<-di2multi(tree)
    
    
   cat(sprintf("Estimating rates of gene %s...",amr))
   
    
    simmap_Q <- make.simmap(tree, xx, model = "ARD", tips = TRUE, pi ="estimated", nsim = 1)
    
    Q<-simmap_Q$Q
    
    simmap_mod <- make.simmap(tree, xx, model = "ARD", tips = TRUE, pi="estimated", nsim = 3, Q=Q)
    
    
    simmap<-summary(simmap_mod)
    
    
    ancstats_all <- as.data.frame(simmap$ace)
    
    ancstats_nodes<-head(ancstats_all, tree$Nnode)
    
    # Loop over all internal nodes (internal nodes are numbered from len(tips)+1 to len(tips) + len(nodes))
    
    ancstats_nodes$node<-length(tree$tip.label):length(tree$edge.length)+1
    
    ancstats_tips<-ancstats_all[-which(rownames(ancstats_all) %in% ancstats_nodes$node),]
    
    ancstats_tips$node<-1:(length(tree$tip.label))
    
    ancstats<-rbind.data.frame(ancstats_tips, ancstats_nodes)
    
    
    cols<-setNames(c("red", "blue"), c("S","R")) ### SOSTITUIRE CON R/S
    
    #tree<-multi2di(tree)
    tree<-di2multi(tree)
    
    
    
    ff<-data.frame(label = names(xx), stat = as.factor(xx) )
    ff<-unique(ff)
    tree2 <- full_join(tree, ff, by="label")
    
    
    node_dt_1 <- data.table(name = rownames(ancstats),
                            pheno = ifelse(ancstats[, "R"] >= ancstats[, "S"], "R", "S"))
    
    
    node_dt_1[, id := 1:.N] 
    node_dt_1[, child1 := Children(tree, id)[1], "id"]  ### take first children of bifurcation of parent node
    node_dt_1[, child2 := Children(tree, id)[2], "id"] ### take second children of bifurcating (no polytomies)
    
    temp <- node_dt_1[, .(id, cpheno1 = pheno)]
    ndt <- merge.data.table(node_dt_1, temp, by.x = "child1", by.y = "id", all.x = TRUE)
    
    temp <- node_dt_1[, .(id, cpheno2 = pheno)]
    ndt <- merge.data.table(ndt, temp, by.x = "child2", by.y = "id", all.x = TRUE)
    
    
    poly_parents <- ndt[((cpheno1 == "R" | cpheno2=="R") & pheno == "S"), id]  
    
    poly_parents_only_R <- ndt[(cpheno1 == cpheno2 & pheno == "S" & cpheno1=="R"), id]  ### polyphiletic group where parents were S (pheno==S is the parent) and one descendants became R
    
    
    
    
    if (length(poly_parents)>0 | length(poly_parents_only_R)>0){
      
      poly_parents<-append(poly_parents,poly_parents_only_R)
      
      
      poly_parents <- poly_parents[order(as.data.frame(nodeHeights(tree))[poly_parents,2], decreasing = TRUE)] ### ordered nodes
      
      
      pheno_nodes <- ndt[pheno == "R", id] ### cases where nodes/leaves are R
      
      pheno_nodes_poly_desc <- pheno_nodes[!pheno_nodes %in% poly_parents]  ### R descendant (nodes/leaves) excluding cases when polyphily occurs (when S became R and parent = S)
      
      
      ########### FILTER OUT THE DESCENDANTS OF POLY PARENTS WHICH CHILDS ARE SINGLETON OR ONLY SENSITIVE
      
      poly_parents_with_R_descendant <- c() 
      
      get_poly_parents_with_R_descendant <- function(tree, poly_parents, pheno_nodes) {
        
        
        for (n in seq_along(poly_parents)) {
          
          
          # Find children of the current poly parent
          childs_poly_gen <- Children(tree, poly_parents[n])
          childs_poly_gen <- childs_poly_gen[childs_poly_gen %in% pheno_nodes]
          
          for (c in childs_poly_gen) {
            # Check for R children
            check_R_child <- Children(tree, c)
            check_R_child <- check_R_child[check_R_child %in% pheno_nodes]
            
            if (length(check_R_child) == 2) {
              # If there are R children, store the parent of the MRCA of those children
              poly_parents_with_R_descendant[n] <- parent(tree, getMRCA(tree, check_R_child))
              
            } 
            
            
            if (length(check_R_child) == 1) {
              
              ### if the children is singleton R and its sibling (node or terminal tips) is also R then add poly_parent to the pp_with_R_desc list
              
              check_R_siblings<-Siblings(tree, check_R_child)
              check_R_siblings <- check_R_siblings[check_R_siblings %in% pheno_nodes]
              
              
              
              if (length(check_R_siblings)>0) {
                
                poly_parents_with_R_descendant[n] <- parent(tree, getMRCA(tree, c(check_R_child, Siblings(tree, check_R_child))))
                
                
                ### if the sibling is S but c R with R desc then add poly_parent to the pp_with_R_desc list
                
              } else if (c %in% seq_len(length(tree$tip.label)) | c %in% pheno_nodes)  {
        
                # If c is a tip, store the parent of that child
                poly_parents_with_R_descendant[n] <- parent(tree, c)
                
                
              } 
              
              
            }
          }
        }
        
        return(poly_parents_with_R_descendant)
        
      }
      
      
      
      poly_parents_with_R_descendant<-get_poly_parents_with_R_descendant(tree, poly_parents, pheno_nodes)
      # 
      # 
      poly_parents<-poly_parents_with_R_descendant[!is.na(poly_parents_with_R_descendant)]
      
      
      if (length(poly_parents)>0) {
        
        
        ll <- list() ### group S nodes becoming R
        k <- 1
        i <- NA
        
        sorted_poly_parents<-sort(poly_parents)
        
        
        
        hier_nodes<-c() ### multiple poly nodes across the same clade
        
        find_clusters<-function(tree, sorted_poly_parents, pheno_nodes_poly_desc, hier_nodes) {
          
          
          for(i in 1:length(sorted_poly_parents)){ ### itera nei casi in cui da S a R
            
            gen1_parent<-Children(tree, sorted_poly_parents[i])
            
            g1<-getDescendants(tree, gen1_parent[1])
            g2<-getDescendants(tree, gen1_parent[2])
            
            
            ### preliminary step to evaluate if any poly parent with other poly nodes across descendants with singletons 
            
            if (length(g2)==1) {
              
              ph2 <- ifelse(g2 %in% pheno_nodes_poly_desc, "R", "S") ### se i figli rientra nei discendenti dei nodi dove S to R
              
              ll[[k]] <- data.table(cluster = i, id = g2, pheno = ph2, poly_parent = sorted_poly_parents[i]) ### create a df with the features of the polyphylies node
              k <- k + 1
              
            }
            
            
            if (length(g1)==1) {
              
              ph1 <- ifelse(g1 %in% pheno_nodes_poly_desc, "R", "S") ### se i figli rientra nei discendenti dei nodi dove S to R
              
              ll[[k]] <- data.table(cluster = i, id = g1, pheno = ph1, poly_parent = sorted_poly_parents[i]) ### create a df with the features of the polyphylies node
              k <- k + 1
              
            }
            
            
            g1_pp<-g1[g1 %in% sorted_poly_parents]
            
            if (length(g1_pp)>0) {
              
              hier_nodes<-append(hier_nodes, g1_pp)
              
            }
            
            
            
            g2_pp<-g2[g2 %in% sorted_poly_parents]
            
            if (length(g2_pp)>0) {
              
              hier_nodes<-append(hier_nodes, g2_pp)
              
            }
            
            
            ph1 <- ifelse(g1 %in% pheno_nodes_poly_desc, "R", "S") ### se i figli rientra nei discendenti dei nodi dove S to R
            ph2 <- ifelse(g2 %in% pheno_nodes_poly_desc, "R", "S") ### se i figli rientra nei discendenti dei nodi dove S to R
            
            
            if (length(g2)>1) {
              
              ll[[k]] <- data.table(cluster = i, id = g2, pheno = ph2, poly_parent = sorted_poly_parents[i]) ### create a df with the features of the polyphylies node
              k <- k + 1
              
            } 
            
            if (length(g2)==1) {
              
              ll[[k]] <- data.table(cluster = i, id = g2, pheno = ph2, poly_parent = sorted_poly_parents[i]) ### create a df with the features of the polyphylies node
              k <- k + 1
              
            }
            
            if (length(g1)>1) {
              
              ll[[k]] <- data.table(cluster = i, id = g1, pheno = ph1, poly_parent = sorted_poly_parents[i]) ### create a df with the features of the polyphylies node
              k <- k + 1
              
            }
            
            if (length(g1)==1) {
              
              ll[[k]] <- data.table(cluster = i, id = g1, pheno = ph1, poly_parent = sorted_poly_parents[i]) ### create a df with the features of the polyphylies node
              k <- k + 1
              
            }
            
            
          } 
          
          return(ll)
          
        }
        
        
        
        ll<-find_clusters(tree, sorted_poly_parents, pheno_nodes_poly_desc, hier_nodes)
        
        
        ll <- rbindlist(ll) ### bind the cluster df
        ll <- ll[pheno == "R"] ### select only cases where the leave/node is R
        
        ll[, cluster := .GRP, by = cluster]  ### performs fast, ordered and unordered, groupings of vectors and data frames (or lists of vectors)
        
        ll<-as.data.frame(ll)
        
        ll<-unique(ll)
        
        
        
        set_cluster_na <- function(df) {
          # Arrange data by poly_parent in decreasing order
          
          df <- df %>%
            arrange(desc(poly_parent))
          
          # Create a temporary variable to hold the previous ids
          prev_ids <- NULL
          
          for (i in 1:nrow(df)) {
            current_id <- df$id[i]
            
            # If the current_id is in prev_ids, set cluster to NA
            if (current_id %in% prev_ids) {
              df$cluster[i] <- NA
            }
            
            # Update prev_ids to add the current_id
            prev_ids <- c(prev_ids, df$id[i])
          }
          
          
          return(df)
        }
        
        
        ll <- set_cluster_na(ll)
        
        ll<-ll[which(!is.na(ll$cluster)),]
        
        ll$node<-ifelse(ll$id>length(1:length(tree$tip.label)),"internal","terminal")
        
        
        
        terminal_single<-ll[which(ll$node=="terminal"),]
        single<-as.data.frame(table(terminal_single$cluster))
        single_clust<-single[which(single$Freq==1),"Var1"]
        
        
        single_clust_id<-terminal_single[which(terminal_single$cluster %in% single_clust),"id"]
        
        polyphyly<-length(unique(ll$cluster))
        
        
        poly_nodes<-pheno_nodes[!pheno_nodes %in% 1:length(tree$tip.label)]
        
        
        ############# FILTER OUT THE POLY NODES WHICH CHILDS ARE SENSITIVE
        
        poly_nodes_with_R_child<-c()
        
        if (length(poly_nodes)>0){
          
          for (n in 1:length(poly_nodes)){
            
            poly_nodes_desc<-getDescendants(tree, poly_nodes[n])
            
            poly_nodes_desc_R<-poly_nodes_desc[which(poly_nodes_desc %in% ll$id)]
            
            if (length(poly_nodes_desc_R)>1){
              
              poly_nodes_with_R_child[n]<-poly_nodes[n]
              
              
            }  
          }
          
          poly_nodes<-poly_nodes_with_R_child[!is.na(poly_nodes_with_R_child)]
          
        }  
        
        
        dist_parent_child<-data.frame(label = tree$edge, stat = tree$edge.length)
        
        tree_label<-1:length(tree$tip.label)
        
        res_nodes_leaves<-pheno_nodes[pheno_nodes %in% tree_label] ### number of R strains as numerator for the emergence rate ratio
        
        
        
        is_integer_zero <- function(obj) {
          is.integer(obj) && length(obj) == 0
        }
        
        
        
        poly_parents_with_terminal_leaves_dist<-c()
        poly_parents_with_terminal_leaves<-c()
        
        singletons<-c()
        sm<-c()
        
        for(k in 1:length(poly_parents)){
          #  k<-8
          child<-Children(tree, poly_parents[k]) ### get child of polyphyly where S became R
          
          ### check if the child+1 of the poly parent has singletons
          
          child_nextgen1<-Children(tree, child)[[1]] 
          child_nextgen2<-Children(tree, child)[[2]] 
          
          singleton1<-child_nextgen1[child_nextgen1 %in% res_nodes_leaves] 
          singleton2<-child_nextgen2[child_nextgen2 %in% res_nodes_leaves] 
          
          
          if (length(singleton1) == 1) {
            
            ### check if siblings have at least one terminal tip childrens which are R, in this case fake singleton, not to be included 
            
            get_siblings1<-Children(tree, Siblings(tree, singleton1))
            check_singleton_siblings1<-get_siblings1[get_siblings1 %in% ll$id] ### direct singleton descendant of poly parent node
            
            if (length(check_singleton_siblings1)>0) {
              
              singletons[k]<-NA
              
            } else if (parent(tree, singleton1) %in% poly_parents) {
              
              ### se il padre del nodo singleton non è un nodo poly allora non va in poly_parent_with_terminal leaves ma è cmq un singleton (inserire distanza?)
              singletons[k]<-singleton1
              
              #poly_parents_with_terminal_leaves_dist[k]<-dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==singleton2),"stat"]
              
              
            } else if (parent(tree, singleton1) %in% poly_parents[k]) {
              
              singletons[k]<-singleton1
              
              poly_parents_with_terminal_leaves_dist[k]<-dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==singleton1),"stat"]
              
              
              
            }
          }
          
          
          
          if (length(singleton2) == 1) {
            
            ### check if siblings have at least one terminal tip childrens which are R, in this case fake singleton, not to be included 
            
            get_siblings2<-Children(tree, Siblings(tree, singleton2))
            check_singleton_siblings2<-get_siblings2[get_siblings2 %in% ll$id] ### direct singleton descendant of poly parent node
            
            if (length(check_singleton_siblings2)>0) {
              
              singletons[k]<-NA
              
            } else if (parent(tree, singleton2) %in% poly_parents) {
              
              ### if parent of singleton not poly then not in poly_parent_with_terminal leaves but it is singleton 
              singletons[k]<-singleton2
              
              #poly_parents_with_terminal_leaves_dist[k]<-dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==singleton2),"stat"]
              
              
            } else if (parent(tree, singleton2) %in% poly_parents[k]) {
              
              singletons[k]<-singleton2
              
              poly_parents_with_terminal_leaves_dist[k]<-dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==singleton2),"stat"]
              
              
              
            }
          }
          
          
          child<-child[child %in% pheno_nodes_poly_desc]
          singleton<-child[child %in% res_nodes_leaves] ### direct singleton descendant of poly parent node
          
          
          child_terminal_tips<-child[child %in% res_nodes_leaves] ### do not consider if a child of a polyphyletic node is a terminal node (especially singleton), otherwise is similar to emergence rate
          
          if (length(child_terminal_tips)>1) {
            
            
            poly_parents_with_terminal_leaves_dist[k]<-max(vapply(child_terminal_tips,function(j){dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==j),"stat"]}, FUN.VALUE = 1))
            
          }
          
          
          if (length(singleton) == 1) {
            
            singletons[k]<-singleton
            #poly_parents_with_terminal_leaves_dist[k]<-singleton
            poly_parents_with_terminal_leaves_dist[k]<-dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==singleton),"stat"]
            
            
          }
          
          
          if (length(child) > 0 & length(child) < 2) {
            
            child<-child[!child %in% res_nodes_leaves]
            
            if (length(child)>0) {
              
              check_singleton_progeny<-Children(tree, child)
              
              ### if the son of the internal node son of the poly parent is singleton not consider (there is the singleton function dist for this)
              check_singleton_progeny<-check_singleton_progeny[check_singleton_progeny %in% single_clust_id]
              
              ### if no singleton as child go on with classic entry rate
              if (is_integer_zero(check_singleton_progeny) | length(check_singleton_progeny)==0 ) {
                
                sm[k]<-dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==child),"stat"]
                
              }
            }
            
          } else if (length(child)==2) {
            
            child<-child[!child %in% res_nodes_leaves]
            # 
            # chil_tips<-child[!child %in% res_nodes_leaves]
            # 
            if (length(child)>0) {
              check_singleton_progeny<-Children(tree, child)
              
              ### if the son of the internal node son of the poly parent is singleton not consider (there is the singleton function dist for this)
              check_singleton_progeny<-check_singleton_progeny[check_singleton_progeny %in% single_clust_id]
              
              ### if no singleton as child go on with classic entry rate
              if (is_integer_zero(check_singleton_progeny) | length(check_singleton_progeny)==0 ) {
                
                test_child_states <- warning(sprintf("Polyphyletic parent node %s has two childs with same opposite states, calculating the max distance", getMRCA(tree, child)), call. = FALSE) ### check if poly parent node has two child with different states instead of one 
                
                ### inserire polyphyly+1 if two R child?
                
                sm[k]<-max(vapply(child,function(j){dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==j),"stat"]}, FUN.VALUE = 1))
                
              }
            }
            
          }
          
        }
        
        singletons<-singletons[!is.na(singletons)]
        
        ###
        poly_parents_with_terminal_leaves_dist<-poly_parents_with_terminal_leaves_dist[!is.na(poly_parents_with_terminal_leaves_dist)]
        
        poly_parents_granchild_singletons<-parent(tree, single_clust_id)[!parent(tree, single_clust_id) %in% poly_parents] ### se i padri dei singleton non sono nei poli parent significa che chi è escluso ha un nonno poly parent con figlio R e nipote singleton R
        
        ### ENTRY RATE INCLUDING SINGLETON
        
        sm_single<-sum(sm,poly_parents_with_terminal_leaves_dist, na.rm = T)
        
        
        df_clust<-as.data.frame(table(ll$cluster))
        colnames(df_clust)<-c("cluster","freq")
        
        #n_singletons<-length(df_clust[which(df_clust$freq==1),"cluster"])
        n_singletons<-length(single_clust_id)
        
        #    singletons<-singletons[!is.na(singletons)]
        singletons_parent<-parent(tree, singletons)
        
        
        
        poly_parents_granchild_singletons<-parent(tree, single_clust_id)[!parent(tree, single_clust_id) %in% poly_parents] ### se i padri dei singleton non sono nei poli parent significa che chi è escluso ha un nonno poly parent con figlio R e nipote singleton R
        
        
        smm <- c()
        #Ne <- data.frame()
        
        for (i in 1:length(poly_parents)){
          #     i<-6
          
          ### if poly parent belong to cluster where other nodes are internal or singleton == next (there is the emergence_rate_singleton to account for this)
          cluster_check<-ll[which(ll$poly_parent==poly_parents[i]),]
          
          cluster_check<-ifelse(cluster_check$id %in% singletons | cluster_check$id %in% ll[ll$node=="internal","id"], NA, poly_parents[i]) ### if poly events has a singleton children == terminal tips -> added to entry rate
          
          cluster_check <- unique(cluster_check[!is.na(cluster_check)])
          
          if (!is.numeric(cluster_check)){
            next  
          }
          
          else {
            
            first_res_poly_child<-Children(tree, poly_parents[i]) ### PICK THE TWO DESCENDANTS, SUPPOSED TO BE ON R AND ONE S
            
            ### CHECK 1
            
            first_res_poly_child<-first_res_poly_child[!first_res_poly_child %in% 1:length(tree$tip.label)] ### check if the node is not a terminal one
            
            ### CHECK 2
            
            first_res_poly_child<-first_res_poly_child[first_res_poly_child %in% poly_nodes_with_R_child] ### check which is resistant
            
            if (length(first_res_poly_child) > 0 & length(first_res_poly_child) < 2) {
              
              first_res_poly_child_desc<-getDescendants(tree, first_res_poly_child) ### get descendants of the resistant child of poly nodes
              
              
              
              ### ADDED FILTER CHECK
              first_res_poly_child_desc<-first_res_poly_child_desc[first_res_poly_child_desc %in% pheno_nodes_poly_desc]
              
              
              smm[i]<- length(first_res_poly_child_desc) / (max(vapply(first_res_poly_child_desc,function(j){nodeheight(tree, node = j)}, FUN.VALUE = 1)-nodeheight(tree, node = first_res_poly_child))/365)
              
              
            }  
            
            if (length(first_res_poly_child)==2) {
              
              
              first_res_poly_child_desc_1<-getDescendants(tree, first_res_poly_child[1]) ### get descendants of the resistant child of poly nodes
              first_res_poly_child_desc_2<-getDescendants(tree, first_res_poly_child[2]) ### get descendants of the resistant child of poly nodes
              
              
              ### ADDED FILTER CHECK
              first_res_poly_child_desc_1<-first_res_poly_child_desc_1[first_res_poly_child_desc_1 %in% pheno_nodes_poly_desc]
              first_res_poly_child_desc_2<-first_res_poly_child_desc_2[first_res_poly_child_desc_2 %in% pheno_nodes_poly_desc]
              
              
              smm_1<-length(first_res_poly_child_desc_1) / (max(vapply(first_res_poly_child_desc_1,function(j){nodeheight(tree, node = j)}, FUN.VALUE = 1)-nodeheight(tree, node = first_res_poly_child[1]))/365)
              
              smm_2<-length(first_res_poly_child_desc_2) / (max(vapply(first_res_poly_child_desc_2,function(j){nodeheight(tree, node = j)}, FUN.VALUE = 1)-nodeheight(tree, node = first_res_poly_child[2]))/365)
         
              smm_merge<-c(smm_1,smm_2)
              #
              smm<-append(smm, smm_merge)
              
              test_child_states <- warning(sprintf("Polyphyletic parent node %s has two childs with same opposite states", getMRCA(tree, first_res_poly_child)), call. = FALSE) ### check if poly parent node has two child with different states instead of one
              #
            }
          }
          
        }
        
        
        emergence_rate_singleton_dist<-c() ### CALCULATE THE DISTANCE BETWEEN PARENT NODE AND TERMINAL R SINGLETON, TO BE ADDED INTO THE SMM VECTOR
        
        
        if (length(poly_parents_granchild_singletons)>0){
          
          for (i in 1:length(poly_parents_granchild_singletons)){
            
            first_res_poly_child_singleton_desc<-getDescendants(tree, poly_parents_granchild_singletons[i]) ### get descendants of the resistant child of poly nodes
            first_res_poly_child_singleton_desc<-first_res_poly_child_singleton_desc[first_res_poly_child_singleton_desc %in% res_nodes_leaves] 
            emergence_rate_singleton_dist[i]<- (max(vapply(first_res_poly_child_singleton_desc,function(j){nodeheight(tree, node = j)}, FUN.VALUE = 1)-nodeheight(tree, node = poly_parents_granchild_singletons[i])))/365
            
          }
        }
        
        
        
        smm_single<-sum(smm,na.rm = T)+sum(emergence_rate_singleton_dist, na.rm = T)
        
        
        ### emergence rate including also the singletons
        # smm_single<-sum(smm,na.rm = T)+sum(singleton_dist, na.rm = T)
        smm_single<-sum(smm,na.rm = T)+sum(emergence_rate_singleton_dist, na.rm = T)
        # 
        
        clust_tips<-length(res_nodes_leaves[res_nodes_leaves %in% ll$id])
        
        n_singletons<-length(singletons)
        
        
        #n_singletons<-length(single_clust_id)
        entry_rate_singletons<-length(single_clust_id)-length(poly_parents_granchild_singletons)
        #
        rates<-data.frame()
        
        rates[amr,"entry_rate_singletons"]<-(polyphyly) / (sum(sm_single, na.rm = T) / 365)
        rates[amr,"entry_rate_no_singletons"]<-(polyphyly-entry_rate_singletons) / (sum(sm, na.rm = T) / 365) ### CHECK HERE IF SINGLETONS ARE ACCOUNTED PROPERLY 
        
        rates[amr,"emergence_rate_singletons"]<- sum(smm_single, na.rm = T)
        #rates[amr,"emergence_rate_no_singletons"]<-(clust_tips-entry_rate_singletons) / (sum(smm, na.rm = T) / 365)
        
        ### RIMUOVERE DA QUI LE SEQUENZE CHE APPARTENGONO A ENTRY RATE!!!
        
        rates[amr,"emergence_rate_no_singletons"]<- sum(smm, na.rm = T) 
        
        rates[amr,"polyphyly"]<-polyphyly
        rates[amr,"polyphyly_no_singleton"]<-polyphyly-n_singletons
        
        
        ### N isolates that harbour the amr
        rates[amr, "n_singletons_all"]<-length(singletons)
        rates[amr, "n_singletons_filtered_by_desc_poly_nodes"]<-n_singletons
        
        rates[amr, "presence"]<-table(xx)[[1]]
        rates[amr, "absence"]<-table(xx)[[2]]
        
        rates[amr, "status"]<-ifelse(table(xx)[[1]]>table(xx)[[2]],"loss","gain")
        
        rates[amr, "S_to_R"]<-Q[2]
        rates[amr, "R_to_S"]<-Q[3]
        
        #write.csv(rates, sprintf("results/3_run/%s_rates.csv", amr))
        
        
        single_clust_id<-terminal_single[which(terminal_single$cluster %in% single_clust),"id"]
        
        
        
        find_nodes_to_remove <- function(tree, node_list) {
          nodes_to_remove <- c()
          
          check_descendants <- function(node) {
            # Get all descendants of the node using ape's descendants function
            descendants_nodes <- getDescendants(tree, node)
            
            # Find the intersection between descendants and the node_list
            matched_descendants <- descendants_nodes[descendants_nodes %in% node_list]
            
            if (length(matched_descendants) > 0) {
              # If multiple matched descendants are found, keep the most recent one
              most_recent_descendant <- matched_descendants[which.min(matched_descendants)]
              
              nodes_to_remove <<- c(nodes_to_remove, node)
            }
          }
          
          for (node in sort(node_list)) {
            check_descendants(node)
          }
          
          return(nodes_to_remove)
        }
  
        
        nmost<-as.data.frame(table(ll$cluster))
        nmost<-as.character(nmost[which(as.numeric(nmost$Freq)>2),"Var1"])
        
        nclust<-unique(ll[which(ll$cluster %in% nmost),"poly_parent"])
        
        
        
        nodes_to_remove<-find_nodes_to_remove(tree, poly_parents)
        
        
        
        
        poly_parents<-poly_parents[!poly_parents %in% nodes_to_remove]
        
        #clust_tips<-length(res_nodes_leaves[res_nodes_leaves %in% ll$id])
        
        n_singletons<-length(singletons)
        
        
        
        
        
        count_lineages_in_interval <- function(tree, intervals) {
          n_lineages <- numeric(length(intervals) - 1)  # Vector to store counts for each interval
          # Iterate over the intervals
          for (i in 1:(length(intervals) - 1)) {
            # Get the tips in the current interval (those whose TMRCA is between the start and end of the interval)
            tips_in_interval <- which(node_depths >= intervals[i] & node_depths <= intervals[i + 1])
            # Count the number of tips in this interval
            
            n_lineages[i] <- length(tips_in_interval)
            
          }
          
          return(n_lineages)
          
        }
        
        
        
        get_lineages_in_interval <- function(tree, intervals) {
          
          n_lineages <- list()
          
          # Iterate over the intervals
          for (i in 1:(length(intervals) - 1)) {
            # Get the tips in the current interval (those whose TMRCA is between the start and end of the interval)
            tips_in_interval <- node_depths_descend_with_trait[which(node_depths_descend_with_trait$node_depths >= intervals[i] & node_depths_descend_with_trait$node_depths <= intervals[i + 1]),"descend_with_trait"]
            
            if (length(tips_in_interval)>0) {
              
              n_lineages[[i]] <- tips_in_interval
              
            }
            
          }
          return(n_lineages)
          
        }
        
        
        df_coal_total<-data.frame()
        df_coal_total_double<-data.frame()
        
        
        
        
        for (p in 1:length(poly_parents)) {
          
          
          
          tt<-extract.clade(tree, poly_parents[p])
          
          tt<-di2multi(tt) ### remove multifurcation 
          
          # 
          entry_rate<-c()
          
          if (tt$Nnode>2) { ### check if the tree has at least 2 internal node to estimate coalescence intervals (singleton excluded by definition)
            
            n_intervals<-5
            
            
            
            #### ENTRY RATE
            
            
            child<-Children(tree, poly_parents[p]) ### get child of polyphyly where S became R
            
            
            child_nextgen1<-child[1]
            child_nextgen2<-child[2]
            
            child_nextgen1_desc<-getDescendants(tree, child_nextgen1) 
            child_nextgen2_desc<-getDescendants(tree, child_nextgen2) 
            
            ### get most recent desc R clade tips
            
            child_nextgen1_trait<-child_nextgen1_desc[child_nextgen1_desc %in%  pheno_nodes]
            child_nextgen2_trait<-child_nextgen2_desc[child_nextgen2_desc %in%  pheno_nodes]
            
            
            
            if (length(child_nextgen1_trait)>0) {
              
              poly_parent_child_1<-nodeheight(tree, child[1]) 
              
              #desc_nextgen1_trait_node_depths<-vapply(child_nextgen1_trait,function(j){nodeheight(tree,j)}, FUN.VALUE = 1)
              
              entry_rate_1 <- abs(poly_parent_child_1 - nodeheight(tree, poly_parents[p]))
              
              entry_rate<-append(entry_rate, entry_rate_1)
              
            }
            
            if (length(child_nextgen2_trait)>0) {
              
              poly_parent_child_2<-nodeheight(tree, child[2]) 
              
              entry_rate_2 <- abs(poly_parent_child_2 - nodeheight(tree, poly_parents[p]))
              
              entry_rate<-append(entry_rate, entry_rate_2)
              
            }
            
            
            entry_rate<-max(entry_rate)
            
            ### TAKE ONLY RES CHILD
            
            first_res_poly_child<-Children(tree, poly_parents[p]) ### PICK THE TWO DESCENDANTS, SUPPOSED TO BE ON R AND ONE S
            
            first_res_poly_child<-first_res_poly_child[!first_res_poly_child %in% 1:length(tree$tip.label)] ### check if the node is not a terminal one
            #
            first_res_poly_child<-first_res_poly_child[first_res_poly_child %in% poly_nodes_with_R_child] ### check which is resistant
            
            if (length(first_res_poly_child) > 0 & length(first_res_poly_child) < 2) {
              
              
              ### if it is both I need to pair both rates in each clade
              
              descend<-getDescendants(tree, first_res_poly_child)
              
              node_depths<-vapply(descend,function(j){nodeheight(tree,j)}, FUN.VALUE = 1)
              
              tt<-extract.clade(tree, first_res_poly_child)
              
              tt<-di2multi(tt) ### remove multifurcation 
              
              descend_tips_trait<-length(tt$tip.label[tt$tip.label %in% names(xx[xx=="R"])])
              
              R_vs_S<-round(descend_tips_trait/length(tt$tip.label)*100,2)
              
              descend_with_trait<-ancstats[which(ancstats$node %in% descend),]
              
              descend_with_trait<-descend_with_trait[which(descend_with_trait$R>descend_with_trait$S),]
              
              descend_with_trait<-descend_with_trait[,"node"] 
              
              
              node_depths<-vapply(descend_with_trait,function(j){nodeheight(tree,j)}, FUN.VALUE = 1)
              
              node_depths_descend_with_trait<-cbind.data.frame(node_depths, descend_with_trait)
              
              
              # Calculate the range of node depths
              min_depth <- min(node_depths)
              max_depth <- max(node_depths)
              
              # Split the range of node depths into intervals
              intervals <- seq(min_depth, max_depth, length.out = n_intervals + 1)
              
              
              lineages_per_interval <- count_lineages_in_interval(tree, intervals)
              
              
              ### get internal and terminal labels for each interval to check singletons
              
              labels_per_interval <- get_lineages_in_interval(tree, intervals)
              
              
              
              time_node_diff<-intervals[2]-intervals[1]
              
              
              
              
              prob_density_coal<-c()  
              
              df_coal<-data.frame()
              
              
              for (f in 1:length(lineages_per_interval)) {
                
                if (lineages_per_interval[f]>1) {
                  
                  
                  teta <- 1 / (lineages_per_interval[f] * time_node_diff) ### coalescente rate approximation
                  
                  prob_density <- exp(-time_node_diff * (choose(lineages_per_interval[f], 2) * teta))*teta ### Not including singletons in contributing to coalescence but only the N
                  
                  prob_density_coal<-append(prob_density_coal, prob_density)
                  
                  #coal_int<-coalescent.intervals(as.ultrametric(tt))$interval.length
                  
                  df_coal[f, "lineages_x_interval"]<-lineages_per_interval[f]
                  df_coal[f, "prob_density"]<-prob_density
                  df_coal[f, "maxnodeheight"]<-max(vapply(labels_per_interval[[f]],function(j){nodeheight(tree,j)}, FUN.VALUE = 1))
                  
                  
                } else {
                  
                  
                  lab_check<-ll[which(ll$id %in% labels_per_interval[[f]]),"node"]
                  #lab_check<-"terminal"
                  
                  
                  if  (length(lab_check>0)) { 
                    
                    if (lab_check=="internal") {
                      
                      teta <- 1 / (1 * time_node_diff) ### coalescente rate approximation
                      
                      prob_density <- exp(-time_node_diff * (choose(lineages_per_interval[f], 2) * teta))*teta ### Not including singletons in contributing to coalescence but only the N
                      
                      prob_density_coal<-append(prob_density_coal, prob_density)
                      
                      #coal_int<-coalescent.intervals(as.ultrametric(tt))$interval.length
                      
                      df_coal[f, "lineages_x_interval"]<-lineages_per_interval[f]
                      df_coal[f, "prob_density"]<-prob_density
                      df_coal[f, "maxnodeheight"]<-max(vapply(labels_per_interval[[f]],function(j){nodeheight(tree,j)}, FUN.VALUE = 1))
                      
                      
                    }
                  }
                }
              }
              
              
              
              df_coal_total[p, "prob_density_coal_rate"]<-prod(df_coal$prob_density, na.rm = T)
              df_coal_total[p, "amr"]<-amr
              df_coal_total[p, "node_lineages"]<-sum(df_coal$lineages_x_interval, na.rm = T)
              df_coal_total[p, "tip_lineages"]<-descend_tips_trait
              
              df_coal_total[p, "poly_parent"]<-poly_parents[p]
              df_coal_total[p, "polyphyly"]<-polyphyly-n_singletons
              df_coal_total[p, "coalescent_interval"]<-time_node_diff
              df_coal_total[p, "poly_parent_nodeheight"]<-nodeheight(tree, poly_parents[p])
              df_coal_total[p, "entry_rate"]<-entry_rate
              df_coal_total[p, "entry_rate"]<-entry_rate
              df_coal_total[p, "R_vs_S_prop"]<-R_vs_S
              df_coal_total[p, "poly_parent_state_prob"]<-ancstats[which(ancstats$node %in% poly_parents[p]),"S"]
              
              
              
              #df_coal_total[paste("node", p, sep = "_"), "ratio"]<-log10(sum(df_coal$lineages)/sum(df_coal$prob_density))
              df_coal_total<-df_coal_total[!is.na(df_coal_total$prob_density_coal_rate),]
              
              write.csv(df_coal_total, sprintf("%s/%s.csv", output_dir, amr))
              
              
            }  
            
            if (length(first_res_poly_child)==2) {
              
              
              for (c in 1:length(first_res_poly_child)) {
                
                descend<-getDescendants(tree, first_res_poly_child[c])
                
                tt<-extract.clade(tree, first_res_poly_child[c])
                tt<-di2multi(tt) ### remove multifurcation 
                
                
                
                descend_tips_trait<-length(tt$tip.label[tt$tip.label %in% names(xx[xx=="R"])])
                
                R_vs_S<-round(descend_tips_trait/length(tt$tip.label)*100,2)
                
                descend_with_trait<-ancstats[which(ancstats$node %in% descend),]
                
                descend_with_trait<-descend_with_trait[which(descend_with_trait$R>descend_with_trait$S),]
                
                descend_with_trait<-descend_with_trait[,"node"] 
                
                
                node_depths<-vapply(descend_with_trait,function(j){nodeheight(tree,j)}, FUN.VALUE = 1)
                
                node_depths_descend_with_trait<-cbind.data.frame(node_depths, descend_with_trait)
                
                
                # Calculate the range of node depths
                min_depth <- min(node_depths)
                max_depth <- max(node_depths)
                
                # Split the range of node depths into intervals
                intervals <- seq(min_depth, max_depth, length.out = n_intervals + 1)
                
                
                lineages_per_interval <- count_lineages_in_interval(tree, intervals)
                
                
                ### get internal and terminal labels for each interval to check singletons
                
                labels_per_interval <- get_lineages_in_interval(tree, intervals)
                
                
                
                time_node_diff<-intervals[2]-intervals[1]
                
                prob_density_coal<-c()  
                
                df_coal_double<-data.frame()
                
                for (f in 1:length(lineages_per_interval)) {
                  
                  if (lineages_per_interval[f]>1) {
                    
                    
                    teta <- 1 / (lineages_per_interval[f] * time_node_diff) ### coalescente rate approximation
                    
                    prob_density <- exp(-time_node_diff * (choose(lineages_per_interval[f], 2) * teta))*teta ### Not including singletons in contributing to coalescence but only the N
                    
                    #prob_density_coal<-append(prob_density_coal, prob_density)
                    
                    #coal_int<-coalescent.intervals(as.ultra
                    df_coal_double[f, "lineages_x_interval"]<-lineages_per_interval[f]
                    df_coal_double[f, "prob_density"]<-prob_density
                    df_coal_double[f, "maxnodeheight"]<-max(vapply(labels_per_interval[[f]],function(j){nodeheight(tree,j)}, FUN.VALUE = 1))
                    
                    
                  } else {
                    
                    
                    lab_check<-ll[which(ll$id %in% labels_per_interval[[f]]),"node"]
                    #lab_check<-"terminal"
                    
                    
                    if  (length(lab_check>0)) { 
                      
                      if (lab_check=="internal") {
                        
                        teta <- 1 / (1 * time_node_diff) ### coalescente rate approximation
                        
                        prob_density <- exp(-time_node_diff * (choose(lineages_per_interval[f], 2) * teta))*teta ### Not including singletons in contributing to coalescence but only the N
                        
                        #prob_density_coal<-append(prob_density_coal, prob_density)
                        
                        #coal_int<-coalescent.intervals(as.ultrametric(tt))$interval.length
                        
                        df_coal_double[f, "lineages_x_interval"]<-lineages_per_interval[f]
                        df_coal_double[f, "prob_density"]<-prob_density
                        df_coal_double[f, "maxnodeheight"]<-max(vapply(labels_per_interval[[f]],function(j){nodeheight(tree,j)}, FUN.VALUE = 1))
                        
                        
                      }
                    }
                  }
                }
                
                
                
                df_coal_total_double[c, "prob_density_coal_rate"]<-prod(df_coal_double$prob_density, na.rm = T)
                df_coal_total_double[c, "amr"]<-amr
                df_coal_total_double[c, "node_lineages"]<-sum(df_coal_double$lineages_x_interval, na.rm = T)
                df_coal_total_double[c, "tip_lineages"]<-descend_tips_trait
                
                df_coal_total_double[c, "poly_parent"]<-poly_parents[p]
                df_coal_total_double[c, "polyphyly"]<-polyphyly-n_singletons
                df_coal_total_double[c, "coalescent_interval"]<-time_node_diff
                df_coal_total_double[c, "poly_parent_nodeheight"]<-nodeheight(tree, poly_parents[p])
                df_coal_total_double[c, "entry_rate"]<-entry_rate
                df_coal_total_double[c, "R_vs_S_prop"]<-R_vs_S
                df_coal_total_double[c, "poly_parent_state_prob"]<-ancstats[which(ancstats$node %in% poly_parents[c]),"S"]
                
                
                
                #df_coal_total_double[paste("node", p, sep = "_"), "ratio"]<-log10(sum(df_coal_double$lineages)/sum(df_coal_double$prob_density))
                df_coal_total_double<-df_coal_total_double[!is.na(df_coal_total_double$prob_density_coal_rate),]
                
                write.csv(df_coal_total_double, sprintf("%s/%s_double.csv", output_dir, amr))
                
                
              }
                
            }      
            
          }    
          
        } # node check
            
      } # poly loop     
    
    } #poly check    
    
  }  #gain
  
  break} #amr pool



