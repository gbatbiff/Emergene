m(list = ls())
library(ggplot2)
library(ggtree)
library(ape)
library(treeio)
library(phytools)
library(phangorn)
library(stringr)
library(data.table)
library(ggplot2)
library(ggtree)
library(phytools)
library(TreeTools)
library(dplyr)
library(treeio)
library(R.utils)




tree_files<-list.files(path = "path_to_trees_folder_", pattern =  "*trees")

#tree_files<-tree_files[1:10]

for (t in tree_files) {
  
  t_sim<-gsub(".trees", "", t)
  
  tree<-read.beast(sprintf("path_to_trees_folder_/%s", t))
  
  ss<-as_tibble(tree)
  
  ss<-ss[which(!is.na(ss$label)),]
  
  pops<-unique(ss$type)
  
  
  df_coal_total<-data.frame()
  
  
  for (pop in pops) {
    
    gene<-ifelse(ss$type==pop,"R","S")
    
    Accession<-ss$label
    
    xx<-setNames(as.factor(gene), Accession)
    
    print (pop)
    
    tree<-read.nexus(sprintf("path_to_trees_folder_/%s",t))
    
    tree<-di2multi(tree)
    
    
    if (length(levels(xx))==1) {
      
      next
    }  
    
    
    
    if (table(xx)[[1]]<table(xx)[[2]]) {
      
      simmap_Q <- make.simmap(tree, xx, model = "ARD", tips = TRUE, pi ="estimated", nsim = 1)
      
      Q<-simmap_Q$Q
      
      simmap_mod <- make.simmap(tree, xx, model = "ARD", tips = TRUE, pi="estimated", nsim = 100, Q=Q)
      
      
      simmap<-summary(simmap_mod)
      
      
      ancstats_all <- as.data.frame(simmap$ace)
      
      ancstats_nodes<-head(ancstats_all, tree$Nnode)
      
      # Loop over all internal nodes (internal nodes are numbered from len(tips)+1 to len(tips) + len(nodes))
      
      ancstats_nodes$node<-length(tree$tip.label):length(tree$edge.length)+1
      
      ancstats_tips<-ancstats_all[-which(rownames(ancstats_all) %in% ancstats_nodes$node),]
      
      ancstats_tips$node<-1:(length(tree$tip.label))
      
      ancstats<-rbind.data.frame(ancstats_tips, ancstats_nodes)
      
      
      cols<-setNames(c("red", "blue"), c("S","R")) 
      
      ff<-data.frame(label = names(xx), stat = as.factor(xx) )
      ff<-unique(ff)
      tree2 <- full_join(tree, ff, by="label")
      
      p <- ggtree(tree2) +
        geom_tippoint(aes(color = stat), size=1) +
        geom_text(aes(label=node), size= 4, nudge_x = -1.5)+
        scale_color_manual(values = cols) +
        theme(legend.position = "right")
      #
      #  pies <- nodepie(ancstats, cols = 1:2)
      #  pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))
      #  p2 <- p + geom_inset(pies, width = .01, height = .01)
      # #
      #ggsave(plot = p2, filename = "01_try/tree_C.svg", device = "svg", height = 60, width = 45, limitsize = F)
      # # # # #
      
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
      
      ### polyphiletic group where parents were S (pheno==S is the parent) and one descendants became R
      poly_parents_only_R <- ndt[(cpheno1 == cpheno2 & pheno == "S" & cpheno1=="R"), id]  
      
      if (length(poly_parents)>0 | length(poly_parents_only_R)>0){
        
        poly_parents<-append(poly_parents,poly_parents_only_R)
        
        
        poly_parents <- poly_parents[order(as.data.frame(nodeHeights(tree))[poly_parents,2], decreasing = TRUE)] ### ordered nodes
        
        
        pheno_nodes <- ndt[pheno == "R", id] ### cases where nodes/leaves are R
        
        ### R descendant (nodes/leaves) excluding cases when polyphily occurs (when S became R and parent = S)
        
        pheno_nodes_poly_desc <- pheno_nodes[!pheno_nodes %in% poly_parents]  
        
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
            
            
            for(i in 1:length(sorted_poly_parents)){ ### iterate cases where from S to R
              
              gen1_parent<-Children(tree, sorted_poly_parents[i])
              
              g1<-getDescendants(tree, gen1_parent[1])
              g2<-getDescendants(tree, gen1_parent[2])
              
              ### preliminary step to evaluate if any poly parent with other poly nodes across descendants with singletons 
              
              if (length(g2)==1) {
                
                ph2 <- ifelse(g2 %in% pheno_nodes_poly_desc, "R", "S") ### checj if childrens belong to parent node descendants
                
                ### create a df with the features of the polyphylies node
            
                ll[[k]] <- data.table(cluster = i, id = g2, pheno = ph2, poly_parent = sorted_poly_parents[i]) 
            
                k <- k + 1
                
              }
              
              
              if (length(g1)==1) {
                
                ph1 <- ifelse(g1 %in% pheno_nodes_poly_desc, "R", "S") 
                
                ll[[k]] <- data.table(cluster = i, id = g1, pheno = ph1, poly_parent = sorted_poly_parents[i]) 
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
              
              
              ph1 <- ifelse(g1 %in% pheno_nodes_poly_desc, "R", "S") 
              ph2 <- ifelse(g2 %in% pheno_nodes_poly_desc, "R", "S") 
              
              if (length(g2)>1) {
                
                ll[[k]] <- data.table(cluster = i, id = g2, pheno = ph2, poly_parent = sorted_poly_parents[i]) 
                
                k <- k + 1
                
              } 
              
              if (length(g2)==1) {
                
                ll[[k]] <- data.table(cluster = i, id = g2, pheno = ph2, poly_parent = sorted_poly_parents[i]) 
                
                k <- k + 1
                
              }
              
              if (length(g1)>1) {
                
                ll[[k]] <- data.table(cluster = i, id = g1, pheno = ph1, poly_parent = sorted_poly_parents[i]) 
                
                k <- k + 1
                
              }
              
              if (length(g1)==1) {
                
                ll[[k]] <- data.table(cluster = i, id = g1, pheno = ph1, poly_parent = sorted_poly_parents[i]) 
                
                k <- k + 1
                
              }
              
              
            } 
            
            return(ll)
            
          }
          
          
          
          ll<-find_clusters(tree, sorted_poly_parents, pheno_nodes_poly_desc, hier_nodes)
          
          
          ll <- rbindlist(ll) ### bind the cluster df
          ll <- ll[pheno == "R"] ### select only cases where the leave/node is R
          
          ll[, cluster := .GRP, by = cluster]  
          
          ll<-as.data.frame(ll)
          
          ll<-unique(ll)
          
          
          
          set_cluster_na <- function(df) {
            
            df <- df %>%
              arrange(desc(poly_parent))
            
            prev_ids <- NULL
            
            for (i in 1:nrow(df)) {
              current_id <- df$id[i]
              
              if (current_id %in% prev_ids) {
                df$cluster[i] <- NA
              }
              
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
          
          ### create a df with the features of the polyphylies node
          
          res_nodes_leaves<-pheno_nodes[pheno_nodes %in% tree_label]
          
          
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
                
                singletons[k]<-singleton1
                
               
                
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
                
              } else if (parent(tree, singleton2) %in% poly_parents[k]) {
                
                singletons[k]<-singleton2
                
                poly_parents_with_terminal_leaves_dist[k]<-dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==singleton2),"stat"]
                
              }
            }
            
            
            child<-child[child %in% pheno_nodes_poly_desc]
            singleton<-child[child %in% res_nodes_leaves] ### direct singleton descendant of poly parent node
            
            child_terminal_tips<-child[child %in% res_nodes_leaves] 
            if (length(child_terminal_tips)>1) {
              
              
              poly_parents_with_terminal_leaves_dist[k]<-max(vapply(child_terminal_tips,function(j){dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==j),"stat"]}, FUN.
                                                                    VALUE = 1))
              
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
                
                check_singleton_progeny<-check_singleton_progeny[check_singleton_progeny %in% single_clust_id]
                
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
                
                if (is_integer_zero(check_singleton_progeny) | length(check_singleton_progeny)==0 ) {
                  
                  test_child_states <- warning(sprintf("Polyphyletic parent node %s has two childs with same opposite states, calculating the max distance", getMRCA(tree, child)), call. = FALSE) ### check 
                
                  if poly parent node has two child with different states instead of one 
                  
                  sm[k]<-max(vapply(child,function(j){dist_parent_child[which(dist_parent_child[,1]==poly_parents[k] & dist_parent_child[,2]==j),"stat"]}, FUN.VALUE = 1))
                  
                }
              }
              
            }
            
          }
          
          singletons<-singletons[!is.na(singletons)]
          
          ###
          poly_parents_with_terminal_leaves_dist<-poly_parents_with_terminal_leaves_dist[!is.na(poly_parents_with_terminal_leaves_dist)]
          
          poly_parents_granchild_singletons<-parent(tree, single_clust_id)[!parent(tree, single_clust_id) %in% poly_parents] ### se i padri dei singleton non sono nei poli parent significa che chi è esclus
          o ha un nonno poly parent con figlio R e nipote singleton R
          
          ### ENTRY RATE INCLUDING SINGLETON
          
          sm_single<-sum(sm,poly_parents_with_terminal_leaves_dist, na.rm = T)
          
          
          df_clust<-as.data.frame(table(ll$cluster))
          colnames(df_clust)<-c("cluster","freq")
          
          #n_singletons<-length(df_clust[which(df_clust$freq==1),"cluster"])
          n_singletons<-length(single_clust_id)
          
          #    singletons<-singletons[!is.na(singletons)]
          singletons_parent<-parent(tree, singletons)
          
          
          
          poly_parents_granchild_singletons<-parent(tree, single_clust_id)[!parent(tree, single_clust_id) %in% poly_parents] ### se i padri dei singleton non sono nei poli parent significa che chi è esclus
          o ha un nonno poly parent con figlio R e nipote singleton R
          
          clust_tips<-length(res_nodes_leaves[res_nodes_leaves %in% ll$id])
          
          n_singletons<-length(singletons)
          
          
          #n_singletons<-length(single_clust_id)
          entry_rate_singletons<-length(single_clust_id)-length(poly_parents_granchild_singletons)
          #
          rates<-data.frame()
          
          
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
              tips_in_interval <- node_depths_descend_with_trait[which(node_depths_descend_with_trait$node_depths >= intervals[i] & node_depths_descend_with_trait$node_depths <= intervals[i + 1]),"descend_
with_trait"]
              
              if (length(tips_in_interval)>0) {
                
                n_lineages[[i]] <- tips_in_interval
                
              }
              
            }
            return(n_lineages)
            
          }
          
          df_coal_total<-data.frame()
          df_coal_total_double<-data.frame()
          
          
          smm <- c()
          
          for (p in 1:length(poly_parents)) {
            
            #  p<-5  
            
            tt<-extract.clade(tree, poly_parents[p])
            
            entry_rate<-c()
            
            if (tt$Nnode>2) { 
              
              n_intervals<-1
              
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
              
              first_res_poly_child<-Children(tree, poly_parents[p]) 
              
              first_res_poly_child<-first_res_poly_child[!first_res_poly_child %in% 1:length(tree$tip.label)] ### check if the node is not a terminal one
              #
              first_res_poly_child<-first_res_poly_child[first_res_poly_child %in% poly_nodes_with_R_child] ### check which is resistant
              
              if (length(first_res_poly_child) > 0 & length(first_res_poly_child) < 2) {
                
                
                
                ### if it is both I need to pair both rates in each clade
                
                descend<-getDescendants(tree, first_res_poly_child)
                
                node_depths<-vapply(descend,function(j){nodeheight(tree,j)}, FUN.VALUE = 1)
                
                tt<-extract.clade(tree, first_res_poly_child)
                
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
                
                lineages_per_interval <- count_lineages_in_interval(tt, intervals)
                
                labels_per_interval <- get_lineages_in_interval(tt, intervals)
                
                time_node_diff<-intervals[2]-intervals[1]
                
                if (lineages_per_interval>1) {
                  
                  
                  
                  cluster_check <- ll[ll$poly_parent == poly_parents[p], ]
                  cluster_check <- ifelse(
                    cluster_check$id %in% singletons | 
                      cluster_check$id %in% ll[ll$node == "internal", "id"],
                    NA,
                    poly_parents[p]
                  )
                  cluster_check <- unique(cluster_check[!is.na(cluster_check)])
                  
                  if (!is.numeric(cluster_check)) next
                  
                  first_res_poly_child <- Children(tree, poly_parents[p])
                  
                  # Exclude terminal tips
                  first_res_poly_child <- first_res_poly_child[!first_res_poly_child %in% seq_along(tree$tip.label)]
                  
                  # Keep only nodes with resistant children
                  first_res_poly_child <- first_res_poly_child[first_res_poly_child %in% poly_nodes_with_R_child]
                  
                  if (length(first_res_poly_child) > 0 & length(first_res_poly_child) < 2) {
                    
                    first_res_poly_child_desc <- getDescendants(tree, first_res_poly_child)
                    first_res_poly_child_desc <- first_res_poly_child_desc[first_res_poly_child_desc %in% pheno_nodes_poly_desc]
                    first_res_poly_child_desc <- first_res_poly_child_desc[first_res_poly_child_desc %in% res_nodes_leaves]
                    
                    if (length(first_res_poly_child_desc) >= 1) {
                      smm[p] <- (1 - (length(first_res_poly_child_desc) / length(tree$tip.label))) / (
                        (max(vapply(first_res_poly_child_desc, function(j) nodeheight(tree, node = j), FUN.VALUE = 1)) - 
                           nodeheight(tree, node = first_res_poly_child)) / 365
                      )
                    }
                  }
                  
                  if (length(first_res_poly_child) == 2) {
                    
                    desc_1 <- getDescendants(tree, first_res_poly_child[1])
                    desc_2 <- getDescendants(tree, first_res_poly_child[2])
                    
                    desc_1 <- desc_1[desc_1 %in% pheno_nodes_poly_desc]
                    desc_2 <- desc_2[desc_2 %in% pheno_nodes_poly_desc]
                    
                    desc_1 <- desc_1[desc_1 %in% res_nodes_leaves]
                    desc_2 <- desc_2[desc_2 %in% res_nodes_leaves]
                    
                    if (length(desc_1) >= 2 | length(desc_2) >= 2) {
                      
                      smm_1 <- length(desc_1) / (
                        (max(vapply(desc_1, function(j) nodeheight(tree, node = j), FUN.VALUE = 1)) - 
                           nodeheight(tree, node = first_res_poly_child[1])) / 365
                      )
                      
                      smm_2 <- length(desc_2) / (
                        (max(vapply(desc_2, function(j) nodeheight(tree, node = j), FUN.VALUE = 1)) - 
                           nodeheight(tree, node = first_res_poly_child[2])) / 365
                      )
                      
                      smm_merge <- c(smm_1, smm_2)
                      smm <- append(smm, smm_merge)
                      
                      warning(sprintf(
                        "Polyphyletic parent node %s has two childs with same opposite states", 
                        getMRCA(tree, first_res_poly_child)
                      ), call. = FALSE)
                    }
                  }
            
                  
                  #if (lineages_per_interval[f]>1) {
                  df_coal_total[p, "pop"]<-pop
                  df_coal_total[p, "node_lineages"]<-lineages_per_interval
                  df_coal_total[p, "tip_lineages"]<-descend_tips_trait
                  df_coal_total[p, "poly_parent"]<-poly_parents[p]
                  df_coal_total[p, "polyphyly"]<-polyphyly-n_singletons
                  df_coal_total[p, "coalescent_interval"]<-time_node_diff
                  df_coal_total[p, "poly_parent_nodeheight"]<-nodeheight(tree, poly_parents[p])
                  df_coal_total[p, "entry_rate"]<-entry_rate
                  df_coal_total[p, "R_vs_S_prop"]<-R_vs_S
                  df_coal_total[p, "poly_parent_state_prob"]<-ancstats[which(ancstats$node %in% poly_parents[p]),"S"]
                  df_coal_total[p, "emergence_rate"]<-smm[p]
                  df_coal_total<-df_coal_total[!is.na(df_coal_total$emergence_rate),]
                  
                  write.csv(df_coal_total, sprintf("path_to_results_folder/%s_%s_sim.csv", pop, t))
                  
                }    
              }
            }
          }
        }      
      }    
    }
    
  } #pop pool
}


