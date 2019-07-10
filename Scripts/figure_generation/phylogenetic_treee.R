library(ape)

tree = read.tree(file = '~/Documents/SINATRA/r_data/phyloT_generated_tree_1560787761_newick.txt')

plot.phylo(tree,type = 'cladogram',cex = 1,show.node.label = FALSE,use.edge.length = TRUE,edge.width = 2,font = 4)
