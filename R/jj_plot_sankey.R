
##use ggsankey instead?
# plot_sankey = function(df, group_df=NULL){
#   #make sankey plot
#   #group df is optional and has same dim as df and contains the group mapping for each label for colouring
#   #df input example (columns are ordered by increasing some parameter), 
#   #each row is one observation where the label is tracked over the increasing parameter:
#   # power03 power04 power05
#   # CD4_1   CD4_2   CD4_3
#   # gdT_1   gdT_2   gdT_3
#   # CD4_1   CD4_2   CD8_3
#   #group_df
#   # power03 power04 power05
#   # CD4   CD4   CD4
#   # gdT   gdT   gdT
#   # CD4   CD4   CD8 
#   
#   library(networkD3)
#   
#   df = as.data.frame(df)
#   #make the labels unique for each column by addint the column number to it
#   for(i in 1:ncol(df)){
#     df[, i] = paste(df[, i], i, sep='_')
#   }
#   #for each column from left to right, make a data.frame of source and target value and count the number of observations for each combination
#   links_list = list()
#   for(i in 1:(ncol(df)-1)){
#     links_list[[i]] = df %>% dplyr::group_by_at(c(i, i+1)) %>% dplyr::summarise(n())
#     colnames(links_list[[i]]) = c('source', 'target', 'value')
#   }
#   links_df = do.call(rbind, links_list)
#   #nodes are required, name them in numbers by the available source and target values
#   nodes = 0:(length( unique(c(links_df$source, links_df$target))) -1)
#   names(nodes) = unique(c(links_df$source, links_df$target))
#   links_df$source = as.integer(from_to(links_df$source, nodes))
#   links_df$target = as.integer(from_to(links_df$target, nodes))
#   #add the column names as node group
#   nodes_df = data.frame(node_name= names(nodes))#, node_group = rep(colnames(df), each=3))
#   
#   if(!is.null(group_df)){
#     group_df = as.data.frame(group_df)
#     mapping_df = cbind(matrix_to_vector(df), matrix_to_vector(group_df))
#     mapping_df = mapping_df[!duplicated(mapping_df), ]
#     nodes_df$group = plyr::mapvalues(nodes_df$node_name, from= mapping_df[,1], to=mapping_df[,2])
#     # #get unique values in each column
#     # unique_mat = apply(df, 2, unique)
#     annotationColors <- msPickSampleColors(nodes_df$group)
#     my_color <- paste0("d3.scaleOrdinal() .domain([\"", 
#                        paste(names(annotationColors), collapse = "\", \""), 
#                        "\"]) .range([\"", paste(annotationColors, collapse = "\", \""), 
#                        "\"])")
#     
#     sankeyNetwork(Links = links_df, Nodes = nodes_df, 
#                   Source = "source", Target = "target", Value = "value", NodeID = "node_name", NodeGroup = "group",
#                   fontSize = 14, nodeWidth = 30, colourScale = my_color )
#   }else{
#     sankeyNetwork(Links = links_df, Nodes = nodes_df, 
#                   Source = "source", Target = "target", Value = "value", NodeID = "node_name", 
#                   fontSize = 14, nodeWidth = 30)
#   }
#   # save the widget
#   # library(htmlwidgets)
#   # saveWidget(p, file=paste0( getwd(), "/HtmlWidget/sankeyColor1.html"))
# }