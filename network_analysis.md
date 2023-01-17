## Networking
```
#load R packages needed
library(cooccur)
library(visNetwork)

#read in data
#data should be species as row names, samples as columns
otu_tbl<-read.delim("otu_table.txt", row.names=1, header=T)

#read in metadata


#calculate cooccurance based on
co <- cooccur(otu_tbl, spp_names = TRUE)
 
#extract edges/nodes
nodes <- data.frame(id = 1:nrow(otu_tbl), label = rownames(otu_tbl), color = "#606482", shadow = TRUE) 
edges <- data.frame(from = co$sp1, to = co$sp2, color = ifelse(co$p_lt <= 0.05, "#B0B2C1", "#3C3F51"), dashes = ifelse(co$p_lt <= 0.05, TRUE, FALSE))

#visualize network
visNetwork(nodes = nodes, edges = edges) %>% visIgraphLayout(layout = "layout_with_kk")




```

https://medium.com/analytics-vidhya/how-to-create-co-occurrence-networks-with-the-r-packages-cooccur-and-visnetwork-f6e1ceb1c523
