Exercise 1

1.3 First part: Examine the PCA plot. Does everything look okay (We wouldn’t ask if it did)?
I would argue that the PCA does not look okay. Fe and LFC-Fe samples seem to cluster near eatch other which is inlikely. This likely means there is a tissue labeling or metadata parsing issue since the two tissues are not transcriptionally similar. The data should be checked and lables should be corrected before downstream k-means analysis.


1.3 What does the PCA plot suggest to you about the data? Why?

When the PCA plot exhibits tight clustering of replicates, it indicates that the data is quite high quality and region-specific expression profiles. Replicates that are from the same tissue will cluster together while distinct tissues are seperated along the prrincipal component axis, suggesting that each midgut region has unique gene expression. Based on the plot, samples align PC1 in a similar order to the physical position along the Drosophila midgut, which means that PC1 is an anterior-posterior gradient in gene expression. On the other hand, PC2 distingusihes tissues that might be involved in metal ion transport or oxidative metabolism.

What feature explains the first principal component (simply saying “tissue” is not sufficient)?

As mentioned in the pervious answer, PC1 is the transcriptional expression change that happens when moving from the anterior digestive regions through the tissues that handle metal ion transport to posterior stem cell rich regions. It basically sorts spatialization and function differention along the midgut anterior-posterior axis.

Exercise 3: Do the categories of enrichments make sense? Why?

I selected cluster 1 and cluster 2 for gene ontology enrichment analysis through PANTHER. For cluster 1, I observed that it was strongly enriched for processes related to lipid metabolism and energy production (fatty acid β-oxidation, lipid transport, carboxylic-acid metabolism). This makes sense because the midgut region like the posterior compartment is knnown to be metabolically active in literature to produce ATP. For cluster 2, I observed it is enriched for protein processing and post-transcriptional modification, probably suggesting that these tissues have high endoplasmic-reticulum activity involved in producing digestive enzymes, pretty well known in the anterior gut region.