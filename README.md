# Misi-n-1

#Importar el countdata
#Para poder importar los datos tuve que descargarlos directamente de la página ya que
#al colocar la URL esta no me daba acceso
countData = read.csv("GSE37704_featurecounts.csv",
                     row.names = 1) %>%
  dplyr::select(-length) %>%
  as.matrix()

#Filtrar los datos quitando las lecturas con 0 y 1
countData = countData[rowSums(countData) > 1,]
head(countData)

#Importar el metadata, de la misma manera la función no me leyó la URL
#así descargué directamente el archivo para que
#pudiese ser leído
colData = read.csv(("GSE37704_metadata.csv"),
                   row.names = 1)
colData

#Configurar DESeqDataSet
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ condition)
#Evaluaión de la expresión diferencial de los genes
#Analiza que genes cambian su expresión diferencial en todaas las réplicas
dds = DESeq(dds)
dds

#Se compara la disminución del gen HOXA1 contra el siRNA control y se reordenarán dependiendo de su p value
#Se utiliza la funión "summary" para conocer cuantos genes estan regulados positiva o negativamente
res = results(dds, contrast = c("condition", "hoxa1_kd", "control_sirna"))
res = res[order(res$pvalue),]
summary(res)

##Instalar y cargar las paqueterías "AnnotationDbi" y "org.Hs.eg.db" 
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

#Con "mapIds" se agregan más columnas a los resultados
res$symbol = mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
res$entrez = mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
res$name = mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "GENENAME",
                     keytype = "ENSEMBL",
                     multiVals = "first")
head(res, 10)


#Pasar al análisis de las vías KEGG y para esto se deben descargar y caragar las paqueterías "gage",
#"pathview" y "gageData"
library(gage)
library(pathview)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

#
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

#
keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)

lapply(keggres, head)

#
keggrespathways = data.frame(id=rownames(keggres$greater),
                             keggres$greater) %>% 
                             tbl_df() %>% 
                             filter(row_number()<=5) %>% 
                             .$id %>% 
                             as.character()
keggrespathways

#
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

#
plot_pathway = function(pid) pathview(gene.data=foldchanges,
                                      pathway.id=pid,
                                      species="hsa",
                                      new.signature=FALSE)
tmp = sapply(keggresids, function(pid)
  pathview(gene.data=foldchanges,
           pathway.id=pid,
           species="hsa"))

#Ontología de los genes
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets = gobpsets, same.dir=TRUE)
lapply(gobpres, head)
