
# Instal Distmap for mapping cells based on known expression patterns 

```{r}
install_github("rajewsky-lab/DistMap")
install.packages('scatterplot3d')
```

```{r}
library(DistMap)
library(scatterplot3d)

```

```{r}
raw.data <- as.matrix(read.table( "~pathway/dge_raw.txt", row.names = 1, check.names = FALSE))

normalized.data <- as.matrix(read.table( "~pathway/dge_normalized.txt.gz", row.names = 1,  check.names = FALSE))

#barcodes<- colnames(normalized.data)
#colnames(raw.data) <- barcodes


insitu.matrix <- as.matrix(read.delim( "~pathway/bdtnp.txt.gz", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))

insitu.matrix2 <- as.matrix(read.table( "~pathway/binarized_bdtnp.csv", header = TRUE, sep = ",", check.names = FALSE))
colnames(insitu.matrix)

geometry <- as.matrix(read.table( "~pathway/geometry.txt.gz", header = TRUE, check.names = FALSE))


namegeo <- c("x", "y", "z")
colnames(geometry) <- namegeo
```


```{r}
dm = new("DistMap",
         raw.data=raw.data,
         data=normalized.data,
         insitu.matrix=insitu.matrix2,
         geometry=geometry)

colnames(insitu.matrix)
colnames(dm@insitu.matrix)


```

```{r}
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)

```

```{r}
VISH <- computeVISH(dm, 'sna', threshold=0.75)
computeGeneGradient(dm, 'PCNA')

```
```{r}

gene1 = "hb"
VISH <- computeVISH(dm, gene1, threshold=0.75)

toPLOT <- as.data.frame(cbind(geometry[,"x"] , geometry[,"y"],geometry[,"z"],  VISH))
toPLOT2 <- as.data.frame(cbind(geometry[,"x"] , geometry[,"z"],geometry[,"y"],  VISH))


tocolor <- rgb(VISH, 0, 0)



scatterplot3d(toPLOT[,1:3],
              color= tocolor,
              pch = 20,
              angle= 45,
              grid=TRUE, box=TRUE,
              main = paste("VISH Embryo, red = ", gene1, sep ="")
              xlab = "A/P",
              ylab = "y",
              zlab = "D/V")

```

