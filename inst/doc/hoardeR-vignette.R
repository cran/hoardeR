## ----eval=FALSE----------------------------------------------------------
#  R> install.packages("hoardeR")

## ----eval=FALSE----------------------------------------------------------
#  R> install.packages("packageName")

## ----eval=FALSE----------------------------------------------------------
#  R> source("https://bioconductor.org/biocLite.R")
#  R> biocLite("Biostrings")

## ----eval=FALSE----------------------------------------------------------
#  R> install.packages("devtools")
#  R> library("devtools")
#  R> install_github("fischuu/hoardeR")

## ---- warning=FALSE, eval=FALSE------------------------------------------
#  R> library(hoardeR)

## ----eval=FALSE----------------------------------------------------------
#  R> head(species)

## ----eval=FALSE----------------------------------------------------------
#         Common.name         Scientific.name Taxon.ID       Ensembl.Assembly       Accession Variation.database Regulation.database Pre.assembly
#  1   Aardvark (Pre)   Orycteropus afer afer  1230840                      -               -                  -                   -      OryAfe1
#  2           Alpaca           Vicugna pacos    30538                vicPac1               -                  -                   -            -
#  3     Amazon molly        Poecilia formosa    48698 Poecilia_formosa-5.1.2 GCA_000485575.1                  -                   -            -
#  4     Anole lizard     Anolis carolinensis    28377              AnoCar2.0 GCA_000090745.1                  -                   -            -
#  5        Armadillo    Dasypus novemcinctus     9361              Dasnov3.0 GCA_000208655.2                  -                   -            -
#  6 Budgerigar (Pre) Melopsittacus undulatus    13146                      -               -                  -                   -    MelUnd6.3

## ----eval=FALSE----------------------------------------------------------
#  R> novelBed <- data.frame(Chr=c(11,18,3),
#                            Start=c(72554673, 62550696, 18148822),
#                            End=c(72555273, 62551296, 18149422),
#                            Gene=c("LOC1", "LOC2", "LOC3"))
#  
#  R> myFasta <- getFastaFromBed(novelBed, species="Bos taurus",
#  +                             fastaFolder="/home/daniel/fasta/")

## ----eval=FALSE----------------------------------------------------------
#  R> myFasta <- getFastaFromBed(novelBed, species="Bos taurus", release = "84",
#  +                             fastaFolder=NULL, version=NULL)

## ----eval=FALSE----------------------------------------------------------
#  Bos_taurus.UMD3.1.dna.chromosome.1.fa.gz
#  Bos_taurus.UMD3.1.dna.chromosome.2.fa.gz
#  Bos_taurus.UMD3.1.dna.chromosome.3.fa.gz
#  ...

## ----eval=FALSE----------------------------------------------------------
#  Ovis_aries.Oar_v4.0.dna.chromosome.1.fa.gz
#  Ovis_aries.Oar_v4.0.dna.chromosome.2.fa.gz
#  Ovis_aries.Oar_v4.0.dna.chromosome.3.fa.gz
#  ...

## ----eval=FALSE----------------------------------------------------------
#  R> myFasta <- getFastaFromBed(novelBed, species="Ovis aries", release = NULL,
#                                fastaFolder="/home/daniel/fasta/", version="Oar_v4.0")

## ----eval=FALSE----------------------------------------------------------
#  R> exportFA(myFasta, file="/home/daniel/myFasta.fa")

## ----eval=FALSE----------------------------------------------------------
#  R> novelFA <- importFA(file="/home/daniel/myFasta.fa")

## ----eval=FALSE----------------------------------------------------------
#  >Chr:Start-End

## ----eval=FALSE----------------------------------------------------------
#  >12:123-456

## ----eval=FALSE----------------------------------------------------------
#  R>  blastSeq(novelFA,
#  +            email="daniel.fischer@luke.fi",
#  +            xmlFolder=file.path(projFolder,"hoardeROut/"),
#  +            logFolder=file.path(projFolder,"hoardeRLog/"),
#  +            keepInMemory=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  Missing: 3
#  Running: 1
#  Finished: 0
#  Avg. Blast Time: 00:00:00
#  Total running time: 00:00:04
#  ---------------------------------------------------------------

## ----eval=FALSE----------------------------------------------------------
#  Run RW99J31C01R : 00:02:23
#  Missing: 1
#  Running: 1
#  Finished: 2
#  Avg. Blast Time: 00:01:10
#  Total running time: 00:02:40
#  ---------------------------------------------------------------

## ----eval=FALSE----------------------------------------------------------
#  R> xmls <- importXML(folder=file.path(projFolder,"hoardeROut/"))

## ----eval=FALSE----------------------------------------------------------
#  R> tableSpecies(xmls, exclude="Bos taurus")

## ----eval=FALSE----------------------------------------------------------
#  R> par(oma=c(5,0,0,0))
#  R> barplot(sort(tableSpecies(xmls, exclude="Bos taurus"), decreasing=TRUE), las=2)

## ----eval=FALSE----------------------------------------------------------
#  R> tableSpecies(xmls, species="Sus scrofa", locations = TRUE)

## ----eval=FALSE----------------------------------------------------------
#                                             Organism hitID hitLen hitChr  hitStart    hitEnd origChr origStart  origEnd
#  28 Sus scrofa breed mixed chromosome 4, Sscrofa10.2   494    644      4 105815870 105816509       3  18148822 18149422

## ----eval=FALSE----------------------------------------------------------
#  R> tableSpecies(xmls)
#  
#  Bos taurus Equus caballus     Sus scrofa     Ovis aries
#           6              1              1              3

## ----eval=FALSE----------------------------------------------------------
#  R> tableSpecies(xmls, species="Sus scrofa")
#  Sus scrofa
#           1

## ----eval=FALSE----------------------------------------------------------
#  R> species[grepl("Sus scrofa", species$Scientific.name),]

## ----eval=FALSE----------------------------------------------------------
#           Common.name Scientific.name Taxon.ID Ensembl.Assembly       Accession Variation.database Regulation.database Pre.assembly
#  57               Pig      Sus scrofa     9823      Sscrofa10.2 GCA_000003025.4                  Y                   Y            -
#  58 Pig FPC_map (Pre)  Sus scrofa map       NA                -               -                  -                   -          MAP

## ----eval=FALSE----------------------------------------------------------
#  R> ssannot <- getAnnotation(species = "Sus scrofa",
#  +                           annotationFolder="/home/daniel/annotation")

## ----eval=FALSE----------------------------------------------------------
#  R> pigHits <- tableSpecies(xmls, species="Sus scrofa", locations = TRUE)
#  R> pigInter <- list()
#  R> for(i in 1:nrow(pigHits)){
#  R>   pigInter[[i]] <- intersectXMLAnnot(pigHits[i,], ssannot)
#  R> }

## ----eval=FALSE----------------------------------------------------------
#  R> pigInter
#  [[1]]
#  Empty data.table (0 rows) of 15 cols: V1,V2,V3,V4,V5,V6...

## ----eval=FALSE----------------------------------------------------------
#  R> pigInter.flank <- list()
#  R> for(i in 1:nrow(pigHits)){
#  R>   pigInter.flank[[i]] <- intersectXMLAnnot(pigHits[i,], ssannot, flanking=100)
#  R> }
#  # This part transforms the list into a data frame and removes those 'hits' that
#  # do not report any intersection
#  R> pigInter.flank <- pigInter.flank[sapply(pigInter.flank,nrow)>0]
#  R> pigInter.flank <- do.call(rbind, pigInter.flank)
#  

## ----eval=FALSE----------------------------------------------------------
#  R> pigInter.flank
#     V1      V2   V3        V4        V5 V6 V7 V8                                                                                                    V9 origChr origStart  origEnd hitChr  hitStart    hitEnd
#  1:  4 ensembl gene 105858080 105858409  .  -  . gene_id "ENSSSCG00000006603"; gene_version "2"; gene_source "ensembl"; gene_biotype "protein_coding";       3  18148822 18149422      4 105815870 105816509

## ----eval=FALSE----------------------------------------------------------
#  R>        plotHit(
#  +             hits=pigInter.flank,
#  +             flanking=100,
#  +             diagonal=0.25 ,
#  +             hitSpecies = "Sus scrofa",
#  +             origSpecies = "Bos taurus",
#  +             fastaFolder = "/home/ejo138/fasta/",
#              # The following options are optional
#  +             window=NULL ,
#  +             which=NULL,
#  +             figureFolder = "/home/daniel/figures/",
#  +             figurePrefix = "pigIS"
#  +             )

## ---- fig.retina = NULL, fig.cap="Similarity between Pig and Cow", echo=FALSE----
knitr::include_graphics("./pigFlanking.png")

## ----eval=FALSE----------------------------------------------------------
#  R> tableSpecies(xmls, species="Ovis aries", locations = TRUE)

## ----eval=FALSE----------------------------------------------------------
#                                                                          Organism hitID hitLen hitChr  hitStart    hitEnd origChr origStart  origEnd
#  8   Ovis aries breed Texel chromosome 3, Oar_v4.0, whole genome shotgun sequence   577    616      3  33956005  33955391      11  72554673 72555273
#  16 Ovis aries breed Texel chromosome 14, Oar_v4.0, whole genome shotgun sequence   553    613     14  59259439  59260042      18  62550696 62551296
#  24  Ovis aries breed Texel chromosome 1, Oar_v4.0, whole genome shotgun sequence   539    604      1 101361981 101362569       3  18148822 18149422

## ---- eval=FALSE---------------------------------------------------------
#  R> oaannot4.0 <- getAnnotation(species = "Ovis aries", release="NCBI", version="Oar_v4.0",
#  +                              type="gff", annotationFolder="/home/daniel/Annotations/")

## ----eval=FALSE----------------------------------------------------------
#  R> oaHits <- tableSpecies(xmlNew, species="Ovis aries", locations=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  R> for(i in 0:24){
#  +    refName <- paste("NC_0",19458+i,".2",sep="")
#  +    oaannot4.0$V1[oaannot4.0$V1==refName] <- i+1
#  +   }
#  R> i <- 25
#  R> refName <- paste("NC_0",19458+i,".2",sep="")
#  R> oaannot4.0$V1[oaannot4.0$V1==refName] <- "X"

## ----eval=FALSE----------------------------------------------------------
#  R> sheepInter4.0 <- list()
#  R> for(i in 1:nrow(oaHits)){
#  +    sheepInter4.0[[i]] <- intersectXMLAnnot(oaHits[i,], oaannot4.0, flanking=50)
#  +  }
#  
#  R> sheepInter4.0 <- sheepInter4.0[sapply(sheepInter4.0,nrow)>0]
#  R> sheepInter4.0 <- do.call(rbind, sheepInter4.0)

## ----eval =FALSE---------------------------------------------------------
#        plotHit(
#              hits=sheepInter4.0,
#              flanking=50,
#              window=NULL ,
#              diagonal=0.25 ,
#              hitSpecies = "Ovis aries",
#              hitSpeciesVersion = "Oar_v4.0",
#              origSpecies = "Bos taurus",
#              fastaFolder = "/home/daniel/fasta/",
#              origAnnot=btannot,
#              hitAnnot=oaannot4.0,
#              coverage=TRUE,
#              bamFolder = "/home/daniel/bams/",
#              which=1:2,
#              figureFolder = "/home/daniel/figures/",
#              figurePrefix = "sheepNOVEL"
#                )

## ---- fig.retina = NULL, fig.cap="Similarity between Sheep and Cow", echo=FALSE----
knitr::include_graphics("./coverageExample.png")

