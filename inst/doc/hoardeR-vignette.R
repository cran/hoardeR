## ----eval=FALSE----------------------------------------------------------
#  install.packages("hoardeR")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("packageName")

## ----eval=FALSE----------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("Biostrings")

## ----eval=FALSE----------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(c('Biostrings', 'GenomicRanges', 'bamsignals', 'IRanges', 'Rsamtools', 'snpStats'))

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  library("devtools")
#  install_github("fischuu/hoardeR")

## ---- warning=FALSE, eval=TRUE-------------------------------------------
library("hoardeR")

## ----eval=FALSE----------------------------------------------------------
#  species[1:6,1:5]

## ----eval=TRUE, echo=FALSE, comment=NA, tidy=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 400)
species[1:6,1:5]

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  findSpecies("Cattle")

## ----eval=TRUE, echo=FALSE, comment=NA, tidy=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
findSpecies("Cattle")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  myFasta
#  
#  $`>11:72554673-72555273`
#  [1] "cccaagaagcaggaatgagagtggcgctttttctgccccaggtaacggtc..."
#  
#  $`>18:62550672-62551296`
#  [1] "aggagatttgcctgcgaaacctctggttctcttagagcttccattcccgt..."
#  
#  Fasta sequences ommited to print: 1

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  summary(myFasta)
#  
#  Summary of fa object
#  ---------------
#  Sequences      : 3
#  Minimum length : 601
#  1st quartile   : 613
#  Median length  : 625
#  Average length : 631
#  3rd quartile   : 646
#  Maximum length : 667

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  myFasta2 <- getFastaFromBed(novelBed, species = "Bos taurus", assembly="UMD_3.1",
#                             fastaFolder = "/home/daniel/fasta/")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  myFasta <- getFastaFromBed(novelBed, species = "Ovis aries",
#                             fastaFolder = "/home/daniel/fasta/",
#                             export = "/home/daniel/fasta",
#                             fileName="seqFasta.fa")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  exportFA(myFasta, file="/home/daniel/myFasta.fa")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  novelFA <- importFA(file="/home/daniel/myFasta.fa")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  >Chr:Start-End

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  >12:123-456

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#    blastSeq(novelFA,
#             email="daniel.fischer@luke.fi",
#             xmlFolder="/home/daniel/results/hoardeR/Proj1-Out",
#             logFolder="/home/daniel/results/hoardeR/Proj1-Log")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Missing: 3
#  Running: 1
#  Finished: 0
#  Avg. Blast Time: NA
#  Total running time: 00:00:04
#  ---------------------------------------------------------------

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Run RW99J31C01R : 00:02:23
#  Missing: 1
#  Running: 1
#  Finished: 2
#  Avg. Blast Time: 00:01:10
#  Total running time: 00:02:40
#  ---------------------------------------------------------------

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  xmls <- importXML(folder=file.path(projFolder,"hoardeROut/"))

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  tableSpecies(xmls, exclude="Bos taurus")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  barplot(tableSpecies(xmls))

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  par(oma=c(5,0,0,0))
#  barplot(sort(tableSpecies(xmls, exclude="Bos taurus"), decreasing=TRUE), las=2)

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  tableSpecies(xmls, species="Sus scrofa", locations = TRUE)

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                              Organism hitID hitLen hitChr  hitStart    hitEnd origChr origStart  origEnd
#  15 Sus scrofa breed mixed chromosome 13, Sscrofa10.2   542    610     13 152539332 152538724       1  62550672 62551296

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  R> tableSpecies(xmls)
#  
#  Bos taurus Equus caballus     Sus scrofa     Ovis aries
#           6              1              1              3

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  R> tableSpecies(xmls, species="Sus scrofa")
#  Sus scrofa
#           1

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  findSpecies("Sus scrofa")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#      Organism.Name Organism.Common.Name Taxid Assembly.Name Assembly.Accession                            Assembly.Submitter    Assembly.Data  ...
#  317    Sus scrofa                  pig  9823   Sscrofa10.2    GCF_000003025.5 The Swine Genome Sequencing Consortium (SGSC) 7 September 2011  ...
#  \normalsize
#  
#  That means that the assembly can be obtained automatically from `hoardeR` for further analysis without using any additional parameters. The command `getAnnotation` downloads automatically the corresponding annotation into a folder specified in `annotationFolder` (here: `/home/daniel/annotation`). Again, if nothing is specified (`annotationFolder=NULL`, default), the working directory is used intead.
#  

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  ssannot <- getAnnotation(species = "Sus scrofa",
#  +                        annotationFolder="/home/daniel/annotation")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  pigHits <- tableSpecies(xmls, species="Sus scrofa", locations = TRUE)
#  pigInter <- intersectXMLAnnot(pigHits, ssannot)
#  

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  pigInter
#  
#  [[1]]
#  Empty data.table (0 rows) of 15 cols: V1,V2,V3,V4,V5,V6...

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  pigInter.flank <- intersectXMLAnnot(pigHits, ssannot, flanking=1000)

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  pigInter.flank
#  
#     V1      V2   V3        V4        V5 V6 V7 V8
#  1: 13 ensembl gene 153238002 153238055  .  -  .
#                                                                                                                     V9 origChr origStart  origEnd hitChr
#  1: gene_id "ENSSSCG00000019624"; gene_version "1"; gene_name "SNORD12"; gene_source "ensembl"; gene_biotype "snoRNA";       1  62550672 62551296     13
#      hitStart    hitEnd
#  1: 152539332 152538724

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
#  +             figurePrefix = "pigInter"
#  +             )

## ---- fig.retina = NULL, fig.cap="Similarity between Pig and Cow", echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::include_graphics("./pigFlanking.png")

## ----eval =FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#        plotHit(
#  +             hits=pigInter.flank,
#  +             flanking=100,
#  +             diagonal=0.25 ,
#  +             hitSpecies = "Sus scrofa",
#  +             origSpecies = "Bos taurus",
#  +             origAnnot=btannot,
#  +             hitAnnot=ssannot,
#  +             fastaFolder = "/home/daniel/fasta/",
#  +             figureFolder = "/home/daniel/figures/",
#  +             figurePrefix = "pigInter",
#  +             coverage=TRUE,
#  +             bamFolder = "/home/daniel/bams/",
#  +             groupIndex = c(1,1,2,2,1,2),
#  +             groupColor = c("blue", "red"))

## ---- fig.retina = NULL, fig.cap="Similarity between Pig and Cow, with added coverage and annotation", echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::include_graphics("./coverageExample.jpg")

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  # Load the package
#    library(hoardeR)
#  
#  # Define the project folder
#    projFolder <- "/home/daniel/hoardeR-Example/"
#  
#  # Set the working directory accordingly
#    setwd(projFolder)
#  
#  # Create some fake coordiantes in bed format
#    novelBed <- data.frame(Chr = c(4,11),
#                           Start = c(104591617,72554673),
#                           End = c(104591916,72555273),
#                           Gene = c("LOC1", "LOC2"))
#  
#  # Print the data
#    novelBed
#  
#  #  Chr     Start       End Gene
#  #1   4 104591671 104591916 LOC1
#  #2  11  72554673  72555273 LOC2
#  
#  # Get the fasta sequence for this
#    myFasta <- getFastaFromBed(novelBed, species="Sus scrofa")
#  
#  #  No directory with fasta files given! Use the working directory:
#  #    /home/daniel/hordeR-Example
#  #  Using species assembly: Sscrofa10.2
#  #  Local file not found! Try to download fasta file:  ssc_ref_Sscrofa10.2_chr4.fa.gz
#  #  Read fasta file:  /home/daniel/hoardeR-Example/ssc_ref_Sscrofa10.2_chr4.fa.gz
#  #  Local file not found! Try to download fasta file:  ssc_ref_Sscrofa10.2_chr11.fa.gz
#  #  Read fasta file:  /home/daniel/hoadeR-Example/ssc_ref_Sscrofa10.2_chr11.fa.gz
#  
#  # Blast the sequences
#    myBlastRes <- blastSeq(myFasta, email="daniel.fischer@luke.fi")
#  
#  # Create/use log folder: hoardeR-10.10.2016@16.43.48/logs
#  # Missing: 2
#  # Running: 1
#  # Finished: 0
#  # Avg. Blast Time: NA
#  # Total running time: 00:00:04
#  # ---------------------------------------------------------------
#  # Missing: 2
#  # Running: 2
#  # Finished: 0
#  # Avg. Blast Time: NA
#  # Total running time: 00:00:08
#  # ---------------------------------------------------------------
#  # Run Z779DC03014 : 00:00:06
#  # Run Z779H9KF014 : 00:00:07
#  # Missing: 2
#  # Running: 2
#  # Finished: 0
#  # Avg. Blast Time: NA
#  # Total running time: 00:00:19
#  # ---------------------------------------------------------------
#  # Run Z779H9KF014 : 00:02:28
#  # Missing: 1
#  # Running: 1
#  # Finished: 1
#  # Avg. Blast Time: 00:02:27
#  # Total running time: 00:02:39
#  # ---------------------------------------------------------------
#  # Missing: 0
#  # Running: 0
#  # Finished: 2
#  # Avg. Blast Time: 00:03:21
#  # Total running time: 00:03:50
#  # ---------------------------------------------------------------
#  
#  # Import the XML
#    xmls <- importXML(folder="/home/daniel/hoardeR-Example/hoardeR-10.10.2016@16.43.48")
#  
#  # Table the species
#    tableSpecies(xmls)
#  
#  #  Bos taurus   Capra hircus Equus caballus     Ovis aries     Sus scrofa
#  #  2              2              1              1              1
#  
#  # This means, we have a couple of good hits in other species: Horse, Pig and Sheep
#  # Lets consider the cow hits first:
#    cowHits <- tableSpecies(xmls, species = "Bos taurus", locations=TRUE)
#    cowHits
#  
#  #                                                                                               Organism hitID hitLen hitChr hitStart   hitEnd origChr origStart
#  # 8 Bos taurus breed Hereford chromosome 3, alternate assembly Btau_5.0.1, whole genome shotgun sequence   280    298      3 16657347 16657644       4 104591617
#  # 9          Bos taurus breed Hereford chromosome 3, Bos_taurus_UMD_3.1.1, whole genome shotgun sequence   280    298      3 16500510 16500807       4 104591617
#  #     origEnd
#  # 8 104591916
#  # 9 104591916
#  
#  # Here we see that there is one match in cow on Chromosome 3 with the Pig Chromosome 4 search area.
#  # Further we can see that one Hit assembly version is UMD_3.1.1 and the other Btau_5.0.1
#  # We check if this assemble is the default in hoardeR
#  
#  findSpecies("Bos taurus")
#  
#  #    Organism.Name Organism.Common.Name Taxid        Assembly.Name Assembly.Accession
#  # 45    Bos taurus               cattle  9913 Bos_taurus_UMD_3.1.1    GCF_000003055.6
#  
#  # The 'Ensembl.Assembly' matches, so we can use the default paramters.
#  
#  # First, get the required annotations
#    ssannot <- getAnnotation(species = "Sus scrofa")
#  #  No assembly version provided, use the default: Sscrofa10.2
#  #  No directory with annotations files given! Use the working directory:
#  #    /home/daniel/hoardeR-Example
#  #  Check if file  /home/daniel/hoardeR-Example/ref_Sscrofa10.2_top_level.gff3.gz  exists ...
#  #  ... file wasn't found. Try to download it from NCBI ftp server.
#  #  ... found!
#  #  Check if file  /home/daniel/hoardeR-Example/chr_accessions_Sscrofa10.2  exists ...
#  #  ... file wasn't found. Try to download it from NCBI ftp server.
#  #  ... found!
#  #  Read 1184334 rows and 1 (of 1) columns from 0.288 GB file in 00:00:03
#  
#    btannot <- getAnnotation(species = "Bos taurus")
#  #  No assembly version provided, use the default: Bos_taurus_UMD_3.1.1
#  #  No directory with annotations files given! Use the working directory:
#  #    /home/hoardeR-Example
#  #  Check if file  /home/hoardeR-Example/ref_Bos_taurus_UMD_3.1.1_top_level.gff3.gz  exists ...
#  #  ... file wasn't found. Try to download it from NCBI ftp server.
#  #  ... found!
#  #  Check if file  /home/daniel/hoardeR-Example/chr_accessions_Bos_taurus_UMD_3.1.1  exists ...
#  #  ... file wasn't found. Try to download it from NCBI ftp server.
#  #  ... found!
#  #  Read 1372833 rows and 1 (of 1) columns from 0.349 GB file in 00:00:04
#  
#  # Now intersect the horse annotation with the matches
#    cattleInter <- intersectXMLAnnot(cowHits, btannot)
#  
#  # Empty data.table (0 rows) of 15 cols: V1,V2,V3,V4,V5,V6...
#  
#  # No matches were found (=no intergenic hits), so we extend the search area with flanking
#  # sites of 20Mb
#  
#  cattleInter.flank <- intersectXMLAnnot(cowHits, btannot, flanking=20)
#  cattleInter.flank
#  
#  #    V1         V2   V3       V4       V5 V6 V7 V8  ...
#  # 1:  3     Gnomon gene 16481126 16481718  .  +  .  ...
#  # 2:  3     Gnomon gene 16498496 16499120  .  -  .  ...
#  # 3:  3 BestRefSeq gene 16503991 16507247  .  -  .  ...
#  # 4:  3     Gnomon gene 16507420 16507918  .  -  .  ...
#  # 5:  3 BestRefSeq gene 16511100 16515182  .  +  .  ...
#  
#  # Here we could find some matches, we are going to visualize it.
#  
#  # For example, we are interested in hit number 3 and we want to plot the anotation tracks
#  # of the hit and the original:
#    plotHit(hits=cattleInter.flank,
#            flanking=20,
#            hitSpecies="Bos taurus",
#            hitAnnot=btannot,
#            origSpecies="Sus scrofa",
#            origAnnot=ssannot,
#            which=3)
#  
#  # If Figure for the 2nd and third should be created and they should be written to the HDD,
#  # the command is as follows
#    plotHit(hits=cattleInter.flank,
#            flanking=20,
#            hitSpecies="Bos taurus",
#            hitAnnot=btannot,
#            origSpecies="Sus scrofa",
#            origAnnot=ssannot,
#            which=c(2,3),
#            figureFolder="/home/hoardeR-Example",
#            figurePrefix="exampleOne")
#  
#  # If in addition a coverage should be plotted, the command is as follows. The additional
#  # option 'bamFolder' is expected to contain sorted bam files that are used for the
#  # coverage plot. Further, the coverage curves can be splitted to several curves.
#  # In that case the index needs to be provided and the corresponding color.
#    plotHit(hits=cattleInter.flank,
#            flanking=20,
#            hitSpecies="Bos taurus",
#            hitAnnot=btannot,
#            origSpecies="Sus scrofa",
#            origAnnot=ssannot,
#            which=c(2,3),
#            figureFolder="/home/daniel/hoardeR-Example",
#            figurePrefix="exampleTwo",
#            bamFolder="/home/daniel/hoardeR-Example",
#            groupIndex=c(1,1,1,2,2,2),
#            groupColor=c("blue","green")
#            )
#  

## ----eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  # Import the XML
#    xmls <- importXML(folder="/home/daniel/hoardeR-Example/hoardeR-10.10.2016@16.43.48")
#  
#  # Table the species
#    tableSpecies(xmls)
#  
#  #  Bos taurus   Capra hircus Equus caballus     Ovis aries     Sus scrofa
#  #  2              2              1              1              1
#  
#  # This means, we have a couple of good hits in other species: Horse, Pig and Sheep
#  # Lets consider the cow hits first:
#    cowHits <- tableSpecies(xmls, species = "Bos taurus", locations=TRUE)
#    cowHits
#  
#  #                                                                                               Organism hitID hitLen hitChr hitStart   hitEnd origChr origStart
#  # 8 Bos taurus breed Hereford chromosome 3, alternate assembly Btau_5.0.1, whole genome shotgun sequence   280    298      3 16657347 16657644       4 104591617
#  # 9          Bos taurus breed Hereford chromosome 3, Bos_taurus_UMD_3.1.1, whole genome shotgun sequence   280    298      3 16500510 16500807       4 104591617
#  #     origEnd
#  # 8 104591916
#  # 9 104591916
#  
#  # Here we see that there is one match in cow on Chromosome 3 with the Pig Chromosome 4
#  # search area. Further we can see that one Hit assembly version is UMD_3.1.1 and the
#  # other Btau_5.0.1. We check if this assemble is the default in hoardeR
#  
#  findSpecies("Bos taurus")
#  
#  #    Organism.Name Organism.Common.Name Taxid        Assembly.Name Assembly.Accession
#  # 45    Bos taurus               cattle  9913 Bos_taurus_UMD_3.1.1    GCF_000003055.6
#  
#  # We are interested in the Btau5.0.1 assembly, so we cannot use the default options for
#  # hoardeR. We see that the hit is on Chromosome 3, so we need the corresponding
#  # chromosome fasta and also the assembly.
#  
#  # Both can be downloaded from NCBI here:
#  # ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/Assembled_chromosomes/seq/bt_alt_Btau_5.0.1_chr3.fa.gz
#  # ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/GFF/alt_Btau_5.0.1_top_level.gff3.gz
#  # ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/Assembled_chromosomes/chr_accessions_Btau_5.0.1
#  
#  # Then we import the annotation:
#  
#  btannot5.0.1 <- getAnnotation(species = "Bos taurus", assembly = "Btau_5.0.1",
#  +                             annotationFolder="/home/ejo138/tmp")
#  
#  # Here we get the error:
#  # check if file  /home/daniel/hoardeR_Example/ref_Btau_5.0.1_top_level.gff3.gz  exists ...
#  # ... file wasn't found. Try to download it from NCBI ftp server.
#  
#  #That means, the downloaded file is wrongly named, so we rename it to
#  # ref_Btau_5.0.1_top_level.gff3.gz
#  
#  # After renaming to the expected filename, the file was imported and used for the
#  # intersection (which provided basically the same results as the UMD3.1.1 assembly)
#  
#  cattleInter.flank <- intersectXMLAnnot(cowHits, btannot5.0.1, flanking=20)
#  cattleInter.flank
#  
#  # Again, we want to visualize the 3 finding, this time we specified directly the fasta
#  # Folder instead of relying on the working directory.
#  
#  plotHit(hits=cattleInter.flank,
#          flanking=20,
#          hitSpecies="Bos taurus",
#          hitSpeciesAssembly = "Btau_5.0.1",
#          hitAnnot=btannot5.0.1,
#          origSpecies="Sus scrofa",
#          origAnnot=ssannot,
#          which=3,
#          fastaFolder="/home/ejo138/tmp")
#  
#  # As the alternative assembly is also in the assembled chromosome folder, hoardeR
#  # is able to automatically adjust the filename and uses the alternative assembly.
#  
#  # In case an entire own genome/assembly is used, then the filenames of both, the
#  # gff and the fasta file need to be adjusted.

