if("TimeCourseVectors" %in% .packages())
   detach("package:TimeCourseVectors", unload=TRUE)
library(RUnit)
library(TimeCourseVectors)
#----------------------------------------------------------------------------------------------------
# prepared in
# ~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/atacPrep.R

f <- system.file(package="TimeCourseVectors", "extdata", "klf1-atac-10k.RData")

tbl.28 <- get(load(f))
checkEquals(dim(tbl.28), c(28, 8))
requiredColumns <- c("chrom", "start", "end", "rep", "day")
checkTrue(all(requiredColumns %in% colnames(tbl.28)))

f2 <- system.file(package="TimeCourseVectors", "extdata", "klf1-atac-200k.RData")
tbl.250 <- get(load(f2))
checkEquals(dim(tbl.250), c(250, 8))
requiredColumns <- c("chrom", "start", "end", "rep", "day")
checkTrue(all(requiredColumns %in% colnames(tbl.250)))

#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    scoreColumn <- "fc"
    tcv <- TimeCourseVectors$new(tbl.28, scoreColumn)
    checkTrue(all(c("R6", "TimeCourseVectors") %in% class(tcv)))

    checkEquals(tcv$getScoreColumn(), scoreColumn)

} # test_ctor
#----------------------------------------------------------------------------------------------------
# small test to start: just 10 regions, 2 of them (4 and 7) are singletons
test_findRegionFamilies.28 <- function()
{
    printf("--- test_findRegionFamilies")
    scoreColumn <- "fc"
    scoreColumn <- "qvalScore"
    tcv <- TimeCourseVectors$new(tbl.28, scoreColumn)
    tcv$findRegionFamilies(minOverlap=100)
    tbl.ov <- tcv$getOverlaps()
    checkEquals(tcv$getFamilyCount(), 4)
    checkEquals(tcv$getFamily(1), c(1,3,6,8,11,14,17,20,24,28))
    checkEquals(tcv$getFamily(2), c(2,5,10,13,16,19,23,27))
    checkEquals(tcv$getFamily(3), c(12,15,18,21,25))
    checkEquals(tcv$getFamily(4), c(22,26))
    checkEquals(tcv$getOrphans(), c(4,7,9))
      #  all of these together should be the same as all rows
    checkEquals(sort(c(tcv$getOrphans(),
                       tcv$getFamily(1),
                       tcv$getFamily(2),
                       tcv$getFamily(3),
                       tcv$getFamily(4))), seq_len(nrow(tbl.28)))
    checkEquals(dim(tcv$getTable()), c(28, 8))
    if(FALSE){  # execute manually to take a look
       showGenomicRegion(igv, "chr19:12,753,089-13,002,509")
       tcv$displayFamily(1)
       tcv$displayFamily(2)
       tcv$displayFamily(3)
       tcv$displayFamily(4)
       removeTracksByName(igv, getTrackNames(igv)[-1])
       tcv$displayByDay()
       }

} # test_findRegionFamilies
#----------------------------------------------------------------------------------------------------
# larger test, spanning 200k
test_findRegionFamilies.250 <- function()
{
    printf("--- test_findRegionFamilies.250")
    scoreColumn <- "fc"
    scoreColumn <- "qvalScore"
    tcv <- TimeCourseVectors$new(tbl.250, scoreColumn)
    checkEquals(dim(tcv$getTable()), c(250, 8))

    tcv$findRegionFamilies(minOverlap=100)
    #tbl.ov <- tcv$getOverlaps()
    checkEquals(tcv$getFamilyCount(), 31)
    orphans <- tcv$getOrphans()
    checkEquals(length(orphans), 43)
    fam.1 <- tcv$getFamily(1)
    checkEquals(fam.1, c(1, 21, 29, 49, 64, 88, 114, 138, 166, 193, 224))
    if(FALSE){
       fam <- 1
       tcv$displayFamily(fam)  # family 1 is 102kb upstream of klf1
       removeTracksByName(igv, getTrackNames(igv)[-1])
       }
} # test_findRegionFamilies.250
#----------------------------------------------------------------------------------------------------
# a reasonable strategy is to have just one value per
# day,  the mean of all reps on that day.
# the first family has multple reps for day 4, 10, 11, 12 and 16,  not 8.
# so the collapse creates 6 days (4,8,10,11,12,16) from 11 rows
test_collapseRepsToDay <- function()
{
    printf("--- test_collapseRepsToDay")

    scoreColumn <- "qvalScore"
    tcv <- TimeCourseVectors$new(tbl.250, scoreColumn)
    tcv$findRegionFamilies(minOverlap=100)
    tcv$getFamilyCount() # 31
    orphans <- tcv$getOrphans()
    checkEquals(length(orphans), 43)
    fam.1 <- tcv$getFamily(1)
    checkEquals(dim(tbl.250[fam.1,]), c(11, 8))
       # combining reps winnows the data down to just 6 days
    tbl.collapse <- tcv$collapseToDay(family=1)
    checkEquals(dim(tbl.collapse), c(6, 2))
    checkEquals(tbl.collapse$day, c(4, 8, 10, 11, 12, 16))
    checkEqualsNumeric(tbl.collapse$score,
                       c(1853.831, 2879.559, 2165.749, 11176.986, 5298.387, 30117.214),
                       tol=1e-2)

} # test_collapseRepsToDay
#----------------------------------------------------------------------------------------------------
# we have two expression data sets at present, with these timepoints
#
# ~/github/TrenaProjectErythropoiesis/prep/import/srm-rna-averaged-final-for-paper/rna.tsv
#  103 TFs by 13 timepoints, nicely matched to srm values on the same days
#  see https://hoodlab.shinyapps.io/tf-srm-rna/
#  c("Day.0", "Day2", "Day4", "Day6", "Day7.5", "Day8", "Day8.5", "Day10", "Day10.5", "Day11", "Day11.5", "Day12", "Day.14")
#
#  ~/github/TrenaProjectErythropoiesis/inst/extdata/expression/mtx.fpkm.tmt.31760x9-interpolated.RData
#  31760 genes by 9 timepoints
#  c("D2", "D7_5", "D8", "D10", "D10_5", "D11", "D11_5", "D12", "D14"))
#
# KLF1 shows similar expression trajectory in both files
#
# the ATAC-seq has these timepoints (days)
#  4  8 10 11 12 16
#
# the task: infer those ATAC-seq values to 2,8,10,11,12,14
test_scoreAtacFamilyAgainstTargetGene <- function()
{
    printf("--- test_scoreAtacFamilyAgainstTargetGene")

      #------------------------------------------------------------------
      # first repeat the collapse of multiple reps to one value per day
      #-----------------------------------------------------------------

    scoreColumn <- "qvalScore"
    #tcv <- TimeCourseVectors$new(tbl.28, scoreColumn)
    tcv <- TimeCourseVectors$new(tbl.250, scoreColumn)
    tcv$findRegionFamilies(minOverlap=100)
    printf("family count: %d", tcv$getFamilyCount())
    orphans <- tcv$getOrphans()
    printf("orphan count: %d", length(orphans))

    f2 <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression/mtx.fpkm.tmt.31760x9-interpolated.RData"
    mtx.2 <- get(load(f2))
    dim(mtx.2)  # 31760 9
    checkEquals(colnames(mtx.2), c("D2", "D7_5", "D8", "D10", "D10_5", "D11", "D11_5", "D12", "D14"))
    targetGene <- "KLF1"
    rna.vec <- as.numeric(mtx.2[targetGene, c("D2", "D8", "D10", "D11", "D12", "D14")])

    for(i in seq_len(tcv$getFamilyCount())){
      fam <- tcv$getFamily(i)
      tbl.collapse <- tcv$collapseToDay(family=i)
      missingDays <- setdiff(c(4,8,10,11,12,16), tbl.collapse$day)
      if(length(missingDays) > 3) next;
      if(length(missingDays) > 0){
         printf("supplying missing days, count: %d ", length(missingDays))
         tbl.faux <- data.frame(day=missingDays, score=0)
         tbl.collapse <- rbind(tbl.faux, tbl.collapse)
         new.order <- order(tbl.collapse$day)
         tbl.collapse <- tbl.collapse[new.order,]
         }
      #browser()
      print(tbl.collapse)
      cor <- cor(rna.vec, tbl.collapse$score, method="spearman")
      #cor <- cor(rna.vec, tbl.collapse$score, method="pearson")
      printf("---- score for family %d: %f", i, cor)
      } # for i

} # test_scoreAtacFamilyAgainstTargetGene
#----------------------------------------------------------------------------------------------------
# at least two sources of brand lab rna measurements:
# ~/github/TrenaProjectErythropoiesis/prep/import/srm-rna-averaged-final-for-paper/
#   rna.tsv
# ~/github/TrenaProjectErythropoiesis/inst/extdata/expression/
#   mtx.fpkm.tmt.31760x9-interpolated.RData
# compare KLF1 expression:
# mtx.1["KLF1",]
#      Day.0      Day2      Day4      Day6    Day7.5      Day8    Day8.5     Day10   Day10.5     Day11   Day11.5     Day12    Day.14
#   12.60498  75.99634 131.98755 182.21960 233.11015 175.58850 202.71625 253.08715 257.27495 285.37190 389.15375 404.56240 564.27105
# mtx.2["KLF1",]
#        D2      D7_5        D8       D10     D10_5       D11     D11_5       D12       D14
#  75.99633 233.11015 175.58850 253.08715 257.27495 285.37190 389.15375 404.56240 564.27105

explore_erythroTfRNA <- function()
{
   printf("--- test_erythroTfRNA")

   f1 <- "~/github/TrenaProjectErythropoiesis/prep/import/srm-rna-averaged-final-for-paper/rna.tsv"
   file.exists(f1)
   mtx.1 <- as.matrix(read.table(f1, sep="\t", as.is=TRUE, header=TRUE, row.names=1, nrow=-1))
   dim(mtx.1)  # 103 13
   checkEquals(colnames(mtx.1), c("Day.0", "Day2", "Day4", "Day6", "Day7.5", "Day8", "Day8.5", "Day10", "Day10.5", "Day11", "Day11.5", "Day12", "Day.14"))
   plot(mtx.1["KLF1",], type="b")

   f2 <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression/mtx.fpkm.tmt.31760x9-interpolated.RData"
   mtx.2 <- get(load(f2))
   dim(mtx.2)  # 31760 9
   colnames(mtx.2)
   checkEquals(colnames(mtx.2), c("D2", "D7_5", "D8", "D10", "D10_5", "D11", "D11_5", "D12", "D14"))
   lines(mtx.2["KLF1",], type="b", col="blue")

} # explore_erythroTfRNA
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_findRegionFamilies.28()
    test_findRegionFamilies.250()
    test_collapseRepsToDay()
    test_scoreAtacFamilyAgainstTargetGene()

} # runTests
#----------------------------------------------------------------------------------------------------
# if(!interactive())
    runTests()
