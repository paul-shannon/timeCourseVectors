if("TimeCourseVectors" %in% .packages())
   detach("package:TimeCourseVectors", unload=TRUE)
library(RUnit)
library(TimeCourseVectors)
#----------------------------------------------------------------------------------------------------
# prepared in
# ~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/atacPrep.R

f <- system.file(package="TimeCourseVectors", "extdata", "klf1-atac-10k.RData")

tbl <- get(load(f))
checkEquals(dim(tbl), c(28, 8))
requiredColumns <- c("chrom", "start", "end", "rep", "day")
checkTrue(all(requiredColumns %in% colnames(tbl)))

f2 <- system.file(package="TimeCourseVectors", "extdata", "klf1-atac-200k.RData")
tbl2 <- get(load(f2))
# checkEquals(dim(tbl), c(28, 8))
requiredColumns <- c("chrom", "start", "end", "rep", "day")
checkTrue(all(requiredColumns %in% colnames(tbl2)))

#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    scoreColumn <- "fc"
    tcv <- TimeCourseVectors$new(tbl, scoreColumn)
    checkTrue(all(c("R6", "TimeCourseVectors") %in% class(tcv)))

    checkEquals(tcv$getScoreColumn(), scoreColumn)

} # test_ctor
#----------------------------------------------------------------------------------------------------
# small test to start: just 10 regions, 2 of them (4 and 7) are singletons
test_findRegionFamilies <- function()
{
    printf("--- test_findRegionFamilies")
    scoreColumn <- "fc"
    scoreColumn <- "qvalScore"
    tbl.sub <- head(tbl, n=10)
    tcv <- TimeCourseVectors$new(tbl, scoreColumn)
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
                       tcv$getFamily(4))), seq_len(nrow(tbl)))

    browser()
    # igv <<- start.igv("KLF1")
    tcv$displayFamily(1)
    tcv$displayFamily(2)
    tcv$displayFamily(3)
    tcv$displayFamily(4)
    tcv$displayByDay()

    #
    #checkEquals(parents, c(1,2,3,5,6,8,10))
    browser()
    xyz <- 99


} # test_findRegionFamilies
#----------------------------------------------------------------------------------------------------
# small test to start: just 10 regions, 2 of them (4 and 7) are singletons
test_findRegionFamilies_200k <- function()
{
    printf("--- test_findRegionFamilies_200k")
    scoreColumn <- "fc"
    scoreColumn <- "qvalScore"
    tcv <- TimeCourseVectors$new(tbl2, scoreColumn)
    tcv$findRegionFamilies(minOverlap=100)
    #tbl.ov <- tcv$getOverlaps()
    checkEquals(tcv$getFamilyCount(), 31)
    for(fam in seq_len(31)){
       printf("    %d: %d", fam, length(tcv$getFamily(fam)))
       }



} # test_findRegionFamilies
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    #test_findRegionFamilies()
    test_findRegionFamilies_200k()

} # runTests
#----------------------------------------------------------------------------------------------------
# if(!interactive())
    runTests()
