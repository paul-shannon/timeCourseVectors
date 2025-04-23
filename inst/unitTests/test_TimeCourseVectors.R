library(RUnit)
library(TimeCourseVectors)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    f <- system.file(package="TimeCourseVectors", "extdata", "klf1-atac-10k.RData")
    tbl <- get(load(f))
    checkEquals(dim(tbl), c(28, 8))
    requiredColumns <- c("chrom", "start", "end", "rep", "day")
    checkTrue(all(requiredColumns %in% colnames(tbl)))
    scoreColumn <- "fc"
    tcv <- TimeCourseVectors$new(tbl, scoreColumn)
    checkTrue(all(c("R6", "TimeCourseVectors") %in% class(tcv)))

    checkEquals(tcv$getScoreColumn(), scoreColumn)

} # test_ctor
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
