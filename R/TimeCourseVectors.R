# TimeCourseVectors.R
#--------------------------------------------------------------------------------
#' @title TimeCourseVectors
#------------------------------------------------------------------------------------------------------------------------
#' @name TimeCourseVectors
#' @rdname TimeCourseVectors
#' @aliases TimeCourseVectors
#----------------------------------------------------------------------------------------------------
#' @import GenomicRanges
#' @import igraph
#----------------------------------------------------------------------------------------------------
#' @description
#' An R6 class to transform experimental matrix data into timecourse vectors.
#'
#' @examples
#'   tcv <- TimeCourseVectors$new(id="abc")
#' @export

TimeCourseVectors = R6Class("TimeCourseVectors",

    #--------------------------------------------------------------------------------
    private = list(tbl=data.frame(),
                   scoreColumn=NULL,
                   tbl.ov=data.frame(),
                   xtab=NULL,
                   families=NULL,
                   parents=c(),
                   orphans=c()
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param tbl data.frame
         #' @param scoreColumn character
         #' @return a new instance of TimeCourseVectors
        initialize = function(tbl, scoreColumn){
            stopifnot(scoreColumn %in% colnames(tbl))
            stopifnot(all(c("chrom", "start", "end") %in% colnames(tbl)))
            private$tbl <- tbl
            private$scoreColumn <- scoreColumn
            },
        #------------------------------------------------------------
        #' @description accessor for the object's scoreColumn field
        #' @return the current value of the id member
        getScoreColumn = function(){
            private$scoreColumn
            },
        #------------------------------------------------------------
        #' @description accessor for the object's data.frame
        #' @return current base table of genomic regions
        getTable = function(){
            private$tbl
            },
        #------------------------------------------------------------
        #' @description find location-sharing genomic regions,
        #' each of which we regard as a family of regions, with
        #' the specified minimum overlap.  singleton families are
        #' included.
        #' @param miniOverlap integer, the number of shared bases
        #' @return the current value of the id member
        findRegionFamilies = function(minOverlap=20){
            ov <- findOverlaps(GRanges(private$tbl),
                               GRanges(private$tbl), type="within")
            tbl.ov <- as.data.frame(ov)
            tbl.ov <- subset(tbl.ov, queryHits != subjectHits)
            g <- graph_from_data_frame(tbl.ov)
               # produces 4 isolated subgraphs
               # extract them.  these are families. display in igv
               #
            all <- sort(unique(with(tbl.ov, c(queryHits, subjectHits))))
            private$orphans <- setdiff(seq_len(nrow(private$tbl)), all)
            private$tbl.ov <- tbl.ov
            xtab <- table(with(tbl.ov, c(queryHits, subjectHits)))
            tbl.dist <- as.data.frame(xtab)
            colnames(tbl.dist) <- c("row", "count")
            tbl.dist <- tbl.dist[rev(order(tbl.dist$count)),]
            rownames(tbl.dist) <- NULL
            membership <- components(g)$membership
            components <- sort(unique(as.integer(components(g)$membership)))
            clusters <- list()
            for(component in components){
               clusters[[component]] <- as.integer(names(membership[membership==component]))
               } # for component
            private$families <- clusters
            }, # findRegionFamilies

        #------------------------------------------------------------
        #' @description accessor for the object's GRanges overlap
        #' @return the overlap table as a data.frame
        getOverlaps = function(){
            private$tbl.ov
            },

        #------------------------------------------------------------
        #' @description accessor for the object's GRanges overlap
        #' @return the overlap table as a data.frame
        getFamilySizes = function(){
            as.integer(lapply(private$clusters, length))
            },


        #------------------------------------------------------------
        #' @description accessor for the parents of all the grouped regions
        #' @return the parents, a list
        getFamilyCount = function(){
            return(length(private$families))
            },

        #------------------------------------------------------------
        #' @description accessor for the parents of all the grouped regions
        #' @return the parents, a list
        getFamily = function(which){
            return(sort(unlist(private$families[which])))
            },

        #------------------------------------------------------------
        #' @description accessor for the parents of all the grouped regions
        #' @return the parents, a list
        getOrphans = function(){
            return(private$orphans)
            },

        #------------------------------------------------------------
        #' @description accessor for the parents of all the grouped regions
        #' @return the parents, a list
        displayFamily = function(which){
            rows <- unlist(private$families[which])
            maxVal <- 1.1 * max(tbl[, private$scoreColumn])
            tbl.fam <- private$tbl[rows,]
            for(i in seq_len(nrow(tbl.fam))){
               trackName <- sprintf("%d-%d-%d-%d", i, tbl$day[i], tbl$rep[i], which)
               tbl.track <- tbl.fam[i, c("chrom", "start", "end", private$scoreColumn)]
               track <- DataFrameQuantitativeTrack(trackName, tbl.track,
                                        color="black", autoscale=FALSE, min=0, max=maxVal)
               displayTrack(igv, track)
               } # for i
            }, # displayFamilies

        #------------------------------------------------------------
        #' @description use igv to display tracks, consolidated by day
        #' @return nothing
        displayByDay = function(){
            tbl <- private$tbl
            maxVal <- 1.1 * max(tbl[, private$scoreColumn])
            days <- sort(unique(tbl$day))
            for(Day in days){
               trackName <- sprintf("day %d", Day)
               tbl.track <- subset(tbl, day==Day)[, c("chrom", "start", "end", private$scoreColumn)]
               track <- DataFrameQuantitativeTrack(trackName, tbl.track,
                                        color="brown", autoscale=FALSE, min=0, max=maxVal)
               displayTrack(igv, track)
               } # for day
            } # displayByDay

        #------------------------------------------------------------
       ) # public

    ) # class
#--------------------------------------------------------------------------------
