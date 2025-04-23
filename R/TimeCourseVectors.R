# TimeCourseVectors.R
#--------------------------------------------------------------------------------
#' @title TimeCourseVectors
#------------------------------------------------------------------------------------------------------------------------
#' @name TimeCourseVectors
#' @rdname TimeCourseVectors
#' @aliases TimeCourseVectors
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
                   scoreColumn=NULL
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
            }
       ) # public

    ) # class
#--------------------------------------------------------------------------------
