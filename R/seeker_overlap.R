#' seeker overlap
#'
#' A generic function to search features that overlap a region defined by the given identifier
#'
#' @param id An Ensembl stable ID in string format
#'
#' @return
#' A list with the features: Band, Gene, Transcript, CDS, Exon, Structural variation, Regulatory, Variation, Motif, ChipSeq.
#'
#' @source
#' https://rest.ensembl.org
#'
#' @author
#' Erick Cuevas-Fernández
#'
#' Heriberto Manuel Rivera
#'
#' @importFrom
#' jsonlite fromJSON
#'
#' @importFrom
#' purrr transpose safely
#'
#' @importFrom
#' furrr future_map
#'
#' @importFrom
#' future plan multiprocess
#'
#' @examples
#' seeker_overlap("ENSG00000157764")
#'
#'
#' @rdname seeker_overlap
#' @export seeker_overlap
seeker_overlap <- function(id){
  server <- "https://rest.ensembl.org"
  if_id <- "/overlap/id/"
  searching <- c(Band = ";feature=band", Gene = ";feature=gene",
                 Transcript = ";feature=transcript",
                 CDS=";feature=cds", Exon =";feature=exon",
                 Structural_variation=";feature=structural_variation",
                 Regulatory=";feature=regulatory", Variation= ";feature=variation",
                 Motif=";feature=motif",
                 ChipSeq=";feature=chipseq")
  id_links <- paste0(server, if_id, id, "?content-type=application/json",
                    searching)



  future::plan("multiprocess")
  contents <- furrr::future_map(id_links, purrr::safely(jsonlite::fromJSON),
                                .progress = FALSE)
  contents_1 <- purrr::transpose(contents)
  contents_request<- contents_1[["result"]]
  names(contents_request) <- names(searching)

  return(contents_request)
  }
