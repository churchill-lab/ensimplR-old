#ENSIMPL_API_URL <- 'http://churchill-lab.jax.org/ensimpl/api'
ENSIMPL_API_URL <- 'http://127.0.0.1:8000/api'


#' Get the version information from ensimpl.
#'
#' @param debug Display some debugging information.
#' @export
#' @examples
#' versions()
versions <- function(debug=FALSE) {
    endpoint <- paste0(ENSIMPL_API_URL, '/versions')

    if (debug) {
        httr::set_config(httr::verbose())
    }

    r <- httr::GET(endpoint)

    httr::reset_config()

    jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))
}


#' Get gene information.
#'
#' @param id The Ensembl identifier.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param version will default to most recent ensimpl version.
#' @param extended TRUE for extra information.
#' @param debug Display some debugging information.
#' @export
#' @examples
#' gene <- getGene('ENSMUSG00000000001')
getGene <- function(id, species=NULL, version=NULL,
                    extended=FALSE, debug=FALSE) {
    endpoint <- paste0(ENSIMPL_API_URL, '/gene')

    params <- list(id = id)

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(version))) {
        params[['version']] <- version
    }

    if (extended) {
        params[['full']] <- 1
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    r <- httr::GET(endpoint, query = params)

    httr::reset_config()

    jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))
}


#' Get batch gene information.
#'
#' @param ids A list of Ensembl identifiers.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param version will default to most recent ensimpl version.
#' @param extended TRUE for extra information.
#' @param debug Display some debugging information.
#' @export
#' @examples
#' genes <- batchGenes(c('ENSMUSG00000000001', 'ENSMUSG00000000002'))
batchGenes <- function(ids, species=NULL, version=NULL,
                       extended=FALSE, debug=FALSE) {
    endpoint <- paste0(ENSIMPL_API_URL, '/genes')

    params <- list('ids[]' = ids)

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(version))) {
        params[['version']] <- version
    }

    if (extended) {
        params[['full']] <- 1
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    r <- httr::POST(endpoint, body = params, encode = 'json')

    httr::reset_config()

    jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))
}


#' Search for gene information.
#'
#' @param term The search term.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param version will default to most recent ensimpl version.
#' @param extended TRUE for extra information.
#' @param debug Display some debugging information.
#' @export
#' @examples
#' results <- searchGenes('Ace')
searchGenes <- function(term, species=NULL, version=NULL,
                        extended=FALSE, debug=FALSE) {
    endpoint <- paste0(ENSIMPL_API_URL, '/search')

    params <- list('term' = term)

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(version))) {
        params[['version']] <- version
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    r <- httr::GET(endpoint, query = params)

    httr::reset_config()

    jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'),
                       simplifyVector = TRUE)
}
