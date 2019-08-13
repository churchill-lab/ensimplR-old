library(dplyr)

ENSIMPL_API_URL <- 'https://churchilllab.jax.org/ensimpl/api'
#ENSIMPL_API_URL <- 'http://127.0.0.1:8000/api'


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

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))

    return(temp$versions %>% as_tibble)
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

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))

    return(temp$gene[[id]])
}


#' Get batch gene information.
#'
#' @param ids A list of Ensembl identifiers.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param version will default to most recent ensimpl version.
#' @param includeAll if TRUE, will give a one to one match
#' @param extended TRUE for extra information. NOT CURRENTLY USED
#' @param debug Display some debugging information.
#' @export
#' @examples
#' genes <- batchGenes(c('ENSMUSG00000000001', 'ENSMUSG00000000002'))
batchGenes <- function(ids, species=NULL, version=NULL, includeAll=TRUE,
                       extended=FALSE, debug=FALSE) {
    endpoint <- paste0(ENSIMPL_API_URL, '/genes')

    # Ensure ids are unique
    ids <- unique(ids)

    params <- list('ids[]' = ids)
    numIds <- length(ids)

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

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)

    ret <- temp$ids %>%
           purrr::map_dfr(~ purrr::map_df(.x, purrr::reduce, stringr::str_c, collapse = ",", sep= "|") ) %>%
           dplyr::rename(gene_id = id,
                         gene_id_version = ensembl_version) %>%
           dplyr::mutate(entrez_id = stringr::str_match(external_ids, '^.*EntrezGene\\|([a-zA-Z0-9]*)(,|$)')[,2],
                         mgi_id = stringr::str_match(external_ids, '^.*MGI\\|(MGI:[a-zA-Z0-9]*)(,|$)')[,2]) %>%
           tibble::add_column(identifier = names(temp$ids))

    if (includeAll) {
        all <- tibble(identifier = ids)
        ret <- dplyr::left_join(all, ret, by = c('identifier' = 'identifier'))
    }

    return(ret)
}




#' Get batch gene information from a file.
#'
#' @param fileName A file name with ensembl ids
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param version will default to most recent ensimpl version.
#' @param includeAll if TRUE, will give a one to one match
#' @param extended TRUE for extra information. NOT CURRENTLY USED
#' @param debug Display some debugging information.
#' @export
#' @examples
#' genes <- batchGenesFile()
batchGenesFile <- function(fileName, header=FALSE,
                           species=NULL, version=NULL, includeAll=TRUE,
                           extended=FALSE, debug=FALSE) {
    tempIds <- read.csv(fileName, header=header, stringsAsFactors = FALSE)
    names(tempIds) <- 'V1'

    return (batchGenes(tempIds$V1, species, version, includeAll,
                       extended, debug))

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
