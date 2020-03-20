library(dplyr)
library(httr)
library(jsonlite)
library(magrittr)
library(purrr)
library(rlang)
library(stringr)
library(tibble)
library(tictoc)

ENSIMPL_API_URL <- 'https://churchilllab.jax.org/ensimpl/api'
ENSIMPL_API_URL <- 'http://127.0.0.1:5000/api'


#' Get all the releases information from ensimpl.
#'
#' @importFrom tictoc tic toc tic.clear
#' @param debug Display some debugging information.
#' @return A tibble consisting of species, release and other columns.
#' @export
#' @examples
#' releases()
releases <- function(debug=FALSE) {
    tic('Releases')

    endpoint <- paste0(ENSIMPL_API_URL, '/releases')

    if (debug) {
        httr::set_config(httr::verbose())
    }

    tic('API')

    r <- httr::GET(endpoint)

    if (debug) {
        toc()
    }

    httr::reset_config()

    if (r$status_code == 404) {
        tic.clear()
        stop(paste0(endpoint, " not found"))
        return(NULL)
    } else if (r$status_code == 500) {
        tic.clear()
        temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)
        stop(temp$message)
        return(NULL)
    }

    tic('Parsing')

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))

    if (debug) {
        toc()
    }

    tic('Converting')

    ret <- temp$releases %>% dplyr::as_tibble()

    if (debug) {
        toc()
        toc()
    }

    tic.clear()

    return(ret)
}


#' Get gene information.
#'
#' @importFrom tictoc tic toc tic.clear
#' @param id The Ensembl identifier.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param release will default to most recent ensimpl release.
#' @param details TRUE for extra information including transcripts.
#' @param debug Display some debugging information.
#' @return A list with vgene information.
#' @export
#' @examples
#' gene <- getGene('ENSMUSG00000000001')
getGene <- function(id, species=NULL, release=NULL,
                    details=FALSE, debug=FALSE) {
    tic('getGene')

    endpoint <- paste0(ENSIMPL_API_URL, '/gene/', id)

    params <- list()

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(release))) {
        params[['release']] <- release
    }

    if (details) {
        params[['details']] <- 1
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    tic('API')

    r <- httr::GET(endpoint, query = params)

    if (debug) {
        toc()
    }

    httr::reset_config()

    if (r$status_code == 404) {
        tic.clear()
        stop(paste0(endpoint, " not found"))
        return(NULL)
    } else if (r$status_code == 500) {
        tic.clear()
        temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)
        stop(temp$message)
        return(NULL)
    }

    tic('Parsing')

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))

    if (debug) {
        toc()
    }

    tic('Converting')

    ret <- temp$gene[[id]]

    if (debug) {
        toc()
        toc()
    }

    tic.clear()

    return(ret)
}


#' Get batch gene information.
#'
#' @importFrom tictoc tic toc tic.clear
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @param ids A list of Ensembl identifiers.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param release will default to most recent ensimpl release.
#' @param details TRUE for extra information.
#' @param debug Display some debugging information.
#' @export
#' @return By default, with details=FALSE, a tibble will be returned.
#'   If details=TRUE, a list will be returned with each element being
#'   the matching "id".
#' @examples
#' genes <- batchGenes(c('ENSMUSG00000000001', 'ENSMUSG00000000002'))
batchGenes <- function(ids=NULL, species=NULL, release=NULL,
                       details=FALSE, debug=FALSE) {

    tic('batchGenes')

    endpoint <- paste0(ENSIMPL_API_URL, '/genes')

    params <- list()

    if (!(is.null(ids))) {
        # Ensure ids are unique
        ids <- unique(ids)
        if (length(ids) == 1) {
            ids <- list(ids)
        }
        params[['ids[]']] <- ids
    }

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(release))) {
        params[['release']] <- release
    }

    if (details) {
        params[['details']] <- 1
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    tic('API')

    r <- httr::POST(endpoint, body = params, encode = 'json')

    if (debug) {
        toc()
    }

    httr::reset_config()

    if (r$status_code == 404) {
        tic.clear()
        stop(paste0(endpoint, " not found"))
        return(NULL)
    } else if (r$status_code == 500) {
        tic.clear()
        temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)
        stop(temp$message)
        return(NULL)
    }

    tic('Parsing')

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)

    if (debug) {
        toc()
    }

    tic('Converting')

    # make sure we have some data
    if (is.null(temp$genes)) {
        return(NULL)
    }

    if (details) {
        ret <- temp$genes
    } else {
        ret <- temp$genes %>%
               purrr::map_dfr(~ purrr::map_df(.x, purrr::reduce, stringr::str_c, collapse = ",", sep= "|") ) %>%
               dplyr::rename(gene_id = .data$id,
                             gene_id_version = .data$ensembl_version) %>%
               tibble::add_column(id = names(temp$genes), .before = 1)

        #if (includeAll) {
        #    all <- tibble(identifier = ids)
        #    ret <- dplyr::left_join(all, ret, by = c('identifier' = 'identifier'))
        #}
    }

    if (debug) {
        toc()
        toc()
    }

    tic.clear()

    return(ret)
}


#' Get batch gene information from a file.
#'
#' @importFrom tictoc tic toc tic.clear
#' @param fileName A file name with ensembl ids
#' @param header TRUE if the file has a header line.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param release will default to most recent ensimpl release.
#' @param details TRUE for extra information.
#' @param debug Display some debugging information.
#' @return By default, with details=FALSE, a tibble will be returned.
#'   If details=TRUE, a list will be returned with each element being
#'   the matching "id".
#' @export
batchGenesFile <- function(fileName, header=FALSE,
                           species=NULL, release=NULL,
                           details=FALSE, debug=FALSE) {
    tic('batchGenesFile')

    tempIds <- utils::read.csv(fileName,
                               header=header,
                               stringsAsFactors = FALSE)
    names(tempIds) <- 'V1'

    if (debug) {
        toc()
    }

    tic.clear()

    return (batchGenes(tempIds$V1, species, release, details, debug))
}


#' Search for gene information.
#'
#' @importFrom tictoc tic toc tic.clear
#' @param term The search term.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param release will default to most recent ensimpl release.
#' @param debug Display some debugging information.
#' @export
#' @examples
#' results <- searchGenes('Ace')
searchGenes <- function(term, species=NULL, release=NULL, debug=FALSE) {

    tic('searchGenes')

    endpoint <- paste0(ENSIMPL_API_URL, '/search')

    params <- list('term' = term)

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(release))) {
        params[['release']] <- release
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    tic('API')

    r <- httr::GET(endpoint, query = params)

    if (debug) {
        toc()
    }

    httr::reset_config()

    if (r$status_code == 404) {
        tic.clear()
        stop(paste0(endpoint, " not found"))
        return(NULL)
    } else if (r$status_code == 500) {
        tic.clear()
        temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)
        stop(temp$message)
        return(NULL)
    }

    tic('Parsing')

    ret <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyVector = TRUE)

    if (debug) {
        toc()
        toc()
    }

    tic.clear()

    return(ret)
}


#' Lookup by ids.
#'
#' @importFrom tictoc tic toc tic.clear
#' @importFrom magrittr %>%
#' @param ids A list of Ensembl identifiers.
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param release will default to most recent ensimpl release.
#' @param source_db A string which defaults to 'Ensembl', but could also
#'   be one of the source_db values from externalDBs() function.
#' @param debug Display some debugging information.
#' @export
#' @return A tibble with the external ids listed.
#' @examples
#' extids <- externalIDs(c('ENSMUSG00000000001', 'ENSMUSG00000000002'))
externalIDs <- function(ids=NULL, species=NULL, release=NULL,
                        source_db='Ensembl', debug=FALSE) {

    tic('externalIDs')

    endpoint <- paste0(ENSIMPL_API_URL, '/external_ids')

    params <- list()

    if (!(is.null(ids))) {
        # Ensure ids are unique
        ids <- unique(ids)
        params[['ids[]']] <- ids
    }

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(release))) {
        params[['release']] <- release
    }

    if (!(is.null(source_db))) {
        params[['source_db']] <- source_db
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    tic('API')

    r <- httr::POST(endpoint, body = params, encode = 'json')

    if (debug) {
        toc()
    }

    httr::reset_config()

    if (r$status_code == 404) {
        tic.clear()
        stop(paste0(endpoint, " not found"))
        return(NULL)
    } else if (r$status_code == 500) {
        tic.clear()
        temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)
        stop(temp$message)
        return(NULL)
    }

    tic('Parsing')

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)

    if (debug) {
        toc()
    }

    tic('Converting')

    if (is.null(temp$ids)) {
        ret <- NULL
    } else {
        ret <- temp$ids
        # To convert to a nice tibble:
        # temp$ids %>%
        #    purrr::map_dfr(~ purrr::map_df(.x, purrr::reduce, stringr::str_c, collapse = ",", sep= "|")) %>%
        #    dplyr::as_tibble() %>%
        #    tibble::add_column(id = names(temp$ids), .before = 1)
    }

    # If the same columns are needed
    #dbs <- c('EntrezGene', 'Uniprot_gn', 'MGI', 'Ensembl_homolog', 'Ensembl')
    #
    #for (db_id in dbs) {
    #    if (db_id %not in% colnames(s) ) {
    #        ret <- ret %>%
    #                  tibble::add_column(!!(db_id) := NA)
    #    }
    #}

    if (debug) {
        toc()
        toc()
    }

    tic.clear()

    return(ret)
}


#' Get the external databases.
#'
#' @importFrom tictoc tic toc tic.clear
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @param species Either 'Mm' or 'Hs', will default to 'Mm'.
#' @param release will default to most recent ensimpl release.
#' @param debug Display some debugging information.
#' @export
#' @return A tibble of external database identifiers.
#' @examples
#' dbs <- externalDBs()
externalDBs <- function(species=NULL, release=NULL, debug=FALSE) {

    tic('externalDBs')

    endpoint <- paste0(ENSIMPL_API_URL, '/external_dbs')

    params <- list()

    if (!(is.null(species))) {
        params[['species']] <- species
    }

    if (!(is.null(release))) {
        params[['release']] <- release
    }

    if (debug) {
        httr::set_config(httr::verbose())
    }

    tic('API')

    r <- httr::GET(endpoint, query = params, encode = 'json')

    if (debug) {
        toc()
    }

    httr::reset_config()

    if (r$status_code == 404) {
        tic.clear()
        stop(paste0(endpoint, " not found"))
        return(NULL)
    } else if (r$status_code == 500) {
        tic.clear()
        temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'), simplifyDataFrame = TRUE)
        stop(temp$message)
        return(NULL)
    }

    tic('Parsing')

    temp <- jsonlite::fromJSON(httr::content(r, type = 'text', encoding = 'UTF-8'))

    if (debug) {
        toc()
    }

    tic('Converting')

    ret <- temp$external_dbs %>%
        dplyr::as_tibble() %>%
        dplyr::select(source_db = .data$external_db_id,
                      source_name = .data$external_db_name)

    if (debug) {
        toc()
        toc()
    }

    tic.clear()

    return(ret)
}
