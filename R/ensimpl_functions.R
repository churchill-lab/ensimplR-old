ENSIMPL_API_URL <- 'http://churchill-lab.jax.org/ensimpl/api'
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

    return(temp$versions)
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
    resultsByIDs <- temp$ids

    ret <- data.frame(identifier = rep(NA, numIds),
                      gene_id = NA,
                      gene_id_version = NA,
                      symbol = NA,
                      name = NA,
                      chromosome = NA,
                      start = NA,
                      end = NA,
                      strand = NA,
                      synonyms = NA,
                      entrez_id = NA,
                      row.names = ids,
                      stringsAsFactors = FALSE)

    for (id in ids) {
        ret[id,]$identifier <- id
        if (id %in% names(resultsByIDs)) {
            ret[id,]$gene_id <- resultsByIDs[[id]]$id
            ret[id,]$gene_id_version <- resultsByIDs[[id]]$ensembl_version
            ret[id,]$symbol <- resultsByIDs[[id]]$symbol
            ret[id,]$name <- resultsByIDs[[id]]$name
            ret[id,]$chromosome <- resultsByIDs[[id]]$chromosome
            ret[id,]$start <- resultsByIDs[[id]]$start
            ret[id,]$end <- resultsByIDs[[id]]$end
            ret[id,]$strand <- resultsByIDs[[id]]$strand
            ret[id,]$synonyms <- paste(resultsByIDs[[id]]$synonyms, collapse = ':')

            extid_idx = which(resultsByIDs[[id]]$external_ids$db=='EntrezGene')
            if (length(extid_idx) != 0) {
                #print(resultsByIDs[[id]]$external_id)
                #print(resultsByIDs[[id]]$external_ids[which(resultsByIDs[[id]]$external_ids$db=='EntrezGene'),]$db_id)
                ret[id,]$entrez_id <- paste(resultsByIDs[[id]]$external_ids[extid_idx,]$db_id, collapse = ',')
            }
        }
    }

    if (!includeAll) {
        ret <- ret[!is.na(ret$gene_id),]
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
    tempIds <- read.csv(fileName, header=header)

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
