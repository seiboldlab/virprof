#' Blacklist of problematic VHDB entries
vhdb_blacklist <- tribble(
    ~vhdb_virus_name, ~vhdb_host_name,
    "IAS virus", "Homo sapiens"
)


#' Load VirusHostDB
load_vhdb <- function(url = NULL, reload = FALSE) {
    global_vhdb <- mget(".virushostdb", env=globalenv(),
                        ifnotfound=list(NULL))[[1]]
    if (!reload && !is.null(global_vhdb)) {
        return (global_vhdb)
    }
    if (is.null(url)) {
        url <- "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv"
    }

    message("Loading VirusHostDB from '", url, "' ...")
    col_spec <- cols(
        `virus tax id` = col_integer(),    # NCBI taxid
        `virus name` = col_character(),
        `virus lineage` = col_character(),
        `refseq id` = col_character(),     # RefSeq Accession
        `KEGG GENOME` = col_character(),   # Dblink
        `KEGG DISEASE` = col_character(),
        DISEASE = col_character(),         # Disease Name
        `host tax id` = col_integer(),     # NCBI taxid
        `host name` = col_character(),
        `host lineage` = col_character(),
        pmid = col_character(),            # PubMed ID
        evidence = col_character(),        # source for host info
        ## if host tax id = 1, type of env sample
        `sample type` = col_character(),
        ## NCBI tax id for env sample (animal, plant, ...)
        `source organism` = col_integer()
    )
    vhdb <- read_tsv(url, col_types = col_spec)
    message("... processing")
    vhdb_filtered <- vhdb %>%
        rename_all(~ gsub(" ", "_", .)) %>%
        rename_all(tolower) %>%
        rename_all(~ paste0("vhdb_", .)) %>%
        mutate(
            vhdb_virus_lineage_name =
                paste(sep="; ",
                      vhdb_virus_lineage,
                      vhdb_virus_name)
        ) %>%
        anti_join(vhdb_blacklist,
                  by=c("vhdb_virus_name", "vhdb_host_name"))

    # Remove annoying spec attribute created by read_tsv
    attr(vhdb_filtered, "spec") <- NULL

    assign(".virushostdb", vhdb_filtered, env=globalenv())
    message("... done")
    vhdb_filtered
}


#' Compute list of hosts for each virus usign VHDB
vhdb_get_host <- function(lineages, taxids=NULL, url=NULL,
                          vhdb=load_vhdb(url=url)) {
    message("> Finding hosts using VirusHostDB for ", length(lineages), " hits")
    result <- tibble(
        lineage=lineages,
        host_names=""
    )

    ## If we have TaxIDs, check those first:
    if (!is.null(taxids)) {
        message(">  matching by NCBI tax ID")
        result$vhdb_virus_tax_id <- taxids
        result %<>%
            mutate(
                rowid=row_number()
            ) %>%
            left_join(select(vhdb, vhdb_virus_tax_id, vhdb_host_name),
                      by="vhdb_virus_tax_id") %>%
            group_by(rowid, lineage) %>%
            summarize(
                host_names = paste(collapse="; ", na.omit(vhdb_host_name)),
                .groups="drop"
            ) %>%
            ungroup() %>%
            select(lineage, host_names)
        message(">    found ", length(which(nzchar(result$host_names))))
    }


    matches <- tibble(lineage=character(), regex=character(), res=character())
    tosearch <- result %>%
        filter(!nzchar(host_names)) %>%
        group_by(lineage) %>%
        summarize(
            regex=paste0("^", gsub("([][$^()])", "\\\\\\1", unique(lineage))),
            .groups="drop"
        ) %>%
        ungroup()

    sum_linname <- function(hlin, hnam, vnam) {
        tibble(lineage=hlin, name=hnam) %>%
            filter(
                !is.na(lineage),
                !is.na(name)
            ) %>%
            mutate(
                domain=sub(";.*", "", lineage),
                kingdom=sub("^([^;]*;[^;]*);.*", "\\1", lineage),
                phylum=sub("^([^;]*;[^;]*;[^;]*);.*", "\\1", lineage)
            ) %>%
            group_by(domain) %>%
            mutate(n=n()/length(hlin)) %>%
            filter(n >= 0.1) %>%
            group_by(kingdom) %>%
            mutate(n=n()/length(hlin)) %>%
            filter(n >= 0.1) %>%
            group_by(phylum) %>%
            mutate(n=n()/length(hlin)) %>%
            filter(n >= 0.1) %>%
            pull(name) %>%
            unique() %>%
            sort() %>%
            paste(collapse=";")
    }

    for(i in seq(25)) {
        if (nrow(tosearch) == 0) {
            break;
        }
        if (i==1) {
            message("  matching by lineage prefix")
        }
        message(sprintf("    iteration %i: %i lineages left",
                        i, nrow(tosearch)))
        ## Match lineages
        found <- vhdb %>%
            regex_inner_join(
                tosearch,
                by=c(vhdb_virus_lineage_name="regex")
            ) %>%
            group_by(lineage, regex) %>%
            summarize(
                #res = paste(collapse=";", unique(vhdb_host_name))
                res = sum_linname(vhdb_host_lineage, vhdb_host_name, vhdb_virus_name),
                .groups="drop"
            ) %>%
            ungroup()
        matches <- rbind(matches, found)

        # Remove new matches from tosearch and drop one level lineage
        tosearch <- tosearch %>%
            anti_join(found, by="lineage") %>%  # remove found
            mutate(regex=sub(";[^;]*$", "", regex))  # drop lineage level

        # Re-use matches we already determined
        found <- matches %>%
            select(regex, res) %>%
            distinct() %>%
            inner_join(
                tosearch,
                by="regex"
            ) %>%
            select(lineage, regex, res)

        if (nrow(found) > 0) {
            message("      re-using ", nrow(found), " previous matches")
            matches <- rbind(matches, found)
            tosearch <- tosearch %>%
                anti_join(found, by="lineage")
        }
    }

    if (nrow(tosearch) > 0) {
        message(">  Failed to find ", nrow(tosearh), " lineages:")
        print(data)
    }

    result <- result %>%
        left_join(matches, by="lineage") %>%
        mutate(
            host_names = ifelse(nzchar(host_names), host_names, res)
        ) %>%
        pull(host_names)
    message("> ... done")
    result
}
