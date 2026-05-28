#!/usr/bin/env Rscript

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
analysis_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "analysis_runtime.R"))
paths <- pic_runtime_paths(analysis_root)

enrichment_out_dir <- paths$enrichment_lib_dir
register_dir <- paths$register_dir
genome_runtime_map_file <- paths$genome_map_file
genome_source_map_file <- file.path(register_dir, "genome_map.tsv")

enrichment_resources <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "WIKIPATHWAYS", "REACTOME", "HDO", "HPO", "MPO")

pick_orgdb_for_resource <- function(resource, species_row) {
  resource <- toupper(as.character(resource))
  if (resource %in% c("HDO", "HPO")) return("org.Hs.eg.db")
  if (resource == "MPO") return("org.Mm.eg.db")
  as.character(species_row$orgdb)
}

normalize_gson_to_symbol <- function(gson_obj, orgdb_pkg) {
  if (!safe_package_available("AnnotationDbi")) {
    return(list(ok = FALSE, reason = "AnnotationDbi not available", gson = NULL))
  }
  if (!safe_package_available(orgdb_pkg)) {
    return(list(ok = FALSE, reason = sprintf("%s not available", orgdb_pkg), gson = NULL))
  }

  gsid2gene <- tryCatch(tibble::as_tibble(gson_obj@gsid2gene), error = function(e) NULL)
  gsid2name <- tryCatch(tibble::as_tibble(gson_obj@gsid2name), error = function(e) NULL)
  if (is.null(gsid2gene) || ncol(gsid2gene) < 2 || !all(c("gsid", "gene") %in% colnames(gsid2gene))) {
    return(list(ok = FALSE, reason = "invalid gsid2gene format", gson = NULL))
  }
  if (is.null(gsid2name) || ncol(gsid2name) < 2 || !all(c("gsid", "name") %in% colnames(gsid2name))) {
    return(list(ok = FALSE, reason = "invalid gsid2name format", gson = NULL))
  }

  genes <- as.character(gsid2gene$gene)
  genes <- trimws(genes)
  genes[is.na(genes)] <- ""
  keytype_raw <- tryCatch(as.character(gson_obj@keytype), error = function(e) character())
  if (length(keytype_raw) == 0 || is.na(keytype_raw[[1]]) || trimws(keytype_raw[[1]]) == "") {
    keytype_raw <- "UNKNOWN"
  } else {
    keytype_raw <- toupper(trimws(keytype_raw[[1]]))
  }
  orgdb_obj <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))

  map_by <- function(keys, keytype) {
    keys <- as.character(keys)
    keys <- trimws(keys)
    keys <- keys[keys != ""]
    if (length(keys) == 0) return(character())
    out <- suppressWarnings(
      suppressMessages(
        AnnotationDbi::mapIds(
          x = orgdb_obj,
          keys = unique(keys),
          column = "SYMBOL",
          keytype = keytype,
          multiVals = "first"
        )
      )
    )
    out
  }

  mapped_symbol <- rep(NA_character_, length(genes))
  if (keytype_raw == "SYMBOL") {
    mapped_symbol <- genes
  } else if (keytype_raw == "ENSEMBL") {
    m <- map_by(genes, "ENSEMBL")
    mapped_symbol <- unname(m[genes])
  } else {
    # ENTREZ/kegg/unknown -> try ENTREZID first.
    genes_entrez <- sub("^[^:]+:", "", genes)
    m_entrez <- map_by(genes_entrez, "ENTREZID")
    mapped_symbol <- unname(m_entrez[genes_entrez])
    # Fallback: some sources may already be SYMBOL.
    na_idx <- which(is.na(mapped_symbol) | mapped_symbol == "")
    if (length(na_idx) > 0) {
      mapped_symbol[na_idx] <- genes[na_idx]
    }
  }

  mapped_symbol <- trimws(as.character(mapped_symbol))
  mapped_symbol[mapped_symbol == ""] <- NA_character_

  gsid2gene_sym <- gsid2gene |>
    dplyr::mutate(gene = mapped_symbol) |>
    dplyr::filter(!is.na(.data$gene), .data$gene != "") |>
    dplyr::distinct(.data$gsid, .data$gene)
  if (nrow(gsid2gene_sym) == 0) {
    return(list(ok = FALSE, reason = "failed to map genes to SYMBOL", gson = NULL))
  }

  as_scalar_chr <- function(x, default = "") {
    x <- as.character(x)
    if (length(x) == 0 || is.na(x[[1]]) || trimws(x[[1]]) == "") return(default)
    trimws(x[[1]])
  }

  species_val <- as_scalar_chr(tryCatch(gson_obj@species, error = function(e) ""), default = "unknown")
  gsname_val <- as_scalar_chr(tryCatch(gson_obj@gsname, error = function(e) ""), default = "unknown")
  version_val <- as_scalar_chr(tryCatch(gson_obj@version, error = function(e) ""), default = "unknown")
  accessed_val <- as_scalar_chr(
    tryCatch(gson_obj@accessed_date, error = function(e) ""),
    default = format(Sys.time(), "%Y-%m-%d", tz = "Asia/Tokyo")
  )

  gson_sym <- gson::gson(
    gsid2gene = gsid2gene_sym,
    gsid2name = gsid2name,
    species = species_val,
    gsname = gsname_val,
    keytype = "SYMBOL",
    version = version_val,
    accessed_date = accessed_val
  )

  list(ok = TRUE, reason = "", gson = gson_sym)
}

parse_dataset_prefix <- function(dataset) {
  x <- tolower(trimws(as.character(dataset)))
  sub("_gene_ensembl$", "", x)
}

normalize_species_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

orgdb_species_table <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    if (!safe_package_available("AnnotationDbi")) {
      cache <<- tibble::tibble(orgdb = character(), species_label = character(), species_key = character(), ensembl_prefix = character())
      return(cache)
    }
    pkgs <- rownames(installed.packages())
    pkgs <- pkgs[grepl("^org\\..+\\.eg\\.db$", pkgs)]
    rows <- lapply(pkgs, function(pkg) {
      if (!safe_package_available(pkg)) return(NULL)
      obj <- tryCatch(get(pkg, envir = asNamespace(pkg)), error = function(e) NULL)
      if (is.null(obj)) return(NULL)
      sp <- tryCatch(as.character(AnnotationDbi::species(obj)), error = function(e) NA_character_)
      sp <- trimws(sp)
      if (is.na(sp) || sp == "") return(NULL)
      parts <- strsplit(tolower(sp), "\\s+")[[1]]
      ens_prefix <- if (length(parts) >= 2) paste0(substr(parts[[1]], 1, 1), parts[[2]]) else ""
      tibble::tibble(
        orgdb = pkg,
        species_label = sp,
        species_key = normalize_species_key(gsub("\\s+", "_", sp)),
        ensembl_prefix = ens_prefix
      )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) {
      cache <<- tibble::tibble(
        orgdb = character(),
        species_label = character(),
        species_key = character(),
        ensembl_prefix = character()
      )
      return(cache)
    }
    cache <<- dplyr::bind_rows(rows) |>
      dplyr::distinct(.data$species_key, .keep_all = TRUE)
    cache
  }
})

species_key_to_kegg <- function(species_key) {
  switch(
    species_key,
    homo_sapiens = "hsa",
    mus_musculus = "mmu",
    rattus_norvegicus = "rno",
    danio_rerio = "dre",
    drosophila_melanogaster = "dme",
    caenorhabditis_elegans = "cel",
    saccharomyces_cerevisiae = "sce",
    gallus_gallus = "gga",
    canis_familiaris = "cfa",
    sus_scrofa = "ssc",
    macaca_mulatta = "mcc",
    ""
  )
}

species_key_to_reactome <- function(species_key) {
  switch(
    species_key,
    homo_sapiens = "human",
    mus_musculus = "mouse",
    rattus_norvegicus = "rat",
    drosophila_melanogaster = "fly",
    danio_rerio = "zebrafish",
    ""
  )
}

species_meta_from_key <- function(species_key) {
  sk <- normalize_species_key(species_key)
  tbl <- orgdb_species_table()
  hit <- tbl |> dplyr::filter(.data$species_key == .env$sk) |> dplyr::slice_head(n = 1)
  if (nrow(hit) == 0) {
    return(list(species_key = sk, species_label = gsub("_", " ", sk, fixed = TRUE), orgdb = NA_character_, kegg = "", reactome = "", wp = ""))
  }
  list(
    species_key = as.character(hit$species_key[[1]]),
    species_label = as.character(hit$species_label[[1]]),
    orgdb = as.character(hit$orgdb[[1]]),
    kegg = species_key_to_kegg(as.character(hit$species_key[[1]])),
    reactome = species_key_to_reactome(as.character(hit$species_key[[1]])),
    wp = as.character(hit$species_label[[1]])
  )
}

dataset_to_species_meta <- function(dataset) {
  prefix <- parse_dataset_prefix(dataset)
  if (is.na(prefix) || prefix == "") return(species_meta_from_key(""))
  tbl <- orgdb_species_table()
  hit <- tbl |>
    dplyr::filter(
      .data$ensembl_prefix == .env$prefix |
        .data$species_key == .env$prefix |
        gsub("_", "", .data$species_key, fixed = TRUE) == .env$prefix
    ) |>
    dplyr::slice_head(n = 1)
  if (nrow(hit) == 0) return(species_meta_from_key(prefix))
  species_meta_from_key(as.character(hit$species_key[[1]]))
}

genome_to_species_meta <- function(genome) {
  g <- tolower(trimws(as.character(genome)))
  alias <- c(
    human = "homo_sapiens", hg18 = "homo_sapiens", hg19 = "homo_sapiens", hg38 = "homo_sapiens", grch36 = "homo_sapiens", grch37 = "homo_sapiens", grch38 = "homo_sapiens",
    mouse = "mus_musculus", mm9 = "mus_musculus", mm10 = "mus_musculus", mm39 = "mus_musculus", grcm37 = "mus_musculus", grcm38 = "mus_musculus", grcm39 = "mus_musculus",
    rat = "rattus_norvegicus", rn4 = "rattus_norvegicus", rn5 = "rattus_norvegicus", rn6 = "rattus_norvegicus", rn7 = "rattus_norvegicus",
    fly = "drosophila_melanogaster", dm3 = "drosophila_melanogaster", dm6 = "drosophila_melanogaster",
    worm = "caenorhabditis_elegans", ce10 = "caenorhabditis_elegans", ce11 = "caenorhabditis_elegans",
    zebrafish = "danio_rerio", danrer10 = "danio_rerio", danrer11 = "danio_rerio",
    yeast = "saccharomyces_cerevisiae", saccer3 = "saccharomyces_cerevisiae",
    chicken = "gallus_gallus", galgal5 = "gallus_gallus", galgal6 = "gallus_gallus",
    dog = "canis_familiaris", canfam3 = "canis_familiaris",
    pig = "sus_scrofa", susscr11 = "sus_scrofa",
    rhesus = "macaca_mulatta", rhemac10 = "macaca_mulatta"
  )
  if (g %in% names(alias)) return(species_meta_from_key(unname(alias[[g]])))
  species_meta_from_key("")
}

safe_package_available <- function(pkg) {
  !is.na(pkg) && pkg != "" && suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE))
}

safe_save_gson <- function(obj, out_file) {
  saveRDS(obj, out_file, compress = "xz")
  out_file
}

quiet_try <- function(expr) {
  suppressWarnings(
    suppressMessages(
      tryCatch(expr, error = function(e) structure(list(message = conditionMessage(e)), class = "pic_error"))
    )
  )
}

is_pic_error <- function(x) inherits(x, "pic_error")

build_species_table <- function(target_genome = NULL) {
  genomes <- suppressMessages(
    readr::read_tsv(
      genome_source_map_file,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )
  )
  lookup_file <- file.path(register_dir, "biomart_lookup.tsv")
  if (!file.exists(lookup_file)) {
    stop(sprintf("biomart lookup was not found: %s", lookup_file), call. = FALSE)
  }
  lookup <- suppressMessages(
    readr::read_tsv(
      lookup_file,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )
  )

  if (!is.null(target_genome) && target_genome != "") {
    genomes <- genomes |>
      dplyr::filter(.data$genome == .env$target_genome)
  }

  dplyr::left_join(genomes, lookup, by = "genome") |>
    dplyr::mutate(
      meta_dataset = purrr::map(.data$biomart_dataset, dataset_to_species_meta),
      meta_genome = purrr::map(.data$genome, genome_to_species_meta),
      orgdb_dataset = purrr::map_chr(.data$meta_dataset, "orgdb"),
      meta = purrr::pmap(
        list(.data$meta_dataset, .data$meta_genome, .data$orgdb_dataset),
        function(m_dataset, m_genome, org_dataset) {
          if (!is.na(org_dataset) && org_dataset != "") m_dataset else m_genome
        }
      )
    ) |>
    dplyr::transmute(
      genome = as.character(.data$genome),
      species_key = purrr::map_chr(.data$meta, "species_key"),
      species_label = purrr::map_chr(.data$meta, "species_label"),
      orgdb = purrr::map_chr(.data$meta, "orgdb"),
      kegg_code = purrr::map_chr(.data$meta, "kegg"),
      reactome_species = purrr::map_chr(.data$meta, "reactome"),
      wp_species = purrr::map_chr(.data$meta, "wp")
    ) |>
    dplyr::distinct(.data$genome, .keep_all = TRUE)
}

build_human_disease_gson <- function(resource) {
  if (!safe_package_available("DOSE")) return(list(ok = FALSE, reason = "DOSE not available", gson = NULL))
  if (!safe_package_available("gson")) return(list(ok = FALSE, reason = "gson not available", gson = NULL))
  if (!safe_package_available("GOSemSim")) return(list(ok = FALSE, reason = "GOSemSim not available", gson = NULL))
  if (!safe_package_available("yulab.utils")) return(list(ok = FALSE, reason = "yulab.utils not available", gson = NULL))
  if (!safe_package_available("R.utils")) return(list(ok = FALSE, reason = "R.utils not available", gson = NULL))

  ontology <- toupper(as.character(resource))
  if (!ontology %in% c("HDO", "HPO", "MPO")) {
    return(list(ok = FALSE, reason = "unknown disease resource", gson = NULL))
  }

  if (ontology == "MPO") {
    run_curl_capture_retry <- function(url, retries = 4L, sleep_sec = 3) {
      curl_bin <- Sys.which("curl")
      if (identical(curl_bin, "")) stop("curl is not available", call. = FALSE)
      last <- NULL
      for (i in seq_len(retries)) {
        out <- tryCatch(
          suppressWarnings(system2(curl_bin, c("-L", "--silent", "--show-error", url), stdout = TRUE, stderr = TRUE)),
          error = function(e) NULL
        )
        rc <- if (is.null(out)) 1L else attr(out, "status")
        if (is.null(rc) || rc == 0) return(paste(out, collapse = "\n"))
        last <- out
        if (i < retries) Sys.sleep(sleep_sec)
      }
      stop(paste(last, collapse = "\n"), call. = FALSE)
    }

    parse_mp_obo <- function(txt) {
      lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
      ids <- character()
      names <- character()
      cur_id <- NULL
      cur_name <- NULL
      in_term <- FALSE
      flush_term <- function() {
        if (!is.null(cur_id) && !is.null(cur_name) && nzchar(cur_id) && nzchar(cur_name)) {
          ids <<- c(ids, cur_id)
          names <<- c(names, cur_name)
        }
      }
      for (ln in lines) {
        if (identical(trimws(ln), "[Term]")) {
          if (in_term) flush_term()
          in_term <- TRUE
          cur_id <- NULL
          cur_name <- NULL
          next
        }
        if (!in_term) next
        if (startsWith(ln, "id: MP:")) {
          cur_id <- trimws(sub("^id:\\s*", "", ln))
        } else if (startsWith(ln, "name:")) {
          cur_name <- trimws(sub("^name:\\s*", "", ln))
        } else if (trimws(ln) == "") {
          flush_term()
          in_term <- FALSE
          cur_id <- NULL
          cur_name <- NULL
        }
      }
      if (in_term) flush_term()
      tibble::tibble(gsid = ids, name = names) |>
        dplyr::filter(!is.na(.data$gsid), .data$gsid != "", !is.na(.data$name), .data$name != "") |>
        dplyr::distinct(.data$gsid, .data$name)
    }

    mp_obo_urls <- c(
      "https://purl.obolibrary.org/obo/mp.obo",
      "https://raw.githubusercontent.com/obophenotype/mammalian-phenotype-ontology/master/mp.obo"
    )
    mgi_urls <- c(
      "https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt",
      "http://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt"
    )
    mgi_coord_urls <- c(
      "https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt",
      "http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt"
    )

    mp_obo_text <- NULL
    for (u in mp_obo_urls) {
      x <- quiet_try(run_curl_capture_retry(u))
      if (!is_pic_error(x) && is.character(x) && nchar(x) > 0) {
        mp_obo_text <- x
        break
      }
    }
    if (is.null(mp_obo_text)) {
      return(list(ok = FALSE, reason = "failed to download mp.obo", gson = NULL))
    }

    mgi_text <- NULL
    for (u in mgi_urls) {
      x <- quiet_try(run_curl_capture_retry(u))
      if (!is_pic_error(x) && is.character(x) && nchar(x) > 0) {
        mgi_text <- x
        break
      }
    }
    if (is.null(mgi_text)) {
      return(list(ok = FALSE, reason = "failed to download MGI_GenePheno.rpt", gson = NULL))
    }
    mgi_coord_text <- NULL
    for (u in mgi_coord_urls) {
      x <- quiet_try(run_curl_capture_retry(u))
      if (!is_pic_error(x) && is.character(x) && nchar(x) > 0) {
        mgi_coord_text <- x
        break
      }
    }
    if (is.null(mgi_coord_text)) {
      return(list(ok = FALSE, reason = "failed to download MGI_Gene_Model_Coord.rpt", gson = NULL))
    }
    gsid2name <- parse_mp_obo(mp_obo_text)
    if (nrow(gsid2name) == 0) {
      return(list(ok = FALSE, reason = "mp.obo parse failed", gson = NULL))
    }

    parse_mgi_pheno <- function(txt, header = TRUE) {
      con <- textConnection(txt)
      on.exit(close(con), add = TRUE)
      base_try <- quiet_try(
        utils::read.delim(
          con,
          header = header,
          sep = "\t",
          quote = "",
          comment.char = "",
          fill = TRUE,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      )
      if (!is_pic_error(base_try) && !is.null(base_try) && nrow(base_try) > 0) {
        return(tibble::as_tibble(base_try, .name_repair = "minimal"))
      }

      readr_try <- quiet_try(
        readr::read_tsv(
          I(txt),
          col_names = header,
          col_types = readr::cols(.default = readr::col_character()),
          progress = FALSE,
          show_col_types = FALSE
        )
      )
      if (!is_pic_error(readr_try) && !is.null(readr_try) && nrow(readr_try) > 0) {
        return(tibble::as_tibble(readr_try, .name_repair = "minimal"))
      }
      NULL
    }

    mgi_tbl <- parse_mgi_pheno(mgi_text, header = FALSE)
    if (is.null(mgi_tbl) || nrow(mgi_tbl) == 0) {
      return(list(ok = FALSE, reason = "failed to parse MGI_GenePheno.rpt", gson = NULL))
    }

    # MGI_GenePheno.rpt sample:
    # col5: MP term id, col7: marker MGI accession id
    if (ncol(mgi_tbl) < 7) {
      return(list(ok = FALSE, reason = "MGI_GenePheno.rpt required columns are missing", gson = NULL))
    }
    colnames(mgi_tbl)[1:7] <- c("genotype", "allele_symbol", "allele_mgi", "background", "gsid", "pmid", "mgi")
    mpo2mgi <- mgi_tbl |>
      dplyr::transmute(
        gsid = as.character(.data$gsid),
        mgi = as.character(.data$mgi)
      ) |>
      dplyr::filter(grepl("^MP:\\d+$", .data$gsid), !is.na(.data$mgi), .data$mgi != "") |>
      dplyr::distinct(.data$gsid, .data$mgi)

    if (nrow(mpo2mgi) == 0) {
      return(list(ok = FALSE, reason = "MPO mpo2mgi is empty", gson = NULL))
    }

    coord_tbl <- parse_mgi_pheno(mgi_coord_text, header = FALSE)
    if (is.null(coord_tbl) || nrow(coord_tbl) == 0) {
      return(list(ok = FALSE, reason = "failed to parse MGI_Gene_Model_Coord.rpt", gson = NULL))
    }
    if (ncol(coord_tbl) < 3) {
      return(list(ok = FALSE, reason = "MGI_Gene_Model_Coord.rpt required columns are missing", gson = NULL))
    }
    colnames(coord_tbl)[1:3] <- c("mgi", "chr", "gene")
    mgi2symbol <- coord_tbl |>
      dplyr::transmute(mgi = as.character(.data$mgi), gene = as.character(.data$gene)) |>
      dplyr::filter(!is.na(.data$mgi), .data$mgi != "", !is.na(.data$gene), .data$gene != "") |>
      dplyr::distinct(.data$mgi, .data$gene)

    gsid2gene <- dplyr::inner_join(mpo2mgi, mgi2symbol, by = "mgi") |>
      dplyr::select("gsid", "gene") |>
      dplyr::distinct(.data$gsid, .data$gene)

    if (nrow(gsid2gene) == 0) {
      return(list(ok = FALSE, reason = "MPO TERM2GENE is empty after mpo2mgi+mgi2symbol merge", gson = NULL))
    }

    gsid2name <- gsid2name |>
      dplyr::filter(.data$gsid %in% unique(gsid2gene$gsid))
    gsid2gene <- gsid2gene |>
      dplyr::filter(.data$gsid %in% unique(gsid2name$gsid))
    if (nrow(gsid2name) == 0 || nrow(gsid2gene) == 0) {
      return(list(ok = FALSE, reason = "MPO TERM2GENE/TERM2NAME overlap is empty", gson = NULL))
    }

    obj <- quiet_try(gson::gson(
      gsid2gene = gsid2gene,
      gsid2name = gsid2name,
      species = "Mus musculus",
      gsname = "MPO",
      keytype = "SYMBOL",
      version = "unknown",
      accessed_date = as.character(Sys.Date())
    ))
    if (is_pic_error(obj)) {
      return(list(ok = FALSE, reason = obj$message, gson = NULL))
    }
    return(list(ok = TRUE, reason = "", gson = obj))
  }

  ont2gene <- quiet_try(DOSE:::get_gene2ont(ontology, output = "data.frame"))
  if (is_pic_error(ont2gene)) {
    return(list(ok = FALSE, reason = ont2gene$message, gson = NULL))
  }
  if (!is.data.frame(ont2gene) || ncol(ont2gene) < 2) {
    return(list(ok = FALSE, reason = sprintf("invalid %s ont2gene format", ontology), gson = NULL))
  }

  gsid2gene <- tibble::as_tibble(ont2gene) |>
    dplyr::transmute(
      gene = as.character(.data[[colnames(ont2gene)[1]]]),
      gsid = as.character(.data[[colnames(ont2gene)[2]]])
    ) |>
    dplyr::select(gsid, gene) |>
    dplyr::filter(!is.na(.data$gsid), .data$gsid != "", !is.na(.data$gene), .data$gene != "") |>
    dplyr::distinct(.data$gsid, .data$gene)

  if (nrow(gsid2gene) == 0) {
    return(list(ok = FALSE, reason = sprintf("%s ont2gene is empty", ontology), gson = NULL))
  }

  termmap <- quiet_try(GOSemSim:::get_onto_data(ontology, table = "term", output = "data.frame"))
  if (is_pic_error(termmap)) {
    return(list(ok = FALSE, reason = termmap$message, gson = NULL))
  }
  if (!is.data.frame(termmap) || ncol(termmap) < 2) {
    return(list(ok = FALSE, reason = sprintf("invalid %s term map format", ontology), gson = NULL))
  }

  gsid2name <- tibble::as_tibble(termmap) |>
    dplyr::transmute(
      gsid = as.character(.data[[colnames(termmap)[1]]]),
      name = as.character(.data[[colnames(termmap)[2]]])
    ) |>
    dplyr::filter(!is.na(.data$gsid), .data$gsid != "", !is.na(.data$name), .data$name != "") |>
    dplyr::distinct(.data$gsid, .data$name) |>
    dplyr::filter(.data$gsid %in% unique(gsid2gene$gsid))

  # Safety: some ontology backends can return TERM2GENE IDs that are not present
  # in the term map (or vice versa). Keep only resolvable IDs.
  gsid2gene <- gsid2gene |>
    dplyr::filter(.data$gsid %in% unique(gsid2name$gsid))

  if (nrow(gsid2name) == 0) {
    return(list(ok = FALSE, reason = sprintf("%s term map is empty", ontology), gson = NULL))
  }

  species_label <- if (resource == "MPO") "Mus musculus" else "Homo sapiens"
  obj <- quiet_try(gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = species_label,
    gsname = ontology,
    keytype = "ENTREZID",
    version = "unknown",
    accessed_date = as.character(Sys.Date())
  ))
  if (is_pic_error(obj)) {
    return(list(ok = FALSE, reason = obj$message, gson = NULL))
  }
  list(ok = TRUE, reason = "", gson = obj)
}

build_gson_for_species <- function(resource, species_row) {
  resource <- toupper(as.character(resource))

  if (resource %in% c("HDO", "HPO")) {
    if (!identical(as.character(species_row$species_key), "homo_sapiens")) {
      return(list(ok = FALSE, reason = "human-only resource", gson = NULL))
    }
    return(build_human_disease_gson(resource))
  }

  if (resource == "MPO") {
    if (!identical(as.character(species_row$species_key), "mus_musculus")) {
      return(list(ok = FALSE, reason = "mouse-only resource", gson = NULL))
    }
    return(build_human_disease_gson(resource))
  }

  if (!safe_package_available("clusterProfiler")) {
    return(list(ok = FALSE, reason = "clusterProfiler not available", gson = NULL))
  }

  if (resource %in% c("GO_BP", "GO_CC", "GO_MF")) {
    if (!safe_package_available(species_row$orgdb)) {
      return(list(ok = FALSE, reason = sprintf("%s not available", as.character(species_row$orgdb)), gson = NULL))
    }
    ont <- sub("^GO_", "", resource)
    out <- quiet_try(clusterProfiler::gson_GO(OrgDb = species_row$orgdb, keytype = "SYMBOL", ont = ont))
    if (is_pic_error(out)) {
      out <- quiet_try(clusterProfiler::gson_GO(OrgDb = species_row$orgdb, keyType = "SYMBOL", ont = ont))
    }
    if (is_pic_error(out)) {
      out <- quiet_try(clusterProfiler::gson_GO(OrgDb = species_row$orgdb, ont = ont))
    }
    if (is_pic_error(out)) return(list(ok = FALSE, reason = out$message, gson = NULL))
    return(list(ok = TRUE, reason = "", gson = out))
  }

  if (resource == "KEGG") {
    if (is.na(species_row$kegg_code) || species_row$kegg_code == "") {
      return(list(ok = FALSE, reason = "KEGG code is missing", gson = NULL))
    }
    out <- quiet_try(clusterProfiler::gson_KEGG(species = species_row$kegg_code, keyType = "kegg"))
    if (is_pic_error(out)) return(list(ok = FALSE, reason = out$message, gson = NULL))
    return(list(ok = TRUE, reason = "", gson = out))
  }

  if (resource == "WIKIPATHWAYS") {
    if (is.na(species_row$wp_species) || species_row$wp_species == "") {
      return(list(ok = FALSE, reason = "WikiPathways species is missing", gson = NULL))
    }
    wp_fun <- get("gson_WP", envir = asNamespace("clusterProfiler"))
    wp_formals <- names(formals(wp_fun))
    wp_args <- list()
    if ("species" %in% wp_formals) wp_args$species <- species_row$wp_species
    if ("organism" %in% wp_formals) wp_args$organism <- species_row$wp_species
    out <- quiet_try(do.call(wp_fun, wp_args))
    if (is_pic_error(out)) return(list(ok = FALSE, reason = out$message, gson = NULL))
    return(list(ok = TRUE, reason = "", gson = out))
  }

  if (resource == "REACTOME") {
    if (!safe_package_available("ReactomePA")) return(list(ok = FALSE, reason = "ReactomePA not available", gson = NULL))
    if (is.na(species_row$reactome_species) || species_row$reactome_species == "") {
      return(list(ok = FALSE, reason = "Reactome species is missing", gson = NULL))
    }
    react_fun <- get("gson_Reactome", envir = asNamespace("ReactomePA"))
    react_formals <- names(formals(react_fun))
    react_args <- list()
    if ("species" %in% react_formals) react_args$species <- species_row$reactome_species
    if ("organism" %in% react_formals) react_args$organism <- species_row$reactome_species
    out <- quiet_try(do.call(react_fun, react_args))
    if (is_pic_error(out)) return(list(ok = FALSE, reason = out$message, gson = NULL))
    return(list(ok = TRUE, reason = "", gson = out))
  }

  list(ok = FALSE, reason = "unsupported resource", gson = NULL)
}

prepare_enrichment_libs <- function(target_genome = NULL, target_resources = NULL) {
  dir.create(enrichment_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(register_dir, recursive = TRUE, showWarnings = FALSE)
  species_tbl <- build_species_table(target_genome = target_genome)
  if (nrow(species_tbl) == 0) {
    stop(sprintf("No genome entry found for --genome %s", as.character(target_genome)), call. = FALSE)
  }
  selected_resources <- if (is.null(target_resources)) {
    enrichment_resources
  } else {
    req <- toupper(trimws(as.character(target_resources)))
    req <- req[req != ""]
    if (length(req) == 0) enrichment_resources else req
  }
  selected_resources <- unique(selected_resources)
  unknown_resources <- setdiff(selected_resources, enrichment_resources)
  if (length(unknown_resources) > 0) {
    stop(
      sprintf(
        "Unknown target_resources: %s (known: %s)",
        paste(unknown_resources, collapse = ","),
        paste(enrichment_resources, collapse = ",")
      ),
      call. = FALSE
    )
  }
  manifest_rows <- list()
  now_jst <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "Asia/Tokyo")
  normalize_updated_at_jst <- function(x, fallback) {
    x <- as.character(x)
    out <- vapply(x, function(v) {
      if (is.na(v) || v == "") return(fallback)
      if (grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}[+-]\\d{4}$", v)) return(v)
      parsed <- suppressWarnings(as.POSIXct(v, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
      if (is.na(parsed)) {
        parsed <- suppressWarnings(as.POSIXct(v, format = "%Y-%m-%dT%H:%M:%S%z", tz = "Asia/Tokyo"))
      }
      if (is.na(parsed)) return(fallback)
      format(parsed, "%Y-%m-%dT%H:%M:%S%z", tz = "Asia/Tokyo")
    }, character(1))
    out
  }

  for (i in seq_len(nrow(species_tbl))) {
    species_row <- species_tbl[i, , drop = FALSE]
    for (resource in selected_resources) {
      message(sprintf("[INFO] Preparing %s for %s (%s)", resource, species_row$genome[[1]], species_row$species_key[[1]]))
      build_res <- build_gson_for_species(resource, species_row[1, ])
      if (!isTRUE(build_res$ok) || is.null(build_res$gson)) {
        message(sprintf("[INFO] Skip %s for %s: %s", resource, species_row$genome[[1]], as.character(build_res$reason)))
        next
      }

      orgdb_pkg <- pick_orgdb_for_resource(resource, species_row[1, ])
      symbol_res <- normalize_gson_to_symbol(build_res$gson, orgdb_pkg = orgdb_pkg)
      if (!isTRUE(symbol_res$ok) || is.null(symbol_res$gson)) {
        message(sprintf("[INFO] Skip %s for %s: %s", resource, species_row$genome[[1]], as.character(symbol_res$reason)))
        next
      }

      out_name <- sprintf("%s__%s.gson", resource, species_row$species_key[[1]])
      out_file <- file.path(enrichment_out_dir, out_name)
      safe_save_gson(symbol_res$gson, out_file)

      manifest_rows[[length(manifest_rows) + 1]] <- tibble::tibble(
        genome = as.character(species_row$genome[[1]]),
        resource = as.character(resource),
        species_key = as.character(species_row$species_key[[1]]),
        updated_at = now_jst
      )
    }
  }

  raw_map_tbl <- if (length(manifest_rows) == 0) {
    tibble::tibble(
      genome = character(),
      resource = character(),
      species_key = character(),
      updated_at = character()
    )
  } else {
    dplyr::bind_rows(manifest_rows) |>
      dplyr::distinct(.data$genome, .data$resource, .keep_all = TRUE)
  }

  empty_wide_map <- function() {
    out <- tibble::tibble(
      genome = character(),
      species_key = character(),
      updated_at = character()
    )
    for (resource in enrichment_resources) out[[resource]] <- logical()
    out
  }

  normalize_wide_map <- function(tbl) {
    if (nrow(tbl) == 0 || !"genome" %in% colnames(tbl)) {
      return(empty_wide_map())
    }
    out <- tbl |>
      dplyr::transmute(
        genome = as.character(.data$genome),
        species_key = if ("species_key" %in% colnames(tbl)) as.character(.data$species_key) else "",
        updated_at = if ("updated_at" %in% colnames(tbl)) as.character(.data$updated_at) else now_jst
      )
    for (resource in enrichment_resources) {
      if (resource %in% colnames(tbl)) {
        vals <- tolower(trimws(as.character(tbl[[resource]])))
        out[[resource]] <- vals %in% c("true", "t", "1", "yes", "y")
      } else {
        out[[resource]] <- FALSE
      }
    }
    out
  }

  existing_raw <- if (file.exists(genome_runtime_map_file)) {
    suppressMessages(
      readr::read_tsv(
        genome_runtime_map_file,
        col_types = readr::cols(.default = readr::col_character()),
        progress = FALSE
      )
    )
  } else {
    tibble::tibble()
  }

  existing_map <- normalize_wide_map(existing_raw)

  current_map <- normalize_wide_map(
    dplyr::left_join(
      species_tbl |>
        dplyr::transmute(
          genome = as.character(.data$genome),
          species_key = as.character(.data$species_key),
          updated_at = now_jst
        ),
      raw_map_tbl |>
        dplyr::transmute(
          genome = as.character(.data$genome),
          resource = as.character(.data$resource),
          available = TRUE
        ) |>
        tidyr::pivot_wider(
          names_from = "resource",
          values_from = "available",
          values_fill = FALSE
        ),
      by = "genome"
    )
  )

  if (!is.null(target_genome) && target_genome != "") {
    target_existing <- existing_map |>
      dplyr::filter(.data$genome == .env$target_genome)
    if (nrow(target_existing) == 0) {
      target_existing <- species_tbl |>
        dplyr::transmute(
          genome = as.character(.data$genome),
          species_key = as.character(.data$species_key),
          updated_at = now_jst
        )
      for (resource in enrichment_resources) target_existing[[resource]] <- FALSE
    }
    current_row <- current_map |>
      dplyr::filter(.data$genome == .env$target_genome)
    if (nrow(current_row) == 0) {
      current_row <- target_existing
    } else {
      for (resource in selected_resources) {
        target_existing[[resource]] <- current_row[[resource]]
      }
      target_existing$species_key <- current_row$species_key
      target_existing$updated_at <- now_jst
      current_row <- target_existing
    }
    existing_map <- existing_map |>
      dplyr::filter(.data$genome != .env$target_genome)
    current_map <- current_row
  }

  map_tbl <- dplyr::bind_rows(existing_map, current_map) |>
    dplyr::distinct(.data$genome, .keep_all = TRUE) |>
    dplyr::mutate(
      species_key = dplyr::if_else(is.na(.data$species_key), "", .data$species_key),
      updated_at = normalize_updated_at_jst(.data$updated_at, now_jst)
    ) |>
    dplyr::arrange(.data$genome) |>
    dplyr::select("genome", dplyr::all_of(enrichment_resources), "species_key", "updated_at")

  readr::write_tsv(map_tbl, genome_runtime_map_file)
  message(sprintf("[INFO] Wrote enrichment map: %s (%d rows)", genome_runtime_map_file, nrow(map_tbl)))
}

if (sys.nframe() == 0) {
  prepare_enrichment_libs()
}
