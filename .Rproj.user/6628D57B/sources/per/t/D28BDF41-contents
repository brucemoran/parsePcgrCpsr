#! R

#' Parse PCGR, CPSR json.gz files for single sample
#' @param PCGR_JSON pcgr file to be parsed
#' @param CPSR_JSON cpsr file to be parsed
#' @return named list, based on VCF_SAMPLE_ID, of PCGR: snv, CNA, TMB, MSI; CPSR: SNV
#' @rdname parse_save_jsons
#' @export

parse_jsons <- function(PCGR_JSON = NULL, CPSR_JSON = NULL){

  ##set flags to be TRUE if we have JSON specified
  pcgr <- cpsr <- FALSE

  ##list to hold what we get
  olist <- list(pcgr_snv = NULL,
                pcgr_cna = NULL,
                pcgr_tmb = NULL,
                pcgr_msi = NULL,
                cpsr_snv = NULL)

  if(!is.null(PCGR_JSON)){
    pcgr <- TRUE
    pcgr_data <- jsonlite::fromJSON(PCGR_JSON)
    olist[["pcgr_snv"]] <- tibble::as_tibble(pcgr_data$content$snv_indel$variant_set$tsv)
    olist[["pcgr_cna"]] <- tibble::as_tibble(pcgr_data$content$cna$variant_set$tsv)
    olist[["pcgr_tmb"]] <- pcgr_data$content$tmb$v_stat$tmb_estimate
    olist[["pcgr_msi"]] <- pcgr_data$content$msi$prediction
  }


  if(!is.null(CPSR_JSON)){
    cpsr <- TRUE
    cpsr_data <- jsonlite::fromJSON(CPSR_JSON)
    olist[["cpsr_snv"]] <- tibble::as_tibble(cpsr_data$content$snv_indel$variant_set$tsv)
  }

  ##if both jsons, check they are from one sample
  if(pcgr && cpsr){
    if(unlist(pcgr_snv$VCF_SAMPLE_ID)[1] != unlist(cpsr_snv$VCF_SAMPLE_ID)[1]){
      print(paste("Sample IDs do not match",
                  c(unlist(pcgr_snv$VCF_SAMPLE_ID)[1],
                    unlist(cpsr_snv$VCF_SAMPLE_ID)[1])), sep = ": ")
    } else {
      return(olist)
    }
    ##if just one json, return result
  } else {
    return(olist)
  }
}

#' Find and parse jsons
#' @param PARSE_DIR character string path where JSONs are a held
#' @param DELIM character delimiter for parsing names from json files (default: '.')
#' @return list of lists, each per sample for the data parsed in parse_jsons
#' @rdname find_parse_jsons
#' @export

find_parse_jsons <- function(PARSE_DIR, DELIM = "\\."){

  ##find names of files available in PARSE_DIR
  print("Recursively searching for JSON files, please be patient...")
  json_paths <- dir(PARSE_DIR, pattern = "json.gz", recursive = TRUE, full.names = TRUE)
  print(paste("Found JSONs:", json_paths), sep = "\n ")
  json_names <- unique(unlist(lapply(json_paths, function(j){
                                strsplit(rev(strsplit(j, "/")[[1]])[1], DELIM)[[1]][1]
                              })
                            ))

  json_file_list <- lapply(json_names, function(f){
                      grep(f, json_paths, value = TRUE)
                    })
  names(json_file_list) <- json_names

  ##parse those files in json_files
  data_list <- lapply(seq_along(json_file_list), function(f){
    print(paste0("Working on: ", names(json_file_list)[f]))
    if(length(json_file_list[[f]])>1){
      cpsr_json <- grep("cpsr", unlist(json_file_list[[f]])[1], value = TRUE)
      pcgr_json <- grep("pcgr", unlist(json_file_list[[f]])[2], value = TRUE)
    } else {
      cpsr_json <- grep("cpsr", unlist(json_file_list[[f]])[1], value = TRUE)
      pcgr_json <- grep("pcgr", unlist(json_file_list[[f]])[1], value = TRUE)
    }
    parse_jsons(PCGR_JSON = pcgr_json,
                CPSR_JSON = cpsr_json)
  })
  names(data_list) <- json_names
  return(data_list)
}

#' Compile multiple sample data into one XLSX workbook
#' @param DATA_LIST list of data returned from find_parse_jsons
#' @param TITLE character string naming output XLSX
#' @param OUTDIR character string path to write file to (default: getwd())
#' @param DELIM character delimiter for parsing names from json files (default: '.')
#' @return none, writes an XLSX workbook with worksheets per data type
#'         where the same GENOMIC_CHANGE is found, compile all samples containing into one line
#' @importFrom magrittr %>%
#' @rdname compile_xlsx
#' @export

compile_xlsx <- function(DATA_LIST, TITLE = "comp_PCGR_CPSR", OUTDIR = "./", DELIM = "\\."){

  ##SNVs
  pcgr_snv <- snv_reduce(DATA_LIST, WHICH = "pcgr_snv")
  cpsr_snv <- snv_reduce(DATA_LIST, WHICH = "cpsr_snv")

  ##TMB, MSI
  msi_tmb <- do.call(rbind, lapply(seq_along(DATA_LIST), function(f){
    msi <- tmb <- "Not estimated"
    if(length(DATA_LIST[[f]]$pcgr_msi)>0){
      msi <- DATA_LIST[[f]]$pcgr_msi
    }
    if(length(DATA_LIST[[f]]$pcgr_tmb)>0){
      tmb <- DATA_LIST[[f]]$pcgr_tmb
    }
    return(tibble::tibble(sampleID = names(DATA_LIST)[f],
                          TMB = tmb, MSI = msi))
  }))

  ##write to workbook
  wb <- openxlsx::createWorkbook(title = TITLE, creator = "parsePcgrCpsr")
  openxlsx::addWorksheet(wb, "Somatic_SNVs")
  openxlsx::addWorksheet(wb, "Germline_SNVs")
  openxlsx::addWorksheet(wb, "MSI_TMB")
  openxlsx::writeData(wb, "Somatic_SNVs", pcgr_snv, colNames = TRUE, borders = "none", borderStyle = "none", rowNames = FALSE)
  openxlsx::writeData(wb, "Germline_SNVs", cpsr_snv, colNames = TRUE, borders = "none", borderStyle = "none", rowNames = FALSE)
  openxlsx::writeData(wb, "MSI_TMB", msi_tmb, colNames = TRUE, borders = "none", borderStyle = "none", rowNames = FALSE)
  openxlsx::saveWorkbook(wb, file = paste0(OUTDIR, "/", TITLE, ".xlsx"), overwrite = TRUE)

}

#' Compile multiple sample data into one XLSX workbook
#' @param DATA_LIST list of data returned from find_parse_jsons
#' @param WHICH character string which list entry to work on: pcgr_snv, cpsr_snv
#' @return tb of reduced SNVs
#' @importFrom magrittr %>%
#' @rdname snv_reduce
#' @export

snv_reduce <- function(DATA_LIST, WHICH){
  snv_tb <- do.call(rbind, lapply(DATA_LIST, function(f){
    return(f[[WHICH]])
  }))
  s2 <- dplyr::group_by(.data = snv_tb, GENOMIC_CHANGE) %>%
    dplyr::add_count() %>%
    dplyr::filter(n > 1)
  s1 <- dplyr::filter(.data = snv_tb, ! GENOMIC_CHANGE %in% unique(unlist(s2$GENOMIC_CHANGE)))
  ##each GENOMIC_CHANGE found in s2 needs to be taken, VCF_SAMPLE_ID compiled and rbind to s1
  s3 <- do.call(rbind, lapply(unique(unlist(s2$GENOMIC_CHANGE)), function(f){
    gcf <- dplyr::filter(.data = s2, GENOMIC_CHANGE %in% !!f)
    samples <- paste(unlist(gcf$VCF_SAMPLE_ID), collapse = ",")
    gcf[1,"VCF_SAMPLE_ID"] <- samples
    gcf[1,]
  }))
  snv_tb <- dplyr::bind_rows(s3, s1)
  return(snv_tb)
}




