

getEstimatesFromIchorCnaRdataFiles <- function(pathToSampleDir, binsize = 10^6){
  
  # pathToSampleDir <- "../runs/miracle_swgs_for_garbage_paper/results/ichorcna/L5730_swgs_ichorCNA_files"
  
  # pathToSampleDir <- pathToSampleDirIchorCNA_DELPHI_colo_swgs[1]
  
  # pathToSampleDir <- "./results/ichorcna_50kb/L9684_ichorCNA_files"
  
  # pathToSampleDir <- "./results/ichorcna_diff/patnr_55_diff_ichorCNA_files/"
  
  sampleFilesIchorCNA <- paste0(pathToSampleDir, "/", list.files(pathToSampleDir))
  
  all_results_path <- sampleFilesIchorCNA |> 
    stringr::str_subset(pattern = "ALL_RESULTS.RData")
  
  load(all_results_path)
  
  standard_rdata_file_path <- sampleFilesIchorCNA |>
    stringr::str_subset(pattern = ".RData") |>
    stringr::str_subset(pattern = "ALL_RESULTS.RData", negate = TRUE)
  
  load(standard_rdata_file_path)
  
  sampleName <- all_results_path |>
    stringr::str_extract(pattern = "L[0-9]{2,5}|EE[0-9]{5}|L[0-9]{4,5}|SE[0-9]{2}-[0-9]{4}|SE[0-9]{2}-[0-9]{4}|I23-1365-[0-9]{2}
                             r|[A-Z]{1}[0-9]{2}|[0-9]{3,4}")
  
  subSample <- case_when(
    stringr::str_detect(all_results_path, "subsample_0.75") ~ 0.75,
    stringr::str_detect(all_results_path, "subsample_0.5") ~ 0.50,
    stringr::str_detect(all_results_path, "subsample_0.25") ~ 0.25,
    stringr::str_detect(all_results_path, "subsample_0.1") ~ 0.1,
    stringr::str_detect(all_results_path, "subsample_0.05") ~ 0.05,
    stringr::str_detect(all_results_path, "subsample_1") ~ 1,
    TRUE ~ 1
  )
  
  segments <- sampleFilesIchorCNA |>
    stringr::str_subset(pattern = ".seg") |> 
    stringr::str_subset(pattern = ".cna.seg|.seg.txt", negate = TRUE) |>
    readr::read_delim(delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
  
  segmentsIchorCNA <- segments |> 
    mutate(
      sampleName = sampleName,
      subSample = subSample,
      seqnames =  paste0("chr", chr)
    ) |>
    select(!c(sample, chr)) |> 
    relocate(sampleName, subSample, seqnames)
  
  segmentsLong <- bind_rows(lapply(1:length(segmentsIchorCNA$bins), function(i){

    #i <- 1

    x <- segmentsIchorCNA[i,]

    segmentsLong <- tibble(
      analysis_ID3 = rep(unique(x$sampleName), x$bins),
      segmentIchorCNA = i,
      start = seq(x$start, x$end, by = binsize),
      seqnames = rep(x$seqnames, x$bins)
    )

    return(segmentsLong)
  }))
  
  x <- lapply(1:length(results), function(i){
    
    #i <- 1
    
    iter <- results[[i]]$results$iter
    
    loglik <- results[[i]]$results$loglik[iter]
    
    return(loglik)
  }) |> unlist()
  
  i <- which(x == max(x))
  
  estimatesIchorCNA <-  results[[i]][[1]][[1]] |>
    #estimatesIchorCNA <-  hmmResults.cor[[1]][[1]] |>
    as_tibble() |> 
    tibble::remove_rownames() |>
    mutate(
      sampleName = sampleName,
      tumorFractionEstimate = 1 - results[[i]]$results$n[[results[[i]]$results$iter]],
      subcloneFractionEstimate = 1 - results[[i]]$results$sp[[results[[i]]$results$iter]],
      ploidyEstimate = results[[i]]$results$phi[[results[[i]]$results$iter]],
      seqnames =  paste0("chr", chr),
      rawReadCount = tumour_copy[[1]]$reads,
      gcContent = tumour_copy[[1]]$gc,
      mappingBias = tumour_copy[[1]]$map,
      sexIchor = gender$gender,
      chrYCovRatio = gender$chrYCovRatio,
      subSample = subSample
    ) |>
    select(!c(sample, chr)) |> 
    left_join(segmentsLong, by = join_by(start, seqnames)) |>
    relocate(sampleName, subSample, tumorFractionEstimate, ploidyEstimate, seqnames)
  
  return(estimatesIchorCNA)
  
}