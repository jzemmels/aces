#' Apply ACES (Automatic Cartridge Evidence Scoring) comparison algorithm to a
#' set of cartridge case scans.
#' @rdname acesAlgorithm
#' @export
#'
#' @param ... either a sequence of x3p objects separated by commas like
#'   'x3p1,x3p2,x3p3,...', a single list of x3p objects separated by commas like
#'   'list(x3p1,x3p2,x3p3,...)', or a vector of x3p file paths like
#'   'c("path/to/file1.x3p", "path/to/file2.x3p", "path/to/file3.x3p")'
#' @param x3p_labels (optional) character vector of x3p labels
#' @param return_registrations whether to return the full scan and cell-based
#'   registrations as tibble columns, each containing nested tibble objects.
#'   This is useful for understanding the feature values, but greatly increases
#'   the size of the returned object.
#' @param standardize_resolutions whether to force resolutions of all scans
#'   equal. If FALSE, then function will throw error when scans aren't same
#'   resolution.
#' @param quiet whether to print internal function messages to the console

comparison_aces <- function(...,
                            x3p_labels = NULL,
                            return_registrations = FALSE,
                            standardize_resolutions = FALSE,
                            quiet = FALSE) {
  input <- list(...)

  p <- progressr::progressor(steps = 1)

  p(message = "Preparing x3p objects for comparison")
  # return a list of labeled x3p objects
  x3pList <- prepare_x3ps(input,
    x3p_labels = x3p_labels,
    standardize_resolutions = standardize_resolutions
  )

  # create a list of all pairwise comparisons (excluding replicates & self
  # comparisons)
  x3pComparisons <- prepare_comparisonNames(x3pList)

  p <- progressr::progressor(steps = nrow(x3pComparisons) * 3)
  progressr::handlers("progress")

  # p(message = "Registering full scans")

  # register full scans to each other
  fullScanRegistrations <- x3pComparisons %>%
    purrr::pmap_dfr(function(...) {
      p(message = paste0("Registering full scans for comparison ", tibble::tibble(...)$comparisonName))
      aces_registration_fullScan(
        x3pList = x3pList,
        x3pComparison = tibble::tibble(...),
        standardize_resolutions = standardize_resolutions,
        quiet = quiet
      )
    })

  # p(message = "Registering cells")
  # then register cells to each other
  cellBasedRegistrations <- fullScanRegistrations %>%
    dplyr::select(comparisonName, direction, cellHeightValues, alignedTargetCell) %>%
    dplyr::group_by(comparisonName, direction) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(~ {
      p(message = paste0("Registering cells for comparison ", unique(.$comparisonName)))
      aces_registration_cellBased(.)
    })

  p(message = "Computing features")
  # compute visual and registration features for full scans
  fullScan_features <- fullScanRegistrations %>%
    dplyr::group_by(comparisonName, direction) %>%
    feature_aLaCarte(
      features = c("visual", "registration"),
      threshold = function(x3p1, x3p2) {
        sd(abs(x3p1$surface.matrix - x3p2$surface.matrix), na.rm = TRUE)
      }
    ) %>%
    dplyr::group_by(comparisonName) %>%
    dplyr::summarize(across(dplyr::where(is.numeric), ~ mean(., na.rm = TRUE))) %>%
    purrr::set_names(paste0("fullScan_", names(.))) %>%
    dplyr::rename(comparisonName = fullScan_comparisonName)

  # compute visual and registration features for cells
  cellBased_features_registVisual <- cellBasedRegistrations %>%
    dplyr::group_by(comparisonName, direction) %>%
    scored::feature_aLaCarte(
      features = c("registration", "visual"), quiet = quiet,
      threshold = function(x3p1, x3p2) {
        sd(abs(x3p1$surface.matrix - x3p2$surface.matrix), na.rm = TRUE)
      }
    ) %>%
    dplyr::group_by(comparisonName) %>%
    dplyr::summarize(across(tidyselect::where(is.numeric), ~ mean(.)))

  # compute density-based features for
  cellBased_features_density <-
    cellBasedRegistrations %>%
    dplyr::select(-c(cellHeightValues, alignedTargetCell)) %>%
    dplyr::arrange(direction) %>%
    dplyr::group_by(comparisonName, direction) %>%
    feature_aLaCarte(
      features = c("density"), quiet = quiet,
      eps = 3, minPts = 3
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      thetaDiff = as.numeric(thetaDiff),
      translationDiff = as.numeric(translationDiff),
      clusterSize = as.numeric(clusterSize)
    ) %>%
    dplyr::group_by(comparisonName) %>%
    dplyr::summarize(across(dplyr::where(is.numeric), ~ mean(.))) %>%
    dplyr::mutate(clusterInd = !is.na(clusterSize))

  # combine cell-based feature data frames
  cellBased_features <-
    cellBased_features_registVisual %>%
    dplyr::left_join(cellBased_features_density, by = "comparisonName") %>%
    purrr::set_names(paste0("cellBased_", names(.), "_4x4")) %>%
    dplyr::select(-c(cellBased_comparisonName_4x4))

  # combine cell-based features with full scan features
  all_features <-
    dplyr::bind_cols(
      fullScan_features,
      cellBased_features
    ) %>%
    dplyr::distinct() %>%
    dplyr::select(
      comparisonName,
      dplyr::contains("fullScan_"), dplyr::contains("cellBased_")
    ) %>%
    # throw away some columns not used in ACES
    dplyr::select(-c(
      contains("fullScan_neighborhoodSizeAve_sd"), contains("fullScan_neighborhoodSizeSD_sd"),
      contains("fullScan_differenceCor_sd"), contains("fullScan_filteredRatio_sd"),
      contains("cellBased_neighborhoodSizeAve_sd"), contains("cellBased_neighborhoodSizeSD_sd"),
      contains("cellBased_differenceCor_sd"),
      contains("fullScan_ccfMean"), contains("cellBased_ccfMean"), contains("cellBased_ccfSD")
    ))

  if (return_registrations) {
    all_features <- all_features %>%
      dplyr::mutate(
        fullScanRegistrations = list(fullScanRegistrations),
        cellBasedRegistrations = list(cellBasedRegistrations)
      )
  }

  return(all_features)

  # TODO:
  # xxx - cell-based comparison based on full scan registration
  # xxx - full scan and cell-based feature calculation
  # xxx - return summarized feature values in a tibble (add metadata too?)
  # NEEDED FEATURES:
  # - should default to scaling all scans to the resolution used in dissertation, but allow users to turn off this setting/set their own resolution
  # xxx - messages about function progress & "quiet" argument
  # - include option to compute similarity score using predict(), with caveat about equal resolutions
}

# function to put x3ps into a standardized format for use in other functions
prepare_x3ps <- function(x3pFiles, x3p_labels = NULL, standardize_resolutions = FALSE) {
  # x3pFiles <- list(...)

  ####### Attempt to load x3p objects into a list

  # input type #1: pass an arbitrary number of x3p objects separated by commas
  # to ... (nothing needs to be done here)
  if (all(purrr::map_chr(x3pFiles, class) == "x3p")) {
    x3pList <- x3pFiles
  }
  # input type #2: pass one list() of x3p objects to ...
  else if (length(x3pFiles) == 1 & all(purrr::map_chr(x3pFiles[[1]], class) == "x3p")) {
    x3pList <- x3pFiles[[1]]
  }
  # input type #3: pass a vector of file paths
  else if (all(purrr::map_chr(x3pFiles[[1]], class) == "character")) {
    stopifnot("Provide a vector of 2 or more file paths for valid comparison" = length(x3pFiles[[1]]) > 1)

    x3pList <- purrr::map(
      x3pFiles[[1]],
      function(fileName) {
        assertthat::assert_that(file.exists(fileName), msg = paste0("Cannot find file ", fileName))

        return(x3ptools::read_x3p(fileName))
      }
    )

    # if no labels provided, guess the x3p labels using the file names
    if (is.null(x3p_labels)) {
      x3p_labels <- stringr::str_remove(basename(x3pFiles[[1]]), "\\.x3p")
    }
  } else {
    stop("Incorrect input type. See ?comparison_aces for more information.")
  }

  # last resort: label each x3p as x3p1, x3p2, etc.
  if (is.null(x3p_labels)) {
    x3p_labels <- paste0("x3p", seq_along(x3pList))
  }
  stopifnot("x3p_labels should be same length as the number of input x3ps." = length(x3p_labels) == length(x3pList))

  ####### Check that x3ps are ready for comparison

  # check that scan resolutions are all equal
  if (!standardize_resolutions) {
    # force the resolutions to be equal
    assertthat::assert_that(are_equal(purrr::map_dbl(x3pList, ~ .$header.info$incrementY)) & are_equal(purrr::map_dbl(x3pList, ~ .$header.info$incrementX)),
      msg = "Scan resolutions are not equal. Either set the resolutions manually yourself or force resolutions to be equal by setting argument 'standardize_resolutions = TRUE'."
    )
  }

  return(purrr::set_names(x3pList, x3p_labels))
}

# create a data frame containing pairwise comparison names
prepare_comparisonNames <- function(x3pList) {
  tidyr::expand_grid(
    reference = names(x3pList),
    target = names(x3pList)
  ) %>%
    dplyr::left_join(
      data.frame(
        scanName = names(x3pList),
        scanInd = 1:length(x3pList)
      ),
      by = c("reference" = "scanName")
    ) %>%
    dplyr::rename(referenceInd = scanInd) %>%
    dplyr::left_join(
      data.frame(
        scanName = names(x3pList),
        scanInd = 1:length(x3pList)
      ),
      by = c("target" = "scanName")
    ) %>%
    dplyr::rename(targetInd = scanInd) %>%
    dplyr::filter(referenceInd < targetInd) %>%
    dplyr::select(-c(referenceInd, targetInd)) %>%
    dplyr::mutate(comparisonName = paste0(reference, "_vs_", target))
}

# checks if elements of a vector are all equal by checking each pair of elements
# adapted from: https://stackoverflow.com/a/27331553/14000041
are_equal <- function(x) {
  # more than one object required
  if (length(x) < 2) stop("More than one object required")

  # matrix of object name pairs
  pairs <- t(combn(x, 2))

  # if only two objects, return all.equal() for them
  if (nrow(pairs) == 1) {
    return(all.equal(pairs[1, 1], pairs[1, 2]))
  }

  eq.fun <- function(z, y) {
    all.eq <- all.equal(z, y)
    name <- paste0(z, " vs. ", y)
    return(list(all.eq, name))
  }

  # list of eq.fun object comparisons
  out <- vector(mode = "list", length = nrow(pairs))

  for (w in 1:nrow(pairs)) {
    eq.list <- eq.fun(pairs[w, 1], pairs[w, 2])
    out[[w]] <- eq.list[[1]]
    names(out)[w] <- eq.list[[2]]
  }

  return(isTRUE(all(as.logical(unlist(out)))))
}

# function to perform the full scan registration using comparison_fullScan
aces_registration_fullScan <- function(x3pList, x3pComparison, standardize_resolutions = FALSE, quiet = FALSE) {
  reference <- x3pList[[which(x3pComparison$reference == names(x3pList))[1]]]
  target <- x3pList[[which(x3pComparison$target == names(x3pList))[1]]]

  if (standardize_resolutions & !isTRUE(all.equal(reference$header.info$incrementX, target$header.info$incrementX))) {
    if (!quiet) message(paste0("Note: resolutions not equal for comparison ", x3pComparison$comparisonName))

    if (reference$header.info$incrementX > target$header.info$incrementX) {
      target <- x3ptools::x3p_interpolate(target, resx = reference$header.info$incrementX)
    } else {
      reference <- x3ptools::x3p_interpolate(reference, resx = target$header.info$incrementX)
    }
  }

  fullScanComparison <-
    scored::comparison_fullScan(reference, target) %>%
    dplyr::group_by(direction) %>%
    dplyr::filter(fft_ccf == max(fft_ccf)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(comparisonName = x3pComparison$comparisonName) %>%
    dplyr::arrange(direction) %>%
    # need to rescale the surface matrix values back to microns:
    dplyr::mutate(
      cellHeightValues = purrr::map(
        cellHeightValues,
        function(dat) {
          if (!all(is.na(dat))) {
            dat$surface.matrix <- dat$surface.matrix * dat$cmcR.info$scaleByVal # convert to micron scale
          }

          return(dat)
        }
      ),
      alignedTargetCell = purrr::map(
        alignedTargetCell,
        function(dat) {
          if (!all(is.na(dat))) {
            dat$surface.matrix <- dat$surface.matrix * dat$cmcR.info$scaleByVal # convert to micron scale
          }

          return(dat)
        }
      )
    )

  return(fullScanComparison)
}

aces_registration_cellBased <- function(fullScanRegist) {
  comparison_cellBased(
    reference = fullScanRegist$cellHeightValues[[1]],
    target = fullScanRegist$alignedTargetCell[[1]],
    thetas = -2:2, numCells = c(4, 4),
    maxMissingProp = .99,
    sideLengthMultiplier = 1.1,
    returnX3Ps = TRUE
  ) %>%
    dplyr::group_by(cellIndex) %>%
    # to save on space, only keep x3ps associated with the max CCF for each cell index
    dplyr::mutate(
      cellHeightValues = ifelse(fft_ccf == max(fft_ccf), cellHeightValues, NA),
      alignedTargetCell = ifelse(fft_ccf == max(fft_ccf), alignedTargetCell, NA)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      comparisonName = unique(fullScanRegist$comparisonName),
      direction = unique(fullScanRegist$direction)
    ) %>%
    dplyr::select(comparisonName, direction, cellIndex, dplyr::everything()) %>%
    # the surface values are normalized (subtract mean & divide by sd) during
    # registration, so we need to convert them back to the original scale
    dplyr::mutate(
      cellHeightValues = purrr::map(
        cellHeightValues,
        function(dat) {
          if (!all(is.na(dat))) {
            dat$surface.matrix <- dat$surface.matrix * dat$cmcR.info$scaleByVal
          }

          return(dat)
        }
      ),
      alignedTargetCell = purrr::map(
        alignedTargetCell,
        function(dat) {
          if (!all(is.na(dat))) {
            dat$surface.matrix <- dat$surface.matrix * dat$cmcR.info$scaleByVal
          }

          return(dat)
        }
      )
    )
}


# function that changes the resolutions of all scans to the same resolution.
# Default to lowest resolution == maximum unit per pixel value
x3p_standardize_resolutions <- function(x3pList, FUN = max) {
  newX <- FUN(purrr::map_dbl(x3pList, ~ .$header.info$incrementX))
  newY <- FUN(purrr::map_dbl(x3pList, ~ .$header.info$incrementY))

  stopifnot("Resolution information missing in one or more x3ps." = is.double(newX) & length(newX) > 0 & is.double(newY) & length(newY) > 0)

  return(purrr::map(x3pList, ~ x3ptools::interpolate_x3p(., resx = newX, resy = newY)))
}
