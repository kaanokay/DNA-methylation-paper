# Permutation test

start.time <- Sys.time()

set.seed(100000)

library(bsseq)
library(data.table)
library(gtools)
library(parallel)
library(locfit)
library(BiocParallel)
library(rtracklayer)
library(sva)

# Function 1

getNullDistribution_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2,
                                              estimate.var, local.correct,
                                              cutoff, stat, maxGap, mc.cores = 1) {
  stopifnot(nrow(idxMatrix1) == nrow(idxMatrix2))
  message(sprintf("[getNullDistribution_BSmooth.tstat] performing %d permutations\n", nrow(idxMatrix1)))
  nullDist <- mclapply(1:nrow(idxMatrix1), function(ii) {
    ptime1 <- proc.time()
    BS.tstat <- BSmooth.tstat(BSseq, estimate.var = estimate.var,
                              group1 = idxMatrix1[ii,],
                              group2 = idxMatrix2[ii,],
                              local.correct = local.correct, maxGap = 10^8,
                              verbose = FALSE, k = 21)
    dmrs0 <- dmrFinder(BS.tstat, stat = stat, cutoff = cutoff, maxGap = 300)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    message(sprintf("[getNullDistribution_BSmooth.tstat] completing permutation %d in %.1f sec\n", ii, stime))
    dmrs0
  }, mc.cores = min(nrow(idxMatrix1), mc.cores), mc.preschedule = FALSE)
  nullDist
}

# Function 2

.isHDF5ArrayBacked <- function(object) {
  if (is(object, "SummarizedExperiment")) {
    return(all(vapply(X = assays(object, withDimnames = FALSE),
                      FUN = .isHDF5ArrayBacked,
                      FUN.VALUE = logical(1L))))
  }
  if (is(object, "DelayedArray")) {
    seed <- object@seed
    seed_classes <- .getSeedClasses(seed)
    is_hdf5_backed <- vapply(unlist(seed_classes, use.names = FALSE),
                             extends, class2 = "HDF5ArraySeed",
                             logical(1L))
    return(any(is_hdf5_backed))
  } else if (is.matrix(object)) {
    FALSE
  } else if (is.null(object)) {
    FALSE
  } else {
    stop("Don't know how to handle object of class ", class(object))
  }
}

# Function 3

makeClusters <- function(hasGRanges, maxGap = 10^8) {
  chrOrder <- as.character(runValue(seqnames(hasGRanges)))
  if(anyDuplicated(chrOrder))
    stop("argument 'hasGRanges' is not properly order")
  grBase <- granges(hasGRanges)
  clusters <- reduce(resize(grBase, width = 2*maxGap + 1, fix = "center"))
  start(clusters) <- pmax(rep(1, length(clusters)), start(clusters))
  clusters.sp <- split(clusters, seqnames(clusters))
  stopifnot(all(sapply(clusters.sp, function(cluster.gr) {
    if(length(cluster.gr) <= 1) return(TRUE)
    all(start(cluster.gr)[-length(cluster.gr)] < end(cluster.gr)[-1])
  }))) # are the clusters ordered within the chromosome? This is probably guranteed
  clusters <- Reduce(c, clusters.sp[chrOrder])
  stopifnot(all(chrOrder == runValue(seqnames(clusters))))
  ov <- findOverlaps(grBase, clusters)
  clusterIdx <- split(as.matrix(ov)[,1], as.matrix(ov)[,2])
  names(clusterIdx) <- NULL
  clusterIdx
}

# Function 4

.rowSds <- function(x, rows = NULL, cols = NULL, ...) {
  row_vars <- .rowVars(x, rows = rows, cols = cols, ...)
  sqrt(row_vars)
}

# Function 5

.rowVars <- function(x, rows = NULL, cols = NULL, ...) {
  if (is(x, "DelayedArray")) {
    if (!is.null(rows)) {
      x <- x[rows, ]
    }
    if (!is.null(cols)) {
      x <- x[, cols]
    }
    row_vars <- rowVars(as.array(x), ...)
  } else {
    row_vars <- rowVars(x, rows = rows, cols = cols, ...)
  }
  row_vars
}

# Function 6

subsetDmrs <- function(nullDist) {
  if (is.null(nullDist) || is(nullDist, "try-error"))
    return(NULL)
  
  # Create an empty list to store the subset DMRs for each NULL distribution
  subsetDmrsList <- list()
  
  # Loop through each NULL distribution in nullDist
  for (i in 1:length(nullDist)) {
    xx <- nullDist[[i]]
    
    # Subset DMRs based on your criteria
    out <- xx[xx[, "n"] >= 3 & abs(xx[, "meanDiff"]) >= 0.05 &
                xx[, "invdensity"] <= 300, ]
    # If there are DMRs that meet the criteria, add them to the list
    if (nrow(out) > 0) {
      subsetDmrsList[[i]] <- out
    }
  }
  # Return the list of subsetted DMRs for each NULL distribution
  return(subsetDmrsList)
}

# Function 7

getFWER <- function(null, type = "blocks") {
  reference <- null[[1]]
  null <- null[-1]
  null <- null[!sapply(null, is.null)]
  better <- sapply(1:nrow(reference), function(ii) {
    # meanDiff <- abs(reference$meanDiff[ii])
    areaStat <- abs(reference$areaStat[ii])
    width <- reference$width[ii]
    n <- reference$n[ii]
    if (type == "blocks") {
      out <- sapply(null, function(nulldist) {
        # any(abs(nulldist$meanDiff) >= meanDiff &
        # nulldist$width >= width)
        any(abs(nulldist$areaStat) >= areaStat &
              nulldist$width >= width)
      })
    }
    if (type == "dmrs") {
      out <- sapply(null, function(nulldist) {
        # any(abs(nulldist$meanDiff) >= meanDiff &
        #     nulldist$n >= n)
        any(abs(nulldist$areaStat) >= areaStat &
              nulldist$n >= n)
      })
    }
    sum(out)
  })
  better
}

# BSseq object

BSseq.obj <- readRDS("smoothed.BSseq.object.no.any.filtering.rds")

filtered_data <- subset(pData(BSseq.obj), !(Samples == "Ctr10.FVBNJ.Batch.25.07.2023" |
                                                Samples == "Ctr10.B6.Batch.25.07.2023"))

BSseq.obj <- BSseq.obj[, rownames(filtered_data)]

# Remove chrX and chrY

BSseq.obj <- subset(BSseq.obj,
                          !seqnames %in% c("chrX", "chrY"))

# Coverage filtering

cov <- getCoverage(BSseq.obj)
keepLoci.ex <- which(rowSums(cov[, BSseq.obj$Groups == "EM.gene"] >= 3) >= 97 &
                       rowSums(cov[, BSseq.obj$Groups == "Control"] >= 3) >= 29)

BSseq.obj <- BSseq.obj[keepLoci.ex, ]

# Upload metadata which has batch info for each sample

metadata <- fread("metadata.csv")

# Assign colors to samples

pData(BSseq.obj)$col <- as.character(metadata$Color)

# Start permutation

group1 = rownames(pData(BSseq.obj)[pData(BSseq.obj)$Groups == 'EM.gene',])
group2 = rownames(pData(BSseq.obj)[pData(BSseq.obj)$Groups == 'Control',])

# Create permutation for group 1

num_permutations <- 100

# --- Set the number of permutations you want to generate end ---

# --- Combine the groups start ---

all_samples <- c(group1, group2)

# --- Combine the groups end ---

# --- Create a matrix to store samples 1000 times permuted start ---

permuted_matrix <- matrix(nrow=num_permutations, ncol=length(group1) + length(group2))

# --- Create a matrix to store samples 1000 times permuted end ---

# --- Perform asymmetric permutations start ---

for (i in 1:num_permutations) {
  permuted_samples <- sample(all_samples)
  permuted_matrix[i, ] <- permuted_samples
}

# --- Perform asymmetric permutations end ---

# --- Separate the permuted samples into case and control groups for each permutation start ---

permuted_groups <- list()
for (i in 1:num_permutations) {
  permuted_group1 <- permuted_matrix[i, 1:length(group1)]
  permuted_group2 <- permuted_matrix[i, (length(group1) + 1):(length(group1) + length(group2))]
  permuted_groups[[i]] <- list(group1=permuted_group1, group2=permuted_group2)
}

# --- Create matrices from permuted results to have two matrices including 1000 permutations for both group1 and group2! start ---

matrix_group1 <- matrix(nrow=num_permutations, ncol=length(group1))
matrix_group2 <- matrix(nrow=num_permutations, ncol=length(group2))

for (i in 1:num_permutations) {
  matrix_group1[i, ] <- permuted_groups[[i]]$group1
  matrix_group2[i, ] <- permuted_groups[[i]]$group2
}

# --- Remove symmetrical permutations from the permutation matrix of group1 start ---

# Initialize a logical vector to mark rows for removal

rows_to_remove <- logical(nrow(matrix_group1))

# Create a list to store unique row combinations
unique_combinations <- list()

# Identify and mark symmetric duplicate rows for removal
for (i in 1:nrow(matrix_group1)) {
  current_row <- matrix_group1[i, ]
  sorted_current_row <- sort(current_row)
  
  # Check if the current row's sorted combination is already seen
  if (paste(sorted_current_row, collapse = "_") %in% unique_combinations) {
    rows_to_remove[i] <- TRUE
  } else {
    unique_combinations[[i]] <- paste(sorted_current_row, collapse = "_")
  }
}

# Keep rows without symmetric duplicates

filtered_matrix_group1 <- matrix_group1[!rows_to_remove, ]

# Check how many permutations were kept

table(rows_to_remove)

original_order_of_group1 <- group1

filtered_matrix_group1 <- rbind(original_order_of_group1, filtered_matrix_group1)
rownames(filtered_matrix_group1) <- NULL


# Check whether all permutations in group1 are unique.


# Function to normalize each row by sorting its elements
normalize_row <- function(row) {
  return(paste(sort(row), collapse = " "))
}

# Apply the normalization function to each row of the matrix
normalized_rows <- apply(filtered_matrix_group1, 1, normalize_row)

# Identify and remove duplicate rows based on the normalized values
unique_indices <- !duplicated(normalized_rows)
unique_df <- filtered_matrix_group1[unique_indices, ]

# Check if the matrix has unique permutations
if (nrow(filtered_matrix_group1) == nrow(unique_df)) {
  print("All rows in the group 1 matrix are unique permutations.")
} else {
  print("The group 1 matrix has duplicate permutations. Duplicates have been removed.")
}

print("Number of unique permutations in group 1 is:")

print(nrow(unique_df))

# Create permutation for group 2

# Initialize a logical vector to mark rows for removal

rows_to_remove2 <- logical(nrow(matrix_group2))

# Create a list to store unique row combinations
unique_combinations <- list()

# Identify and mark symmetric duplicate rows for removal
for (i in 1:nrow(matrix_group2)) {
  current_row <- matrix_group2[i, ]
  sorted_current_row <- sort(current_row)
  
  # Check if the current row's sorted combination is already seen
  if (paste(sorted_current_row, collapse = "_") %in% unique_combinations) {
    rows_to_remove2[i] <- TRUE
  } else {
    unique_combinations[[i]] <- paste(sorted_current_row, collapse = "_")
  }
}

# Keep rows without symmetric duplicates
filtered_matrix_group2 <- matrix_group2[!rows_to_remove2, ]

# --- Add original order of groups to first row of permutation matrix start ---

original_order_of_group2 <- group2
filtered_matrix_group2 <- rbind(original_order_of_group2, filtered_matrix_group2)
rownames(filtered_matrix_group2) <- NULL

# Check whether all permutations in group2 are unique.

# Apply the normalization function to each row of the matrix
normalized_rows <- apply(filtered_matrix_group2, 1, normalize_row)

# Identify and remove duplicate rows based on the normalized values
unique_indices <- !duplicated(normalized_rows)
unique_df <- filtered_matrix_group2[unique_indices, ]

# Check if the matrix has unique permutations
if (nrow(filtered_matrix_group2) == nrow(unique_df)) {
  print("All rows in the group 2 matrix are unique permutations.")
} else {
  print("The group 2 matrix has duplicate permutations. Duplicates have been removed.")
}

print("Number of unique permutations in group 2 is:")

print(nrow(unique_df))


# Obtain NULL DMRs


NULL.DMRs <- getNullDistribution_BSmooth.tstat(BSseq = BSseq.obj, idxMatrix1 = filtered_matrix_group1,
                                                                  idxMatrix2 = filtered_matrix_group2,
                                                                  estimate.var = "same",
                                                                  local.correct = T,
                                                                  cutoff = c(-4.6, 4.6),
                                                                  stat = "tstat",
                                                                  mc.cores = 20,
                                                                  maxGap = 10^8)


# t-test statistics


t.test <- BSmooth.tstat(BSseq.obj,
                               group1 = group1,
                               group2 = group2,
                               estimate.var = "same",
                               local.correct = TRUE,
                               verbose = TRUE, k = 21, mc.cores = 24, maxGap = 10^8)

# Obtain putative DMRs                           

Putative.DMRs <- dmrFinder(t.test, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat", verbose = T)

# Check whether putative and NULL DMRs for original comparison (first permutation) are the same or not

if (identical(NULL.DMRs[[1]], Putative.DMRs)) {
  print("YES, first permutation in NULL DMRs and putative DMRs are the same")
} else {
  print("NO, first permutation in NULL DMRs and putative DMRs are not the same")
}                             

# Subset NULL DMRs

subset.NULL.DMRs <- subsetDmrs(NULL.DMRs)

# Subset putative DMRs

subset.putative.DMRs <- subset(Putative.DMRs, n >= 3 & abs(meanDiff) >= 0.05 & invdensity <= 300)

# Test whether subset NULL and putative DMRs are the same or not

if (identical(subset.NULL.DMRs[[1]], subset.putative.DMRs)) {
  print("YES, first permutation in NULL DMRs and putative DMRs are the same")
} else {
  print("NO, first permutation in NULL DMRs and putative DMRs are not the same")
}

# Obtain FWER from permutation where chrX and chrY should be omitted from all NULL DMRs except first permutation which is putative DMRs

FWER <- getFWER(subset.NULL.DMRs, type = "dmrs")

# Set cutoff for FWER for statistical significance; 0.05 * (number of permutations)

FWER.cutoff <- 0.05*100

# Get genome-wide significant DMRs

print("Number of genome-wide significant DMRs:")

table(FWER < FWER.cutoff)

FWER <- as.data.frame(FWER)

colnames(FWER) <- "FWER.integer.(number.of.permutations)"

FWER$genome.wide.significant.based.on.FWER.integer <- ifelse(FWER$`FWER.integer.(number.of.permutations)` < FWER.cutoff, "Yes", "No")

FWER$FWER.significance <- (FWER$FWER + 1)/(100 + 1)

# Add FWER to putative DMRs

putative.DMRs.with.FWER <- cbind(subset.putative.DMRs, FWER)

# plot all DMRs

pdf(file = "all.putative.DMRs.pdf", width = 10, height = 5)
plotManyRegions(BSseq.obj, putative.DMRs.with.FWER[], extend = 5000, addRegions = putative.DMRs.with.FWER)
dev.off()

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken