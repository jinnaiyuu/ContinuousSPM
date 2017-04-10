## --------------------------------------- ##
## ---------- run three methods ---------- ##
## --------------------------------------- ##
run.script <- function(N = 200, n = 20, r0 = 0.5, data.name = "material", base.dir = "./", 
                       alpha = 0.05, bin.type = 1, 
                       run.prop = TRUE, run.bonf = TRUE, run.bin = TRUE, genonly = FALSE) {
  results <- matrix(rep(NA, 3 * 5), ncol = 3)
  rownames(results) <- c("#testable patterns", "Runtime", "Final freq. thr.", "Precision", "Recall")
  colnames(results) <- c("Proposal", "Bonferroni", "Median-bin")
#  real.data.list <- c("ctg", "faults", "ionosphere", "segment", "waveform", "wdbc")
#  real.data.list <- c("material") # YuuJinnai
  
  file.path.data <- paste0(base.dir, data.name, ".data")
  file.path.class <- paste0(base.dir, data.name, ".class")

  cat(" Using material data")
#  data.name <- "material"
  
  # if (data.name == "synth") {
  #   ## make a synthetic dataset with (N, n, r0)
  #   cat("  > Generate synthetic data\n")
  #   data.name <- paste0(data.name, "_N=", N, "_n=", n, "_r0=", r0, "_alpha=", alpha)
  #   generate.data(file.path.data, N, n)
  #   generate.class(file.path.class, N, r0)
  # } else if (data.name %in% real.data.list) {
  #   command <- paste0("real.data <- ", data.name, "(base.dir)")
  #   eval(parse(text = command))
  #   write.table(real.data$data, file = file.path.data, sep = ",", row.names = FALSE, col.names = FALSE)
  #   write.table(real.data$class, file = file.path.class, sep = ",", row.names = FALSE, col.names = FALSE)
  # } else {
  #   cat("ERROR: please specify a valid data name from:\n")
  #   cat("{synth, ")
  #   cat(real.data.list, sep = ", ")
  #   cat("}\n")
  #   return()
  # }
  
  if (genonly){
    return()
  }

  ## perform the proposed method
  if (run.prop) {
    cat("  > Perform the proposed mehtod\n")
    file.path.prop.sig <- paste0(base.dir, data.name, "_prop.sig")
    file.path.prop.stat <- paste0(base.dir, data.name, "_prop.stat")
    command <- paste("../cc/pmc -i", file.path.data, "-c", file.path.class, 
                     "-o", file.path.prop.sig, "-t", file.path.prop.stat, 
                     "-a", alpha)
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    res <- read.stat(file.path.prop.stat)
    results[1, 1] <- res[4]
    results[2, 1] <- res[7] + res[9]
    results[3, 1] <- res[6]
    cat("    Compute precision and recall\n")
    res <- compute.pr(file.path.prop.sig, n)
    results[4, 1] <- res[1]
    results[5, 1] <- res[2]
  }

  ## perform the Bonferroni correction method
  if (run.bonf) {
    cat("  > Perform the Bonferroni correction mehtod\n")
    file.path.bonf.sig <- paste0(base.dir, data.name, "_bonf.sig")
    file.path.bonf.stat <- paste0(base.dir, data.name, "_bonf.stat")
    command <- paste("../cc/pmc -i", file.path.data, "-c", file.path.class, "-o", file.path.bonf.sig, "-t", file.path.bonf.stat, "-a", alpha, "-b")
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    res <- read.stat(file.path.bonf.stat)
    results[1, 2] <- res[4]
    results[2, 2] <- res[7] + res[9]
    results[3, 2] <- 0
    ## res <- compute.pr(file.path.bonf.sig, n)
    ## we do not compute precision and recall as the FWER contoroll level is different
  }

  ## perform the median-binarization method
  if (run.bin) {
    if (bin.type == 2) {
      bin.method <- "interordinal"
    } else if (bin.type == 3) {
      bin.method <- "interval"
    } else {
      bin.method <- "median"
    }
    cat("  > Perform", bin.method, "discretization + significant itemset mining\n")
    cat("    Binarize data\n")
    file.path.bin.data <- paste0(base.dir, data.name, "_bin_", bin.method, ".data")
    file.path.bin.header <- paste0(base.dir, data.name, "_bin_", bin.method, ".header")
    cat("    Perform significant itemset mining to count the number of testable patterns\n")
    command <- paste("python ../python/binarize.py", file.path.data, file.path.bin.data, file.path.bin.header, bin.type)
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    enum.path = "../c/lcm_lamp_fisher/fim_closed"
    comp.path = "../c/lcm_comp_pvalues_fisher/fim_closed"
    file.path.lamp.stat <- paste0(base.dir, data.name, "_bin_", bin.method, "_lamp")
    cat("    Enumerate significant patterns\n")
    command <- paste(enum.path, file.path.lamp.stat, alpha, file.path.class, file.path.bin.data)
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    res <- readLines(paste0(file.path.lamp.stat, "_results.txt"))
    x <- unlist(strsplit(res[2], " "))
    sig.th <- x[length(x)]
    x <- unlist(strsplit(res[3], " "))
    lcm.th <- x[length(x)]
    x <- unlist(strsplit(res[6], " "))
    results[1, 3] <- as.numeric(x[length(x)])
    results[3, 3] <- as.numeric(lcm.th)
    file.path.lamp.sig <- paste0(base.dir, data.name, "_bin_", bin.method, "_lamp_sig")
    command <- paste(comp.path, file.path.lamp.sig, sig.th, lcm.th, file.path.class, file.path.bin.data)
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    file.path.lamp.sig.sig <- paste0(base.dir, data.name, "_bin_", bin.method, "_lamp_sig_sig_itemsets.txt")
    cat("    Compute precision and recall\n")
    res <- compute.pr(file.path.lamp.sig.sig, n)
    results[4, 3] <- res[1]
    results[5, 3] <- res[2]
  }

  res1 <- readLines(paste0(file.path.lamp.stat, "_timing.txt"))
  res2 <- readLines(paste0(file.path.lamp.stat, "_sig_timing.txt"))
  x <- unlist(strsplit(res1[5], " "))
  results[2, 3] <- as.numeric(x[length(x) - 1])
  x <- unlist(strsplit(res2[5], " "))
  results[2, 3] <- results[2, 3] + as.numeric(x[length(x) - 1])
  results
}



## ----------------------------------------------- ##
## ---------- synthetic data generation ---------- ##
## ----------------------------------------------- ##
generate.data <- function(file.path, N = 1000, n = 100, seed = 1) {
  set.seed(seed)
  X <- matrix(runif(N * n), ncol = n)
  X[, 1] <- sort(X[, 1], decreasing = TRUE)
  if (n > 10) {
    for (i in 2:(n * 0.2)) {
      X[, i] <- X[, 1] + rnorm(100, 0, 0.2)
    }
  }
  write.table(X, file = file.path, row.names = FALSE, col.names = FALSE, sep = ",")
}
generate.class <- function(file.path, N = 1000, ratio = 0.5) {
  if (ratio > 0.5) stop("ratio is too large!")
  N.small <- N * ratio
  cl <- c(rep(0, N.small), rep(1, N - N.small))
  write.table(cl, file = file.path, row.names = FALSE, col.names = FALSE)
}

## ---------------------------------------------- ##
## ---------- read results from a file ---------- ##
## ---------------------------------------------- ##
read.stat <- function(file.path.stat) {
  stat.list <- c("N", "n", "N.small", "num.testable", "corrected.alpha", "freq.threshold", "runtime.enum", "num.significant", "runtime.sig")
  res <- numeric(length(stat.list))
  names(res) <- stat.list

  if (file.access(file.path.stat) != 0) {
    stop(paste0("The file ", file.path.stat, " does not exist!"))
  }

  X <- readLines(file.path.stat)
  idx.lines <- c(3, 5, 4, 7, 8, 9, 10, 12, 13)
  for (i in 1:length(idx.lines)) {
    x <- unlist(strsplit(X[idx.lines[i]], " "))
    if (i == 7 || i == 9) {
      res[i] <- as.numeric(x[length(x) - 1])
    } else {
      res[i] <- as.numeric(x[length(x)])
    }
  }
  res
}

## --------------------------------------------------- ##
## ---------- compute precision and recall  ---------- ##
## --------------------------------------------------- ##
compute.pr <- function(file.path.sig, n) {
  ft <- 0:((n * 0.2) - 1)
  ff <- (n * 0.2):n
  true.positive <- 0

  if (file.access(file.path.sig) != 0) {
    stop(paste0("The file ", file.path.sig, " does not exist!"))
  }

  X <- readLines(file.path.sig)

  for (i in 1:length(X)) {
    x <- unlist(strsplit(X[i], " "))
    x <- as.numeric(x[1:(length(x) - 1)])
    true.positive <- true.positive + prod(x %in% ft)
  }
  if (length(X) == 0) true.positive <- 0
  precision <- true.positive / length(X)
  recall <- true.positive / 2^(n * 0.2)
  res <- c(precision = precision, recall = recall)
  res
}


## ----------------------------------------------------------- ##
## ---------- functions for preprocessing real data ---------- ##
## ----------------------------------------------------------- ##
imputation <- function(X) {
  for (j in 1:(ncol(X) - 1))
    X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
  X
}
check.class.minor <- function(cl) {
  if (sum(cl) < length(cl) / 2)
    cl <- ifelse(cl == 0, 1, 0)
  cl
}

ctg <- function(base.dir) {
  data.name <- "ctg"
  file.path.original <- paste0(base.dir, data.name, "_original.data")
  if (file.access(file.path.original) != 0) {
    cat(paste0("The file ", file.path.original, " does not exist!\n"))
    cat("Please download from the original file from:\n")
    cat("[http://archive.ics.uci.edu/ml/machine-learning-databases/00193/CTG.xls]\n")
    cat("and convert it to a csv file, then save the csv file as \"", file.path.original, "\"\n", sep = "")
    stop()
  }

  X <- read.table(file.path.original, sep = ",", header = TRUE)
  data <- X[, 1:(ncol(X) - 1)]
  data <- imputation(data)
  cl <- X[, ncol(X)]
  cl <- unlist(cl)
  cl[cl == 2] <- 0
  cl[cl == 3] <- 0
  cl <- check.class.minor(cl)
  res <- list(data = data, class = cl)
  res
}

faults <- function(base.dir) {
  data.name <- "faults"
  file.path.original <- paste0(base.dir, data.name, "_original.data")
  url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/00198/Faults.NNA"
  cat("  > download the original data of", data.name, "\n")
  download.file(url, file.path.original)

  X <- read.table(file.path.original, sep = "\t")
  data <- X[, 1:27]
  data <- imputation(data)
  cl <- X[, 28] ## use Pastry
  cl <- check.class.minor(cl)
  res <- list(data = data, class = cl)
  res
}

ionosphere <- function(base.dir) {
  data.name <- "ionosphere"
  file.path.original <- paste0(base.dir, data.name, "_original.data")
  url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/ionosphere.data"
  cat("  > download the original data of", data.name, "\n")
  download.file(url, file.path.original)

  X <- read.table(file.path.original, comment.char = "@", sep = ",", na.strings = "?")
  X <- imputation(X)
  data <- X[, 1:(ncol(X) - 1)]
  cl <- X[, ncol(X)]
  levels(cl) <- c(0, 1)
  cl <- as.numeric(as.character(cl))
  cl <- check.class.minor(cl)
  res <- list(data = data, class = cl)
  res
}

segment <- function(base.dir) {
  data.name <- "segment"
  file.path.original <- paste0(base.dir, data.name, "_original.data")
  url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/statlog/segment/segment.dat"
  cat("  > download the original data of", data.name, "\n")
  download.file(url, file.path.original)

  X <- read.table(file.path.original, sep = " ")
  data <- X[, 1:(ncol(X) - 1)]
  data <- imputation(data)
  cl <- X[, ncol(X)]
  cl[cl == 1] <- 0 ## use "1" and others
  cl[cl > 1] <- 1 ## use "1" and others
  cl <- check.class.minor(cl)
  res <- list(data = data, class = cl)
  res
}

waveform <- function(base.dir) {
  data.name <- "waveform"
  file.path.original <- paste0(base.dir, data.name, "_original.data")
  file.path.original.compressed <- paste0(file.path.original, ".Z")
  url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/waveform/waveform.data.Z"
  cat("  > download the original data of", data.name, "\n")
  download.file(url, file.path.original.compressed)
  cat("  > uncompress", data.name, "\n")
  system(paste("uncompress", file.path.original.compressed))

  X <- read.table(file.path.original, comment.char = "@", sep = ",", na.strings = "?")
  X <- imputation(X)
  data <- X[, 1:(ncol(X) - 1)]
  cl <- X[, ncol(X)]

  data <- data[cl != 1,]
  cl <- cl[cl != 1]
  cl[cl == 2] <- 1

  cl <- check.class.minor(cl)
  res <- list(data = data, class = cl)
  res
}

wdbc <- function(base.dir) {
  data.name <- "wdbc"
  file.path.original <- paste0(base.dir, data.name, "_original.data")
  url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
  cat("  > download the original data of", data.name, "\n")
  download.file(url, file.path.original)

  X <- read.table(file.path.original, comment.char = "@", sep = ",", na.strings = "?")
  data <- X[, 3:ncol(X)]
  data <- imputation(data)
  cl <- X[, 2]
  levels(cl) <- c(0, 1)
  cl <- as.numeric(as.character(cl))
  cl <- check.class.minor(cl)
  res <- list(data = data, class = cl)
  res
}
