#!/usr/bin/env Rscript
#
# Estimate strain frequencies given a marker vcf and a sample vcf
# based on functions from figs.Rmd in Mazzucco et al. (2020) supplement
# using only quadratic programming (for speed reasons)
#
# Rupert Mazzucco, 2024
# Modified by Prerna Goel, 2026

version <- "0.0.9001"

library("docopt")

docstr <- "
  Usage: predef_win_haploFreq.R -i <samples.vcf> -m <markers.vcf> -o <out.rds> -p <pre_win.rds> [options]
    predef_win_haploFreq.R -h
    predef_win_haploFreq.R --version

  Calculate marker line frequencies from observed allele frequencies in
  predefined windows for samples in <samples.vcf> using the marker frequencies
  given in <markers.vcf>. The latter must be normalized, i.e., contain only
  one ALT allele per line, and alleles must be a subset of the alleles in <samples.vcf>.
  A typical workflow might look like this:

  $ bcftools merge --no-index --merge both -Ov pileup1.bcf pileup2.bcf ... |\\
  +   post-merging.awk |\\
  +   bcftools filter -i 'TYPE = \"snp\" & REF != \"N\"' -Ob - > samples.bcf

  where the post-merging.awk script annotates the allele frequency and
  coverage depth tags needed here. The marker file can be extracted,
  filtered and normalized from there, e.g.,

  $ bcftools view --samples line1,line2 samples.bcf |\\
  +   bcftools norm -m- -Ou - |\\
  +   bcftools filter -i 'MAX(XF) - MIN(XF) > 0.9' -Ob > markers.bcf

  Options:
    --af-tag <AF>                 VCF format tag to read as allele frequency [default: XF]
    --dp-tag <DP>                 VCF format tag to read as coverage depth [default: SAD]
    -i --in <samples.vcf>         SNP frequencies in a number of samples
    -m --markers <markers.vcf>	  SNP frequencies in the marker lines. The VCF must be normalized
                                 (i.e., only one ALT allele per line)
    -n --numCores <number>        Number of cores to use in outer loop [default: 1]
    -o --out <out.rds>            Serialized output data table
    -p --pre <pre_win.rds>        Predefined windows
    -h --help                     Print this help text
    --version                     Print the script version number
"

args <- docopt(docstr, version = version, strict = TRUE)

argAfTag <- args[["--af-tag"]]
argDpTag <- args[["--dp-tag"]]
argVcf <- args[["--in"]]
argDrp <- "none"
argMks <- args[["--markers"]]
argNcr <- as.integer(args[["--numCores"]])
argRds <- args[["--out"]]
argPre<- args[["--pre"]]

suppressPackageStartupMessages({
  library("data.table")
  library("matrixStats")
  library("foreach")
  library("doParallel")
  library("MASS")
  library("Matrix")
  library("quadprog")
  library("dplyr")
  library("tidyr")
  library("fuzzyjoin")
})

# solve overdetermined AX = B for unknown strain frequencies X 
# possibly with weights W, i.e., 
# in a least-squares sense using quadratic programming
solve.qp <- function(A, B, W = NULL) {
  # make sure we are working with matrices
  A <- as.matrix(A)
  B <- as.matrix(B)
  # use equal weights if none provided
  W <- (if (is.null(W)) {
    matrix(1, nrow = nrow(B), ncol = ncol(B))
  } else {
    as.matrix(W)
  })
  # make sure weights are normalized
  W <- scale(W, center = FALSE, scale =  colSums(W))
  # sanity checks
  stopifnot(
    nrow(A) >= ncol(A),
    nrow(B) == nrow(A),
    all(dim(W) == dim(B)),
    all(W >= 0)
  )
  
  # Our constraints are implemented as one equality constraint, sum(x)=1, and
  # m=ncol(A) inequality constraints, all(x>=0); expressed in the (m+1) by m
  # constraints matrix t(C) = [1..1;1 0..0; 0 1 0..0; ...; 0..0 1] and rhs
  # y = (1 0..0), with meq=1 indicating that the first is an equality constraint
  C <- cbind(rep(1, ncol(A)), diag(ncol(A)))
  y <- c(1, rep(0, ncol(A)))
  # handle columns of B independently
  foreach(j = seq_len(ncol(B)), .combine = cbind, .multicombine = TRUE) %do% {
    # Minimizing the weighted sum of squared residuals t(r)Wr with r=Ax-b and W=diag(w) can be cast
    # in the canonical form of a quadratic programming problem
    # x = argmin{ t(x)Dx/2-t(d)x } with D=t(A)WA and d=t(A)Wb,
    V <- diag(W[, j])
    D <- t(A) %*% V %*% A
    # not sure why, but apparently D can in obscure cases be not positive-definite
    # (pd means that t(x) D x > 0 for all x)
    # in this case, we use Matrix::nearPD to substitute a pd matrix, emit a warning, and hope for the best
    # rather than catch this (could check if all eigenvalues are positive?), wrap in try() below
    d <- as.vector(t(A) %*% V %*% B[, j, drop = FALSE])
    r <- try(quadprog::solve.QP(D, d, C, y, meq = 1L)[["solution"]])
    if (class(r) == "try-error") {
      warning("quadprog: D probably not positive-definite, forcing nearest ...")
      D <- Matrix::nearPD(D, base.matrix = TRUE)[["mat"]]
      quadprog::solve.QP(D, d, C, y, meq = 1L)[["solution"]]
    } else {
      r
    }
  }
}

# estimate marker frequencies
estimate.frequency <- function(chr, win, sd, lines, samples) {
  require(matrixStats)
  # design matrix (ALT freqs in lines)
  A <- as.matrix(sd[, ..lines])
  # rownames(A) <- .SD[,paste(CHROM,POS,sep=":")]
  # colnames(A) <- lines
  # observed frequencies
  B <- as.matrix(sd[, ..samples])
  # rownames(B) <- rownames(A)
  # colnames(B) <- samples
  # use coverage as weights
  # (the "scaling" below multiplies them with the sum of line coverage depths)
  W <- t(scale(
    t(as.matrix(sd[, paste("i", ..samples, sep = ".")])),
    center = FALSE,
    scale = 1 / rowSums(sd[, paste("i", ..lines, sep = ".")])
  ))
  # disallow actual 0 weights (resulting from occasional DP==0 in samples)
  small.val <- sqrt(.Machine$double.eps)
  W[W < small.val] <- small.val
  # rownames(W) <- rownames(B)
  # colnames(W) <- colnames(B)
  # solve with quadratic programming
  X <- solve.qp(A, B, W)
  # fix numeric zeros
  if (any(X < 0)) {
    warning("Fixing numerical 0s ...")
    X[X < 0] <- 0
    # and rescale to re-enforce constraints
    X <- scale(X, center = FALSE, scale = colSums(X))
  }
  # construct output data table
  mpos <- as.integer(sd[, median(POS)])
  wname <- sd[, paste0(chr, ":", paste(range(POS), collapse = "-"))]
  setcolorder(
    setnames(
      data.table(
        sample = samples,
        t(X)
      )[
        , window := ..wname
      ][
        , pos := ..mpos
      ],
      paste0("V", seq_along(lines)),
      lines
    ),
    c("window", "pos")
  )
}

# main script starts here

message("Reading marker VCF ...")
colNames <- unname(unlist(fread(cmd = paste("/opt/homebrew/bin/bcftools view -h",
                                            argMks,
                                            "| tail -1 | cut -f1,2,4,5,10- | sed 's/^#//'"),
                                header = FALSE)))

if (length(lines <- colNames[-(1:4)]) < 2L) {
  stop("Procedure makes no sense with only one base line in ", argMks, ", exiting.")
}

dt.af.mks <- fread(
  cmd = paste0("/opt/homebrew/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%", argAfTag, "]\n' ", argMks),
  header = FALSE, col.names = colNames,
  key = c("CHROM", "POS", "REF", "ALT")
)

if (nrow(dt.af.mks) == 0L) {
  stop(argMks, " appears to be empty, nothing to be done, exiting.")
}
dt.dp.mks <- fread(
  cmd = paste0("/opt/homebrew/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%", argDpTag, "]\n' ", argMks),
  header = FALSE, col.names = colNames,
  key = c("CHROM", "POS", "REF", "ALT")
)

message("Reading input VCF ...")
colNames <- unname(unlist(fread(cmd = paste("/opt/homebrew/bin/bcftools view -h",
                                            argVcf,
                                            "| tail -1 | cut -f1,2,4,5,10- | sed 's/^#//'"),
                                header = FALSE)))

dropCols <- NULL
if (argDrp != "none") {
  dropCols <- 4L + as.integer(scan(textConnection(argDrp), sep = ","))
  if (is.integer(dropCols)) {
    colNames <- colNames[-dropCols]
  } else {
    message("Cannot understand argument --drop ", argDrp, ", ignoring ...")
  }
}

dt.af.smp <- fread(
  cmd = paste0("/opt/homebrew/bin/bcftools norm -m- ", argVcf, " | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%", argAfTag, "]\n' -"),
  header = FALSE, drop = dropCols, col.names = colNames,
  key = c("CHROM", "POS", "REF", "ALT")
)
dt.dp.smp <- fread(
  cmd = paste0("/opt/homebrew/bin/bcftools norm -m- ", argVcf, " | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%", argDpTag, "]\n' -"),
  header = FALSE, drop = dropCols, col.names = colNames,
  key = c("CHROM", "POS", "REF", "ALT")
)
samples <- colNames[-(1:4)]

message("Merging AF and DP info ...")
aux <- dt.af.smp[dt.af.mks][
  dt.dp.smp[dt.dp.mks]]

message("Removing sites that are missing in one or more samples ...")
keep <- !matrixStats::rowAnys(is.na(aux))
# remove them
message("Keeping ", sum(keep), " of ", nrow(aux), " sites, i.e., ~", round(100*sum(keep)/nrow(aux),2), "%")
aux <- aux[keep]

#Loading window information
pre_win<- readRDS(file=argPre)

window<- unique(pre_win$new_win)
window<- as.data.frame(window)

#Extracting the positions of the start and end of the window
window<- window%>%separate_wider_delim(window, ":", names=c("chr", "pos"))
window<- window%>%separate_wider_delim(pos, "-", names=c("start", "end"))

window$start<- as.numeric(window$start)
window$end<- as.numeric(window$end)

#Dividing the aux into windows 
aux$POS<- as.numeric(aux$POS)
aux_new <- aux |> 
  fuzzy_inner_join(
    window,
    by = c("CHROM"="chr","POS" = "start", "POS" = "end"),
    match_fun = list(`==`,`>=`, `<=`))

aux_new$win<- paste(aux_new$start, aux_new$end, sep="-")

# convert to data.table
aux_new <-  as.data.table(aux_new, key = c("CHROM"))

#Estimating Frequencies
dt.freq <- aux_new[, estimate.frequency(chr, win, .SD, lines, samples), by = .(chr = CHROM, win)]
saveRDS(dt.freq, file = argRds)
