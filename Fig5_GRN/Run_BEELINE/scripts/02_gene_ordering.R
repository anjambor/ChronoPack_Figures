#!/usr/bin/env Rscript
# BEELINE gene ordering (GeneOrdering.csv: VGAMpValue, Variance) for one r9 single-cell context (BEELINE_CONTEXT,
# default K562; prep files in output/prep/<context>/, or output/prep/<context>-cpfilt/ with BEELINE_GENE_FILTER=chronopack).
#
# BEELINE paper (Pratapa et al. 2020, Methods): "We used the general additive model implemented in the gam R
# package to compute the variance and the P value of this variance." The exact call is not published; we use the
# Slingshot workflow's gam test (BEELINE computes pseudotime with Slingshot):
#   gam(z ~ lo(t)), P = summary(fit)[4][[1]][1, 5]  (parametric ANOVA, lo(t) row, Pr(>F))
# with z = log1p expression of a gene across cells and t = time tag (1-12). Variance = var(z). The
# nonparametric-effect and overall-F P values are also written (not used). Bonferroni correction is applied later by
# generateExpInputs.py (-c).
#
# Input:  output/prep/<context>/ExpressionData.full.csv, PseudoTime.csv, gene_filter.csv
# Output: output/prep/<context>/GeneOrdering.csv (genes passing the gene filter, sorted by P)
suppressPackageStartupMessages({
  library(gam)
  library(data.table)
  library(parallel)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
context <- Sys.getenv("BEELINE_CONTEXT", "K562")
gene_filter <- Sys.getenv("BEELINE_GENE_FILTER", "beeline")
stopifnot(gene_filter %in% c("beeline", "chronopack"))
prep <- file.path(dirname(script_dir), "output", "prep", paste0(context, if (gene_filter == "chronopack") "-cpfilt" else ""))
stopifnot(dir.exists(prep))
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

expr <- fread(file.path(prep, "ExpressionData.full.csv"), header = TRUE)
genes <- expr[[1]]
X <- as.matrix(expr[, -1, with = FALSE])
rownames(X) <- genes
pt <- read.csv(file.path(prep, "PseudoTime.csv"), row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(X), rownames(pt)))
t <- pt$PseudoTime

filt <- read.csv(file.path(prep, "gene_filter.csv"))
use <- intersect(filt$gene[filt$kept == "True" | filt$kept == TRUE], genes)
cat(sprintf("%d genes x %d cells; testing %d genes on %d workers\n", nrow(X), ncol(X), length(use), workers))

# Check the summary() layout once: element 4 must be the parametric ANOVA with lo(t) in row 1 and Pr(>F) in col 5.
s <- summary(gam(z ~ lo(t), data = data.frame(z = X[use[1], ], t = t)))
stopifnot(names(s)[4] == "parametric.anova", rownames(s[[4]])[1] == "lo(t)", colnames(s[[4]])[5] == "Pr(>F)")

# Also recorded (not used for selection): the nonparametric-effect P (the other Slingshot-vignette variant,
# summary(fit)[3][[1]][2, 3]) and an overall F test of lo(t) against an intercept-only model.
test_gene <- function(g) {
  d <- data.frame(z = X[g, ], t = t)
  tryCatch({
    fit <- suppressWarnings(gam(z ~ lo(t), data = d))
    sm <- summary(fit)
    overall <- anova(gam(z ~ 1, data = d), fit, test = "F")
    c(p_param = sm[4][[1]][1, 5], p_npar = sm$anova["lo(t)", "Pr(F)"], p_overall = overall[2, "Pr(>F)"],
      variance = var(d$z))
  }, error = function(e) c(p_param = NA, p_npar = NA, p_overall = NA, variance = var(d$z)))
}
start <- Sys.time()
res <- mclapply(use, test_gene, mc.cores = workers)
bad <- !vapply(res, is.numeric, logical(1))
if (any(bad)) stop(sprintf("%d genes returned non-numeric results, e.g. %s", sum(bad), use[which(bad)[1]]))
res <- do.call(rbind, res)
rownames(res) <- use
cat(sprintf("GAM tests done in %.1f min; %d genes with failed fits (NA)\n",
            as.numeric(difftime(Sys.time(), start, units = "mins")), sum(is.na(res[, "p_param"]))))

out <- data.frame(VGAMpValue = res[, "p_param"], Variance = res[, "variance"], NparPValue = res[, "p_npar"],
                  OverallPValue = res[, "p_overall"], row.names = use)
out <- out[!is.na(out$VGAMpValue), ]
out <- out[order(out$VGAMpValue, -out$Variance), ]
write.csv(out, file.path(prep, "GeneOrdering.csv"), quote = FALSE)
bf <- 0.01 / nrow(out)
cat(sprintf("Bonferroni cutoff %.3g: %d of %d genes significant (nonparametric P: %d; overall F: %d)\n", bf,
            sum(out$VGAMpValue < bf), nrow(out), sum(out$NparPValue < bf, na.rm = TRUE), sum(out$OverallPValue < bf, na.rm = TRUE)))
cat("GENE_ORDERING_OK\n")
