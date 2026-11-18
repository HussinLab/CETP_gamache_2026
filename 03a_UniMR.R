# =============================================================================
# 03a_UniMR.R
# Univariable MR (IVW, MR-Egger, weighted median): CETP expression, splicing
# (Jct3up) and sQTL CETP203 (Jct1up) -> outcomes from Data/MR and OpenGWAS.
#
# Requirement: OpenGWAS token in the OPENGWAS_JWT environment variable
#   Sys.setenv(OPENGWAS_JWT = "...")
# Run from the scripts folder (paths are relative: ../../Data).
# =============================================================================
source("03_MR_functions.R")

# --- Parameters ---------------------------------------------------------------
tissues      <- c("Spleen", "Thyroid")
exposure_ids <- c("CETP", "Splicing", "sQTLCETP203")
f_thresholds <- c(10, 15)
ld_max       <- 0.95     # |r| >= ld_max -> keep only one SNP of the pair
min_snps     <- 3
outcome_names <- outcomes$name

# --- Reference ----------------------------------------------------------------
ref   <- load_reference()
exons <- load_cetp_exons()

# --- One analysis = one tissue x one exposure x one F threshold -----------------------
run_univariate <- function(tissue, exposure, f_min, qtl) {
    out_dir <- dir_create(results_dir, "MR", exposure)
    tab_dir <- dir_create(out_dir, "MR_Table")
    pdf_dir <- dir_create(out_dir, "Nuage", paste0("Fstat", f_min))
    stem    <- sprintf("MR.%s.%s.SD.LD%s.Fstat%s", exposure, tissue, ld_max, f_min)

    pdf_open <- FALSE
    on.exit(if (pdf_open) dev.off())
    all_res <- list()
    log     <- list()

    for (outcome in outcome_names) {
        message("[", tissue, " / ", exposure, " / F>=", f_min, "] ", outcome)

        out <- get_outcome(outcome, ref)
        if (is.null(out)) {
            log[[outcome]] <- data.frame(Outcome = outcome, Status = "outcome has no data", nSNPs = NA)
            next
        }

        prep <- prepare_mr_data(setNames(list(qtl), exposure), out, ref, f_min, ld_max, min_snps)
        if (!prep$ok) {
            message("  -> skipped: ", prep$reason)
            log[[outcome]] <- data.frame(Outcome = outcome, Status = prep$reason, nSNPs = prep$n)
            next
        }

        res <- run_uni_mr(prep, exposure, outcome, label = paste0(exposure, "-", tissue))
        all_res[[outcome]] <- res
        log[[outcome]] <- data.frame(Outcome = outcome, Status = "ok", nSNPs = prep$n)

        # 4 panels: MR, p-value comparison, exposure locus, outcome locus
        if (!pdf_open) {
            pdf(file.path(pdf_dir, paste0(stem, ".pdf")), width = 14, height = 10)
            pdf_open <- TRUE
        }
        used <- prep$data$SNP
        print(
            (plot_scatter(prep$data, exposure, res, outcome) |
             plot_pvalue_comparison(prep$all, exposure, outcome, ref)) /
            (plot_locus(prep$all, used, paste0("pval.", exposure), exons, ref, exposure) |
             plot_locus(prep$all, used, "pval.outcome", exons, ref, outcome))
        )
    }

    log <- do.call(rbind, log)
    write.table(log, file.path(tab_dir, paste0(stem, ".Log.txt")),
                quote = FALSE, row.names = FALSE, sep = "\t")

    if (length(all_res) == 0) return(invisible(NULL))
    res <- do.call(rbind, all_res)
    write.table(res, file.path(tab_dir, paste0(stem, ".txt")),
                quote = FALSE, row.names = FALSE, sep = "\t")
    message("Written: ", file.path(tab_dir, paste0(stem, ".txt")))
    invisible(res)
}

# --- Main loop ----------------------------------------------------------
for (tissue in tissues) {
    for (exposure in exposure_ids) {
        qtl <- load_exposure(exposure, tissue, ref)
        if (is.null(qtl)) next
        for (f_min in f_thresholds) {
            run_univariate(tissue, exposure, f_min, qtl)
        }
    }
}
