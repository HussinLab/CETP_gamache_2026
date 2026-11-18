# =============================================================================
# 03b_MVMR.Splicing-CETP203.R
# Multivariable MR: splicing (Jct3up) + CETP expression + sQTL CETP203 (Jct1up)
# estimated jointly (IVW, MR-Egger, median) for each outcome.
#
# `first` = exposure that orients MR-Egger (first argument of mr_mvinput).
# The data (harmonisation, F-stat filter, LD) are prepared once per
# tissue and outcome, then reused for the three orders.
#
# Requirement: OpenGWAS token in OPENGWAS_JWT (Sys.setenv(OPENGWAS_JWT = "...")).
# =============================================================================
source("03_MR_functions.R")

# --- Parameters ---------------------------------------------------------------
tissues  <- c("Spleen", "Thyroid")
f_min    <- 10
ld_max   <- 0.95
set.seed(2026)           # the multivariable median uses a bootstrap

# Exposure order for each `first` (the first one orients MR-Egger)
orders <- list(
    Splicing = c("Splicing", "CETP", "sQTLCETP203"),
    CETP     = c("CETP", "Splicing", "sQTLCETP203"),
    s203     = c("sQTLCETP203", "Splicing", "CETP")
)
exposure_ids  <- orders[[1]]
min_snps      <- length(exposure_ids) + 2     # enough degrees of freedom for IVW and Egger
outcome_names <- outcomes$name

# --- Reference ----------------------------------------------------------------
ref <- load_reference()

out_dir  <- dir_create(results_dir, "MVMR/Splicing/CETP_CETP203")
pdf_dir  <- dir_create(out_dir, "Nuages", paste0("Fstat", f_min))
stem_for <- function(tissue) sprintf("MVMR.Splicing_%s.CETP_CETP203_%s.SD", tissue, tissue)

for (tissue in tissues) {

    # Exposures for this tissue
    exposures <- setNames(lapply(exposure_ids, load_exposure, tissue = tissue, ref = ref), exposure_ids)
    if (any(vapply(exposures, is.null, TRUE))) {
        message("[", tissue, "] one of the exposures is missing, tissue skipped")
        next
    }

    # --- Preparation (once per outcome) ---------------------------------
    prepared <- list()
    log <- list()
    for (outcome in outcome_names) {
        message("[", tissue, "] preparing: ", outcome)
        out <- get_outcome(outcome, ref)
        if (is.null(out)) {
            log[[outcome]] <- data.frame(Outcome = outcome, Status = "outcome has no data", nSNPs = NA)
            next
        }
        prep <- prepare_mr_data(exposures, out, ref, f_min, ld_max, min_snps)
        if (!prep$ok) {
            message("  -> skipped: ", prep$reason)
            log[[outcome]] <- data.frame(Outcome = outcome, Status = prep$reason, nSNPs = prep$n)
            next
        }
        prepared[[outcome]] <- prep
        log[[outcome]] <- data.frame(Outcome = outcome, Status = "ok", nSNPs = prep$n)
    }
    write.table(do.call(rbind, log),
                file.path(out_dir, sprintf("%s.LD%s.Fstat%s.Log.txt", stem_for(tissue), ld_max, f_min)),
                quote = FALSE, row.names = FALSE, sep = "\t")
    if (length(prepared) == 0) next

    # Retained instruments per outcome (PSI betas x 100, F per exposure)
    instruments <- do.call(rbind, lapply(names(prepared), function(o) {
        cbind(Outcome = o, prepared[[o]]$data)
    }))
    write.table(instruments,
                file.path(out_dir, sprintf("MVMR.Instruments.%s.LD%s.Fstat%s.txt", tissue, ld_max, f_min)),
                quote = FALSE, row.names = FALSE, sep = "\t")
    message("[", tissue, "] SNPs used (union): ", length(unique(instruments$SNP)))

    # --- MVMR for each order ---------------------------------------------
    for (first in names(orders)) {
        order_ids <- orders[[first]]
        pdf(file.path(pdf_dir, sprintf("%s.First%s.LD%s.Fstat%s.pdf", stem_for(tissue), first, ld_max, f_min)),
            width = 18, height = 6)

        all_res <- list()
        for (outcome in names(prepared)) {
            message("[", tissue, " / first=", first, "] MVMR : ", outcome)
            prep <- prepared[[outcome]]
            res  <- run_mv_mr(prep, order_ids, outcome, label = tissue)
            res  <- cbind(res[, 1, drop = FALSE], First = first, res[, -1])
            all_res[[outcome]] <- res

            plots <- lapply(order_ids, function(id) plot_scatter(prep$data, id, res, outcome))
            print(wrap_plots(plots, nrow = 1))
        }
        dev.off()

        all_res <- do.call(rbind, all_res)
        out_file <- file.path(out_dir, sprintf("%s.First%s.LD%s.Fstat%s.txt", stem_for(tissue), first, ld_max, f_min))
        write.table(all_res, out_file, quote = FALSE, row.names = FALSE, sep = "\t")
        message("Written: ", out_file)
    }
}
