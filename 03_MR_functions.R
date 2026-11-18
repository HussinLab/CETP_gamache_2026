# =============================================================================
# 03_MR_functions.R
# Shared functions and parameters for 03a_UniMR.R and 03b_MVMR.Splicing-CETP203.R
#
# Contents:
#   1. Parameters (paths, outcomes, exposures)
#   2. Reference: SNP map (Hg19/Hg38/rs) + GTEx LD matrix
#   3. Exposures (QTL + F-stat)
#   4. Outcomes (OpenGWAS or local files)
#   5. Harmonisation, F-stat / LD filters
#   6. MR results and figures
#
# The OpenGWAS token is not stored in the scripts:
#   Sys.setenv(OPENGWAS_JWT = "...")   or   export OPENGWAS_JWT=... in the shell
# =============================================================================

suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(MendelianRandomization)
    library(ieugwasr)
    library(ggplot2)
    library(patchwork)
    library(data.table)
})

freadf <- function(...) as.data.frame(fread(...))

# =============================================================================
# 1. Parameters
# =============================================================================
token <- Sys.getenv("OPENGWAS_JWT")

data_dir    <- "../../Data"
results_dir <- "../../Results"
bfile       <- file.path(data_dir, "vcf/GTEx.Hg38.ADCY9_CETP.MAF01.NoIndel.White.005")

# --- Outcomes ----------------------------------------------------------------
# source = "ieu"   -> extracted with extract_outcome_data(id)
# source = "local" -> read from Data/MR (see read_local_outcome)
outcomes <- data.frame(
    name   = c("HDL cholesterol", "LDL cholesterol", "TG",
               "Trunk fat-free mass", "Whole body fat-free mass", "Height",
               "Basal metabolic rate", "TSH", "Birth weight",
               "Birth weight of first child", "FEV1", "FVC", "CAD",
               "Early age-related macular degeneration", "Hypothyroidism",
               "Hyperthyroidism", "Miscarriage stillbirth"),
    source = c("ieu", "ieu", "ieu",
               "ieu", "ieu", "local",
               "ieu", "ieu", "ieu",
               "ieu", "ieu", "ieu", "local",
               "ieu", "local",
               "local", "local"),
    id     = c("ieu-b-109", "ieu-b-110", "ieu-b-111",
               "ukb-b-17409", "ukb-b-13354", NA,
               "ukb-b-16446", "prot-a-530", "ukb-b-13378",
               "ukb-b-3357", "ukb-b-19657", "ukb-b-7953", NA,
               "ebi-a-GCST010723", NA,
               NA, NA),
    stringsAsFactors = FALSE
)

# --- Expositions -------------------------------------------------------------
# type     : "splicing" (Splicing/Screening) or "expression" (Expression/GTEx_ScreeningeQTL/CETP)
# junction : junction kept from the splicing file
# flip     : Jct1up = 1 - Jct1down  ->  beta up = -1 * beta down
# psi      : PSI betas are multiplied by 100 (effect per 1% PSI)
# short    : label used in the MVMR tables
exposure_specs <- list(
    CETP        = list(type = "expression", junction = NA,      flip = FALSE, psi = FALSE, short = "CETP"),
    Splicing    = list(type = "splicing",   junction = "Jct3up", flip = FALSE, psi = TRUE,  short = "Splicing"),
    sQTLCETP203 = list(type = "splicing",   junction = "Jct1up", flip = TRUE,  psi = TRUE,  short = "sQTL203")
)

# =============================================================================
# 2. Reference: SNP map + LD matrix
# =============================================================================
load_reference <- function() {
    map <- freadf(file.path(data_dir, "ForLiftOverHg19_Hg38.txt"))
    bim <- freadf(paste0(bfile, ".bim"))[, c("V4", "V2", "V6", "V5")]
    colnames(bim) <- c("Hg38", "ID_GTEx", "other_allele.exposure", "effect_allele.exposure")

    map <- merge(map, bim, by = "Hg38")
    map <- map[, c("Hg19", "ID_GTEx", "Hg38", "ID", "other_allele.exposure", "effect_allele.exposure")]
    colnames(map)[c(1, 3)] <- c("POS_Hg19", "POS_Hg38")
    map <- map[!duplicated(map$ID_GTEx), ]

    corr <- ld_matrix(map$ID_GTEx, with_alleles = FALSE, bfile = bfile,
                      plink_bin = genetics.binaRies::get_plink_binary())

    # ld_matrix can drop or reorder SNPs: rename by matching the GTEx IDs,
    # never by position.
    rs <- map$ID[match(rownames(corr), map$ID_GTEx)]
    if (anyNA(rs) || anyDuplicated(rs)) stop("LD matrix IDs not found in the SNP map")
    rownames(corr) <- colnames(corr) <- rs

    list(map = map, corr = corr, mid = sum(range(map$POS_Hg38)) / 2)
}

load_cetp_exons <- function() {
    exons <- read.table(file.path(data_dir, "CoordinateCETP.txt"), header = TRUE, fill = TRUE)
    exons[exons$Region == "Exon", ]
}

# =============================================================================
# 3. Exposures (QTL + F-stat), SNPs as rsIDs
# =============================================================================
load_exposure <- function(id, tissue, ref) {
    spec <- exposure_specs[[id]]
    if (spec$type == "splicing") {
        f_qtl   <- sprintf("%s/Splicing/Screening/GTEx.%s.Splicing.MaxPEER.LM.ForMR.txt", results_dir, tissue)
        f_fstat <- sprintf("%s/Splicing/Screening/GTEx.%s.Splicing.Fstat.CompareFullModel.txt", results_dir, tissue)
    } else {
        f_qtl   <- sprintf("%s/Expression/GTEx_ScreeningeQTL/CETP/GTEx.%s.CETP.MaxPEER.LM.ForMR.SD.txt", results_dir, tissue)
        f_fstat <- sprintf("%s/Expression/GTEx_ScreeningeQTL/CETP/GTEx.%s.CETP.Fstat.SD.CompareFullModel.txt", results_dir, tissue)
    }
    if (!all(file.exists(f_qtl, f_fstat))) {
        message("  [", id, " / ", tissue, "] QTL or F-stat file missing")
        return(NULL)
    }

    qtl   <- read.table(f_qtl, header = TRUE)
    fstat <- read.table(f_fstat, sep = "\t", header = TRUE)
    if (spec$type == "splicing") {
        qtl   <- qtl[qtl[, 2] == spec$junction, ]
        fstat <- fstat[fstat$Transcript == spec$junction, ]
    }
    if (nrow(qtl) == 0 || nrow(fstat) == 0) {
        message("  [", id, " / ", tissue, "] no data for ", id)
        return(NULL)
    }

    if (spec$flip) qtl$beta.exposure <- -qtl$beta.exposure
    if (spec$psi) {
        qtl$beta.exposure <- qtl$beta.exposure * 100
        qtl$se.exposure   <- qtl$se.exposure * 100
    }

    # GTEx name -> rsID, by Hg38 position
    pos2rs <- ref$map[, c("POS_Hg38", "ID")]
    qtl <- merge(qtl, pos2rs, by.x = "Position", by.y = "POS_Hg38", sort = FALSE)
    colnames(qtl)[colnames(qtl) == "ID"] <- "SNP"
    qtl$exposure    <- id
    qtl$id.exposure <- id
    qtl$effect_allele.exposure <- as.character(qtl$effect_allele.exposure)
    qtl$other_allele.exposure  <- as.character(qtl$other_allele.exposure)

    fstat <- merge(fstat, pos2rs, by.x = "Position", by.y = "POS_Hg38")
    colnames(fstat)[colnames(fstat) == "ID"] <- "SNP"

    list(qtl = qtl, fstat = fstat[, c("SNP", "F.stat")])
}

# =============================================================================
# 4. Outcomes
# =============================================================================
std_outcome <- function(d, name) {
    d$effect_allele.outcome <- toupper(as.character(d$effect_allele.outcome))
    d$other_allele.outcome  <- toupper(as.character(d$other_allele.outcome))
    d$id.outcome <- name
    d$outcome    <- name
    d
}

# Local outcomes (chr16). Columns are renamed by position when the file has
# no usable names, as in the original scripts.
read_local_outcome <- function(name, ref) {
    dir <- file.path(data_dir, "MR")

    if (name == "CAD") {
        d <- freadf(file.path(dir, "CAD_META.Chr16"), header = TRUE)
        colnames(d) <- c("MarkerName", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
                         "eaf.se.outcome", "min.freq.outcome", "max.freq.outcome", "beta.outcome",
                         "se.outcome", "pval.outcome", "Direction", "HetISq", "HetChiSq", "HetDf",
                         "HetPVal", "SNP", "CHR", "BP")

    } else if (name %in% c("Hypothyroidism", "Hyperthyroidism")) {
        f <- if (name == "Hypothyroidism") "Hypothyroidism.GWASatlas.Chr16.txt" else "Hyperthyroidism.GWASaltas.Chr16.txt"
        d <- freadf(file.path(dir, f), header = TRUE, sep = ",")
        colnames(d)[1:9] <- c("MarkerName", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
                              "beta.outcome", "se.outcome", "pval.outcome", "N", "I2")
        # Hg19 positions are in MarkerName (chr:pos:...) -> rsID by position
        d$BP <- as.integer(vapply(strsplit(d$MarkerName, ":", fixed = TRUE), `[`, "", 2))
        d$SNP <- ref$map$ID[match(d$BP, ref$map$POS_Hg19)]
        d <- d[!is.na(d$SNP), ]

    } else if (name == "Height") {
        d <- freadf(file.path(dir, "Height.UKB2.Chr16.txt"))
        d <- d[, c("SNPID_UKB", "A1", "A2", "BETA", "SE", "P", "MAF_UKB")]
        colnames(d) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome",
                         "se.outcome", "pval.outcome", "eaf.outcome")

    } else if (name == "Miscarriage stillbirth") {
        d <- freadf(file.path(dir, "PregnancyComplication.Chr16"))
        d <- d[, c("pos", "rsids", "alt", "ref", "af", "beta", "sebeta", "pval")]
        colnames(d) <- c("BP", "SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
                         "beta.outcome", "se.outcome", "pval.outcome")

    } else {
        stop("Local outcome not defined: ", name)
    }

    std_outcome(d, name)
}

# Full outcome for the region (all SNPs in the map), cached:
# one OpenGWAS call per outcome, whatever the tissue, exposure or F threshold.
# To refresh, delete the Results/MR/Temp_outcomes folder.
get_outcome <- function(name, ref) {
    row <- outcomes[outcomes$name == name, ]
    if (nrow(row) != 1) stop("Unknown outcome: ", name)
    if (row$source == "local") return(read_local_outcome(name, ref))

    cache_dir <- file.path(results_dir, "MR/Temp_outcomes")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache <- file.path(cache_dir, paste0(row$id, ".rds"))
    if (file.exists(cache)) return(readRDS(cache))

    if (token == "") stop("OPENGWAS_JWT is empty: run Sys.setenv(OPENGWAS_JWT = '...') first")
    snps <- unique(ref$map$ID[grepl("^rs", ref$map$ID)])
    out <- extract_outcome_data(snps = snps, outcomes = row$id, opengwas_jwt = token)
    if (is.null(out) || nrow(out) == 0) return(NULL)
    saveRDS(out, cache)
    out
}

# =============================================================================
# 5. Harmonisation and filters
# =============================================================================
acgt <- c("A", "C", "G", "T")

# Harmonise the outcome to the exposure; drop palindromic / ambiguous SNPs and indels
harmonise_clean <- function(exp_df, out_df) {
    h <- harmonise_data(exp_df, out_df)
    h <- h[h$mr_keep & !h$palindromic & !h$ambiguous & !h$remove, ]
    h$effect_allele.exposure <- as.character(h$effect_allele.exposure)
    h$effect_allele.outcome  <- as.character(h$effect_allele.outcome)
    ok <- h$effect_allele.exposure %in% acgt & h$other_allele.exposure %in% acgt &
          h$effect_allele.exposure == h$effect_allele.outcome
    if (any(!ok)) message("  ", sum(!ok), " SNP(s) removed (indel, or effect allele differs after harmonisation)")
    h <- h[ok, ]
    h[!duplicated(h$SNP), ]
}

# Reject an analysis and record the reason for the log
skip <- function(reason, n = NA) list(ok = FALSE, reason = reason, n = n)

# Prepare the data for one analysis (univariable or multivariable).
#   exposures : named list of load_exposure() results; the first one is used to harmonise the outcome
#   Returns data.frame `all` (all SNPs), `data` (retained instruments),
#   `corr` (LD of the instruments, same order as `data`) and F-stats per SNP.
prepare_mr_data <- function(exposures, outcome_raw, ref, f_min, ld_max, min_snps) {
    ids <- names(exposures)

    # 1. Outcome harmonised to the first exposure
    h <- harmonise_clean(exposures[[1]]$qtl, outcome_raw)
    if (nrow(h) == 0) return(skip("no SNP left after harmonisation"))
    dat <- data.frame(SNP = h$SNP, effect_allele = h$effect_allele.exposure,
                      beta.outcome = h$beta.outcome, se.outcome = h$se.outcome,
                      pval.outcome = h$pval.outcome + 1e-300, stringsAsFactors = FALSE)
    dat[[paste0("beta.", ids[1])]] <- h$beta.exposure
    dat[[paste0("se.", ids[1])]]   <- h$se.exposure
    dat[[paste0("pval.", ids[1])]] <- h$pval.exposure

    # 2. Other exposures: same GTEx SNPs, but the effect allele is still checked
    for (id in ids[-1]) {
        q <- exposures[[id]]$qtl
        q <- data.frame(SNP = q$SNP, ea = q$effect_allele.exposure,
                        b = q$beta.exposure, s = q$se.exposure, p = q$pval.exposure,
                        stringsAsFactors = FALSE)
        q <- q[!duplicated(q$SNP), ]
        dat <- merge(dat, q, by = "SNP", sort = FALSE)
        bad <- dat$ea != dat$effect_allele
        if (any(bad)) message("  ", sum(bad), " SNP(s) removed: effect allele differs for ", id)
        dat <- dat[!bad, ]
        names(dat)[names(dat) %in% c("b", "s", "p")] <- paste0(c("beta.", "se.", "pval."), id)
        dat$ea <- NULL
    }

    # 3. F-stat: a SNP passes if F >= f_min for at least one exposure
    f_tab <- data.frame(SNP = dat$SNP, stringsAsFactors = FALSE)
    for (id in ids) {
        f_tab[[paste0("F.", id)]] <- exposures[[id]]$fstat$F.stat[match(f_tab$SNP, exposures[[id]]$fstat$SNP)]
    }
    f_cols <- paste0("F.", ids)
    f_tab$F.max <- apply(f_tab[, f_cols, drop = FALSE], 1, function(x) if (all(is.na(x))) -Inf else max(x, na.rm = TRUE))
    dat <- merge(dat, f_tab, by = "SNP", sort = FALSE)

    pass <- dat$SNP[dat$F.max >= f_min]
    if (length(pass) == 0) return(skip(paste0("no SNP with F >= ", f_min)))

    # Strongest SNP of each exposure (protected during LD pruning)
    strongest <- unlist(lapply(ids, function(id) {
        f <- exposures[[id]]$fstat
        f$SNP[which.max(f$F.stat)]
    }))

    # 4. Instruments present in the LD matrix
    keep <- dat$SNP %in% pass & dat$SNP %in% rownames(ref$corr)
    if (sum(keep) < min_snps) return(skip(paste0("fewer than ", min_snps, " SNPs with F >= ", f_min), sum(keep)))
    inst <- dat[keep, ]

    # 5. LD pruning
    corr <- ref$corr[inst$SNP, inst$SNP, drop = FALSE]
    corr <- ld_prune(corr, setNames(inst$F.max, inst$SNP), ld_max, protect = strongest)
    if (nrow(corr) < min_snps) return(skip(paste0("fewer than ", min_snps, " SNPs after LD pruning"), nrow(corr)))
    inst <- inst[match(rownames(corr), inst$SNP), ]

    list(ok = TRUE, all = dat, data = inst, corr = corr, n_pass = length(pass), n = nrow(inst))
}

# While pairs with |r| >= ld_max remain, drop the weakest SNP (lowest F)
# among the SNPs in LD, sparing `protect` unless only protected SNPs are left.
ld_prune <- function(corr, f, ld_max, protect = character(0)) {
    repeat {
        in_ld <- rownames(corr)[rowSums(abs(corr) >= ld_max, na.rm = TRUE) > 1]
        if (length(in_ld) == 0) break
        cand <- setdiff(in_ld, protect)
        if (length(cand) == 0) cand <- in_ld
        fc <- f[cand]
        fc[is.na(fc)] <- -Inf
        drop <- cand[which.min(fc)]
        keep <- rownames(corr) != drop
        corr <- corr[keep, keep, drop = FALSE]
    }
    corr
}

# =============================================================================
# 6. MR results
# =============================================================================
get_slot <- function(x, name) if (name %in% methods::slotNames(x)) methods::slot(x, name) else NA

# One row per exposure and method (same columns as the previous tables)
extract_results <- function(res, label, outcome) {
    test  <- class(res)[1]
    egger <- test %in% c("Egger", "MVEgger")
    data.frame(
        Tissue                = label,
        Outcome               = outcome,
        Exposure              = get_slot(res, "Exposure"),
        Test                  = test,
        Estimate              = get_slot(res, "Estimate"),
        SD                    = get_slot(res, if (egger) "StdError.Est" else "StdError"),
        Pvalue                = get_slot(res, if (egger) "Pvalue.Est" else "Pvalue"),
        `Estimate-Intercept`  = if (egger) get_slot(res, "Intercept") else NA,
        `SD-Intercept`        = if (egger) get_slot(res, "StdError.Int") else NA,
        `Pvalue-Intercept`    = if (egger) get_slot(res, "Pvalue.Int") else NA,
        nSNPs                 = get_slot(res, "SNPs"),
        `Pvalue-Hetero`       = if (test %in% c("IVW", "MVIVW", "Egger", "MVEgger")) get_slot(res, "Heter.Stat")[2] else NA,
        check.names = FALSE, stringsAsFactors = FALSE
    )
}

run_uni_mr <- function(prep, exposure, outcome, label) {
    d <- prep$data
    obj <- mr_input(bx = d[[paste0("beta.", exposure)]], bxse = d[[paste0("se.", exposure)]],
                    by = d$beta.outcome, byse = d$se.outcome,
                    corr = prep$corr, exposure = exposure, outcome = outcome)
    do.call(rbind, lapply(list(mr_ivw(obj), mr_egger(obj), mr_median(obj)),
                          extract_results, label = label, outcome = outcome))
}

# `order` : exposure IDs; the first one orients MR-Egger
run_mv_mr <- function(prep, order, outcome, label) {
    d <- prep$data
    shorts <- vapply(order, function(i) exposure_specs[[i]]$short, "")
    obj <- mr_mvinput(bx = as.matrix(d[, paste0("beta.", order)]), bxse = as.matrix(d[, paste0("se.", order)]),
                      by = d$beta.outcome, byse = d$se.outcome,
                      corr = prep$corr, exposure = unname(shorts), outcome = outcome)
    do.call(rbind, lapply(list(mr_mvivw(obj), mr_mvegger(obj), mr_mvmedian(obj)),
                          extract_results, label = label, outcome = outcome))
}

# =============================================================================
# 7. Figures
# =============================================================================
col_method <- c(IVW = "#ad8aff", MVIVW = "#ad8aff",
                Egger = "#e15512", MVEgger = "#e15512",
                WeightedMedian = "#0097ab", MVMedian = "#0097ab")
col_blue <- "#2a78d6"; col_orange <- "#eb6834"; col_gray <- "#c3c2b7"

fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p >= 0.01, formatC(p, format = "f", digits = 2),
                                                   formatC(p, format = "e", digits = 2)))

mr_title <- function(r, outcome, exposure) {
    lines <- c(outcome, exposure)
    for (i in seq_len(nrow(r))) {
        lines <- c(lines, sprintf("%s : b = %.3f, p = %s", r$Test[i], r$Estimate[i], fmt_p(r$Pvalue[i])))
        if (r$Test[i] %in% c("Egger", "MVEgger"))
            lines <- c(lines, sprintf("Egger intercept : p = %s", fmt_p(r[["Pvalue-Intercept"]][i])))
    }
    paste(lines, collapse = "\n")
}

# Scatter of exposure beta vs outcome beta + IVW, Egger and median lines
plot_scatter <- function(d, exposure, res, outcome) {
    r <- res[res$Exposure == exposure | res$Exposure == exposure_specs[[exposure]]$short, ]
    df <- data.frame(x = d[[paste0("beta.", exposure)]], xse = d[[paste0("se.", exposure)]],
                     y = d$beta.outcome, yse = d$se.outcome)
    p <- ggplot(df, aes(x, y)) +
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
        geom_errorbar(aes(xmin = x - 1.96 * xse, xmax = x + 1.96 * xse), color = col_gray) +
        geom_errorbar(aes(ymin = y - 1.96 * yse, ymax = y + 1.96 * yse), color = col_gray) +
        geom_point() +
        labs(x = paste("Effect on", exposure), y = paste("Effect on", outcome)) +
        theme_bw()
    for (i in seq_len(nrow(r))) {
        inter <- if (is.na(r[["Estimate-Intercept"]][i])) 0 else r[["Estimate-Intercept"]][i]
        p <- p + geom_abline(intercept = inter, slope = r$Estimate[i], colour = col_method[[r$Test[i]]], linewidth = 1.2)
    }
    p + ggtitle(mr_title(r, outcome, exposure)) + theme(plot.title = element_text(size = 9))
}

# -log10(p) exposure vs outcome (CETP half of the region)
plot_pvalue_comparison <- function(all, exposure, outcome, ref) {
    all$POS <- ref$map$POS_Hg38[match(all$SNP, ref$map$ID)]
    all <- all[!is.na(all$POS) & all$POS > ref$mid, ]
    all$top <- ifelse(seq_len(nrow(all)) == which.min(all[[paste0("pval.", exposure)]]), "Max", "Other")
    ggplot(all, aes(-log10(pval.outcome), -log10(.data[[paste0("pval.", exposure)]]), colour = top)) +
        geom_point() +
        scale_colour_manual(values = c(Max = col_orange, Other = "black"), name = NULL) +
        labs(x = paste("-log10(p)", outcome), y = paste("-log10(p)", exposure)) +
        theme_bw() + ggtitle(outcome)
}

# CETP locus with exons; SNPs used as instruments are coloured
plot_locus <- function(all, used, what, exons, ref, title) {
    all$POS  <- ref$map$POS_Hg38[match(all$SNP, ref$map$ID)]
    all <- all[!is.na(all$POS) & all$POS > ref$mid, ]
    all$Used <- ifelse(all$SNP %in% used, "Yes", "No")
    all$y <- -log10(all[[what]])
    ggplot(all, aes(POS, y, colour = Used)) +
        geom_point() +
        scale_colour_manual(values = c(No = col_gray, Yes = col_blue), name = "Instrument") +
        geom_segment(x = min(exons$start), xend = max(exons$end), y = -0.1, yend = -0.1, colour = col_blue, inherit.aes = FALSE) +
        geom_rect(data = exons, aes(xmin = start, xmax = end, ymin = -0.4, ymax = 0.2),
                  inherit.aes = FALSE, fill = "black", colour = "black") +
        labs(x = "Position Hg38", y = "-log10(p)") +
        theme_bw() + ggtitle(title)
}

dir_create <- function(...) {
    p <- file.path(...)
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    p
}
