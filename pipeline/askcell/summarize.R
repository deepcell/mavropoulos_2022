library(data.table)
library(ggplot2)
library(argparse)
library(splitstackshape)
library(RColorBrewer)
library(quadprog)

options(width = 600)
MAX_SCORE <- 100

genotyping <- function(dt, ploidy) {
  dt <- dt[dt$DEPTH > 50]
  if (nrow(dt) == 0) {
    stop("[ERROR]: all the sites with coverage less than 50")
  }

  j <- 0
  N_prev <- 10

  # initial estimates of N (effective genome equivalents sampled) and e (assay error rate)
  N <- 100
  e <- 0.01
  max_iter <- 1000
  N_estimates <- c(N)

  # iteratively estimate N and the genotypes
  iter <- 0
  while (N_prev != N & iter < max_iter) {
    # get likelihoods for all allowed genotypes
    iter <- iter + 1
    for (ploidy in c(2, 3, 4, 5, 6)) {
      allowed_genotypes <- seq(0, ploidy, length.out = (ploidy + 1)) / ploidy
      print("##################################")
      print(paste("ploidy: ", ploidy, sep = ""))
      print(allowed_genotypes)
      likelihoods <- c()
      for (genotype in allowed_genotypes) {
        if (genotype == 0) {
          l <- dbinom(round(dt$AD / dt$DEPTH * N, 0), N, e)
        } else if (genotype == 1) {
          l <- dbinom(round(dt$RD / dt$DEPTH * N, 0), N, e)
        } else {
          l <- dbinom(as.integer(round(dt$AF * N, 0)), N, genotype)
        }
        likelihoods <- cbind(likelihoods, l)
      }
      max_likelihood <- apply(likelihoods, 1, max)

      next_likelihood <- apply(likelihoods, 1, function(x) sort(x, decreasing = T)[2])
      llr <- 10 * log10(max_likelihood / next_likelihood)

      llr[llr > MAX_SCORE] <- MAX_SCORE
      llr <- round(llr, 0)
      for (j in 1:nrow(dt)) {
        dt$G[j] <- allowed_genotypes[which(likelihoods[j, ] == max_likelihood[j])]
      }
      dt$score <- llr
      # re-estimate e and N
      e <- mean(c(dt$AF[dt$G == 0], 1 - dt$AF[dt$G == 1]))
      N_prev <- N
      N <- round(1 / with(dt[dt$G * (1 - dt$G) > 0, ], mean((AF - G)^2 / (G * (1 - G)))), 0)
      N_estimates <- c(N_estimates, N)
    }
  }
  dt
}

sz_genotyping <- function(dt, ploidy) {
  dt <- dt[dt$DEPTH > 50]
  if (nrow(dt) == 0) {
    stop("[ERROR]: all the sites with coverage less than 50")
  }

  e <- mean(c(dt[AF < 0.05]$AF, 1 - dt[AF > 0.95]$AF))
  allowed_genotypes <- c(0, 1/6, 1/5, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 1)
  ploidies <- c(2, 6, 5, 4, 3, 5, 2, 5, 3, 4, 5, 6, 2)
  max_llik <- -9999999
  N <- 100
  N_prev <- 10
  max_iter <- 1000
  iter <- 0
  while (N_prev != N & iter < max_iter) {
    iter <- iter + 1
    likelihoods <- c()
    for (genotype in allowed_genotypes) {
      if (genotype == 0) {
        l <- dbinom(round(dt$AD / dt$DEPTH * N, 0), N, e)
      } else if (genotype == 1) {
        l <- dbinom(round(dt$RD / dt$DEPTH * N, 0), N, e)
      } else {
        l <- dbinom(as.integer(round(dt$AF * N, 0)), N, genotype)
      }
      likelihoods <- cbind(likelihoods, l)
    }
    max_likelihood <- apply(likelihoods, 1, max)
    llik <- sum(log10(max_likelihood))
    if (llik > max_llik) {
      max_llik <- llik
      # print(head(likelihoods))
      max_prob_idx_per_site <- apply(likelihoods, 1, function(x) which(x == max(x)))
      # genotype_est[, i] <- allowed_genotypes[max_prob_idx_per_site]
      dt[, "InferredGenotype" := allowed_genotypes[max_prob_idx_per_site]]
      dt[, "Ploidy" := ploidies[max_prob_idx_per_site]]
      # print(max_llik)
      # print(head(dt))
    }

    # re-estimate e and N
    e <- mean(c(dt[InferredGenotype == 0]$AF, 1 - dt[InferredGenotype == 1]$AF))
    N_prev <- N
    N <- round(1 / with(
      dt[InferredGenotype * (1 - InferredGenotype) > 0, ],
      mean((AF - InferredGenotype)^2 / (InferredGenotype * (1 - InferredGenotype)))
    ), 0)
  }
  dt
}

estimateMultiMixProp <- function(X, Y) {
  Rinv <- solve(chol(t(X) %*% X))
  C <- cbind(rep(1, ncol(X)), diag(ncol(X)))
  b <- c(1, rep(0, ncol(X)))
  d <- t(Y) %*% X
  coef <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)$solution

  coef
}

summarizeSNPs <- function(args) {
  variants_dir <- ""
  if (args$workflow == "tgs") {
    variants_dir <- file.path(args$work_dir, "mutations/")
    # samples <- list.files(variants_dir, "*.mutation.tsv", full.names=F)
  } else {
    variants_dir <- file.path(args$work_dir, "germlines/")
    # samples <- list.files(variants_dir, "*.snp.tsv", full.names=F)
  }
  if (nchar(variants_dir) == 0) {
    stop("Did not find the variants directory")
  }
  if (!dir.exists(variants_dir)) {
    stop(paste(variants_dir, " does not exist", sep = ""))
  }

  bulk_pure_samples <- c()
  meta_dt <- fread(args$metasheet)
  if (nrow(meta_dt) > nrow(unique(meta_dt, by = "SampleName"))) {
    stop("Duplicated samples exist")
  }
  meta_dt[, "SampleName" := gsub("_", "-", SampleName)]
  agg_variants_dt <- data.table()
  agg_dropout_dt <- data.table()
  variants_metrics_dt <- data.table()
  sorted_samples <- c()
  for (i in seq_len(nrow(meta_dt))) {
    sample <- meta_dt[i, ]$SampleName
    if (length(sample) == 0) {
      next
    }
    type <- meta_dt[i, ]$Type
    spikeinsamples <- meta_dt[i, ]$SpikedinSamples
    backgroundsample <- meta_dt[i, ]$Background
    if (type == "B" && spikeinsamples == ".") {
      bulk_pure_samples <- c(bulk_pure_samples, sample)
    }
    if (type != "B") {
      bulk_pure_samples <- c(bulk_pure_samples, backgroundsample)
      bulk_pure_samples <- c(bulk_pure_samples, unlist(strsplit(spikeinsamples, ",")))
    }
    sample_variants <- ""
    if (args$workflow == "tgs") {
      sample_variants <- list.files(variants_dir, paste(sample, ".*.mutation.tsv", sep = ""), full.names = T)
    } else {
      sample_variants <- list.files(variants_dir, paste(sample, ".*.snp.tsv", sep = ""), full.names = T)
    }
    if (length(sample_variants) == 0) {
      stop("Fail to find the sample variants")
    }
    sample_variants_dt <- fread(sample_variants)
    sample_variants_dt <- cSplit(sample_variants_dt,
      splitCols = c("ALT", "AD", "AF"), sep = ",",
      direction = "long", type.convert = F
    )
    sample_variants_dt <- sample_variants_dt[!is.na(ALT)]
    sample_variants_dt[, ":="(DEPTH = as.numeric(DEPTH), RD = as.numeric(RD), AD = as.numeric(AD), AF = as.double(AF))]
    sample_variants_dt[, "Site" := paste(CHR, POSITION, REF, sep = ":")]
    sample_variants_dt <- sample_variants_dt[sample_variants_dt[, .I[AF == max(AF)], by = Site]$V1]
    n_sites <- nrow(unique(sample_variants_dt, by = "Site"))
    n_nonref_calls <- nrow(sample_variants_dt[ALT != "."])
    n_snp_cov_dropout <- nrow(sample_variants_dt[DEPTH == 0])
    n_snp_lt_50_cov <- nrow(sample_variants_dt[DEPTH <= 50 & DEPTH > 0])
    median_cov <- median(unique(sample_variants_dt, by = "Site")$DEPTH)
    min_cov <- min(unique(sample_variants_dt, by = "Site")$DEPTH)
    max_cov <- max(unique(sample_variants_dt, by = "Site")$DEPTH)
    sample_variants_dt[, "VariantID" := paste(CHR, POSITION, REF, ALT, sep = ":")]
    sample_variants_dt[, "Site" := NULL]
    if (spikeinsamples != ".") {
      sorted_samples <- c(sorted_samples, sample)
      # the reason I do not do the mixture analysis here is because
      # the new bulk needed might not be in the db yet
      # need to integrate
      contam_frac <- -1
      num_contam_gt_1pct <- -1
      num_contam_gt_10pct <- -1
      sample_variants_dt[, "InferredGenotype" := -1.0]
      sample_variants_dt[, "Ploidy" := -1]
    } else {
      sample_variants_dt <- sz_genotyping(sample_variants_dt)
      sample_variants_dt <- sample_variants_dt[order(InferredGenotype)]
      homo_variants_dt <- sample_variants_dt[InferredGenotype == 1 | InferredGenotype == 0 & DEPTH > 50]
      homo_variants_dt[, "AF" := ifelse(InferredGenotype == 1, 1 - AF, AF)]
      contam_frac <- 2 * mean(homo_variants_dt$AF)
      num_contam_gt_1pct <- nrow(homo_variants_dt[homo_variants_dt$AF > 0.01])
      num_contam_gt_10pct <- nrow(homo_variants_dt[homo_variants_dt$AF > 0.1])
    }
    variants_metrics_dt <- rbind(
      variants_metrics_dt,
      data.table(
        SampleName = sample,
        NumSites = n_sites,
        NumALTCalls = n_nonref_calls,
        NumLT50Cov = n_snp_lt_50_cov,
        NumDropOut = n_snp_cov_dropout,
        MedianSiteCov = median_cov,
        MinSiteCov = min_cov,
        MaxSiteCov = max_cov,
        ContamFrac = contam_frac,
        NumContamFracGT1pct = num_contam_gt_1pct,
        NumContamFracGT10pct = num_contam_gt_10pct
      )
    )

    sample_variants_dt[, "LocusIndex" := seq(1, nrow(sample_variants_dt))]
    sample_variants_dt[, "Sample" := sample]
    agg_variants_dt <- rbind(agg_variants_dt, sample_variants_dt)
  }
  out_plot_dir <- file.path(args$work_dir, "plots")
  if (!dir.exists(out_plot_dir)) {
    dir.create(out_plot_dir)
  }
  out_gt_plot <- gsub(".tsv", ".genotypes.pdf", basename(args$metasheet))
  out_gt_plot <- file.path(out_plot_dir, out_gt_plot)
  if (nrow(agg_variants_dt[DEPTH > 50 & InferredGenotype != -1]) > 0) {
    plot_dt <- copy(agg_variants_dt[DEPTH > 50 & InferredGenotype != -1]) # this can be bad when the data is huge
    plot_dt[, "InferredGenotype" := ifelse((InferredGenotype == 0 & AF >= 0.01) |
      (InferredGenotype == 1 & (1 - AF) >= 0.01),
    "Contam",
    ifelse(InferredGenotype == -1.0, "NA",
      as.character(InferredGenotype)
    )
    )]
    if (nrow(plot_dt[DEPTH > 50 & InferredGenotype != "NA"]) > 0) {
      m <- ggplot(plot_dt[DEPTH > 50 & InferredGenotype != "NA"], aes(
        x = LocusIndex, y = AF,
        color = as.factor(Ploidy)
        )) +
        geom_point(size = 2, alpha = 0.9) +
        facet_wrap(~Sample, ncol = 3) +
        theme_bw() +
        labs(x = "SNPs", color = "Copy Number State") +
        theme(
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom"
        ) +
        scale_color_brewer(palette = "Set2")
      ggsave(out_gt_plot, m, width = 10, height = 10)
    }
  }
  agg_variants_dt[, "LocusIndex" := NULL]
  out_variants_long <- gsub(".tsv", ".variants.long.tsv", basename(args$metasheet))
  out_variants_long <- file.path(variants_dir, out_variants_long)
  fwrite(agg_variants_dt, out_variants_long, sep = "\t", row.names = F, quote = F)
  # update bulk pure cell line snp DB if some new bulk samples are yet to be included
  bulk_pure_samples <- unique(bulk_pure_samples)
  if (length(bulk_pure_samples) > 0) {
    bulk_pure_snp_db <- data.table()
    if (!file.exists(args$bulk_db)) {
      saveRDS(bulk_pure_snp_db, args$bulk_db)
    } else {
      # if it exists, read the database
      bulk_pure_snp_db <- as.data.table(readRDS(args$bulk_db))
    }
    if (nrow(bulk_pure_snp_db) == 0) {
      bulk_pure_snp_db <- rbind(bulk_pure_snp_db, agg_variants_dt[Sample %in% bulk_pure_samples])
      saveRDS(bulk_pure_snp_db, args$bulk_db)
    } else {
      samples_to_insert <- bulk_pure_samples[which(!bulk_pure_samples %in% bulk_pure_snp_db$Sample)]
      if (length(samples_to_insert) > 0) {
        print("Inserting bulk samples to the existing db")
        print(samples_to_insert)
        bulk_pure_snp_db <- rbind(bulk_pure_snp_db, agg_variants_dt[Sample %in% samples_to_insert])
        saveRDS(bulk_pure_snp_db, args$bulk_db)
      } else {
        print("found all required bulk samples")
      }
    }
  }
  sorted_meta_dt <- meta_dt[SampleName %in% sorted_samples]
  mixture_data_dt <- data.table()
  mixture_res_dt <- data.table()
  if (nrow(sorted_meta_dt) > 0 && args$workflow == "sid") {
    # 	# now lets do the mixture analysis
    for (i in seq_len(nrow(sorted_meta_dt))) {
      sample <- sorted_meta_dt[i, ]$SampleName
      spikedinsamples <- unlist(strsplit(sorted_meta_dt[i, ]$SpikedinSamples, split = ","))
      background <- sorted_meta_dt[i, ]$Background
      sample_variants <- list.files(variants_dir, paste(sample, ".*.snp.tsv", sep = ""), full.names = T)
      sample_variants_dt <- fread(sample_variants)
      sample_variants_dt <- cSplit(sample_variants_dt,
        splitCols = c("ALT", "AD", "AF"), sep = ",",
        direction = "long", type.convert = F
      )
      sample_variants_dt <- sample_variants_dt[!is.na(ALT)]
      sample_variants_dt[, ":="(DEPTH = as.numeric(DEPTH), RD = as.numeric(RD), AD = as.numeric(AD), AF = as.double(AF))]
      sample_variants_dt <- sample_variants_dt[DEPTH > 0]
      sample_variants_dt[, "VariantID" := paste(CHR, POSITION, REF, ALT, sep = ":")]
      sample_variants_dt[, "Site" := paste(CHR, POSITION, REF, sep = ":")]
      sample_variants_dt <- sample_variants_dt[sample_variants_dt[, .I[AF == max(AF)], by = Site]$V1]
      if (length(spikedinsamples) == 1) {
        print(paste("Estimating purity for sample  ", sample, sep = ""))
        # julie: background <- gsub("_", "-", background)
        # julie: spikedinsamples <- gsub("_", "-", spikedinsamples)
        background_snp_db <- bulk_pure_snp_db[Sample == background]
        background_snp_db[, "VariantID" := paste(CHR, POSITION, REF, ALT, sep = ":")]
        background_snp_db[, "ALT" := ifelse(ALT == ".", REF, ALT)]
        background_snp_db <- background_snp_db[, .(VariantID, ALT, AF, InferredGenotype)]
        names(background_snp_db) <- c("VariantID", "BaseALT", "BaseAF", "Base_G")

        spikein_snp_db <- bulk_pure_snp_db[Sample == spikedinsamples]
        spikein_snp_db[, "VariantID" := paste(CHR, POSITION, REF, ALT, sep = ":")]
        spikein_ploidy <- mean(spikein_snp_db$Ploidy)
        spikein_snp_db[, "ALT" := ifelse(ALT == ".", REF, ALT)]
        spikein_snp_db <- spikein_snp_db[, .(VariantID, ALT, AF, InferredGenotype)]
        names(spikein_snp_db) <- c("VariantID", "SpikeALT", "SpikeAF", "Spike_G")

        n_alt_site <- nrow(sample_variants_dt[REF == ALT])
        print(paste("Removing ", n_alt_site, " non-alt sites", sep = ""))
        sample_variants_dt <- sample_variants_dt[REF != ALT]
        sample_variants_dt[, "ALT" := ifelse(ALT == ".", REF, ALT)]
        sample_variants_dt <- merge(sample_variants_dt, background_snp_db, by = c("VariantID"), all.x = T)
        sample_variants_dt <- merge(sample_variants_dt, spikein_snp_db, by = c("VariantID"), all.x = T)
        # sample_variants_dt[, "BaseNA" := is.na(BaseALT)]
        # sample_variants_dt[, "SpikeNA" := is.na(SpikeALT)]
        # print(sample_variants_dt[ is.na(BaseALT) &  is.na(SpikeALT)])
        n_non_informative <- nrow(sample_variants_dt[is.na(BaseALT) & is.na(SpikeALT)])
        print(paste("Removing ", n_non_informative, " non-informative sites", sep = ""))
        sample_variants_dt <- sample_variants_dt[!is.na(BaseALT) | !is.na(SpikeALT)]

        sample_variants_dt[, "BaseAF" := ifelse(is.na(BaseALT), 0.0, BaseAF)]
        sample_variants_dt[, "Base_G" := ifelse(is.na(BaseALT), 0.0, Base_G)]
        sample_variants_dt[, "BaseALT" := ifelse(is.na(BaseALT), REF, BaseALT)]
        sample_variants_dt[, "SpikeAF" := ifelse(is.na(SpikeALT), 0.0, SpikeAF)]
        sample_variants_dt[, "Spike_G" := ifelse(is.na(SpikeALT), 0.0, Spike_G)]
        sample_variants_dt[, "SpikeALT" := ifelse(is.na(SpikeALT), REF, SpikeALT)]
        sample_variants_dt[, "NumALT" := length(unique(c(REF, ALT, BaseALT, SpikeALT))),
          by = seq_len(nrow(sample_variants_dt))
        ]

        # estimate error initially ad hoc and later from genotypes in iterative part
        e <- mean(c(sample_variants_dt[AF < 0.05]$AF, 1 - sample_variants_dt[AF > 0.95]$AF))
        mean_af_het <- mean(c(sample_variants_dt[Base_G == 0.5 & Spike_G == 0]$AF, 1 - sample_variants_dt[Base_G == 0.5 & Spike_G == 1]$AF), na.rm = T)
        max_var <- var(c(sample_variants_dt[Base_G == 0.5 & Spike_G == 0]$AF, 1 - sample_variants_dt[Base_G == 0.5 & Spike_G == 1]$AF), na.rm = T)
        N <- round(mean_af_het * (1 - mean_af_het) / max_var, 0)

        # data <- sample_variants_dt[Base_G != Spike_G & Spike_G != "Contam", ]
        data <- sample_variants_dt[Base_G != Spike_G, ]

        fraction <- seq(0.00001, 0.99998, 0.00001)
        ml_frac <- ml_ll <- NA
        for (frac in fraction) {
          af_expected <- with(data, as.double(Base_G) * (1 - frac) + as.double(Spike_G) * frac)
          logp <- dbinom(round(data$AF * N, 0), N, af_expected, log = T)
          sumlogp <- sum(logp)
          if (is.na(ml_frac)) {
            ml_frac <- frac
            ml_ll <- sumlogp
          } else {
            if (sumlogp > ml_ll) {
              ml_frac <- frac
              ml_ll <- sumlogp
            }
          }
        }
        # use only the cn=2 sites
        simp_frac <- median(c(
          2 * sample_variants_dt[Base_G == 0 & Spike_G == 0.5]$AF,
          2 * (1 - sample_variants_dt[Base_G == 1 & Spike_G == 0.5]$AF),
          sample_variants_dt[Base_G == 0 & Spike_G == 1]$AF,
          1 - sample_variants_dt[Base_G == 1 & Spike_G == 0]$AF,
          1 - 2 * sample_variants_dt[Base_G == 0.5 & Spike_G == 0]$AF,
          2 * sample_variants_dt[Base_G == 0.5 & Spike_G == 1]$AF - 1
        ))
        contam_frac <- 2 * mean(c(
          sample_variants_dt[Base_G == 0 & Spike_G == 0]$AF,
          1 - sample_variants_dt[Base_G == 1 & Spike_G == 1]$AF
        ))
        num_contam_gt_1pct <- nrow(sample_variants_dt[Base_G == 0 & Spike_G == 0 & AF > 0.01, ]) +
          nrow(sample_variants_dt[Base_G == 1 & Spike_G == 1 & AF < 0.99, ])
        num_contam_gt_10pct <- nrow(sample_variants_dt[Base_G == 0 & Spike_G == 0 & AF > 0.1, ]) +
          nrow(sample_variants_dt[Base_G == 1 & Spike_G == 1 & AF < 0.9, ])
        sample_variants_dt[, "PurityEstMLE" := ml_frac]
        mixture_data_dt <- rbind(mixture_data_dt, sample_variants_dt)
        variants_metrics_dt[SampleName == sample, "PurityEstMLE"] <- ml_frac
        variants_metrics_dt[SampleName == sample, "PurityEstSimple"] <- simp_frac
        variants_metrics_dt[SampleName == sample, "ContamFrac"] <- contam_frac
        variants_metrics_dt[SampleName == sample, "NumContamFracGT1pct"] <- num_contam_gt_1pct
        variants_metrics_dt[SampleName == sample, "NumContamFracGT10pct"] <- num_contam_gt_10pct
      } else {
        # multi-mixture prop estimation using quadratic programming
        if (sample %in% unique(bulk_pure_snp_db$Sample)) {
          stop("[ERROR] Some bulk pure spiked in samples are not in the database")
        }
        # julie: background <- gsub("_", "-", background)
        # julie: spikedinsamples <- gsub("_", "-", spikedinsamples)
        mixture_snp_dt <- bulk_pure_snp_db[Sample %in% c(spikedinsamples, background)]
        mixture_snp_dt[, "ALT" := ifelse(InferredGenotype == 0, ".", ALT)]
        mixture_snp_dt[, "AD" := ifelse(InferredGenotype == 0, 0, AD)]
        mixture_snp_dt[, "AF" := ifelse(InferredGenotype == 0, 0.0, AF)]
        mixture_snp_dt[, "VariantID" := paste(CHR, POSITION, REF, ALT, sep = ":")]
        mixture_snp_dt[, "Site" := paste(CHR, POSITION, REF, sep = ":")]
        mixture_snp_dt <- unique(mixture_snp_dt, by = c("VariantID", "Sample"))
        mixture_snp_dt <- rbind(mixture_snp_dt[, ":="(InferredGenotype = NULL, Ploidy = NULL)], sample_variants_dt)
        mixture_snp_dt[, "n" := length(unique(ALT)), by = list(CHR, POSITION, REF)]
        mixture_snp_dt <- mixture_snp_dt[n <= 2]
        mixture_snp_wide_dt <- dcast(mixture_snp_dt, Site ~ Sample, value.var = "AF")
        for (j in seq_len(ncol(mixture_snp_wide_dt))) {
          set(mixture_snp_wide_dt, which(is.na(mixture_snp_wide_dt[[j]])), j, 1e-8)
        }
        print(sample)
        xcols <- c(spikedinsamples, background)
        X <- as.matrix(mixture_snp_wide_dt[, ..xcols])
        X[X == 0] <- 1e-8
        Y <- as.matrix(mixture_snp_wide_dt[, ..sample])
        Y[Y == 0] <- 1e-8
        prop <- estimateMultiMixProp(X, Y)
        prop <- round(prop * 100, 2)
        prop <- paste(prop, "%", sep = "", collapse = ",")
        tmp_dt <- data.table(SampleName = sample, prop = prop)
        tmp_dt[, colnames(X) := tstrsplit(prop, ",", fixed = T)]
        tmp_dt[, "prop" := NULL]
        mixture_res_dt <- rbind(mixture_res_dt, tmp_dt)
      }
    }
    if (nrow(mixture_data_dt) > 0) {
      out_mixture_data <- gsub(".tsv", ".variants.mixture.est.tsv", basename(args$metasheet))
      out_mixture_data <- file.path(variants_dir, out_mixture_data)
	  fwrite(mixture_data_dt, out_mixture_data, sep="\t", row.names=F, quote=F)
    }
    if (nrow(mixture_res_dt) > 0) {
      variants_metrics_dt <- merge(variants_metrics_dt, mixture_res_dt, by = "SampleName")
      for (j in seq_len(ncol(variants_metrics_dt))) {
        set(variants_metrics_dt, which(is.na(variants_metrics_dt[[j]])), j, "NA")
      }
    }
  }
  variants_metrics_dt
}

summarizeMetrics <- function(args) {
  metrics_dir <- file.path(args$work_dir, "metrics/")
  if (!dir.exists(metrics_dir)) {
    stop("[ERROR]: Cannot find the metrics folder")
  }

  meta_dt <- fread(args$metasheet)
  agg_frags_dt <- data.table()
  agg_cov_dt <- data.table()
  agg_metrics_dt <- data.table()
  for (i in seq_len(nrow(meta_dt))) {
    sample <- meta_dt[i, ]$SampleName
    sample <- gsub("_", "-", sample)
    if (length(sample) == 0) {
      next
    }
    sample_metrics <- list.files(file.path(metrics_dir, sample), "*.metrics.csv", full.names = T)
    if (length(sample_metrics) == 0) {
      stop(paste("[ERROR]: Cannot find metrics result for sample ", sample, sep = ""))
    }
    sample_metrics_dt <- fread(sample_metrics)
    names(sample_metrics_dt)[which(names(sample_metrics_dt) == "sample")] <- "SampleName"
    agg_metrics_dt <- rbind(agg_metrics_dt, sample_metrics_dt)

    sample_frags <- list.files(file.path(metrics_dir, sample), "*.frags.csv", full.names = T)
    sample_frags_dt <- fread(sample_frags)
    # sample 100k frags not all of them
    if (nrow(sample_frags_dt) < 100000) {
      agg_frags_dt <- rbind(agg_frags_dt, sample_frags_dt)
    } else {
      agg_frags_dt <- rbind(agg_frags_dt, sample_frags_dt[sample(.N, 100000)])
    }

    sample_cov <- list.files(file.path(metrics_dir, sample), "*.targets.*sv", full.names = T)
    sample_cov_dt <- fread(sample_cov)
    sample_cov_dt[, "sample" := sample]
    agg_cov_dt <- rbind(agg_cov_dt, sample_cov_dt)
  }

  # out_metrics <- file.path(variants_dir, "aggregate.metrics.csv")
  # fwrite(agg_metrics_dt, out_metrics, sep=",", row.names=F, quote=F)
  out_cov_long <- file.path(metrics_dir, "aggregate.targets.cov.long.csv")
  fwrite(agg_cov_dt, out_cov_long, sep = ",", row.names = F, quote = F)

  plot_dir <- file.path(args$work_dir, "plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  frag_hist_plot <- file.path(plot_dir, "aggregated.frags.hist.pdf")
  agg_frags_dt[, "aln_type" := ifelse(aln_type == 1, "Proper",
    ifelse(aln_type == 0, "NonProperSameChr", "NonProperDiffChr")
  )]
  m <- ggplot(agg_frags_dt[frag_len <= 5000], aes(x = frag_len, fill = aln_type)) +
    geom_histogram(binwidth = 25) +
    facet_wrap(~sample, ncol = 3) +
    theme_bw() +
    xlab("Fragment Length")
  ggsave(frag_hist_plot, width = 10)

  frag_box_plot <- file.path(plot_dir, "aggregated.frags.box.pdf")
  m <- ggplot(agg_frags_dt[frag_len <= 5000], aes(x = sample, y = frag_len, color = aln_type)) +
    geom_boxplot() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 65, vjust = 1.0, hjust = 1),
      legend.position = "top"
    ) +
    ylab("Fragment Length")
  ggsave(frag_box_plot, width = 10)

  gc_by_cov_plot <- file.path(plot_dir, "aggregated.cov.gc_by_coverage.pdf")
  m <- ggplot(agg_cov_dt, aes(x = gc, y = coverage)) +
    geom_point(size = 2, alpha = 0.5) +
    facet_wrap(~sample, ncol = 3) +
    theme_bw()
  ggsave(gc_by_cov_plot, width = 10)

  gc_by_ncov_plot <- file.path(plot_dir, "aggregated.cov.gc_by_normcoverage.pdf")
  m <- ggplot(agg_cov_dt, aes(x = gc, y = norm_coverage)) +
    geom_point(size = 2, alpha = 0.5) +
    facet_wrap(~sample, ncol = 3) +
    theme_bw()
  ggsave(gc_by_ncov_plot, width = 10)

  len_by_cov_plot <- file.path(plot_dir, "aggregated.cov.len_by_coverage.pdf")
  m <- ggplot(agg_cov_dt, aes(x = length, y = coverage)) +
    geom_point(size = 2, alpha = 0.5) +
    facet_wrap(~sample, ncol = 3) +
    theme_bw()
  ggsave(len_by_cov_plot, width = 10)

  # agg_cov_dt[, c("buffer", "pbs", "input", "sid", "rep") := tstrsplit(sample, '-', fixed=T)]
  # agg_cov_dt[, "group" := paste(buffer, input, sep='_')]
  targets_plot <- file.path(plot_dir, "aggregated.cov.targets.pdf")
  m <- ggplot(agg_cov_dt, aes(x = amplicon, y = norm_coverage, color = sample)) +
    geom_point(size = 2, alpha = 0.5) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  ggsave(targets_plot, width = 10)

  agg_cov_wide_dt <- dcast(agg_cov_dt[, .(amplicon, norm_coverage, sample)], amplicon ~ sample, value.var = "norm_coverage")
  agg_cov_wide <- file.path(metrics_dir, "aggregate.targets.cov.wide.csv")
  # chr_order <- c(paste("chr", seq(1,22), sep=""), "chrX")
  fwrite(agg_cov_wide_dt, agg_cov_wide, sep = ",", row.names = F, quote = F)

  if (args$workflow == "sid") {
    amplicons <- agg_cov_wide_dt$amplicon
    agg_cov_wide_mat <- as.matrix(agg_cov_wide_dt[, 2:ncol(agg_cov_wide_dt)])
    rownames(agg_cov_wide_mat) <- amplicons
    agg_cov_wide_mat[is.na(agg_cov_wide_mat)] <- 0.0
    agg_cov_wide_mat <- agg_cov_wide_mat + 1e-6
    agg_cov_wide_mat[!is.na(agg_cov_wide_mat)] <- agg_cov_wide_mat[!is.na(agg_cov_wide_mat)] + 1e-06
    res_pca <- prcomp(t(agg_cov_wide_mat))
    variance <- round(res_pca$sdev^2 / sum(res_pca$sdev^2), 3) * 100

    res_pca_df <- data.frame(
      sample = rownames(res_pca$x),
      PC1 = res_pca$x[, 1], PC2 = res_pca$x[, 2]
    )

    pca_plot <- file.path(plot_dir, "aggregated.cov.pca.pdf")
    m <- ggplot(res_pca_df, aes(PC1, PC2, label = sample, color = sample)) +
      geom_text() +
      theme_bw() +
      xlab(paste("PC1 (", variance[1], "%)", sep = "")) +
      ylab(paste("PC2 (", variance[2], "%)", sep = ""))
    ggsave(pca_plot, m, width = 10)
  }

  agg_metrics_dt
}

parser <- ArgumentParser()
parser$add_argument("--sheet",
  dest = "metasheet",
  help = "specify the metasheet"
)
parser$add_argument("--wdir",
  dest = "work_dir",
  help = "specify the working directory"
)
parser$add_argument("--workflow",
  dest = "workflow",
  choices = c("sid", "tgs", "wga", "rna"),
  help = "specify the workflow"
)
parser$add_argument("--bulk_db",
  dest = "bulk_db",
  help = "specify the bulk pure sample AF database"
)
parser$add_argument("--genes",
  dest = "genes",
  help = "specify the gencode annotation at gene level in TSV format"
)

args <- parser$parse_args()

if (!dir.exists(args$work_dir)) {
  stop("[ERROR]: Cannot find the given working directory")
}

if (!file.exists(args$metasheet)) {
  stop("[ERROR]: Cannot find the given metasheet")
}

if (args$workflow == "sid") {
  agg_variants_metrics_dt <- summarizeSNPs(args)
  agg_metrics_dt <- summarizeMetrics(args)

  meta_dt <- fread(args$metasheet)
  # julie: meta_dt[, "SampleName" := gsub("_", "-", SampleName)]
  out_cols <- c(
    names(meta_dt)[1:(ncol(meta_dt) - 1)],
    names(agg_metrics_dt)[2:ncol(agg_metrics_dt)],
    names(agg_variants_metrics_dt)[2:ncol(agg_variants_metrics_dt)],
    "Note"
  )
  meta_dt <- merge(meta_dt, agg_metrics_dt, by = "SampleName")
  meta_dt <- merge(meta_dt, agg_variants_metrics_dt, by = "SampleName")
  meta_dt <- meta_dt[, ..out_cols]
  print(meta_dt)

  out_summary <- gsub(".tsv", ".summary.tsv", basename(args$metasheet))
  out_summary <- file.path(file.path(args$work_dir, "metrics"), out_summary)
  fwrite(meta_dt, out_summary, sep = "\t", row.names = F, quote = F)
} else if (args$workflow == "wga") {
  print("Doing wga stuff")
} else if (args$workflow == "tgs") {
  print("Doing tgs stuff")

agg_metrics_dt <- summarizeMetrics(args)

  meta_dt <- fread(args$metasheet)
  # julie: meta_dt[, "SampleName" := gsub("_", "-", SampleName)]

  meta_dt <- merge(meta_dt, agg_metrics_dt, by = "SampleName")
  print(head(meta_dt))

  out_summary <- gsub(".tsv", ".summary.tsv", basename(args$metasheet))
  out_summary <- file.path(file.path(args$work_dir, "metrics"), out_summary)
  fwrite(meta_dt, out_summary, sep = "\t", row.names = F, quote = F)
} else if (args$workflow == "rna") {
  if (is.null(args$genes)) {
    stop("[ERROR]: you must specify a gencode annotation when running rna workflow")
  }
  genes_dt <- fread(args$genes)
  genes_dt[, c("gene_id", "gene_type", "gene_name", "anno_type") := tstrsplit(V4, split = ",", fixed = T)]
  setkey(genes_dt, gene_id)

  aln_dir <- file.path(args$work_dir, "alignments")
  if (!dir.exists(aln_dir)) {
    stop("[ERROR]: Cannot find the alignment folder")
  }
  metrics_dir <- file.path(args$work_dir, "metrics")
  if (!dir.exists(metrics_dir)) {
    stop("[ERROR]: Cannot find the metrics folder for RNA results")
  }
  metric_files <- list.files(metrics_dir, "*.metrics.csv", full.names = T)
  if (length(metric_files) == 0) {
    stop("[ERROR]: Found no metric files")
  }
  agg_metrics_dt <- data.table()
  for (i in seq_len(length(metric_files))) {
    metrics_dt <- fread(metric_files[i])
    sample <- metrics_dt$SampleName
    sample_cnt <- list.files(file.path(aln_dir, sample), "ReadsPerGene.out.tab", full.names = T)
    if (length(sample_cnt) == 0 || !file.exists(sample_cnt) || is.na(sample_cnt)) {
      stop(paste("[ERROR]: Cannot find the gene count table for sample ", sample, sep = ""))
    }
    sample_cnt_dt <- fread(sample_cnt, skip = 4)
    names(sample_cnt_dt) <- c("gene_id", "unstrand_cnt", "fwd_cnt", "rev_cnt")
    setkey(sample_cnt_dt, gene_id)
    sample_cnt_dt <- sample_cnt_dt[genes_dt]
    sample_cnt_dt[, "tpm" := fwd_cnt / V5] # V5: gene length
    tot_transcript <- sum(sample_cnt_dt$tpm)
    size_factor <- sum(sample_cnt_dt$tpm) / 1000000
    sample_cnt_dt[, "tpm" := tpm / size_factor]
    n_gene_detected <- nrow(sample_cnt_dt[tpm >= 0.1])
    n_mtRNA <- sum(sample_cnt_dt[V1 == "chrM"]$fwd_cnt)
    pct_mtRNA <- n_mtRNA / metrics_dt$NumFrags * 100
    metrics_dt[, "PCT_mtRNA" := pct_mtRNA]
    metrics_dt[, "NumGeneDectected" := n_gene_detected]
    agg_metrics_dt <- rbind(agg_metrics_dt, metrics_dt)
  }
  meta_dt <- fread(args$metasheet)
  # julie: meta_dt[, "SampleName" := gsub("_", "-", SampleName)]
  meta_dt <- merge(meta_dt, agg_metrics_dt, by = "SampleName")

  out_meta <- file.path(metrics_dir, gsub(".tsv", ".summary.tsv", basename(args$metasheet)))
  fwrite(meta_dt, out_meta, sep = "\t", row.names = F, quote = F)

  # aggregate the gene expression profile
  gene_cnt_files <- list.files(aln_dir, "ReadsPerGene.out.tab", recursive = T, full.name = T)
  if (length(gene_cnt_files) == 0) {
    stop("[ERROR]: Found no ReadsPerGene.out.tab files")
  }
  agg_cnt_dt <- data.table()
  for (i in seq_len(length(gene_cnt_files))) {
    sample <- basename(dirname(gene_cnt_files[i]))
    # julie: sample <- gsub("-", "_", sample)
    cnt_dt <- fread(gene_cnt_files[i], skip = 4)
    setkey(cnt_dt, V1)
    if (nrow(agg_cnt_dt) == 0) {
      agg_cnt_dt <- cnt_dt[, c("V1", "V3")]
    } else {
      agg_cnt_dt <- agg_cnt_dt[cnt_dt[, c("V1", "V3")]]
    }
    names(agg_cnt_dt)[length(names(agg_cnt_dt))] <- sample
    setkey(agg_cnt_dt, V1)
  }
  agg_cnt_df <- as.data.frame(agg_cnt_dt[, 2:length(names(agg_cnt_dt))])
  rownames(agg_cnt_df) <- agg_cnt_dt$V1
  expr_dir <- file.path(args$work_dir, "expression")
  if (!dir.exists(expr_dir)) {
    dir.create(expr_dir)
  }
  out_expr <- file.path(expr_dir, "aggregate.gene.counts.tsv")
  fwrite(agg_cnt_df, out_expr, sep = "\t", row.names = T, quote = F)
}
