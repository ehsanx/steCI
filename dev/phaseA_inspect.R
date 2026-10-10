# Step 0: inspect the stored artefacts so the measurements below mirror
# the pipeline's own basis exactly. Read-only.
R <- "./simulation/output/results"
cat("=== metrics_ci.csv ===\n")
ci <- read.csv(file.path(R, "metrics_ci.csv"), stringsAsFactors = FALSE)
cat("dim:", dim(ci), "\n"); cat("names:", paste(names(ci), collapse=", "), "\n")
cat("methods:", paste(sort(unique(ci$method)), collapse=", "), "\n")
cat("estimands:", paste(sort(unique(ci$estimand)), collapse=", "), "\n")
cat("tau vals:", paste(sort(unique(ci$tau)), collapse=", "), "\n\n")

cat("=== metrics_empirical.csv ===\n")
em <- read.csv(file.path(R, "metrics_empirical.csv"), stringsAsFactors = FALSE)
cat("dim:", dim(em), "\n"); cat("names:", paste(names(em), collapse=", "), "\n\n")

cat("=== truth_cache.rds ===\n")
tc <- readRDS("./simulation/output/truth_cache.rds")
cat("class:", class(tc), " names:", paste(names(tc), collapse=", "), "\n")
str(tc, max.level = 2)

cat("\n=== replicates_all.rds (91MB) ===\n")
rp <- readRDS(file.path(R, "replicates_all.rds"))
cat("class:", class(rp), " dim:", dim(rp), "\n")
cat("names:", paste(names(rp), collapse=", "), "\n")
if ("method" %in% names(rp)) cat("methods:", paste(sort(unique(rp$method)), collapse=", "), "\n")
if ("estimand" %in% names(rp)) cat("estimands:", paste(sort(unique(rp$estimand)), collapse=", "), "\n")
cat("--- head ---\n"); print(utils::head(as.data.frame(rp), 4))
cat("\n=== metrics.R truth lookup ===\n")
mm <- readLines("./simulation/metrics.R")
print(grep("truth|tc\\[|alpha", mm, value = TRUE)[1:12])
