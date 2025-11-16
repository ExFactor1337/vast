#!/usr/bin/env Rscript

# Load packages, suppressing package startup messages
suppressPackageStartupMessages({
  library(getopt)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer) 
})

# --- Utility Function for Chromosome Standardization ---
# Converts "chr1", "chrX", "1" to "1", "X", "1" for consistent mapping.
standardize_chrom <- function(chrom_vector) {
  # data.table::fifelse is highly efficient for this conditional replacement
  # Use sub() to remove the leading 'chr' if it exists.
  return(sub("^chr", "", chrom_vector))
}
# --------------------------------------------------------

# --- Test/Debugging Configuration ---
# Uncomment the line below for quick testing/debugging outside of a full CLI run:
# opt <- list(ARGS = character(0), info_data = "test_infos_variants.tsv", format_data = "test_formats_variants.tsv", target = "AF", source = "INFO", reference_coords = "test_genome_hg19.tsv", filter_config = "test_filter_config.tsv", info_key = "test_infos.tsv", format_key = "test_formats.tsv")
# opt <- list(ARGS = character(0), info_data = "ccle.small_infos_variants.tsv", format_data = "ccle.small_formats_variants.tsv", target = "DP", source = "FORMAT", reference_coords = "ccle.small_genome_hg19.tsv", filter_config = "ccle.small_filter_config.tsv", info_key = "ccle.small_infos.tsv", format_key = "ccle.small_formats.tsv")


# --- Define Command Line Options ---
opt_spec = matrix(c(
  'help',             'h', 0, "logical", "Show this help message and exit.",
  'info_data',        'i', 1, "character", "REQUIRED. Path to the INFO data tsv.",
  'target',           't', 1, "character", "REQUIRED. Name of the column representing the target variable.",
  'source',           'o', 1, "character", "REQUIRED. Source data table for the target variable (must be 'INFO' or 'FORMAT').",
  'reference_coords', 'r', 1, "character", "REQUIRED. Path to the reference genome coordinates file.",
  'format_data',      's', 1, "character", "OPTIONAL. Path to FORMAT data tsv. REQUIRED if --source is 'FORMAT'.",
  'is_categorical',   'c', 0, "logical", "OPTIONAL. Flag to treat the target variable as categorical (default: FALSE).",
  'filter_config',    'f', 1, "character", "OPTIONAL. Path to the configuration file for data filtering (e.g., TSV).",
  'transform_config', 'x', 1, "character", "OPTIONAL. Path to the configuration file for data transformations.",
  'info_key',         'k', 1, "character", "REQUIRED. Path to a file listing selected INFO keys for analysis.",
  'format_key',       'y', 1, "character", "OPTIONAL. Path to a file listing selected FORMAT keys for analysis."
), byrow=TRUE, ncol=5)

# --- Parse Arguments ---
opt = getopt(opt_spec)

# --- Handle Help Message ---
if ( !is.null(opt$help) ) {
  cat(getopt(opt_spec, usage=TRUE))
  q(status=0)
}

# --- Check Required Arguments ---
required_args <- c('info_data', 'target', 'reference_coords', 'info_key', 'source')
missing_args <- setdiff(required_args, names(opt))

if (length(missing_args) > 0) {
  stop(
    paste(
      "Missing required arguments:",
      paste(missing_args, collapse=", "),
      "\nUse --help for usage details."
    )
  )
}

# --- Validate Source and Dependencies ---

# Validate that the source is a known value
valid_sources <- c("INFO", "FORMAT")
if (!opt$source %in% valid_sources) {
  stop(sprintf(
    "Invalid value for --source: '%s'. Must be one of: %s",
    opt$source,
    paste(valid_sources, collapse = ", ")
  ))
}

# Conditional requirement: If source is FORMAT, then format_data must be provided.
if (opt$source == "FORMAT" && is.null(opt$format_data)) {
  stop("Error: --format_data (-s) is REQUIRED when --source is set to 'FORMAT'.")
}


# --- Apply Defaults for Optional Arguments ---

if (is.null(opt$is_categorical)) {
  opt$is_categorical = FALSE
}

# Set paths to NULL if not provided.
if (is.null(opt$format_data)) {
  opt$format_data = NULL
}
if (is.null(opt$filter_config)) {
  opt$filter_config = NULL
}
if (is.null(opt$transform_config)) {
  opt$transform_config = NULL
}
if (is.null(opt$format_key)) {
  opt$format_key = NULL
}

# --- File Existence Checks ---

# Function to check file existence and stop if not found
check_file_exists <- function(path, arg_name) {
  if (!is.null(path) && !file.exists(path)) {
    stop(sprintf("Error: File specified by --%s not found at path: %s", arg_name, path))
  }
}

# Check REQUIRED file paths
check_file_exists(opt$info_data, 'info_data')
check_file_exists(opt$reference_coords, 'reference_coords')
check_file_exists(opt$info_key, 'info_key')

# Check OPTIONAL file paths (only if they were provided/not NULL)
check_file_exists(opt$format_data, 'format_data')
check_file_exists(opt$filter_config, 'filter_config')
check_file_exists(opt$transform_config, 'transform_config')
check_file_exists(opt$format_key, 'format_key')

# --- Metadata loading ---
info_key <- data.table::fread(opt$info_key)
format_key <- data.table::fread(opt$format_key)
target_desc <- ifelse(opt$source == "INFO", info_key$DESCRIPTION[which(info_key$ID == opt$target)],
                      format_key$DESCRIPTION[which(format_key$ID == opt$target)])
reference_version <- gsub("^.*genome_|.tsv$","",opt$reference_coords) #relies on vast naming convention


# --- Data Filtering Function ---

# Applies filtering rules defined in the config to the INFO and FORMAT data tables.
# Returns a list containing the filtered INFO and FORMAT data tables.
apply_filters <- function(info_dt, format_dt, filter_config_path) {
  
  cat("\n--- Applying Data Filters ---\n")
  
  # 1. Load Filter Configuration (Assumed to be TSV with header)
  if (is.null(filter_config_path)) {
    cat("No filter config path provided. Skipping data filtering.\n")
    return(list(info_dt = info_dt, format_dt = format_dt))
  }
  
  filter_config_dt <- data.table::fread(filter_config_path)
  
  if (nrow(filter_config_dt) == 0) {
    cat("Filter config is empty. Skipping data filtering.\n")
    return(list(info_dt = info_dt, format_dt = format_dt))
  }
  
  # Ensure mandatory columns are present
  required_config_cols <- c("source", "metric", "operation", "threshold")
  if (!all(required_config_cols %in% names(filter_config_dt))) {
    stop("Filter config file must contain columns: source, metric, operation, threshold.")
  }
  
  # Initialize VPK sets. Start with all VPKs from the existing data.
  # Note: VPK must be a character for reliable set intersection.
  passing_info_vpks_set <- as.character(info_dt[["VPK"]])
  passing_format_vpks_set <- if (!is.null(format_dt)) as.character(format_dt[["VPK"]]) else character(0)
  
  initial_info_rows <- nrow(info_dt)
  initial_format_rows <- if (!is.null(format_dt)) nrow(format_dt) else 0
  
  # Iterate over rules and apply filters
  filter_rules = c()
  for (i in 1:nrow(filter_config_dt)) {
    rule <- filter_config_dt[i]
    source <- as.character(rule$source)
    metric <- as.character(rule$metric)
    operation <- as.character(rule$operation)
    threshold <- as.character(rule$threshold)
    
    rule = sprintf("Rule %d: %s %s %s %s\n", i, source, metric, operation, threshold)
    cat(rule)
    filter_rules[i] <- rule
    
    dt_to_filter <- if (source == "INFO") info_dt else format_dt
    
    if (is.null(dt_to_filter)) {
      cat(sprintf("Warning: Skipping rule, %s data table is NULL.\n", source))
      next
    }
    if (!metric %in% names(dt_to_filter)) {
      cat(sprintf("Warning: Metric '%s' not found in %s data. Skipping rule.\n", metric, source))
      next
    }
    
    # --- Dynamic Filtering Logic ---
    
    # Check if the threshold can be safely parsed as numeric
    is_numeric_op <- tryCatch({
      as.numeric(threshold)
      TRUE
    }, warning = function(w) FALSE, error = function(e) FALSE)
    
    
    if (is_numeric_op) {
      # Ensure the column is numeric for comparison
      # Using data.table syntax to change column type in place
      if (!is.numeric(dt_to_filter[[metric]])) {
        dt_to_filter[, (metric) := as.numeric(get(metric))]
      }
      # Numeric comparison expression
      condition_str <- sprintf("%s %s %s", metric, operation, threshold)
    } else {
      # String comparison expression (wrapping threshold in single quotes)
      condition_str <- sprintf("%s %s '%s'", metric, operation, threshold)
    }
    
    # Execute the filter within data.table's context to get passing VPKs
    # Use a local temporary variable to store the result of the subset
    filtered_vpks <- dt_to_filter[eval(parse(text = condition_str)), VPK]
    
    # Update the overall passing set for the source (intersect with previous passes)
    if (source == "INFO") {
      passing_info_vpks_set <- intersect(passing_info_vpks_set, filtered_vpks)
    } else {
      passing_format_vpks_set <- intersect(passing_format_vpks_set, filtered_vpks)
    }
  }
  
  # Final Filtering of the Data Tables
  
  # VPK must pass ALL INFO rules AND ALL FORMAT rules to be kept.
  final_vpks_to_keep <- intersect(passing_info_vpks_set, passing_format_vpks_set)
  
  # Filter INFO data table
  info_dt_filtered <- info_dt[VPK %in% final_vpks_to_keep]
  
  # Filter FORMAT data table
  format_dt_filtered <- NULL
  if (!is.null(format_dt)) {
    format_dt_filtered <- format_dt[VPK %in% final_vpks_to_keep]
  }
  
  cat(sprintf("\nInitial INFO rows: %d. Final INFO rows: %d.\n", initial_info_rows, nrow(info_dt_filtered)))
  cat(sprintf("Initial FORMAT rows: %d. Final FORMAT rows: %d.\n", initial_format_rows, if (!is.null(format_dt_filtered)) nrow(format_dt_filtered) else 0))
  
  return(list(info_dt = info_dt_filtered, format_dt = format_dt_filtered, filter_rules = filter_rules))
}

# --- Plotting Functions ---

# Prepares the data.table for plotting, ensuring all required columns (target, ABS_POS, CHROM) are present.
prepare_plotting_data <- function(primary_dt, info_dt, target_col, source, is_categorical) {
  cat("\n--- Preparing Plotting Data ---\n")
  
  plotting_dt <- NULL
  
  # 1. Select essential columns and ensure CHROM is present
  if (source == "INFO") {
    # INFO table (wide) already contains VPK, ABS_POS, CHROM
    cols_to_select <- c("VPK", "ABS_POS", "CHROM", target_col)
    
    if (!target_col %in% names(primary_dt)) {
      stop(sprintf("Error: Target column '%s' not found in INFO data.", target_col))
    }
    
    plotting_dt <- primary_dt[, ..cols_to_select]
    
  } else if (source == "FORMAT") {
    # FORMAT table (long) needs CHROM from the filtered INFO table (info_dt)
    cols_to_select <- c("VPK", "ABS_POS", "SAMPLE_NAME", target_col)
    
    if (!target_col %in% names(primary_dt)) {
      stop(sprintf("Error: Target column '%s' not found in FORMAT data.", target_col))
    }
    if (!"SAMPLE_NAME" %in% names(primary_dt)) {
      stop("Error: 'SAMPLE_NAME' column is missing from FORMAT data.")
    }
    
    # 1.1 Select FORMAT data columns
    plotting_dt_base <- primary_dt[, ..cols_to_select]
    
    # 1.2 Get VPK and CHROM from filtered INFO data (must use the VPKs that passed all filters)
    chrom_lookup <- info_dt[, c("VPK", "CHROM")]
    
    # 1.3 Merge to add CHROM to the plotting data
    plotting_dt <- merge(plotting_dt_base, chrom_lookup, by = "VPK", all.x = TRUE)
    
    # Remove any rows where CHROM lookup failed (if any VPK was inconsistent, though unlikely)
    plotting_dt <- plotting_dt[!is.na(CHROM)]
  }
  
  # 2. Convert relevant columns to correct types
  
  # Convert ABS_POS to numeric
  plotting_dt[, ABS_POS := as.numeric(ABS_POS)]
  
  # Handle the target column type
  if (!is_categorical) {
    # Coerce to numeric, suppressing warnings if non-numeric values become NA
    plotting_dt[, (target_col) := suppressWarnings(as.numeric(get(target_col)))]
  } else {
    # Coerce to factor for discrete coloring/plotting
    plotting_dt[, (target_col) := as.factor(get(target_col))]
  }
  
  # 3. Apply Chromosome Standardization
  plotting_dt[, CHROM := standardize_chrom(CHROM)]
  
  return(plotting_dt)
}

plot_genomic_data <- function(plotting_dt, 
                              ref_coords_dt, 
                              target_col, 
                              source, 
                              is_categorical, 
                              reference_version,
                              target_desc, 
                              rules) {
  cat("\n--- Generating Genomic Plot ---\n")
  
  # Corrected column renaming: The reference file has 5 columns
  required_names <- c("CHROM", "start", "stop", "length", "midpoint")
  if (ncol(ref_coords_dt) >= 5) {
    # Rename the first 5 columns to ensure consistency and correct types
    setnames(ref_coords_dt, 1:5, required_names)
  } else {
    stop(sprintf("Error: Reference coordinates table must have at least 5 columns. Found %d.", ncol(ref_coords_dt)))
  }
  
  # Convert relevant columns to numeric
  ref_coords_dt[, start := as.numeric(start)]
  ref_coords_dt[, stop := as.numeric(stop)]
  ref_coords_dt[, midpoint := as.numeric(midpoint)]
  
  # Apply Chromosome Standardization to Reference Coords
  ref_coords_dt[, CHROM := standardize_chrom(CHROM)]
  
  # Filter reference coords to only include chromosomes present in the filtered data
  chroms_in_data <- unique(plotting_dt$CHROM)
  ref_coords_dt <- ref_coords_dt[CHROM %in% chroms_in_data]
  
  # Calculate Axis Breaks and Labels
  axis_breaks <- ref_coords_dt$midpoint
  axis_labels <- ref_coords_dt$CHROM
  
  # Determine X-axis limits
  x_min <- min(ref_coords_dt$start)
  x_max <- max(ref_coords_dt$stop)
  
  # Prepare Rectangle Data for Shaded Chromosome Backgrounds
  ref_coords_dt[, Index := 1:.N]
  ref_coords_dt[, AlternatingColor := factor(Index %% 2 == 0, levels = c(FALSE, TRUE))]
  rect_dt <- ref_coords_dt[, .(
    CHROM,
    xmin = start,
    xmax = stop,
    AlternatingColor
  )]
  
  #  Add Chromosome Ordering for Point Coloring
  plotting_dt[, CHROM_ORDERED := factor(CHROM, levels = ref_coords_dt$CHROM)]
  y_label <- target_col
  
  # --- DYNAMIC PLOTTING LOGIC: Determine Facetting ---
  num_samples <- if (source == "FORMAT") length(unique(plotting_dt$SAMPLE_NAME)) else 1
  facet_by_sample <- source == "FORMAT" && num_samples > 1
  
  # Coloring is always by CHROM_ORDERED in the facetting approach
  aes_color <- aes(x = ABS_POS, y = get(target_col), color = CHROM_ORDERED)
  
  # Manual color palette for CHROM coloring
  chrom_color_palette <- rep(RColorBrewer::brewer.pal(8, "Dark2"), length.out = length(unique(plotting_dt$CHROM_ORDERED)))
  
  p <- ggplot(plotting_dt, aes_color) +
    geom_rect(data = rect_dt,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = AlternatingColor), 
              inherit.aes = FALSE, 
              alpha = 0.4) + 
    geom_point(size = 1.5, shape = 20, alpha = 0.8) +
    scale_x_continuous(
      name = sprintf("Genomic Position (%s)", reference_version),
      breaks = axis_breaks,
      labels = axis_labels,
      expand = c(0.01, 0.01),
      limits = c(x_min, x_max)
    ) +
    labs(y = y_label) +
    
    # Apply color scales
    scale_color_manual(values = chrom_color_palette) + 
    scale_fill_manual(values = c("FALSE" = "grey90", "TRUE" = "grey80")) + 
    { if (!is_categorical && is.numeric(plotting_dt[[target_col]])) geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") else NULL } +
    
    # --- FACETTING LOGIC ---
    # Changed scales = "free_y" to scales = "fixed" to force consistent Y-axes.
    { if (facet_by_sample) facet_grid(SAMPLE_NAME ~ ., scales = "fixed") else NULL } +
    
    # General plot theme and cleanup
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0, size = 10, color = "grey40"),
      # Legends: CHROM legend is suppressed here. Fill legend is suppressed in guides below.
      legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "grey30"),
      axis.title.x = element_text(vjust = -0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linetype = "dotted"),
      strip.text.y = element_text(angle = 0, face = "bold") # Format facet labels (Sample Names)
    ) + 
    # Remove both fill and color legends
    guides(fill = "none", color = "none") 
  
  
  filter_rules_string <- paste(gsub("Rule.*: |\\n$", "", rules), collapse = " && ")
  plot_title <- sprintf("Genomic Distribution Plot: %s from %s", target_col, source)
  plot_subtitle <- sprintf("%s\nFilters: %s", target_desc, filter_rules_string)
  p <- p + labs(title = plot_title, subtitle = plot_subtitle)
  
  # Save the plot to a file
  # Plot height remains dynamic based on the number of samples
  plot_height <- if (facet_by_sample) max(6, 3 * num_samples) else 6
  plot_filename <- sprintf("vast_plot_%s_%s.png", source, target_col)
  ggsave(plot_filename, plot = p, width = 12, height = plot_height, units = "in")
  
  cat(sprintf("Plot saved successfully to %s\n", plot_filename))
}

# --- 10. Main Function ---

main <- function(opt) {
  
  # --- Configuration Summary ---
  cat("--- VAST Plotter Configuration ---\n")
  cat(sprintf("INFO Data Path: %s\n", opt$info_data))
  cat(sprintf("FORMAT Data Path: %s\n", if (is.null(opt$format_data)) "NULL" else opt$format_data))
  cat(sprintf("Target Variable: %s\n", opt$target))
  cat(sprintf("Target Data Source: %s\n", opt$source))
  cat(sprintf("Reference Coords Path: %s\n", opt$reference_coords))
  cat(sprintf("Is Target Categorical: %s\n", opt$is_categorical))
  cat(sprintf("Filter Config Path: %s\n", if (is.null(opt$filter_config)) "NULL" else opt$filter_config))
  cat(sprintf("Transform Config Path: %s\n", if (is.null(opt$transform_config)) "NULL" else opt$transform_config))
  cat(sprintf("INFO Key File Path: %s (REQUIRED)\n", opt$info_key))
  cat(sprintf("FORMAT Key File Path: %s\n", if (is.null(opt$format_key)) "NULL" else opt$format_key))
  cat("-----------------------------------\n")
  
  cat("\n--- Loading Data and Configuration Files ---\n")
  
  # 1. Load Data Tables using data.table::fread (assumes TSV format)
  
  # Required data files
  info_data_dt <- data.table::fread(opt$info_data)
  # Ensure VPK is treated as character for consistent joining/filtering
  info_data_dt[, VPK := as.character(VPK)]
  cat(sprintf("Loaded INFO Data: %d rows, %d columns\n", nrow(info_data_dt), ncol(info_data_dt)))
  
  reference_coords_dt <- data.table::fread(opt$reference_coords)
  cat(sprintf("Loaded Reference Coords: %d rows, %d columns\n", nrow(reference_coords_dt), ncol(reference_coords_dt)))
  
  # Optional/Conditionally Required FORMAT data file
  format_data_dt <- NULL
  if (!is.null(opt$format_data)) {
    format_data_dt <- data.table::fread(opt$format_data)
    # Ensure VPK is treated as character
    format_data_dt[, VPK := as.character(VPK)]
    cat(sprintf("Loaded FORMAT Data: %d rows, %d columns\n", nrow(format_data_dt), ncol(format_data_dt)))
  }
  
  cat("\nData loading complete. Starting filtering.\n")
  
  # --- 4. Filtering Step ---
  filtered_data <- apply_filters(
    info_dt = info_data_dt,
    format_dt = format_data_dt,
    filter_config_path = opt$filter_config
  )
  
  # Update the data tables with filtered results
  info_data_dt <- filtered_data$info_dt
  format_data_dt <- filtered_data$format_dt
  
  # 5. Select the primary data table
  primary_data_dt <- if (opt$source == "INFO") info_data_dt else format_data_dt
  
  if (is.null(primary_data_dt) || nrow(primary_data_dt) == 0) {
    stop("\nProcessing aborted: The primary data source is empty after filtering or was not loaded.")
  }
  
  # --- 6. Prepare Plotting Data ---
  plotting_dt <- prepare_plotting_data(
    primary_dt = primary_data_dt,
    info_dt = info_data_dt, # Pass the filtered INFO data for CHROM lookup
    target_col = opt$target,
    source = opt$source,
    is_categorical = opt$is_categorical
  )
  
  # --- 7. Apply transformations (if needed, placeholder) ---
  # if (!is.null(transform_config_content)) {
  #  plotting_dt <- apply_transforms(plotting_dt, transform_config_content)
  # }
  
  # --- 8. Generate Plot ---
  plot_genomic_data(
    plotting_dt = plotting_dt,
    ref_coords_dt = reference_coords_dt,
    target_col = opt$target,
    source = opt$source,
    is_categorical = opt$is_categorical,
    target_desc = target_desc,
    reference_version = reference_version,
    rules = filtered_data$filter_rules
  )
  
  cat("\nData processing and plotting pipeline complete.\n")
}

# Execute the main function with the parsed options
main(opt)