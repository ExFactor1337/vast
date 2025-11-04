use std::{path::PathBuf, fs::File};
use clap::Parser;
use std::io::{self, BufReader, BufRead, Write};
use std::collections::{HashSet, HashMap};
// -----------------------------------------------------------
// 1. PUBLIC TYPE ALIAS AND STRUCTS
// -----------------------------------------------------------
// Helper function to convert a standard error into the boxed trait object required by VastResult
fn box_err<E: std::error::Error + Send + Sync + 'static>(e: E) -> Box<dyn std::error::Error> {
    Box::new(e)
}

// Custom Result type
pub type VastResult<T> = Result<T, Box<dyn std::error::Error>>;

// CLI Argument Struct
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to the VCF file to be processed
    pub path: PathBuf,

    #[arg(long, default_value_t = false)]
    pub print_dfs: bool,
}

/// Holds all parsed content from a VCF file, separating metadata from data.
pub struct VcfContent {
    // Header metadata lines starting with '##'
    pub contigs: Vec<String>,
    pub filters: Vec<String>,
    pub formats: Vec<String>,
    pub infos: Vec<String>,
    pub misc: Vec<String>,

    // The single column header line (e.g., #CHROM POS ID...)
    pub column_header: String,

    // The VCF data/variant lines (everything after the column_header)
    pub variants: Vec<String>,
}

pub fn read_and_split_vcf(path: &PathBuf) -> VastResult<VcfContent> {
    // 1. Open the file efficiently
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut contigs = Vec::new();
    let mut filters = Vec::new();
    let mut formats = Vec::new();
    let mut infos = Vec::new();
    let mut misc = Vec::new();
    let mut column_header = String::new();
    let mut variants = Vec::new();

    let mut header_complete = false;

    // 2. Read and categorize lines
    for line_result in reader.lines() {
        let line = line_result?;

        if line.starts_with("##") {
            if header_complete {
                // FIX 1: Manually create io::Error and box it
                let err = io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Found '##' header line after variant data began."
                );
                return Err(box_err(err));
            }

            // Sort header lines based on prefix (simple and direct)
            if line.starts_with("##contig=<") {
                contigs.push(line);
            } else if line.starts_with("##FILTER=<") {
                filters.push(line);
            } else if line.starts_with("##FORMAT=<") {
                formats.push(line);
            } else if line.starts_with("##INFO=<") {
                infos.push(line);
            } else {
                misc.push(line);
            }
        } else if line.starts_with("#") {
            if header_complete {
                // FIX 2: Manually create io::Error and box it
                let err = io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Found second column header line."
                );
                return Err(box_err(err));
            }
            // This is the column header line
            column_header = line;
            header_complete = true;
        } else if !line.trim().is_empty() {
            // All non-empty lines that don't start with '#' are variant data
            variants.push(line);
        }
    }
    // 3. Simple validation (ensure column header was found)
    if column_header.is_empty() {
        // FIX 3: Manually create io::Error and box it
        let err = io::Error::new(
            io::ErrorKind::InvalidData,
            "VCF file is missing the required column header line ('#CHROM...')."
        );
        return Err(box_err(err));
    }

    // 4. Return the structured content
    Ok(VcfContent {
        contigs,
        filters,
        formats,
        infos,
        misc,
        column_header,
        variants,
    })
}

// Helper function to inspect and print the contents of the VcfContent struct
fn inspect_vcf_content(content: &VcfContent) {
    println!("\nVCF Content Inspection:");

    // Header Counts
    println!("  - Column Header: {}", content.column_header);
    println!("  - Contig lines: {}", content.contigs.len());
    println!("  - Filter lines: {}", content.filters.len());
    println!("  - Format lines: {}", content.formats.len());
    println!("  - Info lines:   {}", content.infos.len());
    println!("  - Misc lines:   {}", content.misc.len());

    // Variant/Data Count
    println!("  - **VARIANT LINES**: {}", content.variants.len());
}

// Helper function to extract values from a single VCF meta-information line
fn parse_meta_line(line: &String) -> Vec<String> {
    // 1. Find the content inside the angled brackets: <...>
    if let Some(start_index) = line.find('<') {
        if let Some(end_index) = line.rfind('>') {
            let content = &line[start_index + 1..end_index];

            let mut values = Vec::new();
            let mut current_pair = String::new();
            let mut in_quotes = false;

            // Simple state machine to handle quoted commas
            for char in content.chars() {
                match char {
                    // Toggle the state when a quote is encountered
                    '"' => {
                        in_quotes = !in_quotes;
                        current_pair.push(char);
                    }
                    // Split only if outside of quotes
                    ',' if !in_quotes => {
                        // Process the complete key=value pair
                        if let Some((_, value)) = current_pair.as_str().split_once('=') {
                            // Trim quotes and push the value
                            values.push(value.trim_matches('"').to_string());
                        } else {
                            // Handle malformed pair before comma
                            values.push("".to_string());
                        }
                        current_pair.clear();
                    }
                    _ => {
                        current_pair.push(char);
                    }
                }
            }
            // Process the final pair after the loop finishes
            if !current_pair.is_empty() {
                if let Some((_, value)) = current_pair.as_str().split_once('=') {
                    values.push(value.trim_matches('"').to_string());
                } else {
                    values.push("".to_string());
                }
            }

            return values;
        }
    }
    // Return empty if the line is not a valid meta-information format
    Vec::new()
}

// The main function to format the header vectors
pub fn format_header_vecs(header_lines: Vec<String>) -> Vec<Vec<String>> {
    let mut table: Vec<Vec<String>> = Vec::new();

    for line in header_lines {
        // Call the helper to extract the values from the current line
        let values = parse_meta_line(&line);

        if !values.is_empty() {
            table.push(values);
        }
    }
    table
}

/// Parses a key-value string (like the INFO or AD fields) into a HashMap.
/// Takes a string and two delimiters (one for pair separation, one for key/value separation).
fn split_variant_data(data: &str, pair_delimiter: char, kv_delimiter: char) -> HashMap<String, String> {
    let mut map = HashMap::new();

    // INFO fields are usually KEY=VALUE;KEY=VALUE
    // FORMAT/Sample fields are KEY:VALUE:KEY:VALUE
    for pair in data.split(pair_delimiter) {
        if let Some((key, value)) = pair.split_once(kv_delimiter) {
            // Remove leading/trailing whitespace and store
            map.insert(key.to_string(), value.trim().to_string());
        } else if !pair.is_empty() {
            // Handle flag fields (e.g., 'NT') where there is no '='
            // Set the value to 'FLAG' as requested
            map.insert(pair.to_string(), "FLAG".to_string());
        }
    }
    map
}

pub type VariantTables = (
    Vec<String>,        // Wide Variant Headers (VPK + Fixed + INFO)
    Vec<Vec<String>>,   // Wide Variant Table Data (INFO)
    Vec<String>,        // Long Sample Headers (VPK + Sample Name + FORMAT keys)
    Vec<Vec<String>>,   // Long Sample Table Data (FORMAT/Sample)
);

// A tuple for returning the column headers and the variant table
pub fn process_variants_to_table(column_header: &str, variants: Vec<String>) -> VastResult<VariantTables> {
    let mut info_data_maps: Vec<HashMap<String, String>> = Vec::new();
    let mut all_info_keys: HashSet<String> = HashSet::new();

    let mut all_format_keys: HashSet<String> = HashSet::new();

    // Changed from variant_ids to vpk_indices
    let mut vpk_indices: Vec<String> = Vec::new();

    // Get Sample names and the first 9 fixed column names
    let cols: Vec<&str> = column_header.trim().split('\t').collect();
    // Sample names start at column index 9
    let sample_names: Vec<&str> = cols.get(9..).unwrap_or(&[]).to_vec();

    // --- First Pass: Parse all rows and collect all unique INFO and FORMAT keys ---
    for (i, variant_line) in variants.iter().enumerate() {
        let fields: Vec<&str> = variant_line.split('\t').collect();

        // Ensure the line has enough fixed fields (CHROM, POS, REF, ALT)
        if fields.len() < 5 {
            continue; // Skip malformed lines
        }

        // Generate VPK (Variant Primary Key) as a string index
        vpk_indices.push((i + 1).to_string());

        // INFO field is at index 7
        if fields.len() > 7 {
            let info_data = fields[7];

            // FIX: If the INFO field is the missing value indicator ('.'), skip parsing it.
            if info_data != "." {
                let info_map = split_variant_data(info_data, ';', '=');
                for key in info_map.keys() {
                    all_info_keys.insert(key.clone());
                }
                info_data_maps.push(info_map);
            } else {
                // If INFO is '.', push an empty map so no '.' key is generated.
                info_data_maps.push(HashMap::new());
            }
        } else {
            // Handle variants without an INFO column (field count too low)
            info_data_maps.push(HashMap::new());
        }

        // FORMAT field is at index 8; Sample data starts at index 9
        if fields.len() > 9 {
            let format_keys: Vec<&str> = fields[8].split(':').collect();

            // Collect all unique FORMAT keys
            for key in format_keys.iter() {
                all_format_keys.insert(key.to_string());
            }
        }
    }

    // --- Second Pass: Build the Wide INFO Table (Variant Level) ---
    let fixed_headers: Vec<String> = vec!["VPK".to_string(),
                                          "CHROM".to_string(),
                                          "POS".to_string(),
                                          "REF".to_string(),
                                          "ALT".to_string(),
                                          "QUAL".to_string(),
                                          "FILTER".to_string()];
    let mut wide_info_headers = fixed_headers.clone();

    let mut info_headers: Vec<String> = all_info_keys.into_iter().collect();
    (&mut info_headers).sort_unstable();
    wide_info_headers.extend(info_headers.clone());

    let mut wide_info_table: Vec<Vec<String>> = Vec::new();

    for i in 0..variants.len() {
        let fields: Vec<&str> = variants[i].split('\t').collect();
        let current_map = &info_data_maps[i];

        // Start the row with fixed columns + VPK
        let mut row: Vec<String> = vec![
            vpk_indices[i].clone(), // VPK
            fields[0].to_string(), // CHROM
            fields[1].to_string(), // POS
            fields[3].to_string(), // REF
            fields[4].to_string(), // ALT
            fields[5].to_string(), // QUAL
            fields[6].to_string(), // FILTER
        ];

        // Append INFO columns, inserting "NA" for missing
        for header in &info_headers {
            match current_map.get(header) {
                Some(value) => row.push(value.to_string()),
                None => row.push("NA".to_string()),
            }
        }
        wide_info_table.push(row);
    }

    // --- Third Pass: Build the Long Sample Table (Sample Level) ---

    // Headers for the long table: VPK, SAMPLE_NAME, and all FORMAT keys
    let mut long_sample_headers: Vec<String> = vec!["VPK".to_string(), "SAMPLE_NAME".to_string()];
    let format_headers_vec: Vec<String> = {
        let mut v: Vec<String> = all_format_keys.into_iter().collect();
        (&mut v).sort_unstable();
        long_sample_headers.extend(v.clone());
        v
    };

    let mut long_sample_table: Vec<Vec<String>> = Vec::new();

    // We must re-iterate over the original variants to process samples.
    // We use the original variant lines from the parameter, which are consumed by `into_iter()`.
    // Since the first pass also iterated over references `&variants`, the VPK indices are correctly aligned.
    for (var_index, variant_line) in variants.into_iter().enumerate() {
        let fields: Vec<&str> = variant_line.split('\t').collect();

        // Only proceed if there is FORMAT data and at least one sample
        if fields.len() >= 10 {
            let format_keys: Vec<&str> = fields[8].split(':').collect();

            // Iterate over all sample columns for this variant
            for (sample_index, sample_name) in sample_names.iter().enumerate() {

                // Sample data starts at index 9
                if let Some(sample_data_str) = fields.get(9 + sample_index) {
                    let format_values: Vec<&str> = sample_data_str.split(':').collect();

                    let mut row: Vec<String> = vec![
                        vpk_indices[var_index].clone(), // VPK (Link to wide table)
                        sample_name.to_string(),        // SAMPLE_NAME
                    ];

                    // Create a map for easy lookup of this sample's data
                    let sample_map: HashMap<&str, &str> = format_keys.iter()
                        .zip(format_values.iter())
                        .map(|(k, v)| (*k, *v))
                        .collect();

                    // Append FORMAT values, inserting "NA" for missing
                    for header in &format_headers_vec {
                        match sample_map.get(header.as_str()) {
                            Some(value) => row.push(value.to_string()),
                            None => row.push("NA".to_string()),
                        }
                    }
                    long_sample_table.push(row);
                }
            }
        }
    }

    Ok((wide_info_headers, wide_info_table, long_sample_headers, long_sample_table))
}

// Helper function to write a table (Vec<Vec<String>>) to a TSV file
fn write_table_to_tsv(
    filename: &str,
    headers: Option<&Vec<String>>,
    data: &Vec<Vec<String>>
) -> VastResult<()> {

    let mut file = File::create(filename)?;

    // Write headers if provided
    if let Some(h) = headers {
        file.write_all(h.join("\t").as_bytes())?;
        file.write_all(b"\n")?;
    }

    // Write data rows
    for row in data.iter() {
        file.write_all(row.join("\t").as_bytes())?;
        file.write_all(b"\n")?;
    }

    println!("-> Successfully wrote {} rows to {}", data.len(), filename);
    Ok(())
}

fn write_raw_lines_to_file(
    filename: &str,
    lines: &Vec<String>
) -> VastResult<()> {

    let mut file = File::create(filename)?;

    // Write raw lines
    for line in lines.iter() {
        // Just write the line content and a newline
        file.write_all(line.as_bytes())?;
        file.write_all(b"\n")?;
    }

    println!("-> Successfully wrote {} raw lines to {}", lines.len(), filename);
    Ok(())
}


//noinspection ALL
pub fn run(cli: Cli) -> VastResult<()> {

    println!("\nStarting processing for file: {}", (&cli.path).display());

    // 1. Extract file stem for naming output files
    let file_stem = cli.path.file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| box_err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Could not extract file stem from path for output file naming."
        )))?;

    // 2. Read and split the VCF file.
    let vcf_content = read_and_split_vcf(&cli.path)?;
    inspect_vcf_content(&vcf_content);

    // --- Header Table Processing ---

    // Process all parsable header metadata into structured tables
    let contig_table = format_header_vecs(vcf_content.contigs);
    let info_table = format_header_vecs(vcf_content.infos);
    let format_table = format_header_vecs(vcf_content.formats);
    let filter_table = format_header_vecs(vcf_content.filters);

    // The misc lines cannot be reliably parsed into a table, so we keep the raw lines
    let raw_misc_lines = vcf_content.misc;

    // --- Variant Processing ---

    println!("\n--- Variant Processing ---");
    println!("Processing {} variant lines...", vcf_content.variants.len());

    let (
        wide_info_headers,      // 0: Wide Variant Headers (VPK + Fixed + INFO)
        wide_info_table,        // 1: Wide Variant Table Data
        long_sample_headers,    // 2: Long Sample Headers (VPK + Sample Name + FORMAT keys)
        long_sample_table,      // 3: Long Sample Table Data
    ) = process_variants_to_table(
        &vcf_content.column_header,
        vcf_content.variants
    )?;

    // --- Output File Writing ---

    println!("\n--- Writing Output Files (.tsv) ---");

    // Standard Header Column Names for inspection/writing (where applicable)
    let id_len_assembly_headers = vec!["ID".to_string(), "LENGTH".to_string(), "ASSEMBLY".to_string()];
    let id_number_type_desc_headers = vec!["ID".to_string(), "NUMBER".to_string(), "TYPE".to_string(), "DESCRIPTION".to_string()];
    let id_desc_headers = vec!["ID".to_string(), "DESCRIPTION".to_string()];

    // 1. contigs -> <vcf name>_contigs.tsv
    write_table_to_tsv(
        &format!("{}_contigs.tsv", file_stem),
        Some(&id_len_assembly_headers),
        &contig_table
    )?;

    // 2. filters -> <vcf name>_filters.tsv
    write_table_to_tsv(
        &format!("{}_filters.tsv", file_stem),
        Some(&id_desc_headers),
        &filter_table
    )?;

    // --- NEW: VCF Metadata Definitions ---

    // 3. infos -> <vcf name>_infos.tsv (NEW)
    write_table_to_tsv(
        &format!("{}_infos.tsv", file_stem),
        Some(&id_number_type_desc_headers),
        &info_table
    )?;

    // 4. formats -> <vcf name>_formats.tsv (NEW)
    write_table_to_tsv(
        &format!("{}_formats.tsv", file_stem),
        Some(&id_number_type_desc_headers),
        &format_table
    )?;

    // --- End NEW: VCF Metadata Definitions ---


    // 5. misc (RAW LINES) -> <vcf name>_misc.tsv
    // Writes raw lines, one per row, in a single unnamed column
    write_raw_lines_to_file(
        &format!("{}_misc.tsv", file_stem),
        &raw_misc_lines
    )?;

    // 6. INFO Variant Table (Wide format) -> <vcf_name>_infos_variants.tsv
    write_table_to_tsv(
        &format!("{}_infos_variants.tsv", file_stem),
        Some(&wide_info_headers),
        &wide_info_table
    )?;

    // 7. FORMAT/Sample Table (Long format) -> <vcf_name>_formats_variants.tsv
    write_table_to_tsv(
        &format!("{}_formats_variants.tsv", file_stem),
        Some(&long_sample_headers),
        &long_sample_table
    )?;

    // --- Inspection of Processed Data (Condensed) ---

    let rows_to_show = 5;

    println!("\n--- Inspection of Processed Metadata ---");

    // Contig Inspection
    println!("\nTabular Contig Data (First 3 rows):");
    if contig_table.get(0).is_some() {
        println!("{}", id_len_assembly_headers.join("\t"));
        println!("--------------------------------------------------");
        for row in contig_table.iter().take(3) {
            println!("{}", row.join("\t"));
        }
    } else {
        println!("(No contig lines found or parsed.)");
    }

    // INFO Inspection
    println!("\nTabular INFO Data (First 3 rows):");
    if info_table.get(0).is_some() {
        println!("{}", id_number_type_desc_headers.join("\t"));
        println!("--------------------------------------------------");
        for row in info_table.iter().take(3) {
            println!("{}", row.join("\t"));
        }
    } else {
        println!("(No INFO lines found or parsed.)");
    }

    // FORMAT Inspection
    println!("\nTabular FORMAT Data (First 3 rows):");
    if format_table.get(0).is_some() {
        println!("{}", id_number_type_desc_headers.join("\t"));
        println!("--------------------------------------------------");
        for row in format_table.iter().take(3) {
            println!("{}", row.join("\t"));
        }
    } else {
        println!("(No FORMAT lines found or parsed.)");
    }

    // FILTER Inspection
    println!("\nTabular FILTER Data (First 3 rows):");
    if filter_table.get(0).is_some() {
        println!("{}", id_desc_headers.join("\t"));
        println!("--------------------------------------------------");
        for row in filter_table.iter().take(3) {
            println!("{}", row.join("\t"));
        }
    } else {
        println!("(No FILTER lines found or parsed.)");
    }

    // MISC Inspection (Now printing raw lines)
    println!("\nRaw MISC Lines (Total lines: {}):", raw_misc_lines.len());
    if raw_misc_lines.get(0).is_some() {
        println!("--------------------------------------------------");
        for line in raw_misc_lines.iter().take(3) {
            println!("{}", line);
        }
    } else {
        println!("(No MISC lines found.)");
    }

    // --- Inspection of WIDE INFO Table (Variant Level) ---
    println!("\nWide INFO Table Inspection (Variant Level):");

    // 1. Print all Headers/Columns
    println!("Total INFO Columns: {}", wide_info_headers.len());
    // Updated header description to reflect VPK
    println!("Headers (Starting with VPK): {}", wide_info_headers.join("\t"));
    println!("--------------------------------------------------");

    // 2. Print the first 5 Rows
    let mut rows_printed = 0;
    for row in wide_info_table.iter().take(rows_to_show) {
        println!("{}", row.join("\t"));
        rows_printed += 1;
    }

    // 3. Print Summary
    let total_rows = wide_info_table.len();
    if total_rows > rows_printed {
        println!("[... {} more rows not shown]", total_rows - rows_printed);
    }
    println!("Total Unique Variants Processed: {}", total_rows);


    // --- Inspection of LONG SAMPLE Table (Genotype Level) ---
    println!("\nLong Sample Table Inspection (Genotype Level):");

    // 1. Print all Headers/Columns
    println!("Total Sample Columns: {}", long_sample_headers.len());
    // Updated header description to reflect VPK
    println!("Headers (Starting with VPK, SAMPLE_NAME): {}", long_sample_headers.join("\t"));
    println!("--------------------------------------------------");

    // 2. Print the first 5 Rows
    rows_printed = 0;
    for row in long_sample_table.iter().take(rows_to_show) {
        println!("{}", row.join("\t"));
        rows_printed += 1;
    }

    // 3. Print Summary
    let total_sample_rows = long_sample_table.len();
    if total_sample_rows > rows_printed {
        println!("[... {} more rows not shown]", total_sample_rows - rows_printed);
    }
    println!("Total Genotype Records Processed: {}", total_sample_rows);

    println!("\nVCF file successfully read, structured, and saved as TSV files.");

    Ok(())
}
