use std::{io::{self, BufReader, BufRead}, fs::File, path::Path, collections::{BTreeMap, BTreeSet}, path::PathBuf};
use polars::prelude::*;
use polars::datatypes::PlSmallStr;
use regex::Regex;
use clap::Parser;


// -----------------------------------------------------------
// 1. PUBLIC TYPE ALIAS AND STRUCTS
// -----------------------------------------------------------

// Custom Result type
pub type VastResult<T> = Result<T, Box<dyn std::error::Error>>;

// CLI Argument Struct
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to the VCF file to be processed
    pub path: PathBuf,

    // Optional: Add a flag to print the dataframes in the run function
    #[arg(long, default_value_t = false)]
    pub print_dfs: bool,
}


// Data Structure for VCF content
pub struct VcfContent {
    /// Lines defining FILTER fields (e.g., ##FILTER=<ID=q10,...>).
    pub filters: Vec<String>,
    /// Lines defining FORMAT fields (e.g., ##FORMAT=<ID=GT,...>).
    pub formats: Vec<String>,
    /// Lines defining INFO fields (e.g., ##INFO=<ID=DP,...>).
    pub infos: Vec<String>,
    /// Lines defining reference contigs (e.g., ##contig=<ID=chr1,...>).
    pub contigs: Vec<String>,
    /// All other ## lines (e.g., ##fileformat, ##source, ##date).
    pub misc: Vec<String>,
    /// Remaining lines, typically the single column header line (#CHROM)
    /// and all subsequent variant records.
    pub data: Vec<String>,
}

// -----------------------------------------------------------
// 2. I/O UTILITY
// -----------------------------------------------------------

/// Reads the entire contents of a file line by line, collecting all lines
/// into a single vector of strings.
pub fn read_all_lines<P: AsRef<Path>>(path: P) -> io::Result<Vec<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let lines: Vec<String> = reader
        .lines()
        .filter_map(|line_result| {
            line_result.ok()
        })
        .collect();

    Ok(lines)
}

// -----------------------------------------------------------
// 3. VCF SPLITTING LOGIC
// -----------------------------------------------------------

/// Splits a vector of raw VCF lines into metadata categories and variant data.
pub fn split_vcf_lines(all_lines: Vec<String>) -> VcfContent {
    let mut filters = Vec::new();
    let mut formats = Vec::new();
    let mut infos = Vec::new();
    let mut contigs = Vec::new();
    let mut misc = Vec::new();
    let mut data = Vec::new();

    for line in all_lines {
        if line.starts_with("##") {
            // Sort the header line based on the exact VCF keyword prefix
            if line.starts_with("##FILTER=<") {
                filters.push(line);
            } else if line.starts_with("##FORMAT=<") {
                formats.push(line);
            } else if line.starts_with("##INFO=<") {
                infos.push(line);
            } else if line.starts_with("##contig=<") {
                contigs.push(line);
            } else {
                // Anything else that starts with ## (fileformat, source, etc.)
                misc.push(line);
            }
        } else {
            // Data or column header line
            data.push(line);
        }
    }

    VcfContent {
        filters,
        formats,
        infos,
        contigs,
        misc,
        data
    }
}

// -----------------------------------------------------------
// 4. DISPLAY UTILITY (re-added)
// -----------------------------------------------------------

/// Prints all header and metadata sections of the VcfContent struct.
pub fn display_vcf_metadata(content: &VcfContent) {
    let print_section = |title: &str, lines: &[String]| {
        println!("\n[ {} ({}) ]", title, lines.len());
        if lines.is_empty() {
            println!("  (None found)");
        } else {
            for line in lines {
                println!("  {}", line);
            }
        }
    };

    println!("--- VCF Metadata Summary ---");

    print_section("INFO Definitions", &content.infos);
    print_section("FILTER Definitions", &content.filters);
    print_section("FORMAT Definitions", &content.formats);
    print_section("CONTIG Definitions", &content.contigs);
    print_section("MISC Header Lines", &content.misc);

    println!("\n[ RAW Data Lines Found: {} ]", content.data.len());
}

/// Helper function to display the generated DataFrames
pub fn display_polars_tables(tables: &BTreeMap<String, DataFrame>) {
    println!("\n--- Polars DataFrame Summary ---");
    for (name, df) in tables {
        println!("\n[ DataFrame: {} ({} rows, {} columns) ]", name, df.height(), df.width());
        println!("{}", df);
    }
}

// -----------------------------------------------------------
// 5. POLARS HELPER LOGIC
// -----------------------------------------------------------
/// Helper function to parse a single header category (e.g., INFO, FORMAT) into a Polars DataFrame.
fn parse_header_category(lines: &[String]) -> Result<DataFrame, PolarsError> {
    let regex_start = Regex::new(r"^##[A-Z_a-z]+=<").expect("Failed to compile start regex");
    let regex_end = Regex::new(r">$").expect("Failed to compile end regex");

    let mut all_rows: Vec<BTreeMap<String, String>> = Vec::new();
    let mut unique_keys: BTreeSet<String> = BTreeSet::new();

    for line in lines {
        // 1. Cleanup line: Remove ##KEYWORD=< and >
        let cleaned_line = regex_end.replace(
            &regex_start.replace(line, ""),
            ""
        ).to_string();

        // 2. Split by comma and process key=value pairs
        let parts = cleaned_line.split(',');
        let mut row: BTreeMap<String, String> = BTreeMap::new();

        for part in parts {
            // Find the first equals sign to split key and value
            if let Some((key, value)) = part.split_once('=') {
                let key_str = key.trim().to_string();
                // Remove surrounding quotes from the value (e.g., from "Read Depth")
                let value_str = value.trim().trim_matches('"').to_string();

                unique_keys.insert(key_str.clone());
                row.insert(key_str, value_str);
            }
        }
        all_rows.push(row);
    }

    // 3. Build Polars Series (and convert to Column)
    let mut ordered_keys: Vec<String> = Vec::new();
    let id_key = "ID".to_string();

    // 3a. If 'ID' exists, push it first
    if unique_keys.contains(&id_key) {
        ordered_keys.push(id_key.clone());
    }

    // 3b. Push all other keys from the BTreeSet (maintaining their sorted order)
    for key in unique_keys.iter() {
        if key != &id_key {
            ordered_keys.push(key.clone());
        }
    }

    // 3c. Build Polars Columns from the new, ordered key list
    let columns: Vec<Column> = ordered_keys.iter()
        .map(|key| {
            let col_data: Vec<Option<String>> = all_rows.iter()
                .map(|row| row.get(key).cloned())
                .collect();

            // Create the Series
            let series = Series::new(
                PlSmallStr::from(key.to_string()),
                col_data);

            // Convert the Series into a Column
            series.into()

        })
        .collect();

    // 4. Create DataFrame
    DataFrame::new(columns)
}

// -----------------------------------------------------------
// 6. MAIN GENERATOR FUNCTION
// -----------------------------------------------------------

/// Generates a Polars DataFrame for each categorized VCF header section.
///
/// # Returns
/// A map of DataFrames, keyed by the VCF section name (e.g., "INFO", "FORMAT").
pub fn generate_polars_tables(
    content: &VcfContent,
) -> Result<BTreeMap<String, DataFrame>, Box<dyn std::error::Error>> {
    let mut tables = BTreeMap::new();

    // Process each categorized header section, propagating errors if DataFrame creation fails
    tables.insert(
        "INFO".to_string(),
        parse_header_category(&content.infos)?
    );
    tables.insert(
        "FORMAT".to_string(),
        parse_header_category(&content.formats)?
    );
    tables.insert(
        "FILTER".to_string(),
        parse_header_category(&content.filters)?
    );
    tables.insert(
        "CONTIG".to_string(),
        parse_header_category(&content.contigs)?
    );

    Ok(tables)
}

pub fn run(cli: Cli) -> VastResult<()> {

    println!("Processing VCF file: {}", cli.path.display());

    // 1. Read all lines from the VCF file
    let all_lines = read_all_lines(&cli.path)?;

    // 2. Split the lines into metadata and data
    let vcf_content = split_vcf_lines(all_lines);

    // 3. Display the raw metadata summary
    display_vcf_metadata(&vcf_content);

    // 4. Generate Polars DataFrames from the VCF metadata
    let polars_tables = generate_polars_tables(&vcf_content)?;


    display_polars_tables(&polars_tables);


    println!("\nSuccessfully parsed VCF metadata into {} Polars DataFrames.", polars_tables.len());

    Ok(())
}
