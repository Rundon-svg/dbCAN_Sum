import re
import argparse
import glob
from collections import defaultdict, Counter
import os

# Predefined gene categories and their order
CATEGORIES = ["AA", "CBM", "CE", "GH", "GT", "PL"]

# Parse results for each tool
def parse_results(result, tool):
    if result == "N":
        return []
    families = []
    if tool == "HMMER" or tool == "DIAMOND":
        parts = result.split('+')
        for part in parts:
            match = re.match(r"([A-Z]+)\d+", part)
            if match:
                families.append(match.group(0))
    elif tool == "dbCAN_sub":
        parts = result.split('+')
        for part in parts:
            match = re.match(r"([A-Z]+)\d+_e\d+", part)
            if match:
                families.append(re.match(r"([A-Z]+)\d+", match.group(0)).group(0))
    return families

# Select the final gene family subclass
def decide_family(hmmer, dbcan, diamond):
    combined = hmmer + dbcan + diamond
    counter = Counter(combined)
    if len(counter) == 1:
        return list(counter.keys())[0]

    # Use frequency count to determine the most common family
    most_common_families = counter.most_common()
    if most_common_families[0][1] != most_common_families[1][1]:
        return most_common_families[0][0]

    # Prioritize HMMER and dbCAN results
    combined_hmmer_dbcan = hmmer + dbcan
    counter_hmmer_dbcan = Counter(combined_hmmer_dbcan)
    most_common_hmmer_dbcan = counter_hmmer_dbcan.most_common()
    if len(most_common_hmmer_dbcan) > 1 and most_common_hmmer_dbcan[0][1] == most_common_hmmer_dbcan[1][1]:
        return dbcan[0] if dbcan else most_common_hmmer_dbcan[0][0]
    else:
        return most_common_hmmer_dbcan[0][0]

# Process gene data from file
def process_file(input_file):
    family_counter = Counter()
    category_counter = Counter({cat: 0 for cat in CATEGORIES})
    with open(input_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            columns = line.strip().split('\t')
            gene_id = columns[header.index("Gene ID")]
            hmmer_result = columns[header.index("HMMER")]
            dbcan_result = columns[header.index("dbCAN_sub")]
            diamond_result = columns[header.index("DIAMOND")]
            tools_count = int(columns[header.index("#ofTools")])

            # Skip data with only one tool
            if tools_count < 2:
                continue

            # Extract family information
            hmmer_families = parse_results(hmmer_result, "HMMER")
            dbcan_families = parse_results(dbcan_result, "dbCAN_sub")
            diamond_families = parse_results(diamond_result, "DIAMOND")

            # Decide the final family based on rules
            final_family = decide_family(hmmer_families, dbcan_families, diamond_families)
            family_counter[final_family] += 1
            main_category = re.match(r"([A-Z]+)", final_family).group(1)
            category_counter[main_category] += 1

    return family_counter, category_counter

# Aggregate results from all input files
def aggregate_results(files):
    # Dictionary to store family and category counts per species
    aggregated_family_counts = defaultdict(lambda: Counter())
    aggregated_category_counts = defaultdict(lambda: Counter({cat: 0 for cat in CATEGORIES}))
    
    # Process each file
    for file_path in files:
        # Use filename (without path or extension) as species identifier
        species_name = os.path.splitext(os.path.basename(file_path))[0]
        family_counter, category_counter = process_file(file_path)
        
        # Update statistics
        for family, count in family_counter.items():
            aggregated_family_counts[family][species_name] = count
        for category, count in category_counter.items():
            aggregated_category_counts[category][species_name] = count

    # Output aggregated results to file
    with open("aggregated_family_category_statistics.txt", 'w') as f:
        # Output family statistics with families as rows, species as columns
        f.write("Family\t" + "\t".join(os.path.splitext(os.path.basename(fp))[0] for fp in files) + "\n")
        sorted_families = sorted(aggregated_family_counts.keys(), key=lambda x: (CATEGORIES.index(re.match(r"([A-Z]+)", x).group(1)), int(re.search(r"\d+", x).group(0))))
        for family in sorted_families:
            counts = [str(aggregated_family_counts[family].get(species, 0)) for species in [os.path.splitext(os.path.basename(fp))[0] for fp in files]]
            f.write(f"{family}\t" + "\t".join(counts) + "\n")
        f.write("\n")

        # Output category statistics with categories as rows, species as columns
        f.write("Category\t" + "\t".join(os.path.splitext(os.path.basename(fp))[0] for fp in files) + "\n")
        for category in CATEGORIES:
            counts = [str(aggregated_category_counts[category].get(species, 0)) for species in [os.path.splitext(os.path.basename(fp))[0] for fp in files]]
            f.write(f"{category}\t" + "\t".join(counts) + "\n")

        # Total gene count for each species
        f.write("\nTotal\t" + "\t".join(str(sum(aggregated_category_counts[cat].get(species, 0) for cat in CATEGORIES)) for species in [os.path.splitext(os.path.basename(fp))[0] for fp in files]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze gene function prediction results for multiple species.")
    parser.add_argument("-in", dest="input_files", nargs='+', required=True, help="List of txt files to analyze, wildcards supported.")
    args = parser.parse_args()

    # Get list of files based on wildcard patterns
    files = []
    for pattern in args.input_files:
        files.extend(glob.glob(pattern))

    # Run analysis and generate summary results
    aggregate_results(files)
