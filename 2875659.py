# Importing necessary modules
import re
import sys

# Checking if the correct number of command line arguments is provided
if len(sys.argv) != 3:
    print("Error: Both SAM file and Genes Summary file not provided as command line arguments.")
    sys.exit(1)  # Exiting the script with an error code if the expected number of arguments is not provided
    
# Extracting command line arguments
sam_file = sys.argv[1]  # SAM file path provided as the first command-line argument
genes_file = sys.argv[2]  # Genes file path provided as the second command-line argument

# Creating a function to process the CIGAR string
def cigar_process(cigar_string, pos, junction_counts):
    # Ensuring the correct type in variables
    assert isinstance(cigar_string, str), "Invalid cigar_string type"
    assert isinstance(pos, int), "Invalid pos type"
    assert isinstance(junction_counts, dict), "Invalid junction_counts type"

    position = pos # pos is the position where the alignment starts

    # Writing the regular expression pattern to extract junction information
    # Only including the Match(M), Deletions(D), and Skipped Region denoting Introns(N)
    for match in re.finditer(r'(\d+)([MDN])', cigar_string):
        # Assigning variables to each group in the pattern
        num = int(match.group(1)) # The number that would be added if the subsequent symbol is relevant
        sym = match.group(2) # The symbol following the digit

        # Checking if the number is an integer and if the symbol is contained in the required list
        assert isinstance(num, int), "Invalid num type"
        assert sym in ['M', 'D', 'N'], "Invalid sym type"

        if sym == 'N': # N denotes the skipped region
            # Updating the start and the end of the junction
            intron_start = position 
            intron_end = position + num
            
            
            # Represents a tuple containing the start and end positions of a potential junction identified in the alignment
            intron_key = (intron_start, intron_end)

            # Checking if junction has already been encountered
            # If yes, incrementing the corresponding count value by 1
            if intron_key in junction_counts:
                junction_counts[intron_key] += 1
            # If not, creating a new value and setting count as 1
            else:
                junction_counts[intron_key] = 1
            # updating the genomic position (position) by adding the length of the current segment (num) in the CIGAR string   
            position += num

        if sym in ['M', 'D']:
            position += num
            
# Checking if the function works
# Sample CIGAR string
sample_cigar = "5M10N15M"

# Initialising junction_counts dictionary
junction_counts = {}

# Calling the cigar function
cigar_process(sample_cigar, 0, junction_counts)

# Assertions to check the function with the sample cigar string
assert len(junction_counts) == 1, "Number of junctions unexpected"
assert (5, 15) in junction_counts, "Junction key not found"
assert junction_counts[(5, 15)] == 1, "Incorrect read count"

# Creating a function to parse the Genes Summary file
def gene_parsing(line):
    
    # Checking line type
    assert isinstance(line, str), "Invalid line type"
    
    # Setting variables for each column
    # gene_loc has the locations of the start and end positions of the gene
    gene_id, gene, gene_loc = line.rstrip().split("\t")
    
    # Ensuring correct variable type
    assert isinstance(gene_id, str), "Invalid gene_id type"
    assert isinstance(gene, str), "Invalid gene type"
    assert isinstance(gene_loc, str), "Invalid gene_loc type"
    
    # Parsing the third column to get the gene start and end and checking the correct number of splits
    gene_loc_parts = gene_loc.split(":")
    assert len(gene_loc_parts) == 2, "Invalid gene_loc format"

    gene_position = gene_loc_parts[1].split("..")
    assert len(gene_position) == 2, "Invalid gene_loc format"

    gene_start = int(gene_position[0].replace(",", ""))
    gene_end = int(gene_position[1].replace(",", "").removesuffix("(+)").removesuffix("(-)"))

    # Checkcing if the start and end positions are integers
    assert isinstance(gene_start, int), "Invalid gene_start type"
    assert isinstance(gene_end, int), "Invalid gene_end type"

    return gene_id, gene_start, gene_end


# Creating a function to check if the intron junctions are within the genes
def junctions_in_range(junction_counts, gene_start, gene_end):
    # Dictionary to store sorted junctions
    junctions_range = {}

    # Iterating through junction_counts
    for junction, read_count in junction_counts.items():
        start, end = junction
        # Checking if the junction is within the specified gene range and making the read counts as the value
        if (gene_start <= start <= gene_end) and (gene_start <= end <= gene_end):
            junctions_range[junction] = read_count

    # Returning the sorted junctions which are within the genes
    return junctions_range

# Dictionary to store junction counts
junction_counts = {}

try:
    # Attempting to open the SAM file
    with open(sam_file) as sam:
        for line in sam:
            # Skipping header lines starting with "@"
            if line.startswith("@"):
                continue
            else:
                line = line.rstrip().split("\t")
                rname = line[2]
                pos = int(line[3])
                cigar_string = line[5]
                nh = line[-1]

                # Checking for the presence of introns in the CIGAR string and reads aligning just once
                if "NH:i:1" in nh and "N" in cigar_string:
                    # Processing the CIGAR string and updating junction_counts
                    cigar_process(cigar_string, pos, junction_counts)

# Handling the case where the SAM file is not found
except FileNotFoundError:
    print(f"Error: File not found - {sam_file}")

# Creating an output file to store junction information
with open("2875659P", "w") as output_file:
    # Writing the header 
    output_file.write(f'Gene_ID\tJunction_Start\tJunction_End\tRead_Counts\n')

    try: 
        # Attempting to open the genes file
        with open(genes_file) as text:
            # Skipping the header line in the genes file
            next(text)

            # Iterating through lines in the genes file
            for line in text:
                # Parsing gene information from the line by calling the function
                gene_id, gene_start, gene_end = gene_parsing(line)

                # Finding junctions for the current gene
                final_junctions = junctions_in_range(junction_counts, gene_start, gene_end)

                # Iterating through the final junctions and (printing information)
                for junction, read_count in final_junctions.items():
                    start, end = junction
                    #print(f"Gene ID: {gene_id}, Junction Start: {start}, Junction End: {end}, READ_COUNT: {read_count}")
                    # Writing the output file
                    output_file.write(f'{gene_id}\t{start}\t{end}\t{read_count}\n')

                # Adding a newline after each gene
                output_file.write('\n')

    # Handling the case where the Gene Summary file is not found
    except FileNotFoundError:
        print(f"Error: File not found - {genes_file}")