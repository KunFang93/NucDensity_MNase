#!/bin/bash
prefix=$1
chromsize=$2
# Loop through chromosome numbers and letters
for chr in {1..19} X Y; do
    # Set the input file name based on the current chromosome
    input_file=${prefix}_chr${chr}.like_wig
    
    # Set the output file name based on the current chromosome
    output_file="foo_chr${chr}.wig"
    
    # Check if the input file exists before attempting to process it
    if [ -f "$input_file" ]; then
        # Use the command with the loop variable
        { 
          echo "variableStep chrom=chr${chr} span=10"; 
          tail -n +4 "$input_file" | awk -F "\t" '{print $1, $2}'; 
        } > "$output_file"
        echo "Processed $input_file into $output_file"
    else
        echo "Input file $input_file does not exist, skipping..."
    fi
done

# Concatenate all .wig files into one
cat foo_chr*.wig > allchr.wig

# Delete the individual chromosome .wig files
rm foo_chr*.wig

# Convert the concatenated .wig to .bigWig
wigToBigWig allchr.wig $chromsize "$prefix"_inps.bigWig
rm allchr.wig

echo "Conversion to bigWig completed."


