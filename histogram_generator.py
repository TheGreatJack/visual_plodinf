"""
Processes a TSV file to calculate frequency distributions from a specific column
and generates a histogram of these frequencies.
"""

import argparse
import csv
import os
import sys
import matplotlib.pyplot as plt

def create_histogram(data, bins, output_path, title):
    """
    Generates and saves a histogram from a list of data.

    Args:
        data (list): A list of floats or integers to be plotted.
        bins (int): The number of bins for the histogram.
        output_path (str): The full path to save the output PNG file.
        title (str): The title for the histogram chart.
    """
    if not data:
        print("Warning: No data to plot. Histogram will not be generated.")
        return

    print(f"\nGenerating histogram with {len(data)} data points and {bins} bins...")

    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, color='skyblue', edgecolor='black')

    plt.title(f'Frequency Distribution: {title}')
    plt.xlabel('Calculated Frequency')
    plt.ylabel('Count')
    plt.grid(axis='y', alpha=0.75)

    try:
        # Ensure the output directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            
        plt.savefig(output_path)
        print(f"Success! Histogram saved to: {output_path}")
    except Exception as e:
        print(f"Error: Could not save the histogram file. {e}", file=sys.stderr)
    finally:
        plt.close()

def main():
    """
    Main function to parse arguments and orchestrate the file processing and plotting.
    """
    parser = argparse.ArgumentParser(
        description="Reads a TSV file, calculates frequencies from the 6th column, "
                    "and outputs a histogram.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('tsv_file',
                        type=str,
                        help="Path to the input TSV file.")
    parser.add_argument('bins',
                        type=int,
                        help="Number of bins for the output histogram.")
    parser.add_argument('output_dir',
                        type=str,
                        help="Path to the directory where the output PNG histogram will be saved.")

    args = parser.parse_args()

    # --- 1. Validate inputs ---
    if not os.path.isfile(args.tsv_file):
        print(f"Error: Input file not found at '{args.tsv_file}'", file=sys.stderr)
        sys.exit(1)
        
    if args.bins <= 0:
        print(f"Error: Number of bins must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    # --- 2. Process the TSV file ---
    frequencies = []
    print(f"Reading data from: {args.tsv_file}")

    try:
        with open(args.tsv_file, mode='r', newline='', encoding='utf-8') as infile:
            reader = csv.reader(infile, delimiter='\t')
            
            for i, row in enumerate(reader):
                row_num = i + 1
                # The 6th column is at index 5
                try:
                    # Get the string from the 6th column
                    numbers_str = row[5]
                    
                    if not numbers_str.strip():
                        print(f"Warning: Skipping row {row_num} due to empty 6th column.")
                        continue
                        
                    # Split the string by commas and convert to integers
                    numbers = [int(n) for n in numbers_str.split(',')]
                    
                    # Calculate the sum
                    total = sum(numbers)
                    
                    if total == 0:
                        print(f"Warning: Skipping row {row_num} because the sum of its integers is 0.")
                        continue
                        
                    # Calculate frequencies and add to the main list
                    for num in numbers:
                        frequency = num / total
                        frequencies.append(frequency)
                        
                except IndexError:
                    print(f"Warning: Skipping row {row_num} as it has fewer than 6 columns.", file=sys.stderr)
                except ValueError:
                    print(f"Warning: Skipping row {row_num} due to non-integer value in the 6th column.", file=sys.stderr)
                except Exception as e:
                    print(f"An unexpected error occurred at row {row_num}: {e}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.tsv_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 3. Generate the output histogram ---

    # Construct the output file name from the input file name
    base_name = os.path.basename(args.tsv_file)
    file_name_no_ext = os.path.splitext(base_name)[0]
    output_filename = f"{file_name_no_ext}.png"
    output_path = os.path.join(args.output_dir, output_filename)

    create_histogram(frequencies, args.bins, output_path, title=file_name_no_ext)


if __name__ == '__main__':
    main()