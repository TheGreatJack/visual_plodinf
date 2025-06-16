import argparse
import sys


def file_reader(file_path,min_reads_per_allele, min_depth,max_depth):

    while True:
        line = file_path.readline()
        if not line:
            break
        line = line.strip().split("\t")
        alleles=line[2]
        allele_counts=line[5]
        
        alleles,allele_counts = allele_depth_parser(alleles,allele_counts)
        
        real_depth,alleles,allele_counts,allele_proportions = allele_QC(alleles,allele_counts,min_reads_per_allele)

        if (real_depth <= max_depth) and (real_depth >= min_depth) and len(alleles) >= 2:
            if len(alleles) == 2:
                mutation_type = transition_transversion_type(alleles)
                prop_diff = biallelic_diff(allele_proportions)
            elif len(alleles) > 2:
                mutation_type = "-"
                prop_diff = "-"

            chrom= line[0]
            pos = line[1]
            allele_counts = [str(counts) for counts in allele_counts]
            allele_proportions = [str(props) for props in allele_proportions]
            allele_number = len(alleles)

            to_print = [chrom,
                        pos,
                        ",".join(alleles),
                        allele_number,
                        real_depth,
                        ",".join(allele_counts),
                        ",".join(allele_proportions),
                        mutation_type,
                        prop_diff]

            print(*to_print,sep = "\t")



def allele_depth_parser(alleles,allele_counts):
    allele_list = alleles.split(",")
    allele_counts_list = allele_counts.split(",")

    allele_list = [allele for allele in allele_list if allele != "<*>"]
    allele_counts_list= [int(ac) for ac in allele_counts_list if ac != "0"]
    #print(allele_list,allele_counts_list)

    return allele_list,allele_counts_list


def allele_QC(alleles,allele_counts,min_reads_per_allele):
    alleles_cleaned = []
    allele_counts_cleaned = []

    for allele, count in zip(alleles,allele_counts):
        if count >= min_reads_per_allele:
            alleles_cleaned.append(allele)
            allele_counts_cleaned.append(count)

    real_depth = sum(allele_counts_cleaned)
    allele_proportions = [count/real_depth for count in allele_counts_cleaned]

    return real_depth, alleles_cleaned,allele_counts_cleaned,allele_proportions


def transition_transversion_type(alleles):
    type_dict = {"transition":[["A","G"],
                               ["C","T"]],
                "transversion":[["A","C"],
                                ["G","T"],
                                ["A","T"],
                                ["G","C"]]}
    
    type_classf = ""
    for transition_type,type_data in type_dict.items():
        for posibility in type_data:
            if all(item in posibility for item in alleles):
                type_classf = transition_type
                break
        if type_classf != "":
            break
    
    return type_classf

def biallelic_diff(allele_proportions):
    return allele_proportions[0] - allele_proportions[1]




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read the output of the following command from Bcftools: 'bcftools mpileup -q 20 --skip-all-unset 3 -a FORMAT/AD'." + 
                                                " The idea is the get the per allele snp statistics. ")
    parser.add_argument("filename", nargs="?", type=argparse.FileType('r'),
                        help="File to read. If omitted, reads from standard input.")
    parser.add_argument('-m','--min-allele-depth', dest='min_reads_per_allele',
                        default=3,help='Minimum number of reads per allele in a site')
    parser.add_argument('-d','--min-depth', dest='min_depth',
                        default=12,help='Minimum real depth per site')
    parser.add_argument('-M','--max-depth', dest='max_depth',
                        default=80,help='Max real depth per site')
    args = parser.parse_args()

    # Use standard input if no filename provided
    if args.filename is None:
        file = sys.stdin
    else:
        file = args.filename

    #print(args.min_reads_per_allele)
    #print(args.min_depth)
#    line_procceser_eff(file)

    min_reads_per_allele = int(args.min_reads_per_allele)
    min_depth = int(args.min_depth)
    max_depth = int(args.max_depth)
    file_reader(file,min_reads_per_allele, min_depth,max_depth)

    # Close the file if it's not standard input
    if file != sys.stdin:
        file.close()

