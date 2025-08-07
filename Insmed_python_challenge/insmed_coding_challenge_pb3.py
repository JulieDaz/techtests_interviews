#-*- coding: utf-8 -*-
# Python 3.7.1

import argparse
import re

def read_file(input_f) :

    """Read the input file and divide it into 2 lists

    Args:
        input_f {str}: file name

    Returns:
        str: str containing the first set of genes coordinates.
        str: str containing the second set of genes coordinates.
    """

    with open(input_f, "r") as gene_list_f :nxjjmr9h7n2z
    
        
        myfile = gene_list_f.readlines()

        # first set contained in the 1st line of the file
        gene_set_1 = myfile[0]
        # second set contained the 2nd line of the file
        gene_set_2 = myfile[1]

    return gene_set_1, gene_set_2


def create_list_of_coord(coord2parse) :

    """Function to parse the string of coordinates to return a list.

    Args:
        coord2parse {str}: long string of coordinates separated by round brackets
    
    Returns:
        list: list with the genes coordinates from the input string.
    """

    # initiate the list that will contain the coordinates
    list_coord_genes = []

    # use of a regex to find the coordinates couples
    coord_regex = re.findall("([0-9]+, [0-9]+)", coord2parse)

    for couple in coord_regex : 

        ## each couple is composed of 2 numbers separated by a ", "
        ## we can split on the "," which returns a list of 2 elements that
        ## we convert to int and add to a list
        couple_split =  couple.split(",")
        list_coord_genes.append([int(couple_split[0]), int(couple_split[1])])

    return list_coord_genes


def compare_ovl(list_coord_genes_1, list_coord_genes_2) :

    """Calculate the overlapping regions between 2 lists of genes, and then 
    calculate the sum of all the overlaps.

    Args:
        list_coord_genes_1 {list}: 1st list of genes coordinates.
        list_coord_genes_2 {list}: 2nd list of genes coordinates.

    Returns:
        list: list of the overlapping regions
        int: sum of all the overlaps
    """

    list_overlap = []
    total_distance_ovl = 0 

    # loop on the 2 lists of genes
    for start1, end1 in list_coord_genes_1 :

        for start2, end2 in list_coord_genes_2 : 

            # check for overlaps 
            if start1 <= end2 and end1 >= start2 : 

                # we need the biggest start to calculate the overlap
                if start1 > start2 :
                    biggest_start = start1
                
                elif start1 < start2 :
                    biggest_start = start2

                else :
                    # start1 == start2, we select one
                    biggest_start = start1

                # and we need the smallest end
                if end1 > end2 :
                    smallest_end = end2

                elif end1 < end2 :
                    smallest_end = end1

                else :
                    # end1 == end2, we select one
                    smallest_end = end1

                list_overlap.append([biggest_start, smallest_end])

                # calculate the distance between the overlaps
                total_distance_ovl += smallest_end - biggest_start

    return list_overlap, total_distance_ovl


def write_output(output_f, list_overlap, total_distance_ovl) :

    """Write the final output in a file.
    Converts the list of overlaps to a string that matches the format of the 
    input.

    Args:
        output_f {str}: path of the output file
        list_overlap {list}: list of the overlapping regions
        total_distance_ovl {int}: value of the sum of the overlaps.

    """

    seq2write = ""
    with open(output_f, "w") as output_ovl :

        output_ovl.write(f"The overlapping regions are as follows:\n")
        for ovl_coord in list_overlap :

            # write the final sequence of overlaps to write in the output file.
            seq2write += f"({ovl_coord[0]}, {ovl_coord[1]}),"

        output_ovl.write(seq2write.strip(","))
        output_ovl.write(f"\nThe total overlap is {total_distance_ovl}.")


#### MAIN

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", help="Input file with the reference seq and the CIGAR operations.")
    parser.add_argument("--output", help="Path of the output file.")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    set1_genes, set2_genes = read_file(input_file)

    coord_genes_1 = create_list_of_coord(set1_genes)
    coord_genes_2 = create_list_of_coord(set2_genes)

    list_ovl, distance_ovl_total = compare_ovl(coord_genes_1, coord_genes_2)

    write_output(output_file, list_ovl, distance_ovl_total)



