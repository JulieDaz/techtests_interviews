#-*- coding: utf-8 -*-
# Python 3.7.1

import argparse

def read_data_and_build_strings(input_f) :

    """Read the file and create the 3 strings

    Args:
        input_f {str}: file name

    Returns:
        str: ref_seq which is the original sequence
        str: the alignment characters composed of "|" and " "
        str: the sequence from the CIGAR operations
    """

    # initiate the variables
    # we will build each sequence as we read the file
    ref_seq = ""
    aln_char = ""
    operation_seq = ""
    # the index is to keep up of where we are in the original sequence, and 
    # will be incremented for 3 out of 4 operations
    index_original_seq = 0

    with open(input_f, "r") as cigar_file :

        for line in cigar_file :

            # check if the line starts with any nt
            # if so it means it's the sequence that we'll use
            if line.startswith(("A", "T", "C", "G")) :
                original_seq = line.strip()

            else :
                # split the line on the comma
                # 2nd element is always the number of nt impacted
                myline = line.strip().split(",")
                nb_nt = int(myline[1])

                if line.startswith("M") :

                    # ref_seq and operation_seq have the same action in this case
                    # we copy the oringal_seq by using the index and the nb_nt
                    ref_seq += original_seq[index_original_seq:index_original_seq+nb_nt]
                    aln_char += nb_nt*"|"
                    operation_seq += original_seq[index_original_seq:index_original_seq+nb_nt]

                    # increment the index
                    index_original_seq += nb_nt

                elif line.startswith("I") :

                    # retrieving the inserted nucleotides
                    new_nt = myline[2]

                    # adding a gap for theref
                    ref_seq += nb_nt*"-"
                    # white space as mismatch
                    aln_char += nb_nt*" "
                    # adding the new nt in the operation seq
                    operation_seq += new_nt

                    # do not increment the index in case of an insertion 
                    # as the original sequence won't change

                elif line.startswith("X") :

                    new_nt = myline[2]

                    # ref seq unchanged in this condition, so we use the original seq
                    ref_seq += original_seq[index_original_seq:index_original_seq+nb_nt]
                    aln_char += nb_nt*" "
                    # adding the new nt instead of the original seq
                    operation_seq += new_nt

                    index_original_seq += nb_nt

                elif line.startswith("D") :

                    # ref seq unchanged in this condition, so we use the original seq
                    ref_seq += original_seq[index_original_seq:index_original_seq+nb_nt]
                    aln_char += nb_nt*" "
                    # gap to add 
                    operation_seq += nb_nt*"-"

                    index_original_seq += nb_nt

    return ref_seq, aln_char, operation_seq


def calculate_ASCII_sum(seq) : 

    """Calculate the ASCII of all the characters of 1 string

    Args:
        seq {str}: characters sequence.

    Returns:
        int: the ASCII sum for the input string
    """

    # initiate the variable at 0
    sumASCII = 0

    for char in seq : 
        # using ord() function that returns the number representing the unicode 
        # code of a specified character
        sumASCII += ord(char)

    return sumASCII


def write_output(output_f, ref_seq, aln_char, operation_seq, ASCII_score) :

    """Write the 3 strings in one output file

    Args:
        output_f {str}: output file name
        ref_seq {str}: reference sequence string
        aln_char {str}: alignment characters string
        operation_seq {str}: operation sequence string
        ASCII_score {int}: ASCII value

    """

    with open(output_f, "w") as output : 

        output.write(f"> The product of the ASCII score between the 3 following strings is {ASCII_score}.\n")
        output.write(f"{ref_seq}\n")
        output.write(f"{aln_char}\n")
        output.write(f"{operation_seq}\n")


#### MAIN

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", help="Input file with the reference seq and the CIGAR operations.")
    parser.add_argument("--output", help="Path of the output file.")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    ref_sequence, aln_char_sequence, operation_sequence = read_data_and_build_strings(input_file)

    # reuse the same function to calculate the ASCII sum for each of the 3 strings
    ref_sequence_ASCII = calculate_ASCII_sum(ref_sequence)
    aln_char_sequence_ASCII = calculate_ASCII_sum(aln_char_sequence)
    operation_sequence_ASCII = calculate_ASCII_sum(operation_sequence)

    # calculate the ASCII score
    total_ASCII_score = ref_sequence_ASCII * aln_char_sequence_ASCII * operation_sequence_ASCII
    
    write_output(output_file, ref_sequence, aln_char_sequence, operation_sequence, total_ASCII_score)


