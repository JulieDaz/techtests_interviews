#-*- coding: utf-8 -*-

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict


def parse_blast_matrix(matrix) :

    """Function to parse the matrix output from the blastp command.

    Returns:
        dict: dictionary containing the ids from assembly2 as a key, and point
        to a list of list that have the ids from assembly1 and the pidentity
        list: a list of the proteins that are duplicated between the two files
    """

    list_duplicates = []
    dict_seq2desambiguate = defaultdict(lambda : [])

    with open(matrix, "r") as blast_matrix: 

        for line in blast_matrix :

            myline = line.strip().split("\t")

            query = myline[0]
            subject = myline[1]
            pident = float(myline[2])

            aln_len = int(myline[3])
            qlen = int(myline[6])
            slen = int(myline[9])

            # calculate the coverages based on the length of the aln
            # and the query/subject length
            covQ = aln_len/qlen
            covS = aln_len/slen

            ## filter on the %identity and the 2 coverages
            if pident >= 95 and covQ >= 0.9 and covS >= 0.9 :
                dict_seq2desambiguate[query].append([subject, pident])
                list_duplicates.append(query)
                list_duplicates.append(subject)

    return dict_seq2desambiguate, list_duplicates


def read_data(assembly1, assembly2) :

    """Function to read the fasta files and save the data into several lists

    Args:
        assembly1 (str): path of the first assembly file 
        assembly2 (str): path of the second assembly file

    Returns:
        list: List of 2 lists
        - list containing all the records from both files
        - list containing all the record ids
    """

    ## creating empty lists
    all_records = [] # all the records info
    all_records_id = [] # only the ids of the records
    assembly2_data = []

    ## read the two files as SeqIO objects
    seq_records_1 = SeqIO.parse(assembly1, "fasta") 
    seq_records_2 = SeqIO.parse(assembly2, "fasta") 

    ## extract the data for assembly1
    for seq_record in seq_records_1 :

        ## removes the "." found at the end of each seq
        ## return an error otherwise
        seq_record.seq = seq_record.seq.strip(".")

        ## put the records into 1 list
        ## and the record ids into another one
        all_records.append(seq_record)
        all_records_id.append(seq_record.id)

    ## same as above
    for seq_record in seq_records_2 :

        ## need the assembly 2 data for task 2
        assembly2_data.append(seq_record)
        all_records.append(seq_record)
        all_records_id.append(seq_record.id)

    return assembly2_data, all_records, all_records_id


def desambiguate_data(dict_seq2desambiguate, list_duplicates, all_records, all_records_id) :

    """Function to disambiguate the files. Uses the dict produced in the first 
    function.
    If one of the key has a list > 1, it means several possible matches, so we
    select the one with the highest %identity

    Returns:
        list: list of all the records, including the unique sequences and the 
        desambiguated ones.
    """

    ## create the empty list
    final_record = []

    for subject, query_pident in dict_seq2desambiguate.items() :

        ## if true, several matches possible
        if len(query_pident) > 1 :

            highest_pident = max(query_pident, key=lambda x: x[1])
            gene_selected = highest_pident[0]

        ## only one match
        else :
            gene_selected = query_pident[0][0]

        ## loop on all the records to retrieve the records matching the ids from
        ## above
        for rec in all_records :

            if subject == rec.id :
                rec_subject = rec
            
            if gene_selected == rec.id :
                rec_query = rec
            
        ## write new sequence name by concatenating names of seq 1 and 2
        ## I choose to use the description, as it contains all the information
        ## so the initial header are kept, though "<unknown description>"
        ## is added for the second half
        ## "||" to separate the 2 names        
        newname = f"{rec_query.description}||{rec_subject.description}"
        ## create new record
        newrecord = SeqRecord(rec_subject.seq, id = newname)

        final_record.append(newrecord)

    ## remove the ids of duplicates sequences from the list that contains all the 
    ## ids from the 2 files
    for seq_dupl in list_duplicates :
        if seq_dupl in all_records_id :
            ## all_records_id only now has the ids of the unique sequences 
            all_records_id.remove(seq_dupl)

    ## retrieve the whole record based on the unique ids of all_records_id
    for ident in all_records_id :
        for rec in all_records :
            if ident == rec.id :
                # print("coucou")
                final_record.append(rec)

    return final_record


def write_fasta_output(final_record, output_f) :

    """Write the sequences into an output in fasta format

    Args:
        final_record (list): List of all the records without duplicates
        output_f (str): Path to the output file
    """

    SeqIO.write(final_record, output_f, "fasta")


##### MAIN

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--fasta1", help="File for 1st assembly, eg: assembly_1.test.prot.fa")
    parser.add_argument("--fasta2", help="File for 2nd assembly, eg: assembly_2.test.prot.fa")
    parser.add_argument("--matrix", help="Matrix from blastp, eg: blastp_output.tsv")
    parser.add_argument("--out", help="Output path, eg: assembly1_2_desambiguate.prot.fa")
    args = parser.parse_args()

    file1 = args.fasta1
    file2 = args.fasta2
    blast_matrix = args.matrix
    output_file = args.out

    ## run functions
    dict_seq2desambiguate, list_duplicates = parse_blast_matrix(blast_matrix)
    assembly2_data, all_records, all_records_id = read_data(file1, file2)
    final_record = desambiguate_data(dict_seq2desambiguate, list_duplicates, all_records, all_records_id)
    write_fasta_output(final_record, output_file)



