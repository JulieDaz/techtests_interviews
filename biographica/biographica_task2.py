#-*- coding: utf-8 -*-

import argparse
from Bio import SeqIO
import uniprotAPI
import biographica_task1


def get_hydrolase(assembly2_data) :

    """ Function to retrieve the gene ids that have the "hydrolase activity, acting on ester bonds"
    in the GO term

    Use of the API found at https://www.uniprot.org/help/id_mapping

    Args:
        assembly2_data (list): List of records from the second assembly file

    Returns:
        list: List of gene names that have the "hydrolase activity, acting on ester bonds" GO term
    """

    list_geneid = []
    list_geneid_hydrolase_ester = []

    ## retrieve the gene ids in assembly2
    for rec in assembly2_data : 

        gene_id = rec.description.split(" ")[1]
        list_geneid.append(gene_id)

    ## use of the API script found on uniprot
    job_id = uniprotAPI.submit_id_mapping(
        from_db="Gene_Name", to_db="UniProtKB", ids=list_geneid
    )
    if uniprotAPI.check_id_mapping_results_ready(job_id):
        link = uniprotAPI.get_id_mapping_results_link(job_id)
        results = uniprotAPI.get_id_mapping_results_search(link)

    for result in results["results"]:

        for ref in result["to"]["uniProtKBCrossReferences"] :

            ## hydrolase activity found in the GO terms
            if ref["database"] == "GO" :
                for keys in ref["properties"] :

                    ## find the function in the properties
                    if "F:hydrolase activity, acting on ester bonds" in keys["value"] :
                        ## put result["from"] which is the gene id in a list
                        list_geneid_hydrolase_ester.append(result["from"])

    return list_geneid_hydrolase_ester


def get_assembly1_hydrolase_genes(list_geneid_hydrolase_ester, desambiguate_f, output_f) :

    """ Function to retrieve the sequences that have the hydrolase function in
    the assembly file 1, by using the desambiguate file created in task1. Write
    the fasta file containing the proteins that have the hydrolase function.

    Args:
        list_geneid_hydrolase_ester (list): List of gene names with the esterase function
        desambiguate_f (str): Path to the file in which to look for gene names
        output_f (str): Path to the output file
    """

    list_genes_hydrolase_assembly1 = []

    for rec in SeqIO.parse(desambiguate_f, "fasta")  :

        rec_description = rec.description

        ## headers of sequences desambiguated have ||, by splitting on it and
        ## checking the length with can check is the sequence is one of the
        ## desambiguated ones we want to check for the hydrolase
        if len(rec_description.split("||")) > 1 :

            ## retrieve the gene id from the ref genome
            geneid_2 = rec_description.split("||")[1].split(" ")[1]
            ## if the id is in list_geneid_hydrolase_ester
            ## write the seqrecord in a new output
            if geneid_2 in list_geneid_hydrolase_ester :

                list_genes_hydrolase_assembly1.append(rec)

    SeqIO.write(list_genes_hydrolase_assembly1, output_f, "fasta")



### main
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--fasta1", help="assembly_1.test.prot.fa")
    parser.add_argument("--fasta2", help="assembly_2.test.prot.fa")
    parser.add_argument("--merged", help="assembly1_2_desambiguate.test.prot.fa")
    parser.add_argument("--out", help="proteins_hydrolase_ester.fa")

    args = parser.parse_args()

    file1 = args.fasta1
    file2 = args.fasta2
    desambiguate_file = args.merged
    output_file = args.out

    assembly2_data, all_records, all_records_id = biographica_task1.read_data(file1, file2)
    list_geneid_hydrolase_ester = get_hydrolase(assembly2_data)
    get_assembly1_hydrolase_genes(list_geneid_hydrolase_ester, desambiguate_file, output_file) 
