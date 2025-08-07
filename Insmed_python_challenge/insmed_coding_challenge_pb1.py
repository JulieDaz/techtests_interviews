#-*- coding: utf-8 -*-
# Python 3.7.1

import argparse

def read_data(input_f) :

    """Function to read the data (cell types and cell sizes) and add it to a dictionary.

    Args:
        input_f {str}: input file name

    Returns:
        dict: dictionary where key=cell type ; and value is a list of cell sizes.
    """

    # initiate the variables
    dict_cell_type_size = {}

    with open(input_f, "r") as cell_type_f :

        for line in cell_type_f :

            myline = line.strip()
            # split the line on the colon and put the data in 2 variables
            cell_type, cell_size = myline.split(":")

            # dictionary with keys = cell type, and value = list of cell sizes
            dict_cell_type_size.setdefault(cell_type, []).append(int(cell_size))

    return dict_cell_type_size


def calculate_diff_growth(dict_cell_type_size) :

    """Calculate the differential growth between the first and last groups, 
    after sorted the dictionary keys in the ascending order, first by alphabetical
    order, then numerical order.

    Args: 
        dict_cell_type_size {dict}: contains the cell_type as a key, and a list 
        of cell_size as a value.

    Returns:
        str: the first cell type name.
        str: the last cell type name.
        str: the differential growth value.
    """

    dict_cell_type = {}

    # create a dictionary with the key = cell type, and value = list of the 
    # alpha and numerical part of the name
    for cell_type, cell_size in dict_cell_type_size.items() :
        alpha_val = cell_type[:2]
        num_val = cell_type[2:]

        dict_cell_type[cell_type] = [alpha_val, int(num_val)]

    # sort the values of the dict on the first item, then the second item
    # returns list of tuples
    # [
    #     key, [alpha_val, int(num_val)],
    #     key, [alpha_val, int(num_val)]
    #     key, [alpha_val, int(num_val)]
    #     key, [alpha_val, int(num_val)]] -> x[1] = list alpha val/num val
    dict_cell_type_sorted = sorted(dict_cell_type.items(), key = lambda x:(x[1][0],x[1][1]))
    
    # retrieve first and last groups by accessing the 1st and last elements of 
    # the dict
    first_group = dict_cell_type_sorted[0][0]
    last_group = dict_cell_type_sorted[-1][0]

    # calculate the differential growth between the first and last groups.
    differential_growth = sum(dict_cell_type_size[first_group]) - sum(dict_cell_type_size[last_group])

    return first_group, last_group, differential_growth


###### MAIN 

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", help="Input file where each line consists of a cell type and cell size.")
    args = parser.parse_args()

    input_file = args.input
    dico_cell_type = read_data(input_file)
    first_gp, last_gp, diff_growth = calculate_diff_growth(dico_cell_type)

    print(f"The differential growth between cell types {first_gp} and {last_gp} is {diff_growth}.")
