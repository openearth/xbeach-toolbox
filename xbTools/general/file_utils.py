"""
Functions that help you output or work with files
"""

# Standard imports
import json
import os

def write_2d_arr_2_file(arr, file_path, decimal_digits=3):
    """
    Writes a 2d array into a file with each row of the 2d array written as a row in the file.

    Inputs:
        file_path      : Name of the file that the data should be written to.
        arr            : 2d array that should be written to the file.
        decimal_digits : Number of decimal digits to format the elements in the array (default is 2).
    """

    # Create a format string based on the number of decimal digits
    format_str = f"{{:.{decimal_digits}f}}"

    # Open the file in write mode
    with open(file_path, 'w') as file:
        # Write each row of the array to the file
        for row in arr:
            formatted_row = ' '.join(format_str.format(value) for value in row)
            file.write(formatted_row + '\n')

def format_header_line(text):
    """
    Used to format the header line the in params.txt file
    
    inputs:
        text: string
            The text should be put in the header
    """
    return f"%%% {text.ljust(71)} %%%"

def format_subsection_header(text, length = 78):
    """
    Format the subsection headers in the text file
    """
    string = f"%%% {text}"

    num_perc_signs = length - len(string)

    perc_signs = "%" * num_perc_signs
    return f"{string} {perc_signs}"

def get_json(folder_path, file_name):
    """
    Reads a json and returns it as a dict
    NOTE: Can be used to look at the json in script so you can see the inputs
    """

    ## Open the json file
    json_file         = open(os.path.join(folder_path, file_name),'r')

    # Read the json file and store it as a python dict
    json_dict = json.loads(json_file.read())

    return json_dict