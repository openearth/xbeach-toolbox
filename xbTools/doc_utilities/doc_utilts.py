"""
Functions to help modify and write documents
"""

# Standard imports
import json
import os

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