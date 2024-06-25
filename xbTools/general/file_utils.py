"""
Functions that help you output or work with files
"""

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


