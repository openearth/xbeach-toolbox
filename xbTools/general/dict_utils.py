# dict_utils.py

def check_keys_exist(dict1, dict2):
    """
    Check if all keys from dict1 are present in dict2.

    Parameters:
    - dict1 (dict): Dictionary to check keys.
    - dict2 (dict): Dictionary to check against.

    Returns:
    - bool: True if all keys from dict1 are in dict2, False otherwise.
    """
    return all(key in dict2 for key in dict1)

def check_strings_in_dict(list_of_strings, dict_strings):
    """
    Check if all strings in list_of_strings are present as values in dict_strings.

    Parameters:
    - list_of_strings (list): List of strings to check.
    - dict_strings (dict): Dictionary with string values.

    Returns:
    - bool: True if all strings in list_of_strings are in dict_strings, False otherwise.
    """
    dict_values = set(dict_strings.values())
    return all(string in dict_values for string in list_of_strings)

# Other utility functions can go here if needed
# def other_utility_function():
#     pass

# Test for the code
if __name__ == "__main__":
    dict1 = {'a': 1, 'b': 2, 'c': 3}
    dict2 = {'a': 10, 'b': 20, 'c': 30, 'd': 40}

    print(check_keys_exist(dict1, dict2))  # Output: True

    dict3 = {'a': 10, 'c': 30}
    print(check_keys_exist(dict1, dict3))  # Output: False
