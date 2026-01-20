import re
import json

def extract_parameters_by_section(file_path):
    section_params = {}
    current_section = None

    with open(file_path, 'r', encoding='utf-8') as file:
        # Skip the first 22 lines (header)
        for _ in range(22):
            next(file)

        for line in file:
            stripped = line.strip()

            # Skip comment lines unless it's a section marker
            if stripped.startswith('!'):
                if not stripped.lower().startswith('! [section]'):
                    continue

            # Handle section headers
            if stripped.lower().startswith('! [section]'):
                current_section = stripped[11:].strip()
                section_params[current_section] = []
                continue

            # Look for lines with "::" indicating parameter declarations
            if "::" in line:
                try:
                    # Extract portion after '::'
                    after_colon = line.split("::", 1)[1]

                    # Determine where the parameter name ends
                    if "=" in after_colon:
                        param_part = after_colon.split("=")[0]
                    elif "!" in after_colon:
                        continue  # Parameter with an unitialized default value
                    else:
                        continue  # Malformed or not a declaration

                    # Clean the parameter name: remove brackets and whitespace
                    param_cleaned = re.sub(r'\(.*?\)', '', param_part).strip()

                    # if current_section and param_cleaned:
                    section_params[current_section].append(param_cleaned)

                except Exception:
                    continue  # Skip any badly formatted lines

    return section_params


# Example usage

file_path = r"c:\checkouts\xbeach_trunk\src\xbeachlibrary\params.def"  # replace with your actual file path
parameters = extract_parameters_by_section(file_path)

# this one doesn't format well
try:
    parameters['Output variables'].remove('pointtypes                                   !  [-]  Point types (0')
except:
    pass



# remove these from the standard key
underscore_list = [
        "nglobalvar",
        "nmeanvar",
        "npointvar",
        "npoints"
    ]
for par in underscore_list:
    parameters['Output variables'].remove(par)
# add the underscore list to the hidden output variables
parameters['_Output variables'] = underscore_list

# remove these from the list as they are not read from params.txt
parameters.pop('Constants, not read in params.txt')
parameters.pop('Variables, not read in params.txt')

# write to file
with open(r'c:\checkouts\python\xbeachToolbox\xbTools\par.json', 'w', encoding='utf-8') as f:
    json.dump(parameters, f, indent=4)  # indent for pretty formatting