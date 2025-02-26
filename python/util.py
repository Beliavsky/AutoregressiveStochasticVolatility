import os
import shutil
import glob
import shutil

def copy_files(src_dir,dest_dir):
    """
    Copy all files from src_dir to dest_dir.
    """
    for root, dirs, files in os.walk(src_dir):
        for f in files:
            shutil.copy(os.path.join(root,f), dest_dir)

def files_in_dir(path,extension):
    """
    Return a list of files with the given extension in the specified directory.
    """
    os.chdir(path)
    result = [i for i in glob.glob('*.{}'.format(extension))]
    return result

def lines_in_f1_not_f2(f1,f2):
    """
    Return lines that are in f2 but not in f1.
    """
    s1 = open(f1,"r").readlines()
    s2 = open(f2,"r").readlines()
    lines = []
    for word in s2:
        if (word not in s1):
            lines.append(word)
    return lines
    
def truncated_string(original_string, sentinel):
    """ return original string up to the occurrence of sentinel """
    index = original_string.find(sentinel)
    if index != -1:
        truncated_string = original_string[:index]
        return truncated_string
    else:
        return original_string

def save_numbered_file(file_path):
    """
    Copy a file to an unused file name with a numbered suffix.
    Args:
        file_path (str): The path to the file to be copied.
    Returns:
        str: The path of the newly copied file.
    """
    base_dir = os.path.dirname(file_path)
    base_name, ext = os.path.splitext(os.path.basename(file_path))
    new_file_path = file_path
    counter = 1
    while os.path.exists(new_file_path):
        new_file_name = f"{base_name}{counter}{ext}"
        new_file_path = os.path.join(base_dir, new_file_name)
        counter += 1
    shutil.copyfile(file_path, new_file_path)
    return new_file_path

def modify_lines(input_string, start_strings, append_char, ignore_case=False,
    ignore_spaces=False):
    """
    Modifies lines of a string that start with any string from a list.

    Args:
    input_string (str): The string to modify.
    start_strings (list of str): The list of strings to check if a line starts with.
    append_char (str): The character to append to matching lines. Defaults to '\n'.
    ignore_case (bool, optional): If True, ignores case when matching start strings. Defaults to False.
    ignore_spaces (bool, optional): If True, ignores leading spaces when matching start strings. Defaults to False.

    Returns:
    str: The modified string.
    """
    # Break the string into lines
    lines = input_string.split('\n')
    # Iterate over lines and append character if the line starts with any string from start_strings
    for i, line in enumerate(lines):
        # Remove leading spaces if ignore_spaces is True
        line_to_check = line.lstrip() if ignore_spaces else line
        # Convert to lower case if ignore_case is True
        line_to_check = line_to_check.lower() if ignore_case else line_to_check
        # Check if line starts with any string from start_strings, considering ignore_case option
        if any(line_to_check.startswith(start_str.lower() if ignore_case else start_str) for start_str in start_strings):
            lines[i] += append_char
    # Join modified lines back into a string
    output_string = '\n'.join(lines)
    return output_string

def replace_leading_spaces(file_path, nspaces=1, rep="", output_file=None):
    """
    Replace leading spaces in each line of a file.

    Parameters:
    file_path (str): The path to the file to be processed.
    nspaces (int): The number of leading spaces to remove from each line. Defaults to 1.
    rep (str): The string to replace the removed spaces with. Defaults to an empty string.
    output_file (str): The path to the output file. If None, the original file is overwritten.

    Notes:
    - Lines with between 1 and `nspaces-1` leading spaces will have those spaces replaced by `rep`.
    - Lines with `nspaces` or more leading spaces will have exactly `nspaces` spaces replaced by `rep`.
    - Lines with no leading spaces are left unchanged.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    if output_file is None:
        output_file = file_path

    with open(output_file, 'w') as file:
        for line in lines:
            leading_spaces = len(line) - len(line.lstrip())
            if 1 <= leading_spaces < nspaces:
                # Replace between 1 and nspaces-1 leading spaces
                line = rep + line[leading_spaces:]
            elif leading_spaces >= nspaces:
                # Replace exactly nspaces leading spaces
                line = rep + line[nspaces:]
            file.write(line)

def powers(n, primes=[2, 3, 5]):
    """
    Partially factors an integer n into powers of the given primes.
    
    Args:
        n (int): The integer to factorize.
        primes (list): A list of primes to use for factorization (default: [2, 3, 5]).
    
    Returns:
        tuple: A dictionary with the prime factors and their exponents, and the remaining factor.
    """
    factors = {}
    remaining = n
    for p in primes:
        count = 0
        while remaining % p == 0:
            remaining //= p
            count += 1
        if count > 0:
            factors[p] = count    
    return factors, remaining

def print_vec(vec, fmt, label=None, end=None):
    """ print elements of a list with format fmt """
    if label is not None:
        print(label, end="")
    print("".join(fmt%x for x in vec))
    if end is not None:
        print(end=end)
