#

""" Package for string/regular-expressions manipulation"""

##@package regexp

import sys;
import string;
import re;

# --
def replace(filename, option, value, section=None):
    """ Function to edit a "key value" configuration file.
    
    Input:
     - filename  str : (text) filename to be edited
     - option    str : Option/key string to search for and change corresponding 'value'
     - value     str : Value to substitute for given 'option'
     - section   str : If 'filename' contains multiple sections, define the one to be looked
    
    Ouput:
     - ? 
    
    ---
    """
    
    section_flag = False
    for line in fileinput.input(filename, inplace=1):
        if (section in line) or (section == None):
            section_flag = True
    
        if (option in line) and (section_flag == True):
            line = option + " " + str(value) + "\n"
            sys.stdout.write(line)
        
    return

# --
def str2lst(word, sep=",", valid_digits="a-zA-Z0-9_\ .", valid_symbols="\-\+\^"):
    """
    Function to read a string with special chars and return without them.

    This function cleans the given 'word' from non-valid ('digits' and 'symbols') 
    characters and splits it if 'sep' character is found. A list is returned with 
    resultant list of strings.

    Character used as "valid ones" are passed through 'valid_digits' and 'valid_symbols';
    These two parameters are concateneted to form one "valid-characters" string.

    * The function workflow first split the string and then removes the non-given chars

    Input:
     - word           str : The word/string to clean and split
     - sep            str : string separator
     - valid_digits   str : valid characters
     - valid_symbols  str : valid characters

    Output:
    - list with given 'word' parts and cleaned

    ---
    """

    valid_chars = digits + symbols

    lista = [ re.sub("[^"+valid_chars+"]", "", i)  for i in string.split(word, sep=sep) ]

    return lista

# ---

# ==========================================================
def line_filter(lines, comments="#", inverse=False):
    """
    Function to return uncomment line entries from a list.

    'lines' is a list of strings, those that do not begin with 'comments' symbol/char 
    are returned (also as a list of strings). If 'inverse' is set to True, just commented 
    lines are returned in the same way.

    Input:
     - lines    [str] : List with string lines
     - comments   str : Symbol used for commented lines
     - inverse   bool : Return commented lines(?)

    Output:
     - list with selected string lines

    ---
    """

    data_lines = [];
    comment_lines = [];

    for line in lines:
    
        line = line.rstrip("\n");
        line = re.sub("^\s*", "", line);
        line = re.sub("\s+", " ", line);

        if re.search("^$", line): continue;

        if re.search("^"+comments, line):
            comment_lines.append( re.sub("^"+comments, "", line) );
            continue;

        data_lines.append(line);

    if inverse:
        return comment_lines;

    return data_lines;
