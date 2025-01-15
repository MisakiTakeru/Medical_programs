import re
from tkinter.filedialog import askopenfilename
"""
Updated on Mon Sep  2 2024

@author: Joachim Normann Larsen
"""

# the regex to find the lines we want from the report
prog = re.compile(r'[\d-]+\s+[a-zA-Z():\/]*\s\d+\s+[a-zA-Z():\/\s.0-9 -]*')

units = ['Date', 'Clearance', 'Normalized Clearance', 'GFR']

def read_report(filename):
    with open (filename) as file:
        lines = file.readlines()
    return lines


def get_data(lines):
    """
    Parameters
    ----------
    lines : String
        Input is the lines of a rtf file.

    Returns
    -------
    all_data : list
        list of all found data in the rtf file.

    """
    all_data = []
    for line in lines:
        data = prog.match(line)
        if data != None:
            split_line = re.split(r'\s\s\s\s', data.group(0))
            split_line = [l.strip() for l in split_line]


            split_line[1] = re.findall(r'[0-9]+', split_line[1])[0]
            find_1 = re.findall(r'\d+$', split_line[2])
# There are currently found 2 different variations for Normalized Clearance.
# One where we have the explanation and ends with the number and one where we
# have explanations but ends with the number, and an area of expected normal (xxx - yyy)
            if find_1 == []:
                split_line[2] = re.findall(r'\d{2,3}(?=\s\(\d+\s-\s\d+\))', split_line[2])[0]
            else:
                split_line[2] = find_1[0]
            split_line[0] = split_line[0].replace('-','')
#            split_line[2] = re.findall(r'\d+\s\([\d\s-]+\)', split_line[2])[0]
            print(dict(zip(units,split_line)))
            res = dict(zip(units,split_line))
            all_data.append(res)
    return all_data


if __name__ == '__main__':
    filename = askopenfilename()
    text = read_report(filename)

    d = get_data(text)



