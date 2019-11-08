"""
Script used to pull data from TASOC ASCII .dat files.
- Loads data from file (name specified with text prompt)
- Skips over commented text
- For each line, checks if any quality flags arises. If not, saves column 1 (time) and column 4 (corrected flux) in
  new file.

The TASOC lightcurve file is expected structured with the following columns (ordered from column number low to high):
Time (BJD-2457000, days)
RAW FLUX (erg/s) (e-/s)
RAW FLUX ERROR (erg/s)
CORRECTED FLUX (PPM)
CORRECTED FLUX ERROR (PPM)
PIXEL QUALITY FLAGS
QUALITY FLAGS (OPTIONAL, SEE PROMPT COMMENT in corrector function)
"""

import os


def correcter(filename):
    f = open(os.path.join('Datafiles', filename), 'r')  # Open in readmode
    open(os.path.join('Datafiles', 'corrected' + filename), 'w').close()  # Clears the write file
    b = 0
    print('Input column number for pixel quality flag (as indicated by datafile).')
    inpt1 = input()
    pqual_col = int(inpt1)
    print('Input general quality flag column number (if not applicable, input 0).')
    inpt2 = input()
    gqual_col = int(inpt2)
    if gqual_col == 0:
        gqual_col = pqual_col

    for line in f:
        columns = line.split()
        if columns[0] == '#-------------------------------------------------':
            if b == 0:
                b += 1
            elif b == 1:
                print('done')
        elif columns[0] == '#':
            pass
        elif columns[pqual_col-1] and columns[gqual_col-1] == '0':
            content = columns[0] + '   ' + columns[3]
            with open(os.path.join('Datafiles', 'converted_' + filename), 'a') as g:
                g.write(content + '\n')
    f.close()


print('What is the filename? (F.ex. file.dat)')
inpt = input()
correcter(inpt)
