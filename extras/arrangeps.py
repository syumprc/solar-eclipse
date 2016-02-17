#!/usr/bin/env python
"""
Name     : arrangeps
Date     : Tue Apr 18 15:28:31 CDT 2000
Author   : kpolk
Usage    : [-o outputfile] [-p #plots/page] [-L Landscape] [-P Portrait] (files)"
Location : /net/spuhler/opt1/local/share/python

Group a list of Landscape Postscript files onto a series of single pages.
"""
sccs_ident = "@(#)arrangeps.py	1.2 04/18/00 SFBR"

import sys, os, string, getopt

PORTRAIT = 0
LANDSCAPE = 1
PAGESIZE = [564, 740]

getpspage = 'grep -v "%%" %s'
getps = 'grep -v "%%" %s | grep -v "showpage"'
plot_template = """gsave
%%%%arrangeps - kpolk
%.2f %.2f translate
%.2f %.2f scale
%.2f rotate
%s
grestore
"""

# Portrait mode arrangement control dictionary
# Format: plots/page : [#cols, #rows, xoffset, yoffset, reduction]
pinfo={ 1 : [1,1,  0, 0, .85],
        2 : [1,2, 40,45, .68],
        3 : [2,2,  0, 0, .47], 4 : [2,2,  0, 0, .47],
        5 : [2,3,  0,35, .44], 6 : [2,3,  0,35, .44],
        7 : [2,4, 50,25, .31], 8 : [2,4, 50,25, .31],
        9 : [3,3, 10, 0, .30],
        10: [3,4, 10,25, .30], 11: [3,4, 10,25, .30], 12: [3,4, 10,25, .30],
        13: [4,4, 15, 5, .23], 14: [4,4, 15, 5, .23],
        15: [4,4, 15, 5, .23], 16: [4,4, 15, 5, .23],
        17: [4,5, 15,20, .23], 18: [4,5, 15,20, .23],
        19: [4,5, 15,20, .23], 20: [4,5, 15,20, .23],
        21: [4,6, 30,25, .20], 22: [4,6, 30,25, .20],
        23: [4,6, 30,25, .20], 24: [4,6, 30,25, .20]}

# Landscape mode arrangement control dictionary
# Format: plots/page : [#cols, #rows, xoffset, yoffset, reduction]
linfo={ 1 : [1,1,  0,10, 1.0],
        2 : [1,2,120, 0, .58],
        3 : [2,2, 10,15, .52], 4 : [2,2, 10,15, .52],
        5 : [2,3, 40, 0, .40], 6 : [2,3, 40, 0, .40],
        7 : [2,4, 60, 8, .31], 8 : [2,4, 60, 8, .31],
        9 : [3,3, 10,20, .36],
        10: [3,4, 25,10, .30], 11: [3,4, 25,10, .30], 12: [3,4, 10,25, .30],
        13: [4,4, 15,25, .26], 14: [4,4, 15,25, .26],
        15: [4,4, 15,25, .26], 16: [4,4, 15,25, .26],
        17: [4,5, 20,15, .24], 18: [4,5, 20,15, .24],
        19: [4,5, 20,15, .24], 20: [4,5, 20,15, .24],
        21: [4,6, 35,20, .20], 22: [4,6, 35,20, .20],
        23: [4,6, 35,20, .20], 24: [4,6, 35,20, .20]}

def runcommand(command):
    pfrom = os.popen(command)
    outputstr = pfrom.read()
    pfrom.close()
    return outputstr

def writefile(fn, data):
    try:
        if fn == '-':
            fd = sys.stdout
            if type(data) == type([]):
                fd.writelines(data)
            elif type(data) == type(''):
                fd.write(data)
        else:
            fd = open(fn, 'w')
            if type(data) == type([]):
                fd.writelines(data)
            elif type(data) == type(''):
                fd.write(data)
            fd.close()
    except IOError:
        print "IOError writing : %s " %fn
        return None

def format (filelist, plots, orient):

    offset = [0,0]
    psfile = "%!PS-Adobe-2.0\n"

    if orient == PORTRAIT :
        cols, rows, offset[0], offset[1], scaling = pinfo[plots]

        width = PAGESIZE[0] / cols
        height = PAGESIZE[1] / rows

        rrows = range(rows)
        rrows.reverse()
        for row in rrows:
            y = height * (row) + offset[1] + height

            for col in range(cols):
                x = width * (col) + offset[0] 

                if len(filelist) > 1:
                    psfile = psfile + plot_template %(x, y, \
                        scaling, scaling, -90, \
                        runcommand(getps %filelist[0]))
                    del filelist[0]
                elif filelist : 
                    psfile = psfile + plot_template %(x, y, \
                        scaling, scaling, -90, \
                        runcommand(getpspage %filelist[0]))
                    del filelist[0]

    elif orient == LANDSCAPE :
        rows, cols, offset[0], offset[1], scaling = linfo[plots]

        width = PAGESIZE[1] / cols
        height = PAGESIZE[0] / rows

        for row in range(rows):
            x = height * (row) + offset[0] 

            for col in range(cols):
                y = width * (col) + offset[1] 

                if len(filelist) > 1:
                    psfile = psfile + plot_template %(x, y, \
                        scaling, scaling, 0, \
                        runcommand(getps %filelist[0]))
                    del filelist[0]
                elif filelist : 
                    psfile = psfile + plot_template %(x, y, \
                        scaling, scaling, 0, \
                        runcommand(getpspage %filelist[0]))
                    del filelist[0]
    return psfile

def usage(msg) :
    sys.stdout = sys.stderr
    if msg: print msg
    print "usage: %s [-o outputfile] [-p #plots/page] [-L Landscape] [-P Portrait] (files)" \
        % os.path.basename(sys.argv[0])
    sys.exit(2)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'p:o:LPh?')
    except getopt.error, msg:
        usage(msg)

    orient = LANDSCAPE
    plots = len(args)
    filename = ''

    for o, a in opts:
        if o in ['-h','-?']:
            usage()
        elif o == '-p' :
            plots = int(a)
        elif o == '-L' :
            orient = LANDSCAPE
        elif o == '-P' :
            orient = PORTRAIT
        elif o == '-o' :
            filename = a

    else:
        if len(args) <= 1:
            usage('')

    if plots > 24: plots = 24

    data =''
    while args:
        data = data +format(args[:plots], plots, orient)
        del(args [:plots])

    if filename:
        writefile(filename, data)
    else:
        sys.stdout.write(data)

if __name__ == '__main__' :
    main()

"""
Example script to change the control parameters for printing 3 plots per page
in Portrait mode:

---------------------------------------------------------
#!/usr/bin/env python
from arrangeps import *

# Fool the height calculation into squeezing plots closer together
PAGESIZE[1] = 560

# Change to 1 column, 3 rows, no x-offset, 200 point y offset and 40%
# reduction
pinfo[3] = [1, 3,  0, 200, .40]

main()
---------------------------------------------------------
"""
