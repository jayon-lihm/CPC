import os
import os.path
import numpy

def execute(com, debug=False):
    if debug==True:
        print (com)
    os.system(com)

def isExist(filename):
    if os.path.exists(filename) and os.path.isfile(filename):
        return True
    else:
        return False

def fileSize(filename):
    return int(os.path.getsize(filename))

def delete(filename):
    if os.path.exists(filename) and os.path.isfile(filename):
        os.unlink(filename)

""" Linear search, returns first index """
def find_first_index(lst, elem):
    ind = -1
    for l in lst:
        if str(elem)==str(l):
            ind=0

    return ind

def read_one_str_col(filename):
    fh = open(filename, "r")
    values = []
    for line in fh:
        values.append(line.strip('\r\n'))
    return values


def txt2gzip(txtfile, debug=True):
    if (isExist(txtfile) and fileSize(txtfile) > 0):
        com = 'gzip -f ' + txtfile
        execute(com, debug)
        if debug==True:
            print(txtfile + " File created ")


def gzip2txt(zipfile, debug=True):
    if (isExist(zipfile) and fileSize(zipfile) > 0):
        com = 'gunzip -f ' + zipfile
        execute(com, debug)
        if debug==True:
            print(zipfile + " File created ")

