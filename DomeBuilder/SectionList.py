import csv

import numpy as np

class SectionList(object):

    __CHS = None

    @staticmethod
    def load_CHS(verbose=True):
        if SectionList.__CHS is None:
            SectionList.__CHS = SectionList('resources/CHS.csv', verbose)
        return SectionList.__CHS
        
    def __init__(self, file_name, verbose):
        self.catalogue = ['name','d-yy','d-zz','t-yy','t-zz','area','weight','Iyy','Qyy','ryy','Izz','Qzz','rzz','Zyy','Zzz','J']
        self.index = {}
        self.fields = len(self.catalogue)

        for c in range(0, len(self.catalogue)):
            self.index[self.catalogue[c]] = c

        self.headers = {}
        self.names = []

        self.count = 0
        self.data = None

        if verbose:
            print("Loading section data from '{f}'...".format(f=file_name))

        with open(file_name, newline='') as csvfile:
            reader = csv.reader(csvfile)
            bHeader = True
            for row in reader:
                if bHeader:
                    bHeader = False
                    for c in range(0, len(self.catalogue)):
                        self.headers[self.catalogue[c]] = row[c]
                else:
                    self.names.append(row[0])
                    values = list(map(float, row[1:self.fields]))
                    if self.count == 0:
                        self.data = np.asarray([[0, *values],])
                    else:
                        self.data = np.append(self.data, [[self.count, *values],], axis=0)
                    self.count = self.count + 1

        if verbose:
            print(' => {c} sections loaded'.format(c=self.count))

    def get(self, section_index, section_property):
        p = None
        if section_index < self.count:
            i = self.index[section_property]
            if i is not None:
                if i == 0:
                    p = self.names[int(self.data[section_index,0])]
                else:
                    p = self.data[section_index,i]
        return p
