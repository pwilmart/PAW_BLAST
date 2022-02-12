"""make_subset_DB_from_list_3.py -
writes a subset FASTA protein database from a list of accessions in a text file
and a FASTA protein databases file.

MIT License

Copyright (c) 2019 Phillip Wilmarth, OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Written by Phil Wilmarth, OHSU, 2016
"""
import os
import sys
import time
import copy
import fasta_lib

# control flags
KEEP_CONTAMS = False

def make_subset_db(list_file_name, database_name, subset_DB_name):
    # open the accession list file and save accessions in dictionary
    print('=======================================================')
    print(' make_subset_DB_from_list.py, Phil Wilmarth, OHSU 2018 ')
    print('=======================================================')
    print('Processing accessions file:', time.ctime())
    IDs = {}
    read = 0
    for acc in open(list_file_name, 'r'):
        # this assumes accessions from PAW pipeline results files
        try:
            acc = acc.split()[0]    # split at whitespace to remove "additional match" numbers
        except IndexError:
            pass
        
        # could add a split on ";" for PD or MaxQuant where multiple accessions can be in one cell
        # the design choice is to take first one only or all of them...
        
        acc = acc.replace('_family', '')     # remove family designations
        if acc == '' or acc == 'Accession' or acc == 'Accessions' or acc.startswith('REV') or acc.startswith('('):
            print('...WARNING:', acc, 'was skipped')
            continue # skip possible header element and/or decoys
        read += 1
        IDs[acc] = True # load a dictionary keyed by accession   
    print('...there were %s protein accessions in file' % read)
    print('...there were %s unique protein accessions in file' % len(IDs))

    # create instances of reader object and protein object, initialize counters
    print('\n...Opening %s and reading proteins' % os.path.split(database_name)[1])
    f = fasta_lib.FastaReader(database_name)
    p = fasta_lib.Protein()
    prot = 0
    sub_prot = 0
    contam = 1000
    subset_DB = open(subset_DB_name, 'w')

    """Results files may parse out a subset of the accession string so simple
    exact matching to protein database accessions may not work. We use a slower
    'in' test that can match accession substrings to longer accessions.
    """
    # read proteins until EOF
    while f.readNextProtein(p):
        prot += 1
        for key in IDs:
            if key in p.accession: # this can have some unexpected, extra matches
                if p.accession.startswith('CONT|'): # make contaminant accession format consistent
                    contam += 1
                    new_cont = 'CONT_%04d|' % (contam,)
                    p.accession = p.accession.replace('CONT|', new_cont)
                    p.parseCONT()
                if p.accession.startswith('CONT_') and not KEEP_CONTAMS:
                    continue
                if IDs[key] == True:
                    # replace value in IDs dictionary with protein object
                    IDs[key] = [(prot, copy.deepcopy(p))]
                else:
                    # extra matches
                    IDs[key].append((prot, copy.deepcopy(p)))

    # we may have some contaminants in the imput list that we need to skip
    # these might not be in the selected FASTA file and would still have True values
    new_IDs = copy.deepcopy(IDs)
    for key in IDs:
        if IDs[key] == True:
            del new_IDs[key]

    # we have all possible matches, try to remove some obvious extra matches
    for key in new_IDs:    
        new_list = [x for x in new_IDs[key] if key == x[1].accession] # maybe removes stuff like '85' in '853', '185', etc.
        if new_list and len(new_IDs[key]) > 1:
            new_IDs[key] = new_list

    # write out the matching proteins, need to collect all tuples in dictionary values
    out_list = []
    for key in new_IDs:
        for tup in new_IDs[key]:
            out_list.append(tup)

    for pout in [x[1] for x in sorted(out_list)]: # put matches in the same order as in the FASTA file
        pout.printProtein(subset_DB)
        sub_prot += 1

    # close files and print summaries
    subset_DB.close()
    print('...%s proteins read' % prot)
    print('...%s proteins written: %s' % (sub_prot, time.ctime()))
    # end

# FASTA databases and identification files can be passed via command line or interactively selected

if __name__ == '__main__':

    # FASTA files from the command line
    if len(sys.argv) == 4:

        list_file_name = sys.argv[1]
        database_name = sys.argv[2]
        subset_DB_name = sys.argv[3]

        # check the list file exists
        if not os.path.exists(list_file_name):
            print('FATAL : invalid identifications list file')
            sys.exit()
        # check the query database exists
        if not os.path.exists(database_name):
            print('FATAL : invalid query database file')
            sys.exit()

    # browse and select the FASTA databases and identification files
    elif len(sys.argv) == 1:
        # get the files, etc.
        list_file_name = fasta_lib.get_file(os.getcwd(), [('Text files', '*.txt')],
                                            'Browse to accession list text file')
        if not list_file_name: sys.exit()
        results_location = os.path.split(list_file_name)[0]
        database_name = fasta_lib.get_file(r'C:\Xcalibur\database', [('FASTA files', '*.fasta')],
                                        'Select the database')
        if not database_name: sys.exit()
        new_name = os.path.split(database_name)[1]
        new_name = os.path.splitext(new_name)[0]
        subset_DB_name = fasta_lib.save_file(results_location, [('FASTA files', '*.fasta')],
                                            default_file=new_name + '_subset.fasta',
                                            title_string='Name of subset database')

        print("saving subset datbase as", subset_DB_name)
        if not subset_DB_name: sys.exit()
        if os.path.splitext(subset_DB_name)[1] == '':
            subset_DB_name += '.fasta'

    else:
        # invalid command line
        print('FATAL: invalid command line argument')
        print('   Usage: make_subset_DB_from_list_3 [identifications_list_file, database_file, subset_output_path]')
        sys.exit()

    # run make subset database
    make_subset_db(list_file_name, database_name, subset_DB_name)
