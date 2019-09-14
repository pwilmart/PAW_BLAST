"""Writes a subset FASTA protein database from a list of accessions in a text file
and a FASTA protein databases file.

Written by Phil Wilmarth, OHSU, 2016
"""
import os
import sys
import time
import copy
import fasta_lib

# control flags
KEEP_CONTAMS = False

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
if not subset_DB_name: sys.exit()
if os.path.splitext(subset_DB_name)[1] == '':
    subset_DB_name += '.fasta'

# open the accession list file and save accessions in dictionary
print('=======================================================')
print(' make_subset_DB_from_list.py, Phil Wilmarth, OHSU 2018 ')
print('=======================================================')
print('Processing accessions file:', time.ctime())
IDs = {}
read = 0
for acc in open(list_file_name, 'r'):
    try:
        acc = acc.split()[0]    # remove additional match numbers
    except IndexError:
        pass
    acc = acc.replace('_family', '')     # remove family designations
    if acc == '' or acc == 'Accession' or acc == 'Accessions' or acc.startswith('REV') or acc.startswith('('):
        print('...WARNING:', acc, 'was skipped')
        continue
    read += 1
    IDs[acc] = True    
print('...there were %s protein accessions in file' % read)
print('...there were %s unique protein accessions in file' % len(IDs))

# create instances of reader object and protein object, initialize counters
print('\n...Opening %s and reading proteins' % os.path.split(database_name)[1])
f = fasta_lib.FastaReader(database_name)
p = fasta_lib.Protein()
prot = 0
sub_prot = 0
contam = 1000
subset_DB = open(os.path.join(results_location, subset_DB_name), 'w')

"""Results files may parse out a subset of the accession string so simple
exact matching to protein database accessions may not work.
"""
# read proteins until EOF
while f.readNextProtein(p):
    prot += 1
    for key in IDs:
        if key in p.accession: # this can have some unexpected, extra matches
            if p.accession.startswith('CONT|') and KEEP_CONTAMS:
                contam += 1
                new_cont = 'CONT_%04d|' % (contam,)
                p.accession = p.accession.replace('CONT|', new_cont)
                p.parseCONT()
            if p.accession.startswith('CONT_') and not KEEP_CONTAMS:
                continue
            if IDs[key] == True:
                IDs[key] = [(prot, copy.deepcopy(p))]
            else:
                IDs[key].append((prot, copy.deepcopy(p)))

# we may have some contaminants in the imput list that we need to skip
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
for pout in [x[1] for x in sorted(out_list)]:
    pout.printProtein(subset_DB)
    sub_prot += 1

# close files and print summaries
subset_DB.close()
print('...%s proteins read' % prot)
print('...%s proteins written: %s' % (sub_prot, time.ctime()))

# end

