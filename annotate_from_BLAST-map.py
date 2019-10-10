"""annotate_from_BLAST_map.py -
adds "query" annotations to UniProt "hit" database.

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
SOFTWARE."""

import os
import sys
import copy
import fasta_lib

# bit scores of 75-100 seem to be above random scores
MIN_BIT = 100.0

# get the two PAW results databases to blast against each other
default = os.getcwd()

print('Select the BLAST mapping file file')
anno_file = fasta_lib.get_file(default, [('TXT files', '*.txt')], 'Select mapping file')
if not anno_file:
    sys.exit()     # cancel button was hit
    
default = os.path.dirname(anno_file)
print('Select the FASTA file')
orig_database = fasta_lib.get_file(default, [('FASTA files', '*.fasta')], 'Select database')
if not orig_database:
    sys.exit()    # cancel button was hit

new_database = orig_database.replace('.fasta', '_fixed.fasta')

# echo database names to console output
print('Mapping file:', os.path.basename(anno_file))
print('Original database:', os.path.basename(orig_database))
print('New database:', os.path.basename(new_database))

class Annotation:
    """Annotation record."""
    def __init__(self, line, col_map, counters, matches):
        self.query_acc = None
        self.query_header = None
        self.query_desc = None
        self.query_gene = None
        self.hit_ident = None
        self.hit_desc = None
        self.hit_gene = None
        self.new_desc = None

        self.parse_line(line, col_map, counters, matches)
        self.make_new_desc()

    def parse_line(self, line, col_map, counters, matches):
        """Parses BLAST mapping line."""
        line = line.rstrip()
        items = line.split('\t')
        self.query_acc = items[col_map['query_acc']]
        self.query_header = items[col_map['query_desc']].lstrip()
        self.query_desc, self.query_gene = self.parse_desc(self.query_header)
        status = items[col_map['match_status']]
        for idx, match in enumerate(matches):
            if status == match:
                counters[idx] += 1
        if status == 'No_match':
            return
        if status == 'Poor_match':
            bit = float(items[col_map['bit_score']])
            if bit < MIN_BIT:
                counters[-1] += 1
                return
        self.hit_ident = items[col_map['hit_acc']].split('|')[-1]
        self.hit_desc, self.hit_gene = self.parse_desc(items[col_map['hit_desc']])
        
    def parse_desc(self, desc):
        """Parses description and gene symbol from UniProt FASTA headers."""
        desc = desc.strip()
        tags = [' OS=', ' OX=', ' GN=', ' PE=', ' SV=']
        indexes = [desc.find(x) for x in tags]
        first = indexes[0]
        if first == -1:
##            print('NO OS= string')
            description = desc
        else:
            description = desc[: first]

        if indexes[2] == -1:
            gene = None
        else:
            gene = desc[indexes[2] + 4:].split()[0]

        return description, gene

    def make_new_desc(self):
        """Combines description and gene symbols into composite FASTA header."""
        new_desc = self.query_desc
        if self.query_gene:
            new_desc += ' GN=%s' % self.query_gene
        if self.hit_ident and self.hit_desc:
            new_desc += ' [%s: %s' % (self.hit_ident, self.hit_desc)
        else:
            new_desc += ' [NA'
        if self.hit_gene:
            new_desc += ' GN=%s]' % self.hit_gene
        else:
            new_desc += ']'

        self.new_desc = new_desc
        

# figure out where top of table is
map_buffer = list(open(anno_file, 'rt').readlines())
for skip, line in enumerate(map_buffer):
    if 'query_number\tquery_acc\t' in line:
        col_map = {header.strip(): index for (index, header) in enumerate(line.rstrip().split('\t'))}
        break

# parse information out of the table and save in object disctionary
maps = {}
matches = ['No_match', 'OK', 'Partial_match_both', 'Partial_match_hit',
           'Partial_match_query', 'Poor_match', 'rejected (low bit_score)']
counters = [0 for dummy in matches] 

for line in map_buffer[skip+1:-5]:
    next_map = Annotation(line, col_map, counters, matches)
    maps[next_map.query_acc] = next_map

# read FASTA file, change description, write new FASTA file
new_db_obj = open(new_database, 'wt')
f = fasta_lib.FastaReader(orig_database)

prot = fasta_lib.Protein()
pcount = 0
while f.readNextProtein(prot):
    pcount += 1
    prot.new_desc = maps[prot.new_acc].new_desc
    prot.printProtein(new_db_obj)
new_db_obj.close()

# print status message
print(pcount, 'proteins were processed')
for idx, match in enumerate(matches):
    print('...%d proteins had %s status' % (counters[idx], match))
print((pcount-counters[0]-counters[-1]), 'proteins were annotated')

# double check that nothing got altered
f_orig = fasta_lib.FastaReader(orig_database)
p_orig = {}
while f_orig.readNextProtein(prot):
    p_orig[prot.new_acc] = copy.deepcopy(prot)
f_new = fasta_lib.FastaReader(new_database)
p_new = {}
while f_new.readNextProtein(prot):
    p_new[prot.new_acc] = copy.deepcopy(prot)

for key in p_orig.keys():
    if p_orig[key].accession != p_new[key].accession:
        print('acc mismatch:', key)
    if p_orig[key].sequence != p_new[key].sequence:
        print('seq mismatch', key)

    
    
