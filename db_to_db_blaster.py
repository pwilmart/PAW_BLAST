"""Local BLAST program - launches local BLAST runs of small protein
databases against other small databases.  Caputures XML output and
parses the information.

Written by Phil Wilmarth, OHSU, 2009-2016.
"""
import os
import sys
from xml.sax.handler import ContentHandler
import math
#
#
#===============================================================
def get_file(default_location, extension_list, title_string=""):
#===============================================================
    """Dialog box to browse to a file.  Returns full file name.

    Usage: full_file_name = get_file(default_location, extension_list, [title]),
        where "default_location" is a starting folder location,
        extension_list is a list of (label, pattern) tuples,
        e.g. extension_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "full_file_name" is the complete name of the selected file.

    Written by Phil Wilmarth, OHSU, 2008.
    """
    import os
    import tkinter
    from tkinter import filedialog
    #
    # set up GUI elements
    #
    root = tkinter.Tk()
    root.withdraw()
    #
    # set default title string if not passed
    #
    if title_string == "":   
        title_string = 'Select a single FILE'
    #
    # create dialog box for file selection
    #
    root.update()
    filename = filedialog.askopenfilename(parent=root, initialdir=default_location,
                                          filetypes=extension_list, title=title_string)
    #
    # return full filename
    #
    return filename   
    #
    # end
    #
#
# data structure for general Blast information (version, parameters, etc.)
#
class Blast_header:
    def __init__(self):
        self.program = ''
        self.version = ''
        self.reference = ''
        self.query_db = ''
        self.hit_db = ''
        self.matrix = ''
        self.expect = 0.0
        self.gap_open = 11
        self.gap_extend = 1
        self.hit_num = 0
        self.hit_len = 0
        self.hit_hsp_len = 0
        self.hit_kappa = 0.0
        self.hit_lambda = 0.0
        self.hit_entropy = 0.0
#
# data structure for query information
#
class Blast_query:
    def __init__(self):
        self.num = 0
        self.id = ''
        self.desc = ''
        self.len = 0
        self.hits = []
        self.status = None
#
# data structure for hit information
#
class Blast_hit:
    def __init__(self):
        self.num = 0
        self.id = ''
        self.desc = ''
        self.accession = 0
        self.len = 0
        self.hit_hsps = []
#
# data structure for hsp information
#
class Blast_hsp:
    def __init__(self):
        self.num = 0
        self.bit_score = 0.0
        self.score = 0.0
        self.evalue = 0.0
        self.query_start = 0
        self.query_end = 0
        self.hit_start = 0
        self.hit_end = 0
        self.query_frame = 0
        self.hit_frame = 0
        self.identity = 0
        self.ident_pc = 0.0
        self.positive = 0
        self.pos_pc = 0.0
        self.gaps = 0
        self.align_len = 0
        self.query_seq = ''
        self.hit_seq = ''
        self.midline = ''
#
# simple content handler object for BLAST xml output
#
class SimpleHandler(ContentHandler):
    """Bare bones content handler object
    """
    def __init__(self):
        #
        # atrributes
        #
        self.inHeader = False
        self.inQuery = False
        self.inHit = False
        self.inHit_hsps = False
        self.inHsp = False
        self.inStats = False
        self.header = Blast_header()
        self.queries = []
        self.buff = ''
        self.get_chars = False
        self.last_name = ''
    
    def startElement(self, name, attrs):
        #
        # set state via flags as tags are encountered
        #
##        print('...start: ', name)
        if name == 'BlastOutput':
            self.inHeader = True
        elif name == 'Iteration':
            self.inQuery = True
            self.queries.append(Blast_query())
        elif name == 'Hit':
            self.inHit = True
            self.queries[-1].hits.append(Blast_hit())
        elif name == 'Hit_hsps':
            self.inHit_hsps = True
            self.inHit = False  # skip unecessay testing
        elif name == 'Hsp':
            self.inHsp = True
            self.queries[-1].hits[-1].hit_hsps.append(Blast_hsp())
        elif name == 'Statistics':
            self.inStats = True
        elif self.inHeader or self.inQuery or self.inHit \
             or self.inHsp or self.inStats:
            if name != self.last_name:  # read new chars when tag changes
                self.buff = ''          # reset character buffer
                self.get_chars = True   # set flag to capture text
                self.last_name = name   # save current tag name
        else:
            pass

    def characters(self, chars):
        #
        # capture text for the desired tags
        #
        if self.get_chars:
            self.buff += chars

    def endElement(self, name):
##        print('...end: ', name)
        #
        # save the desired XML elements into data structures
        #
        if name == 'Parameters':
            self.inHeader = False
        elif name == 'Iteration':
            self.inQuery = False
        elif name == 'Hits':
            self.inHit = False
        elif name == 'Hit_hsps':
            self.inHit_hsps = False
        elif name == 'Hsp':
            self.inHsp = False
        elif name == 'Statistics':
            self.inStats = False
        #
        # header elements
        #
        elif name == 'BlastOutput_program':
            self.header.program = self.buff
        elif name == 'BlastOutput_version':
            self.header.version = self.buff
        elif name == 'BlastOutput_reference':
            self.header.reference = self.buff
        elif name == 'Parameters_matrix':
            self.header.matrix = self.buff
        elif name == 'Parameters_expect':
            self.header.expect = float(self.buff)
        elif name == 'Parameters_gap-open':
            self.header.gap_open = float(self.buff)
        elif name == 'Parameters_gap-extend':
            self.header.gap_extend = float(self.buff)        
        #
        # Query elements
        #
        elif name == 'Iteration_iter-num':
            self.queries[-1].num = int(self.buff)
        elif name == 'Iteration_query-ID':
            self.queries[-1].id = self.buff
        elif name == 'Iteration_query-def':
            self.queries[-1].desc = self.buff
        elif name == 'Iteration_query-len':
            self.queries[-1].len = int(self.buff)
        #
        # Hit elements
        #
        elif name == 'Hit_num':
            self.queries[-1].hits[-1].num = int(self.buff)
        elif name == 'Hit_id':
            self.queries[-1].hits[-1].id = self.buff
        elif name == 'Hit_def':
            self.queries[-1].hits[-1].desc = self.buff
        elif name == 'Hit_accession':
            self.queries[-1].hits[-1].accession = self.buff
        elif name == 'Hit_len':
            self.queries[-1].hits[-1].len = int(self.buff)
        #
        # Hit hsp elements
        #
        elif name == 'Hsp_num':
            self.queries[-1].hits[-1].hit_hsps[-1].num = int(self.buff)
        elif name == 'Hsp_bit-score':
            self.queries[-1].hits[-1].hit_hsps[-1].bit_score = float(self.buff)
        elif name == 'Hsp_score':
            self.queries[-1].hits[-1].hit_hsps[-1].score = float(self.buff)
        elif name == 'Hsp_evalue':
            self.queries[-1].hits[-1].hit_hsps[-1].evalue = float(self.buff)
        elif name == 'Hsp_query-from':
            self.queries[-1].hits[-1].hit_hsps[-1].query_start = int(self.buff)
        elif name == 'Hsp_query-to':
            self.queries[-1].hits[-1].hit_hsps[-1].query_end = int(self.buff)
        elif name == 'Hsp_hit-from':
            self.queries[-1].hits[-1].hit_hsps[-1].hit_start = int(self.buff)
        elif name == 'Hsp_hit-to':
            self.queries[-1].hits[-1].hit_hsps[-1].hit_end = int(self.buff)
        elif name == 'Hsp_query-frame':
            self.queries[-1].hits[-1].hit_hsps[-1].query_frame = int(self.buff)
        elif name == 'Hsp_hit-frame':
            self.queries[-1].hits[-1].hit_hsps[-1].hit_frame = int(self.buff)
        elif name == 'Hsp_identity':
            self.queries[-1].hits[-1].hit_hsps[-1].identity = int(self.buff)
        elif name == 'Hsp_positive':
            self.queries[-1].hits[-1].hit_hsps[-1].positive = int(self.buff)
        elif name == 'Hsp_gaps':
            self.queries[-1].hits[-1].hit_hsps[-1].gaps = int(self.buff)
        elif name == 'Hsp_align-len':
            self.queries[-1].hits[-1].hit_hsps[-1].align_len = int(self.buff)
        elif name == 'Hsp_qseq':
            self.queries[-1].hits[-1].hit_hsps[-1].query_seq = self.buff
        elif name == 'Hsp_hseq':
            self.queries[-1].hits[-1].hit_hsps[-1].hit_seq = self.buff
        elif name == 'Hsp_midline':
            self.queries[-1].hits[-1].hit_hsps[-1].midline = self.buff
        #
        elif name == 'Statistics_db-num':
            self.header.hit_num = int(self.buff)
        elif name == 'Statistics_db-len':
            self.header.hit_len = int(self.buff)
        elif name == 'Statistics_hsp-len':
            self.header.hit_hsp_len = int(self.buff)
        elif name == 'Statistics_kappa':
            self.header.hit_kappa = float(self.buff)
        elif name == 'Statistics_lambda':
            self.header.hit_lambda = float(self.buff)
        elif name == 'Statistics_entropy':
            self.header.hit_entropy= float(self.buff)            
        #
        else:
            pass
    #            
    # end SimpleHandler class
    #

class Blast_results:
    """Place holder class for parsed BLAST xml results.
    """
    def __init__(self):
        self.header = None
        self.queries = None

    def snoop(self):
        print('Blast program:', self.header.program)
        print('program version:', self.header.version)
        print('Blast reference:', self.header.reference)
        print('query database:', self.header.query_db)
        print('hit database:', self.header.hit_db)
        print('substitution matrix:', self.header.matrix)
        print('output e-value cutoff:', self.header.expect)
        print('gap open penalty:', self.header.gap_open)
        print('gap extend penalty:', self.header.gap_extend)
        print('number of proteins in hit DB:', self.header.hit_num)
        print('number of amino acids in hit DB:', self.header.hit_len)
        print('hit hsp length:', self.header.hit_hsp_len)
        print('kappa statistic:', self.header.hit_kappa)
        print('lambda statistic:', self.header.hit_lambda)
        print('entropy statistic:', self.header.hit_entropy)
        print('number of query calls:', len(self.queries))

    def print_top_hits(self, out=None):
        print('\nProtein to Protein Blast Summary:\n', file=out)
        counter = 0
        for query in self.queries:
            print('\n%s> %s' % (query.num, query.desc[:]), file=out)
            try:
                h = query.hits[0]
                hsp = h.hit_hsps[0]
                print('... %s' % (h.desc[:],), file=out)
                print('...... ident:%d/%d pos:%d/%d hit:%d align:%d bit:%0.1f' %
                      (hsp.identity, query.len, hsp.positive, query.len,
                       h.len, hsp.align_len, hsp.bit_score), file=out)
                print('...... %s' % (query.status), file=out)
                if query.status == 'Poor_match' or query.status == 'Partial_match':
                    counter += 1
            except IndexError:
                counter += 1
                print('...... No Match', file=out)
        return(counter)
    
    def print_top_hits_tabs(self, out=None):
        print('\nTab-delimited Protein to Protein Blast Summary:\n', file=out)
        print('Query database:', self.header.query_db, file=out)
        print('Hit database:', self.header.hit_db, file=out)
        print('\nquery_number\tquery_acc\tquery_desc\thit_acc\thit_desc\tblast_scores\tstatus', file=out)
        counter = 0
        for query in self.queries:
            query_acc = query.desc.split()[0]
            query_desc = query.desc[len(query_acc):]
            try:
                h = query.hits[0]
                hsp = h.hit_hsps[0]
                hit_acc = h.desc.split()[0]
                hit_desc = h.desc[len(hit_acc)+1:]
                blast_scores = ('ident:%d/%d pos:%d/%d query:%d hit:%d align:%d bit:%0.1f' %
                                (hsp.identity, query.len, hsp.positive, query.len, query.len,
                                 h.len, hsp.align_len, hsp.bit_score))
                string = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (query.num, query_acc, query_desc,
                                                         hit_acc, hit_desc, blast_scores, query.status)
                print(string, file=out)
                if query.status == 'Poor_match' or query.status == 'Partial_match':
                    counter += 1
            except IndexError:
                counter += 1
                string = '%s\t%s\t%s\tNo_match\tNA\tNA\tNo_match' % (query.num, query_acc, query_desc)
                print(string, file=out)
        return(counter)
               
    def test_top_hits(self, score, cutoff):
        for query in self.queries:            
            try:
                h = query.hits[0]
                hsp = h.hit_hsps[0]
                query.status = 'OK'
                #
                # eliminate fragments
                #
                if (hsp.hit_end - hsp.hit_start)/float(h.len) < 0.5:
                    query.status = 'Partial_match_hit'
                if (hsp.query_end - hsp.query_start)/float(query.len) < 0.5:
                    if query.status == 'OK':
                        query.status = 'Partial_match_query'
                    else:
                        query.status = 'Partial_match_both'
                if query.status != 'OK':
                    continue
                if score == 'i':
                    value = hsp.ident_pc
                elif score == 'p':
                    value = hsp.pos_pc
                elif score == 'b':
                    value = hsp.bit_score
                else:
                    pass
                if value < cutoff:
##                    print(value, cutoff)
                    query.status = 'Poor_match'
            except IndexError:
                query.status = 'No_Blast_match'
        return

    def calc_percentage(self):
        for query in self.queries:
            try:
                hsp = query.hits[0].hit_hsps[0]
                self.ident_pc = 100.0*float(hsp.identity)/float(query.len)
                self.pos_pc = 100.0*float(hsp.positive)/float(query.len)
            except IndexError:
                pass

    def calc_ave_match(self, score):
        quant = []
        for query in self.queries:
            try:
                hsp = query.hits[0].hit_hsps[0]
                if score == 'i':
                    hsp.ident_pc = 100.0*float(hsp.identity)/float(query.len)
                    quant.append(100.0*float(hsp.identity)/float(query.len))
                elif score == 'p':
                    hsp.pos_pc = 100.0*float(hsp.positive)/float(query.len)
                    quant.append(100.0*float(hsp.positive)/float(query.len))
                elif score == 's':
                    quant.append(hsp.score)
                elif score == 'b':
                    quant.append(hsp.bit_score)
                elif score == 'e':
                    try:
                        quant.append(math.log10(hsp.evalue))
                    except ValueError:
                        pass
            except IndexError:
                pass
        #
        mean = sum(quant) / float(len(quant))
        try:
            stdev = sum([(float(x)-mean)**2 for x in quant]) / (float(len(quant))-1)
            stdev = math.sqrt(stdev)
        except ZeroDivisionError:
            stdev = 0.0
        lookup = {'i':'identity', 'p':'positive', 's':'score', \
                  'b':'bit_score', 'e':'log10_evalue'}
##        print('\n%s had a mean of %0.3f and stdev of %0.3f\n' % (lookup[score], mean, stdev))
        return(mean, stdev)
    #
    # end Blast_results class
    #

def better_match_check(results):
    """Flags matches if another protein is a better match

    written by Phil Wilmarth, OHSU, 2009.
    """
    #
    # build dictionary of best bit scores and accessions for the hit sequences
    #
    best_score = {}
    for q in results.queries:
        try:
            h = q.hits[0] # top hit
            h_acc = h.desc.split()[0]
#            hit_desc = h.desc[len(hit_acc):]            
            value = (h.hit_hsps[0].bit_score, h_acc)
            if h_acc in best_score:
                if best_score[h_acc][0] < value[0]:
                    best_score[h_acc] = value
            else:
                best_score[h_acc] = value
        except IndexError:
            pass
    #
    for q in results.queries:
        try:
            h = q.hits[0]
            h_acc = h.desc.split()[0]
            if best_score[h_acc][0] > h.hit_hsps[0].bit_score:
                q.status = 'Better_match to %s' % best_score[h_acc][1]
        except IndexError:
            pass
#        
#==============================================================        
# MAIN program to call local copy of blastp
# user supplies a local BLAST database path and
# a query sequence in FASTA format (supplied via stdin stream).
# XML output is captured, parsed, and saved in content handler.
#==============================================================        
#
import subprocess
import xml.sax
import copy
#
# print program information
#
print('\n=====================================================================')
print(' program "db_to_db_blaster.py", v1.1, Phil Wilmarth, OHSU, 2011, 2017')
print('=====================================================================')
#
# test platform and set the BLAST program path
#
if sys.platform == 'win32':
    blast_path = r'C:\Program Files\NCBI\blast-2.2.30+\bin'
else:
    blast_path = r'/usr/local/ncbi/blast/bin'
if not os.path.exists(blast_path):
    print('WARNING: BLAST path is not set correctly for this computer')
    print('...BLAST path was set to:', blast_path)
    print('...Aborting program. Please change "blast_path" and re-launch')
    sys.exit()
#
# get the two PAW results databases to blast against each other
#
if os.path.exists(r'C:\Xcalibur\database'):
    default = r'C:\Xcalibur\database'
else:
    default = os.getcwd()
#
print('Select first FASTA file')
first_db = get_file(default, [('FASTA files', '*.fasta')],\
                                'Select first database')
if not first_db: sys.exit()     # cancel button was hit
default = os.path.dirname(first_db)
print('Select the second FASTA file')
second_db = get_file(default, [('FASTA files', '*.fasta')],\
                                'Select second database')
if not second_db: sys.exit()    # cancel button was hit
#
# echo database names to console output
#
print('Query database:', os.path.basename(first_db))
print('Hit database:', os.path.basename(second_db))
print('Results files will be in:', os.path.dirname(first_db))
#
# create a data structure to hold the results
#
results = Blast_results()
#
# make the BLAST databases if they don't exist
#
for db in [first_db, second_db]:
    make_db = False
    for ext in ['.phr', '.pin', '.psq']:
        if not os.path.exists(db+ext):
            make_db = True
    if make_db:
        command = [os.path.join(blast_path, 'makeblastdb'),
                   '-in', db, '-dbtype', 'prot']
        print('Command line:', command)
        p = subprocess.call(command)
        print('Formatting %s as BLAST database' % (os.path.basename(db),))
#
# launch the local BLAST run
#
i = 0
query = first_db
hit = second_db
#
# make the Blast command line and launch subprocess
#
out_name = os.path.join(os.path.dirname(query),
                        os.path.basename(query)+'_vs_'+
                        os.path.basename(hit)+'.xml')
out_name = out_name.replace('.fasta', '')
if os.path.exists(out_name):
    print('BLAST XML results file exists, skipping BLAST run')
else:
    print('Starting local BLAST run (may take some time)...')
    command = [os.path.join(blast_path, 'blastp'), '-query', query, \
               '-db', hit, '-evalue', '10.0', '-outfmt', '5', '-out', out_name]
    print('Command line:', ' '.join(command))
##    blastp = subprocess.Popen(command,
##                              stdin=subprocess.PIPE,
##                              stdout=subprocess.PIPE,
##                              stderr=subprocess.PIPE)
    blastp = subprocess.Popen(command)
    blastp.wait()
    #
    # close output streams
    #
##    blastp.stdin.close()
##    blastp.stdout.close()
##    blastp.stderr.close()
###########################################################################    
###
### convert the temp.xml file to properly named file and strips bad chars
###
##temp_name = os.path.join(os.path.dirname(query), 'temp.xml')
##os.rename(out_name, temp_name)
##out_obj = open(out_name, 'w')
##for line in open(temp_name, 'r'):
##    line = line.rstrip()
##    print >>out_obj, line.replace(chr(255), '*')
##out_obj.close()
##os.remove(temp_name)
############################################################################
#
# set up to parse the XML output.  Content handler holds the BLAST results
#
print('Starting XML results file parsing (may take a few minutes)...')
out_obj = open(out_name, 'r')
ch = SimpleHandler()
ch.header.query_db = query
ch.header.hit_db = hit
saxparser = xml.sax.make_parser()
saxparser.setContentHandler(ch)
saxparser.setFeature(xml.sax.handler.feature_validation, 0)
saxparser.setFeature(xml.sax.handler.feature_namespaces, 0)
saxparser.setFeature(xml.sax.handler.feature_external_pes, 0)
saxparser.setFeature(xml.sax.handler.feature_external_ges, 0)
saxparser.parse(out_name)
#
# save the results before the next iteration
#
results.header = copy.deepcopy(ch.header)
results.queries = copy.deepcopy(ch.queries)
better_match_check(results)
#
# do something with the results next
#
result_file = out_name.replace('.xml', '.txt')
out = open(result_file, 'w')
score = 'i'
mean, stdev = results.calc_ave_match(score)
cutoff = mean - 3*stdev
if (mean-cutoff) > (0.9*mean):
    cutoff = 0.9*mean
results.test_top_hits(score, cutoff)
total = results.print_top_hits_tabs(out)
print('\n', total, 'proteins had no or poor matches', file=out)
print('Identity scores had a mean of %0.2f%% (%0.2f%%)' % (mean, stdev), file=out)
print('Identity cutoff for OK was %0.2f%%' % cutoff, file=out)
print(len(results.queries), 'proteins processed.', file=out)
out.close()
print(len(results.queries), 'proteins processed.')
#
# end
#
