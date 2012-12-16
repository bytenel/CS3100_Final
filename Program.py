from Bio.Seq import *
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Alphabet import generic_protein

from DNAdfa import *
from DotDFA import *
from blast import *

"""Yes, this uses Ganesh's dfa.py file from assignment7! It was the easiest way to ensure that the dfa portion worked"""

def file_to_list(file,format,alphabet=0):
    recordList = list()

    if(alphabet==0):
        for record in SeqIO.parse("Source/" + file, format):
            print record.id
            print repr(record.seq)
            print len(record)
            recordList.append(record.seq)
    else:
        for record in SeqIO.parse("Source/" + file, format, alphabet):
            print record.id
            print repr(record.seq)
            print len(record)
            recordList.append(record.seq)
    return recordList

def print_annotations(record):
    for i in record.annotations:
        print str(i) + (": ") + str(rec.annotations[i])
        print ("")

#to implement later:
def generate_random_sequence(type, organism):

    start_codons = CodonTable.unambiguous_dna_by_id[organism].start_codons
    stop_codons = CodonTable.unambiguous_dna_by_id[organism].stop_codons

    #generate a random dna sequence
    if str(type.upper())=="DNA":
        print("dna")
    
    #generate a random protein sequence
    if str(type.upper())=="PROTEIN":
        print("protein")

    #generate a random rna sequence
    if str(type.upper())=="RNA":
        print("protein")

if __name__ == "__main__":
    print("Welcome to Carlos and Ben's Biopython Environment!")
    