from Bio import Blast
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

#compares the generated sequence to the nonredundant blast database
def check_blast_dna(strand):
    strand = str(strand)
    "opens the blast database and gets the results for the strand that we have inputted"
    result_handle = NCBIWWW.qblast("blastn", "nr", strand)
    "saves the results to the file blast_output and closes the file"
    save_file = open("blast_output.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

    result_handle = open("blast_output.xml")
    blast_record = NCBIXML.read(result_handle)
    if blast_record.alignments != None:
        for alignment in blast_record.alignments: 
            for hsp in alignment.hsps:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')
    else:
         print("Blast does not recognize this sequence!")


def check_blast_protein(strand):
    strand = str(strand)
    "opens the blast database and gets the results for the strand that we have inputted"
    result_handle = NCBIWWW.qblast("blastp", "nr", strand)
    "saves the results to the file blast_output and closes the file"
    save_file = open("blast_output.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

    result_handle = open("blast_output.xml")
    blast_record = NCBIXML.read(result_handle)
    if blast_record.alignments != None:
        for alignment in blast_record.alignments: 
            for hsp in alignment.hsps:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')
    else:
         print("Blast does not recognize this sequence!")


def test_blast():
    result_handle = open("blast_output.xml")
    blast_record = NCBIXML.read(result_handle)
    if blast_record.alignments != None:
        for alignment in blast_record.alignments: 
            for hsp in alignment.hsps:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')
    else:
         print("Blast does not recognize this sequence!")