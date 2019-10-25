#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from geopy.geocoders import Nominatim
from Bio.Alphabet import generic_dna, generic_protein

# Module: Gene Features Extract
button = ""

# Convert File to FASTA ##################################################
def converttofasta():
    while True:
        # Module Name
        print("Module 1: Convert File to FASTA")
        print("")
        # Parsing File
        print("Step 1: Parsing FIle, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        print("Continue...")
        print("")
        # Writing File
        print("Step 2: Writing File, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_convert_db = SeqIO.convert(os.path.expanduser("~/extractor/local/import/sequence.gb.txt"), "genbank", os.path.expanduser("~/extractor/local/import/sequence.fasta"), "fasta")
        print("Convert to FATSA : %i records" % handle_convert_db)
        print("")
        # Module Complete
        print("Module 1: Convert File to FASTA Complete")
        print("")
        break

# Convert File to GenBank ##################################################
def converttogenbank():
    while True:
        # Module Name
        print("Module 2: Convert File to GenBank")
        print("")
        # Parsing File
        print("Step 1: Parsing FIle, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        print("Continue...")
        print("")
        # Writing File
        print("Step 2: Writing File, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_convert_db = SeqIO.convert(os.path.expanduser("~/extractor/local/import/sequence.gb.txt"), "genbank", os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        print("Convert to GenBank : %i records" % handle_convert_db)
        print("")
        # Module Complete
        print("Module 2: Convert File to GenBank Complete")
        print("")
        break

# Convert File FASTA to GenBank ##################################################
def convertfastatogenbank():
    while True:
        # Module Name
        print("Module 3: Convert File FASTA to GenBank")
        print("")
        # Parsing File
        print("Step 1: Parsing FIle, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        print("Continue...")
        print("")
        # Writing File
        print("Step 2: Writing File, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #handle_convert_db = SeqIO.convert("/Users/chalermpong/extractor/local/import/sequence.gb.txt", "genbank", "/Users/chalermpong/extractor/local/import/sequence.gbk", "genbank")
        input_handle = open(os.path.expanduser("~/extractor/local/import/sequence.fasta"), "rU")
        output_handle = open(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "w")
        sequences = list(SeqIO.parse(input_handle, "fasta"))
        #asign generic_dna or generic_protein
        for seq in sequences:
            seq.seq.alphabet = generic_dna
        count = SeqIO.write(sequences, output_handle, "genbank")
        output_handle.close()
        input_handle.close()
        print("Convert FASTA to GenBank : %i records" % count)
        print("")
        # Module Complete
        print("Module 3: Convert File FASTA to GenBank Complete")
        print("")
        break

# All Convert File ##################################################
def allconvert():
    while True:
        # Module Name
        print("Module 4: All Convert File")
        print("")
        # Parsing File
        print("Step 1: Parsing FIle, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        print("Continue...")
        print("")
        # Writing File
        print("Step 2: Writing File, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_convert_db = SeqIO.convert(os.path.expanduser("~/extractor/local/import/sequence.gb.txt"), "genbank", os.path.expanduser("~/extractor/local/import/sequence.fasta"), "fasta")
        handle_convert_db = SeqIO.convert(os.path.expanduser("~/extractor/local/import/sequence.gb.txt"), "genbank", os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        print("Convert to FATSA : %i records" % handle_convert_db)
        print("Convert to GenBank : %i records" % handle_convert_db)
        print("")
        # Module Complete
        print("Module 4: Convert File to GenBank Complete")
        print("")
        break

# Module Main ##################################################
def mainconvert():
    while True:
        print("")
        print("\tMenu: Convert File")
        print("")
        print("\t\tModule 1. To FASTA \t\t\t: Button[1]")
        print("\t\tModule 2. TO GenBank \t\t\t: Button[2]")
        print("\t\tModule 3. FASTA TO GenBank \t\t: Button[3]")
        print("\t\tModule 4. All Convert \t\t\t: Button[4]")
        print("\t\tModule 5. Main Menu \t\t\t: Button[5]")        
        print("\t\tModule 6. Exit. \t\t\t: Button[6]")
        print("")
        try:
            button = int(input("\t\t\tPlease push your Button. \t: "))
        except ValueError:
            print("")
            print("\t\t\t\t\t\t\t: No Module.")
            print("")
            continue
        else:
            if button == 1:
                print("")
                converttofasta()
            elif button == 2:
                print("")
                converttogenbank()
            elif button == 3:
                print("")
                convertfastatogenbank()
            elif button == 4:
                print("")
                allconvert()
            elif button == 5:
                print("")
                import main
                main.main()
            elif button == 6:
                print("")
                bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                for i in range(100):
                    time.sleep(0.1)
                    bar.update(i)
                print(" Exit.")
                print("")
                sys.exit(0)
            else:
                print("")
                print("\t\t\t\t\t\t\t: No Module.")
                print("")

mainconvert()
