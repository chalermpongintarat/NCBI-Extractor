#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from geopy.geocoders import Nominatim
from Bio.Alphabet import generic_dna, generic_protein

# Module: MISC Features Extract
button = ""

# Note Extract ##################################################
def note():
    while True:
        # Module Name
        print("Module 1: Note Extract")
        print("")
        # Parsing Note
        print("Step 1: Parsing Note, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        note_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "misc_feature":
                        note = feature.qualifiers.get('note')
                        accession_list_with_none.append(record.id)
                        note_list_with_none.append(note)
        note_list_with_none_list = [['None'] if note_value is None else note_value for note_value in note_list_with_none]
        note_list_with_none_list_none = []
        for note_none in note_list_with_none_list:
            note_list_with_none_list_none.append('\n'.join(note_none))
        # Zipping Note
        note_zip = zip(accession_list_with_none,note_list_with_none_list_none)
        print(note_zip)
        note_list = []
        for accession_note in note_zip:
            note_list.append("\t".join(accession_note))
        print("")
        # Writing Note
        print("Step 2: Writing Note, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        note_string = '\n'.join(str(note_value) for note_value in note_list)
        print(note_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/misc_feature_note.txt"), "w")
        handle_outputfile.write(str(note_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 1: Note Extract Complete")
        print("")
        break

# Module Main ##################################################
def mainmiscfeature():
    while True:
        print("")
        print("\tMenu: Extract MISC Features")
        print("")
        print("\t\tModule 1. Note \t\t\t\t: Button[1]")
        print("\t\tModule 2. Main Menu \t\t\t: Button[2]")         
        print("\t\tModule 3. Exit. \t\t\t: Button[3]")
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
                note()
            elif button == 2:
                print("")
                import main
                main.main()
            elif button == 3:
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

mainmiscfeature()
