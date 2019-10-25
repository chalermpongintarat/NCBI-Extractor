#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from geopy.geocoders import Nominatim
from Bio.Alphabet import generic_dna, generic_protein

# Module: Extract Features
button = ""

def main():
    while True:
        #handle_parse_db = list(SeqIO.parse("/Users/chalermpong/Desktop/Bioinformatics/Bio_Info_Learning/06_Weeks/05_DB Grouping/Local/sequence.gbk", "genbank"))
        #db_length = len(handle_parse_db)
        print("")
        print("====================================================================")
        print(" Chalermpong Intarat")
        print(" NBT, NSTDA")
        print(" chalermpong.int@biotec.or.th")
        print(" Oct 22, 2019")
        print("====================================================================")
        print("")
        print("Main Menu: Feature Extractor")
        print("")
        print("\tModule 1. Convert File\t\t\t: Button[1]")
        print("\tModule 2. Extract Description\t\t: Button[2]")
        print("\tModule 3. Extract Source\t\t: Button[3]")
        print("\tModule 4. Extract Gene\t\t\t: Button[4]")
        print("\tModule 5. Extract CDS\t\t\t: Button[5]")
        print("\tModule 6. Extract MISC Features\t\t: Button[6]")
        print("\tModule 7. Exit. \t\t\t: Button[7]")
        print("")
        #print("\t\tDB Download... " + str(db_length)) + " Sequences"
        #print("")
        try:
            button = int(input("\t\tPlease push your Button. \t: "))
        except ValueError:
            print("")
            print("\t\t\t\t\t\t\t: No Module.")
            print("")
            continue
        else:
            if button == 1:
                print("")
                import convert
                convert.mainconvert()
            elif button == 2:
                print("")
                import description
                description.maindescription()
            elif button == 3:
                print("")
                import source
                source.mainsource()
            elif button == 4:
                print("")
                import gene
                gene.maingene()
            elif button == 5:
                print("")
                import cds
                cds.maincds()
            elif button == 6:
                print("")
                import misc_feature
                misc_feature.mainmiscfeature()
            elif button == 7:
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

if __name__ == "__main__":
    main()
