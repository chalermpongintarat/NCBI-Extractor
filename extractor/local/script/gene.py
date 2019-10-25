#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from geopy.geocoders import Nominatim
from Bio.Alphabet import generic_dna, generic_protein

# Module: Gene Features Extract
button = ""

# Exception Extract ##################################################
def exception():
    while True:
        # Module Name
        print("Module 1: Exception Extract")
        print("")
        # Parsing Exception
        print("Step 1: Parsing Exception, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        exception_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        exception = feature.qualifiers.get('exception')
                        accession_list_with_none.append(record.id)
                        exception_list_with_none.append(exception)
        exception_list_with_none_list = [['None'] if exception_value is None else exception_value for exception_value in exception_list_with_none]
        exception_list_with_none_list_none = []
        for exception_none in exception_list_with_none_list:
            exception_list_with_none_list_none.append('\n'.join(exception_none))
        # Zipping Exception
        exception_zip = zip(accession_list_with_none,exception_list_with_none_list_none)
        print(exception_zip)
        exception_list = []
        for accession_exception in exception_zip:
            exception_list.append("\t".join(accession_exception))
        print("")
        # Writing Exception
        print("Step 2: Writing Exception, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        exception_string = '\n'.join(str(exception_value) for exception_value in exception_list)
        print(exception_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_exception.txt"), "w")
        handle_outputfile.write(str(exception_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 1: Exception Extract Complete")
        print("")
        break

# Gene Extract ##################################################
def gene():
    while True:
        # Module Name
        print("Module 2: Gene Extract")
        print("")
        # Parsing Gene
        print("Step 1: Parsing Gene, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        gene_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        gene = feature.qualifiers.get('gene')
                        accession_list_with_none.append(record.id)
                        gene_list_with_none.append(gene)
        gene_list_with_none_list = [['None'] if gene_value is None else gene_value for gene_value in gene_list_with_none]
        gene_list_with_none_list_none = []
        for gene_none in gene_list_with_none_list:
            gene_list_with_none_list_none.append('\n'.join(gene_none))
        # Zipping Gene
        gene_zip = zip(accession_list_with_none,gene_list_with_none_list_none)
        print(gene_zip)
        gene_list = []
        for accession_gene in gene_zip:
            gene_list.append("\t".join(accession_gene))
        print("")
        # Writing Gene
        print("Step 2: Writing Gene, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        gene_string = '\n'.join(str(gene_value) for gene_value in gene_list)
        print(gene_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_gene.txt"), "w")
        handle_outputfile.write(str(gene_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 2: Gene Extract Complete")
        print("")
        break

# Gene DESC Extract ##################################################
def genedesc():
    while True:
        # Module Name
        print("Module 3: Gene DESC Extract")
        print("")
        # Parsing Gene DESC
        print("Step 1: Parsing Gene DESC, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        gene_desc_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        gene_desc = feature.qualifiers.get('gene_desc')
                        accession_list_with_none.append(record.id)
                        gene_desc_list_with_none.append(gene_desc)
        gene_desc_list_with_none_list = [['None'] if gene_desc_value is None else gene_desc_value for gene_desc_value in gene_desc_list_with_none]
        gene_desc_list_with_none_list_none = []
        for gene_desc_none in gene_desc_list_with_none_list:
            gene_desc_list_with_none_list_none.append('\n'.join(gene_desc_none))
        # Zipping Gene DESC
        gene_desc_zip = zip(accession_list_with_none,gene_desc_list_with_none_list_none)
        print(gene_desc_zip)
        gene_desc_list = []
        for accession_gene_desc in gene_desc_zip:
            gene_desc_list.append("\t".join(accession_gene_desc))
        print("")
        # Writing Gene DESC
        print("Step 2: Writing Gene DESC, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        gene_desc_string = '\n'.join(str(gene_desc_value) for gene_desc_value in gene_desc_list)
        print(gene_desc_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_gene_desc.txt"), "w")
        handle_outputfile.write(str(gene_desc_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 3: Gene DESC Extract Complete")
        print("")
        break

# Locus Tag Extract ##################################################
def locustag():
    while True:
        # Module Name
        print("Module 4: Locus Tag Extract")
        print("")
        # Parsing Locus Tag
        print("Step 1: Parsing Locus Tag, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        locus_tag_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        locus_tag = feature.qualifiers.get('locus_tag')
                        accession_list_with_none.append(record.id)
                        locus_tag_list_with_none.append(locus_tag)
        locus_tag_list_with_none_list = [['None'] if locus_tag_value is None else locus_tag_value for locus_tag_value in locus_tag_list_with_none]
        locus_tag_list_with_none_list_none = []
        for locus_tag_none in locus_tag_list_with_none_list:
            locus_tag_list_with_none_list_none.append('\n'.join(locus_tag_none))
        # Zipping Locus Tag
        locus_tag_zip = zip(accession_list_with_none,locus_tag_list_with_none_list_none)
        print(locus_tag_zip)
        locus_tag_list = []
        for accession_locus_tag in locus_tag_zip:
            locus_tag_list.append("\t".join(accession_locus_tag))
        print("")
        # Writing Locus Tag
        print("Step 2: Writing Locus Tag, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        locus_tag_string = '\n'.join(str(locus_tag_value) for locus_tag_value in locus_tag_list)
        print(locus_tag_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_locus_tag.txt"), "w")
        handle_outputfile.write(str(locus_tag_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 4: Locus Tag Extract Complete")
        print("")
        break

# Note Extract ##################################################
def note():
    while True:
        # Module Name
        print("Module 5: Note Extract")
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
                    if feature.type == "gene":
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
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_note.txt"), "w")
        handle_outputfile.write(str(note_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 5: Note Extract Complete")
        print("")
        break

# Pseudo Extract ##################################################
def pseudo():
    while True:
        # Module Name
        print("Module 6: Pseudo Extract")
        print("")
        # Parsing Pseudo
        print("Step 1: Parsing Pseudo, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        pseudo_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        pseudo = feature.qualifiers.get('pseudo')
                        accession_list_with_none.append(record.id)
                        pseudo_list_with_none.append(pseudo)
        pseudo_list_with_none_list = [['None'] if pseudo_value is None else pseudo_value for pseudo_value in pseudo_list_with_none]
        pseudo_list_with_none_list_none = []
        for pseudo_none in pseudo_list_with_none_list:
            pseudo_list_with_none_list_none.append('\n'.join(pseudo_none))
        # Zipping Pseudo
        pseudo_zip = zip(accession_list_with_none,pseudo_list_with_none_list_none)
        print(pseudo_zip)
        pseudo_list = []
        for accession_pseudo in pseudo_zip:
            pseudo_list.append("\t".join(accession_pseudo))
        print("")
        # Writing Pseudo
        print("Step 2: Writing Pseudo, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        pseudo_string = '\n'.join(str(pseudo_value) for pseudo_value in pseudo_list)
        print(pseudo_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_pseudo.txt"), "w")
        handle_outputfile.write(str(pseudo_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 6: Pseudo Extract Complete")
        print("")
        break

# Trans Splicing Extract ##################################################
def transsplicing():
    while True:
        # Module Name
        print("Module 7: Trans Splicing Extract")
        print("")
        # Parsing Trans Splicing
        print("Step 1: Parsing Trans Splicing, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        trans_splicing_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        trans_splicing = feature.qualifiers.get('trans_splicing')
                        accession_list_with_none.append(record.id)
                        trans_splicing_list_with_none.append(trans_splicing)
        trans_splicing_list_with_none_list = [['None'] if trans_splicing_value is None else trans_splicing_value for trans_splicing_value in trans_splicing_list_with_none]
        trans_splicing_list_with_none_list_none = []
        for trans_splicing_none in trans_splicing_list_with_none_list:
            trans_splicing_list_with_none_list_none.append('\n'.join(trans_splicing_none))
        # Zipping Trans Splicing
        trans_splicing_zip = zip(accession_list_with_none,trans_splicing_list_with_none_list_none)
        print(trans_splicing_zip)
        trans_splicing_list = []
        for accession_trans_splicing in trans_splicing_zip:
            trans_splicing_list.append("\t".join(accession_trans_splicing))
        print("")
        # Writing Trans Splicing
        print("Step 2: Writing Trans Splicing, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        trans_splicing_string = '\n'.join(str(trans_splicing_value) for trans_splicing_value in trans_splicing_list)
        print(trans_splicing_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_trans_splicing.txt"), "w")
        handle_outputfile.write(str(trans_splicing_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 7: Trans Splicing Extract Complete")
        print("")
        break

# All Extract ##################################################
def allextract():
    while True:
        # Module Name
        print("Module 8: All Extract")
        print("")
        # Parsing All
        print("Step 1: Parsing All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        # Exception
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        exception_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        exception = feature.qualifiers.get('exception')
                        accession_list_with_none.append(record.id)
                        exception_list_with_none.append(exception)
        exception_list_with_none_list = [['None'] if exception_value is None else exception_value for exception_value in exception_list_with_none]
        exception_list = []
        for exception_none in exception_list_with_none_list:
            exception_list.append('\n'.join(exception_none))
        # Gene
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        gene_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        gene = feature.qualifiers.get('gene')
                        accession_list_with_none.append(record.id)
                        gene_list_with_none.append(gene)
        gene_list_with_none_list = [['None'] if gene_value is None else gene_value for gene_value in gene_list_with_none]
        gene_list = []
        for gene_none in gene_list_with_none_list:
            gene_list.append('\n'.join(gene_none))
        # Gene DESC
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        gene_desc_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        gene_desc = feature.qualifiers.get('gene_desc')
                        accession_list_with_none.append(record.id)
                        gene_desc_list_with_none.append(gene_desc)
        gene_desc_list_with_none_list = [['None'] if gene_desc_value is None else gene_desc_value for gene_desc_value in gene_desc_list_with_none]
        gene_desc_list = []
        for gene_desc_none in gene_desc_list_with_none_list:
            gene_desc_list.append('\n'.join(gene_desc_none))
        # Locus Tag
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        locus_tag_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        locus_tag = feature.qualifiers.get('locus_tag')
                        accession_list_with_none.append(record.id)
                        locus_tag_list_with_none.append(locus_tag)
        locus_tag_list_with_none_list = [['None'] if locus_tag_value is None else locus_tag_value for locus_tag_value in locus_tag_list_with_none]
        locus_tag_list = []
        for locus_tag_none in locus_tag_list_with_none_list:
            locus_tag_list.append('\n'.join(locus_tag_none))
        # Note
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        note_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        note = feature.qualifiers.get('note')
                        note_list_with_none.append(note)
        note_list_with_none_list = [['None'] if note_value is None else note_value for note_value in note_list_with_none]
        note_list = []
        for note_none in note_list_with_none_list:
            note_list.append('\n'.join(note_none))
        # Pseudo
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        pseudo_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        pseudo = feature.qualifiers.get('pseudo')
                        accession_list_with_none.append(record.id)
                        pseudo_list_with_none.append(pseudo)
        pseudo_list_with_none_list = [['None'] if pseudo_value is None else pseudo_value for pseudo_value in pseudo_list_with_none]
        pseudo_list = []
        for pseudo_none in pseudo_list_with_none_list:
            pseudo_list.append('\n'.join(pseudo_none))
        # Trans Splicing
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        trans_splicing_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "gene":
                        trans_splicing = feature.qualifiers.get('trans_splicing')
                        accession_list_with_none.append(record.id)
                        trans_splicing_list_with_none.append(trans_splicing)
        trans_splicing_list_with_none_list = [['None'] if trans_splicing_value is None else trans_splicing_value for trans_splicing_value in trans_splicing_list_with_none]
        trans_splicing_list = []
        for trans_splicing_none in trans_splicing_list_with_none_list:
            trans_splicing_list.append('\n'.join(trans_splicing_none))
        # Zipping All Extract
        all_extract_zip = zip(accession_list_with_none,exception_list,gene_list,gene_desc_list,locus_tag_list,note_list,pseudo_list,trans_splicing_list)
        print(all_extract_zip)
        all_extract_list = []
        for all_feature in all_extract_zip:
            all_extract_list.append("\t".join(all_feature))
        print("")
        # Writing All
        print("Step 2: Writing All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        all_string = '\n'.join(str(all_value) for all_value in all_extract_list)
        print(all_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/gene_all_extract.txt"), "w")
        handle_outputfile.write(str(all_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 8: All Extract Complete")
        print("")
        break

# Module Main ##################################################
def maingene():
    while True:
        print("")
        print("\tMenu: Extract Gene")
        print("")
        print("\t\tModule 1.  Exception \t\t\t: Button[1]")
        print("\t\tModule 2.  Gene \t\t\t: Button[2]")
        print("\t\tModule 3.  Gene DESC \t\t\t: Button[3]")
        print("\t\tModule 4.  Locus Tag \t\t\t: Button[4]")
        print("\t\tModule 5.  Note \t\t\t: Button[5]")
        print("\t\tModule 6.  Pseudo \t\t\t: Button[6]")
        print("\t\tModule 7.  Trans Splicing \t\t: Button[7]")
        print("\t\tModule 8.  All Extract \t\t\t: Button[8]")
        print("\t\tModule 9.  Main Menu \t\t\t: Button[9]")
        print("\t\tModule 10. Exit. \t\t\t: Button[10]")
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
                exception()
            elif button == 2:
                print("")
                gene()
            elif button == 3:
                print("")
                genedesc()
            elif button == 4:
                print("")
                locustag()
            elif button == 5:
                print("")
                note()
            elif button == 6:
                print("")
                pseudo()
            elif button == 7:
                print("")
                transsplicing()
            elif button == 8:
                print("")
                allextract()
            elif button == 9:
                print("")
                import main
                main.main()
            elif button == 10:
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

maingene()
