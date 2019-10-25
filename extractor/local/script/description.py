#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from geopy.geocoders import Nominatim
from Bio.Alphabet import generic_dna, generic_protein

# Module: Description Features Extract
button = ""

# Accession Number Extract ##################################################
def accessionnumber():
    while True:
        # Module Name
        print("Module 1: Accession Number Extract")
        print("")
        # Parsing Accession Number
        print("Step 1: Parsing Accession Number, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_number_list_with_none = []
        for record in handle_parse_db:
            accession_number_list_with_none.append(record.id)
        accession_number_list_with_none_list = [['None'] if accession_number_value is None else accession_number_value for accession_number_value in accession_number_list_with_none]
        accession_number_list = []
        for accession_number_none in accession_number_list_with_none_list:
            accession_number_list.append(accession_number_none)
        print(accession_number_list)
        print("")
        # Writing Accession Number
        print("Step 2: Writing Accession Number, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        accession_number_string = '\n'.join(str(accession_number_value) for accession_number_value in accession_number_list)
        print(accession_number_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_accession_number.txt"), "w")
        handle_outputfile.write(accession_number_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 1: Accession Number Extract Complete")
        print("")
        break

# Annotations Accessions Extract ##################################################
def annotationsaccessions():
    while True:
        # Module Name
        print("Module 2: Annotations Accessions Extract")
        print("")
        # Parsing Annotations Accessions
        print("Step 1: Parsing Annotations Accessions, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_accessions_list_with_none = []
        for record in handle_parse_db:
            annotations_accessions_list_with_none.append(record.annotations["accessions"])
        annotations_accessions_list_with_none_list = [['None'] if annotations_accessions_value is None else annotations_accessions_value for annotations_accessions_value in annotations_accessions_list_with_none]
        annotations_accessions_list = []
        for annotations_accessions_none in annotations_accessions_list_with_none_list:
            annotations_accessions_list.append(annotations_accessions_none)
        print(annotations_accessions_list)
        print("")
        # Writing Annotations Accessions
        print("Step 2: Writing Annotations Accessions, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_accessions_string = '\n'.join(str(annotations_accessions_value) for annotations_accessions_value in annotations_accessions_list)
        print(annotations_accessions_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_accessions.txt"), "w")
        handle_outputfile.write(annotations_accessions_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 2: Annotations Accessions Extract Complete")
        print("")
        break

# Annotations Authors Extract ##################################################
def annotationsauthors():
    while True:
        # Module Name
        print("Module 3: Annotations Authors Extract")
        print("")
        # Parsing Annotations Authors
        print("Step 1: Parsing Annotations Authors, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_authors_list_with_none = []
        for record in handle_parse_db:
            annotations_authors_list_with_none.append(record.annotations["references"][0].authors)
        annotations_authors_list_with_none_list = [['None'] if annotations_authors_value is None else annotations_authors_value for annotations_authors_value in annotations_authors_list_with_none]
        annotations_authors_list = []
        for annotations_authors_none in annotations_authors_list_with_none_list:
            annotations_authors_list.append(annotations_authors_none)
        print(annotations_authors_list)
        print("")
        # Writing Annotations Authors
        print("Step 2: Writing Annotations Authors, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_authors_string = '\n'.join(str(annotations_authors_value) for annotations_authors_value in annotations_authors_list)
        print(annotations_authors_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_authors.txt"), "w")
        handle_outputfile.write(annotations_authors_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 3: Annotations Authors Extract Complete")
        print("")
        break

# Annotations Comment Extract ##################################################
def annotationscomment():
    while True:
        # Module Name
        print("Module 4: Annotations Comment Extract")
        print("")
        # Parsing Annotations Comment
        print("Step 1: Parsing Annotations Comment, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_comment_list_with_none = []
        for record in handle_parse_db:
            annotations_comment_list_with_none.append(record.annotations["references"][0].comment)
        annotations_comment_list_with_none_list = [['None'] if annotations_comment_value is None else annotations_comment_value for annotations_comment_value in annotations_comment_list_with_none]
        annotations_comment_list = []
        for annotations_comment_none in annotations_comment_list_with_none_list:
            annotations_comment_list.append(annotations_comment_none)
        print(annotations_comment_list)
        print("")
        # Writing Annotations Comment
        print("Step 2: Writing Annotations Comment, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_comment_string = '\n'.join(str(annotations_comment_value) for annotations_comment_value in annotations_comment_list)
        print(annotations_comment_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_comment.txt"), "w")
        handle_outputfile.write(annotations_comment_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 4: Annotations Comment Extract Complete")
        print("")
        break

# Annotations Consrtm Extract ##################################################
def annotationsconsrtm():
    while True:
        # Module Name
        print("Module 5: Annotations Consrtm Extract")
        print("")
        # Parsing Annotations Consrtm
        print("Step 1: Parsing Annotations Consrtm, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_consrtm_list_with_none = []
        for record in handle_parse_db:
            annotations_consrtm_list_with_none.append(record.annotations["references"][0].consrtm)
        annotations_consrtm_list_with_none_list = [['None'] if annotations_consrtm_value is None else annotations_consrtm_value for annotations_consrtm_value in annotations_consrtm_list_with_none]
        annotations_consrtm_list = []
        for annotations_consrtm_none in annotations_consrtm_list_with_none_list:
            annotations_consrtm_list.append(annotations_consrtm_none)
        print(annotations_consrtm_list)
        print("")
        # Writing Annotations Consrtm
        print("Step 2: Writing Annotations Consrtm, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_consrtm_string = '\n'.join(str(annotations_consrtm_value) for annotations_consrtm_value in annotations_consrtm_list)
        print(annotations_consrtm_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_consrtm.txt"), "w")
        handle_outputfile.write(annotations_consrtm_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 5: Annotations Consrtm Extract Complete")
        print("")
        break

# Annotations Date Extract ##################################################
def annotationsdate():
    while True:
        # Module Name
        print("Module 6: Annotations Date Extract")
        print("")
        # Parsing Annotations Date
        print("Step 1: Parsing Annotations Date, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_date_list_with_none = []
        for record in handle_parse_db:
            annotations_date_list_with_none.append(record.annotations["date"])
        annotations_date_list_with_none_list = [['None'] if annotations_date_value is None else annotations_date_value for annotations_date_value in annotations_date_list_with_none]
        annotations_date_list = []
        for annotations_date_none in annotations_date_list_with_none_list:
            annotations_date_list.append(annotations_date_none)
        print(annotations_date_list)
        print("")
        # Writing Annotations Date
        print("Step 2: Writing Annotations Date, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_date_string = '\n'.join(str(annotations_date_value) for annotations_date_value in annotations_date_list)
        print(annotations_date_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_date.txt"), "w")
        handle_outputfile.write(annotations_date_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 6: Annotations Date Extract Complete")
        print("")
        break

# Annotations Data File Division Extract ##################################################
def annotationsdatafiledivision():
    while True:
        # Module Name
        print("Module 7: Annotations Data File Division Extract")
        print("")
        # Parsing Annotations Data File Division
        print("Step 1: Parsing Annotations Data File Division, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_data_file_division_list_with_none = []
        for record in handle_parse_db:
            annotations_data_file_division_list_with_none.append(record.annotations["data_file_division"])
        annotations_data_file_division_list_with_none_list = [['None'] if annotations_data_file_division_value is None else annotations_data_file_division_value for annotations_data_file_division_value in annotations_data_file_division_list_with_none]
        annotations_data_file_division_list = []
        for annotations_data_file_division_none in annotations_data_file_division_list_with_none_list:
            annotations_data_file_division_list.append(annotations_data_file_division_none)
        print(annotations_data_file_division_list)
        print("")
        # Writing Annotations Data File Division
        print("Step 2: Writing Annotations Data File Division, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_data_file_division_string = '\n'.join(str(annotations_data_file_division_value) for annotations_data_file_division_value in annotations_data_file_division_list)
        print(annotations_data_file_division_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_data_file_division.txt"), "w")
        handle_outputfile.write(annotations_data_file_division_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 7: Annotations Data File Division Extract Complete")
        print("")
        break

# Annotations Evidence Extract ##################################################
def annotationsevidence():
    while True:
        # Module Name
        print("Module 8: Annotations Evidence Extract")
        print("")
        # Parsing Annotations Evidence
        print("Step 1: Parsing Annotations Evidence, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_evidence_list_with_none = []
        for record in handle_parse_db:
            annotations_evidence_list_with_none.append(record.annotations["evidence"])
        annotations_evidence_list_with_none_list = [['None'] if annotations_evidence_value is None else annotations_evidence_value for annotations_evidence_value in annotations_evidence_list_with_none]
        annotations_evidence_list = []
        for annotations_evidence_none in annotations_evidence_list_with_none_list:
            annotations_evidence_list.append(annotations_evidence_none)
        print(annotations_evidence_list)
        print("")
        # Writing Annotations Evidence
        print("Step 2: Writing Annotations Evidence, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_evidence_string = '\n'.join(str(annotations_evidence_value) for annotations_evidence_value in annotations_evidence_list)
        print(annotations_evidence_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_evidence.txt"), "w")
        handle_outputfile.write(annotations_evidence_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 8: Annotations Evidence Extract Complete")
        print("")
        break

# Annotations GI Extract ##################################################
def annotationsgi():
    while True:
        # Module Name
        print("Module 9: Annotations GI Extract")
        print("")
        # Parsing Annotations GI
        print("Step 1: Parsing Annotations GI, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_gi_list_with_none = []
        for record in handle_parse_db:
            annotations_gi_list_with_none.append(record.annotations["gi"])
        annotations_gi_list_with_none_list = [['None'] if annotations_gi_value is None else annotations_gi_value for annotations_gi_value in annotations_gi_list_with_none]
        annotations_gi_list = []
        for annotations_gi_none in annotations_gi_list_with_none_list:
            annotations_gi_list.append(annotations_gi_none)
        print(annotations_gi_list)
        print("")
        # Writing Annotations GI
        print("Step 2: Writing Annotations GI, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_gi_string = '\n'.join(str(annotations_gi_value) for annotations_gi_value in annotations_gi_list)
        print(annotations_gi_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_gi.txt"), "w")
        handle_outputfile.write(annotations_gi_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 9: Annotations GI Extract Complete")
        print("")
        break

# Annotations Journal Extract ##################################################
def annotationsjournal():
    while True:
        # Module Name
        print("Module 10: Annotations Journal Extract")
        print("")
        # Parsing Annotations Journal
        print("Step 1: Parsing Annotations Journal, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_journal_list_with_none = []
        for record in handle_parse_db:
            annotations_journal_list_with_none.append(record.annotations["references"][0].journal)
        annotations_journal_list_with_none_list = [['None'] if annotations_journal_value is None else annotations_journal_value for annotations_journal_value in annotations_journal_list_with_none]
        annotations_journal_list = []
        for annotations_journal_none in annotations_journal_list_with_none_list:
            annotations_journal_list.append(annotations_journal_none)
        print(annotations_journal_list)
        print("")
        # Writing Annotations Journal
        print("Step 2: Writing Annotations Journal, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_journal_string = '\n'.join(str(annotations_journal_value) for annotations_journal_value in annotations_journal_list)
        print(annotations_journal_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_journal.txt"), "w")
        handle_outputfile.write(annotations_journal_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 10: Annotations Journal Extract Complete")
        print("")
        break

# Annotations Keywords Extract ##################################################
def annotationskeywords():
    while True:
        # Module Name
        print("Module 11: Annotations Keywords Extract")
        print("")
        # Parsing Annotations Keywords
        print("Step 1: Parsing Annotations Keywords, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_keywords_list_with_none = []
        for record in handle_parse_db:
            annotations_keywords_list_with_none.append(record.annotations["keywords"])
        annotations_keywords_list_with_none_list = [['None'] if annotations_keywords_value is None else annotations_keywords_value for annotations_keywords_value in annotations_keywords_list_with_none]
        annotations_keywords_list = []
        for annotations_keywords_none in annotations_keywords_list_with_none_list:
            annotations_keywords_list.append(annotations_keywords_none)
        print(annotations_keywords_list)
        print("")
        # Writing Annotations Keywords
        print("Step 2: Writing Annotations Keywords, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_keywords_string = '\n'.join(str(annotations_keywords_value) for annotations_keywords_value in annotations_keywords_list)
        print(annotations_keywords_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_keywords.txt"), "w")
        handle_outputfile.write(annotations_keywords_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 11: Annotations Keywords Extract Complete")
        print("")
        break

# Annotations Location Extract ##################################################
def annotationslocation():
    while True:
        # Module Name
        print("Module 12: Annotations Location Extract")
        print("")
        # Parsing Annotations Location
        print("Step 1: Parsing Annotations Location, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_location_list_with_none = []
        for record in handle_parse_db:
            annotations_location_list_with_none.append(record.annotations["references"][0].location)
        annotations_location_list_with_none_list = [['None'] if annotations_location_value is None else annotations_location_value for annotations_location_value in annotations_location_list_with_none]
        annotations_location_list = []
        for annotations_location_none in annotations_location_list_with_none_list:
            annotations_location_list.append(annotations_location_none)
        print(annotations_location_list)
        print("")
        # Writing Annotations Location
        print("Step 2: Writing Annotations Location, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_location_string = '\n'.join(str(annotations_location_value) for annotations_location_value in annotations_location_list)
        print(annotations_location_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_location.txt"), "w")
        handle_outputfile.write(annotations_location_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 12: Annotations Location Extract Complete")
        print("")
        break

# Annotations Medline ID Extract ##################################################
def annotationsmedlineid():
    while True:
        # Module Name
        print("Module 13: Annotations Medline ID Extract")
        print("")
        # Parsing Annotations Medline ID
        print("Step 1: Parsing Annotations Medline ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_medlineid_list_with_none = []
        for record in handle_parse_db:
            annotations_medlineid_list_with_none.append(record.annotations["references"][0].medline_id)
        annotations_medlineid_list_with_none_list = [['None'] if annotations_medlineid_value is None else annotations_medlineid_value for annotations_medlineid_value in annotations_medlineid_list_with_none]
        annotations_medlineid_list = []
        for annotations_medlineid_none in annotations_medlineid_list_with_none_list:
            annotations_medlineid_list.append(annotations_medlineid_none)
        print(annotations_medlineid_list)
        print("")
        # Writing Annotations Medline ID
        print("Step 2: Writing Annotations Medline ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_medlineid_string = '\n'.join(str(annotations_medlineid_value) for annotations_medlineid_value in annotations_medlineid_list)
        print(annotations_medlineid_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_medlineid.txt"), "w")
        handle_outputfile.write(annotations_medlineid_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 13: Annotations Medline ID Extract Complete")
        print("")
        break
# Annotations Molecule Type Extract ##################################################
def annotationsmoleculetype():
    while True:
        # Module Name
        print("Module 14: Annotations Molecule Type Extract")
        print("")
        # Parsing Annotations Molecule Type
        print("Step 1: Parsing Annotations Molecule Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_molecule_type_list_with_none = []
        for record in handle_parse_db:
            annotations_molecule_type_list_with_none.append(record.annotations["molecule_type"])
        annotations_molecule_type_list_with_none_list = [['None'] if annotations_molecule_type_value is None else annotations_molecule_type_value for annotations_molecule_type_value in annotations_molecule_type_list_with_none]
        annotations_molecule_type_list = []
        for annotations_molecule_type_none in annotations_molecule_type_list_with_none_list:
            annotations_molecule_type_list.append(annotations_molecule_type_none)
        print(annotations_molecule_type_list)
        print("")
        # Writing Annotations Molecule Type
        print("Step 2: Writing Annotations Molecule Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_molecule_type_string = '\n'.join(str(annotations_molecule_type_value) for annotations_molecule_type_value in annotations_molecule_type_list)
        print(annotations_molecule_type_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_molecule_type.txt"), "w")
        handle_outputfile.write(annotations_molecule_type_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 14: Annotations Molecule Type Extract Complete")
        print("")
        break

# Annotations Organism Extract ##################################################
def annotationsorganism():
    while True:
        # Module Name
        print("Module 15: Annotations Organism Extract")
        print("")
        # Parsing Annotations Organism
        print("Step 1: Parsing Annotations Organism, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_organism_list_with_none = []
        for record in handle_parse_db:
            annotations_organism_list_with_none.append(record.annotations["organism"])
        annotations_organism_list_with_none_list = [['None'] if annotations_organism_value is None else annotations_organism_value for annotations_organism_value in annotations_organism_list_with_none]
        annotations_organism_list = []
        for annotations_organism_none in annotations_organism_list_with_none_list:
            annotations_organism_list.append(annotations_organism_none)
        print(annotations_organism_list)
        print("")
        # Writing Annotations Organism
        print("Step 2: Writing Annotations Organism, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_organism_string = '\n'.join(str(annotations_organism_value) for annotations_organism_value in annotations_organism_list)
        print(annotations_organism_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_organism.txt"), "w")
        handle_outputfile.write(annotations_organism_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 15: Annotations Organism Extract Complete")
        print("")
        break

# Annotations References Lenght Extract ##################################################
def annotationsreferenceslength():
    while True:
        # Module Name
        print("Module 16: Annotations References Length Extract")
        print("")
        # Parsing Annotations References Length
        print("Step 1: Parsing Annotations References Length, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        #accession_list_with_none = []
        annotations_references_length_list_with_none = []
        for record in handle_parse_db:
        #    accession_list_with_none.append(record.id)
            annotations_references_length_list_with_none.append(len(record.annotations["references"]))
        annotations_references_length_list_with_none_list = [['None'] if annotations_references_length_value is None else annotations_references_length_value for annotations_references_length_value in annotations_references_length_list_with_none]
        #annotations_references_length_list = []
        annotations_references_length_list_with_none_list_none = []
        for annotations_references_length_none in annotations_references_length_list_with_none_list:
            annotations_references_length_list_with_none_list_none.append(annotations_references_length_none)
        #annotations_references_length_zip = zip(accession_list_with_none,annotations_references_length_list_with_none_list_none)
        #print(annotations_references_length_zip)
        print(annotations_references_length_list_with_none_list_none)
        #annotations_references_length_list = []
        #for accession_annotations_references_length in annotations_references_length_zip:
        #    annotations_references_length_list.append(accession_annotations_references_length)
        print("")
        # Writing Annotations References Length
        print("Step 2: Writing Annotations References Length, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_references_length_string = '\n'.join(str(annotations_references_length_value) for annotations_references_length_value in annotations_references_length_list_with_none_list_none)
        print(annotations_references_length_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_references_length.txt"), "w")
        handle_outputfile.write(str(annotations_references_length_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 16: Annotations References Length Extract Complete")
        print("")
        break

# Annotations Ref.[0] Authors Extract ##################################################
def annotationsref0authors():
    while True:
        # Module Name
        print("Module 17: Annotations Ref.[0] Authors Extract")
        print("")
        # Parsing Annotations Ref.[0] Authors
        print("Step 1: Parsing Annotations Ref.[0] Authors, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_ref0_authors_list_with_none = []
        for record in handle_parse_db:
            annotations_ref0_authors_list_with_none.append(record.annotations["references"][0].authors)
        annotations_ref0_authors_list_with_none_list = [['None'] if annotations_ref0_authors_value is None else annotations_ref0_authors_value for annotations_ref0_authors_value in annotations_ref0_authors_list_with_none]
        annotations_ref0_authors_list = []
        for annotations_ref0_authors_none in annotations_ref0_authors_list_with_none_list:
            annotations_ref0_authors_list.append(annotations_ref0_authors_none)
        print(annotations_ref0_authors_list)
        print("")
        # Writing Annotations Ref.[0] Authors
        print("Step 2: Writing Annotations Ref.[0] Authors, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_ref0_authors_string = '\n'.join(str(annotations_ref0_authors_value) for annotations_ref0_authors_value in annotations_ref0_authors_list)
        print(annotations_ref0_authors_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_ref0_authors.txt"), "w")
        handle_outputfile.write(annotations_ref0_authors_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 17: Annotations Ref.[0] Authors Extract Complete")
        print("")
        break

# Annotations Ref.[0] Title Extract ##################################################
def annotationsref0title():
    while True:
        # Module Name
        print("Module 18: Annotations Ref.[0] Title Extract")
        print("")
        # Parsing Annotations Ref.[0] Title
        print("Step 1: Parsing Annotations Ref.[0] Title, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_ref0_title_list_with_none = []
        for record in handle_parse_db:
            annotations_ref0_title_list_with_none.append(record.annotations["references"][0].title)
        annotations_ref0_title_list_with_none_list = [['None'] if annotations_ref0_title_value is None else annotations_ref0_title_value for annotations_ref0_title_value in annotations_ref0_title_list_with_none]
        annotations_ref0_title_list = []
        for annotations_ref0_title_none in annotations_ref0_title_list_with_none_list:
            annotations_ref0_title_list.append(annotations_ref0_title_none)
        print(annotations_ref0_title_list)
        print("")
        # Writing Annotations Ref.[0] Title
        print("Step 2: Writing Annotations Ref.[0] Title, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_ref0_title_string = '\n'.join(str(annotations_ref0_title_value) for annotations_ref0_title_value in annotations_ref0_title_list)
        print(annotations_ref0_title_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_ref0_title.txt"), "w")
        handle_outputfile.write(annotations_ref0_title_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 18: Annotations Ref.[0] Title Extract Complete")
        print("")
        break

# Annotations Ref.[0] Journal Extract ##################################################
def annotationsref0journal():
    while True:
        # Module Name
        print("Module 19: Annotations Ref.[0] Journal Extract")
        print("")
        # Parsing Annotations Ref.[0] Journal
        print("Step 1: Parsing Annotations Ref.[0] Journal, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_ref0_journal_list_with_none = []
        for record in handle_parse_db:
            annotations_ref0_journal_list_with_none.append(record.annotations["references"][0].journal)
        annotations_ref0_journal_list_with_none_list = [['None'] if annotations_ref0_journal_value is None else annotations_ref0_journal_value for annotations_ref0_journal_value in annotations_ref0_journal_list_with_none]
        annotations_ref0_journal_list = []
        for annotations_ref0_journal_none in annotations_ref0_journal_list_with_none_list:
            annotations_ref0_journal_list.append(annotations_ref0_journal_none)
        print(annotations_ref0_journal_list)
        print("")
        # Writing Annotations Ref.[0] Journal
        print("Step 2: Writing Annotations Ref.[0] Journal, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_ref0_journal_string = '\n'.join(str(annotations_ref0_journal_value) for annotations_ref0_journal_value in annotations_ref0_journal_list)
        print(annotations_ref0_journal_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_ref0_journal.txt"), "w")
        handle_outputfile.write(annotations_ref0_journal_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 19: Annotations Ref.[0] Journal Extract Complete")
        print("")
        break

# Annotations Ref.[0] Pubmed ID Extract ##################################################
def annotationsref0pubmedid():
    while True:
        # Module Name
        print("Module 20: Annotations Ref.[0] Pubmed ID Extract")
        print("")
        # Parsing Annotations Ref.[0] Pubmed ID
        print("Step 1: Parsing Annotations Ref.[0] Pubmed ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_ref0_pubmedid_list_with_none = []
        for record in handle_parse_db:
            annotations_ref0_pubmedid_list_with_none.append(record.annotations["references"][0].pubmed_id)
        annotations_ref0_pubmedid_list_with_none_list = [['None'] if annotations_ref0_pubmedid_value is None else annotations_ref0_pubmedid_value for annotations_ref0_pubmedid_value in annotations_ref0_pubmedid_list_with_none]
        annotations_ref0_pubmedid_list = []
        for annotations_ref0_pubmedid_none in annotations_ref0_pubmedid_list_with_none_list:
            annotations_ref0_pubmedid_list.append(annotations_ref0_pubmedid_none)
        print(annotations_ref0_pubmedid_list)
        print("")
        # Writing Annotations Ref.[0] Pubmed ID
        print("Step 2: Writing Annotations Ref.[0] Pubmed ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_ref0_pubmedid_string = '\n'.join(str(annotations_ref0_pubmedid_value) for annotations_ref0_pubmedid_value in annotations_ref0_pubmedid_list)
        print(annotations_ref0_pubmedid_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_ref0_pubmedid.txt"), "w")
        handle_outputfile.write(annotations_ref0_pubmedid_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 20: Annotations Ref.[0] Pubmed ID Extract Complete")
        print("")
        break

# Annotations Ref.[All] Authors Extract ##################################################
def annotationsrefallauthors():
    while True:
        # Module Name
        print("Module 21: Annotations Ref.[All] Authors Extract")
        print("")
        # Parsing Annotations Ref.[all] Authors
        print("Step 1: Parsing Annotations Ref.[All] Authors, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        annotations_refall_authors_list_with_none = []
        for record in handle_parse_db:
            #annotations_refall_authors_list_with_none.append(record.annotations["references"][all].authors)
            #print(type(record.annotations["references"]))
            #print(len(record.annotations["references"]))
            for refall in record.annotations["references"]:
                #print(refall.journal)
                accession_list_with_none.append(record.id)
                annotations_refall_authors_list_with_none.append(refall.authors)
        #print(annotations_refall_authors_list_with_none)
        annotations_refall_authors_list_with_none_list = [['None'] if annotations_refall_authors_value is None else annotations_refall_authors_value for annotations_refall_authors_value in annotations_refall_authors_list_with_none]
        annotations_refall_authors_list = []
        for annotations_refall_authors_none in annotations_refall_authors_list_with_none_list:
            annotations_refall_authors_list.append(annotations_refall_authors_none)
        journal_zip = zip(accession_list_with_none,annotations_refall_authors_list)
        print(journal_zip)
        journal_list = []
        for accession_journal in journal_zip:
            journal_list.append("\t".join(accession_journal))
        #print(annotations_refall_authors_list)
        print("")
        # Writing Annotations Ref.[All] Authors
        print("Step 2: Writing Annotations Ref.[All] Authors, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #annotations_refall_authors_string = '\n'.join(str(annotations_refall_authors_value) for annotations_refall_authors_value in annotations_refall_authors_list)
        annotations_refall_authors_string = '\n'.join(str(accession_journal_value) for accession_journal_value in journal_list)
        print(annotations_refall_authors_string)
        handle_outputfile = open(os.path.expanduser("~/export/description_annotations_refall_authors.txt"), "w")
        handle_outputfile.write(annotations_refall_authors_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 21: Annotations Ref.[All] Authors Extract Complete")
        print("")
        break

# Annotations Ref.[All] Title Extract ##################################################
def annotationsrefalltitle():
    while True:
        # Module Name
        print("Module 22: Annotations Ref.[All] Title Extract")
        print("")
        # Parsing Annotations Ref.[all] Title
        print("Step 1: Parsing Annotations Ref.[All] Title, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        annotations_refall_title_list_with_none = []
        for record in handle_parse_db:
            #annotations_refall_title_list_with_none.append(record.annotations["references"][all].title)
            #print(type(record.annotations["references"]))
            #print(len(record.annotations["references"]))
            for refall in record.annotations["references"]:
                #print(refall.title)
                accession_list_with_none.append(record.id)
                annotations_refall_title_list_with_none.append(refall.title)
        #print(annotations_refall_title_list_with_none)
        annotations_refall_title_list_with_none_list = [['None'] if annotations_refall_title_value is None else annotations_refall_title_value for annotations_refall_title_value in annotations_refall_title_list_with_none]
        annotations_refall_title_list = []
        for annotations_refall_title_none in annotations_refall_title_list_with_none_list:
            annotations_refall_title_list.append(annotations_refall_title_none)
        title_zip = zip(accession_list_with_none,annotations_refall_title_list)
        print(title_zip)
        title_list = []
        for accession_title in title_zip:
            title_list.append("\t".join(accession_title))
        #print(annotations_refall_title_list)
        print("")
        # Writing Annotations Ref.[All] Title
        print("Step 2: Writing Annotations Ref.[All] Title, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #annotations_refall_title_string = '\n'.join(str(annotations_refall_title_value) for annotations_refall_title_value in annotations_refall_title_list)
        annotations_refall_title_string = '\n'.join(str(accession_title_value) for accession_title_value in title_list)
        print(annotations_refall_title_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_refall_title.txt"), "w")
        handle_outputfile.write(annotations_refall_title_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 22: Annotations Ref.[All] Title Extract Complete")
        print("")
        break

# Annotations Ref.[All] Journal Extract ##################################################
def annotationsrefalljournal():
    while True:
        # Module Name
        print("Module 23: Annotations Ref.[All] Journal Extract")
        print("")
        # Parsing Annotations Ref.[All] Journal
        print("Step 1: Parsing Annotations Ref.[All] Journal, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        annotations_refall_journal_list_with_none = []
        for record in handle_parse_db:
            #annotations_refall_journal_list_with_none.append(record.annotations["references"][all].journal)
            #print(type(record.annotations["references"]))
            #print(len(record.annotations["references"]))
            for refall in record.annotations["references"]:
                #print(refall.journal)
                accession_list_with_none.append(record.id)
                annotations_refall_journal_list_with_none.append(refall.journal)
        #print(annotations_refall_journal_list_with_none)
        annotations_refall_journal_list_with_none_list = [['None'] if annotations_refall_journal_value is None else annotations_refall_journal_value for annotations_refall_journal_value in annotations_refall_journal_list_with_none]
        annotations_refall_journal_list = []
        for annotations_refall_journal_none in annotations_refall_journal_list_with_none_list:
            annotations_refall_journal_list.append(annotations_refall_journal_none)
        journal_zip = zip(accession_list_with_none,annotations_refall_journal_list)
        print(journal_zip)
        journal_list = []
        for accession_journal in journal_zip:
            journal_list.append("\t".join(accession_journal))
        #print(annotations_refall_journal_list)
        print("")
        # Writing Annotations Ref.[All] Journal
        print("Step 2: Writing Annotations Ref.[All] Journal, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #annotations_refall_journal_string = '\n'.join(str(annotations_refall_journal_value) for annotations_refall_journal_value in annotations_refall_journal_list)
        annotations_refall_journal_string = '\n'.join(str(accession_journal_value) for accession_journal_value in journal_list)
        print(annotations_refall_journal_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_refall_journal.txt"), "w")
        handle_outputfile.write(annotations_refall_journal_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 23: Annotations Ref.[All] Journal Extract Complete")
        print("")
        break

# Annotations Ref.[All] Pubmed ID Extract ##################################################
def annotationsrefallpubmedid():
    while True:
        # Module Name
        print("Module 24: Annotations Ref.[All] Pubmed ID Extract")
        print("")
        # Parsing Annotations Ref.[All] Pubmed ID
        print("Step 1: Parsing Annotations Ref.[All] Pubmed ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        annotations_refall_pubmedid_list_with_none = []
        for record in handle_parse_db:
            #annotations_refall_pubmedid_list_with_none.append(record.annotations["references"][all].pubmed_id)
            #print(type(record.annotations["references"]))
            #print(len(record.annotations["references"]))
            for refall in record.annotations["references"]:
                #print(refall.pubmed_id)
                accession_list_with_none.append(record.id)
                annotations_refall_pubmedid_list_with_none.append(refall.pubmed_id)
        #print(annotations_refall_pubmedid_list_with_none)
        annotations_refall_pubmedid_list_with_none_list = [['None'] if annotations_refall_pubmedid_value is None else annotations_refall_pubmedid_value for annotations_refall_pubmedid_value in annotations_refall_pubmedid_list_with_none]
        annotations_refall_pubmedid_list = []
        for annotations_refall_pubmedid_none in annotations_refall_pubmedid_list_with_none_list:
            annotations_refall_pubmedid_list.append(annotations_refall_pubmedid_none)
        pubmedid_zip = zip(accession_list_with_none,annotations_refall_pubmedid_list)
        print(pubmedid_zip)
        pubmedid_list = []
        for accession_pubmedid in pubmedid_zip:
            pubmedid_list.append("\t".join(accession_pubmedid))
        #print(annotations_refall_pubmedid_list)
        print("")
        # Writing Annotations Ref.[All] Pubmed ID
        print("Step 2: Writing Annotations Ref.[All] Pubmed ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #annotations_refall_pubmedid_string = '\n'.join(str(annotations_refall_pubmedid_value) for annotations_refall_pubmedid_value in annotations_refall_pubmedid_list)
        annotations_refall_pubmedid_string = '\n'.join(str(accession_pubmedid_value) for accession_pubmedid_value in pubmedid_list)
        print(annotations_refall_pubmedid_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_refall_pubmedid.txt"), "w")
        handle_outputfile.write(annotations_refall_pubmedid_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 24: Annotations Ref.[All] Pubmed ID Extract Complete")
        print("")
        break

# Annotations Sequence Version Extract ##################################################
def annotationssequenceversion():
    while True:
        # Module Name
        print("Module 25: Annotations Sequence Version Extract")
        print("")
        # Parsing Annotations Sequence Version
        print("Step 1: Parsing Annotations Sequence Version, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_sequence_version_list_with_none = []
        for record in handle_parse_db:
            annotations_sequence_version_list_with_none.append(record.annotations["sequence_version"])
        annotations_sequence_version_list_with_none_list = [['None'] if annotations_sequence_version_value is None else annotations_sequence_version_value for annotations_sequence_version_value in annotations_sequence_version_list_with_none]
        annotations_sequence_version_list = []
        for annotations_sequence_version_none in annotations_sequence_version_list_with_none_list:
            annotations_sequence_version_list.append(annotations_sequence_version_none)
        print(annotations_sequence_version_list)
        print("")
        # Writing Annotations Sequence Version
        print("Step 2: Writing Annotations Sequence Version, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_sequence_version_string = '\n'.join(str(annotations_sequence_version_value) for annotations_sequence_version_value in annotations_sequence_version_list)
        print(annotations_sequence_version_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_sequence_version.txt"), "w")
        handle_outputfile.write(annotations_sequence_version_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 25: Annotations Sequence Version Extract Complete")
        print("")
        break

# Annotations Source Extract ##################################################
def annotationssource():
    while True:
        # Module Name
        print("Module 26: Annotations Source Extract")
        print("")
        # Parsing Annotations Source
        print("Step 1: Parsing Annotations Source, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_source_list_with_none = []
        for record in handle_parse_db:
            annotations_source_list_with_none.append(record.annotations["source"])
        annotations_source_list_with_none_list = [['None'] if annotations_source_value is None else annotations_source_value for annotations_source_value in annotations_source_list_with_none]
        annotations_source_list = []
        for annotations_source_none in annotations_source_list_with_none_list:
            annotations_source_list.append(annotations_source_none)
        print(annotations_source_list)
        print("")
        # Writing Annotations Source
        print("Step 2: Writing Annotations Source, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_source_string = '\n'.join(str(annotations_source_value) for annotations_source_value in annotations_source_list)
        print(annotations_source_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_source.txt"), "w")
        handle_outputfile.write(annotations_source_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 26: Annotations Source Extract Complete")
        print("")
        break

# Annotations Taxonomy Extract ##################################################
def annotationstaxonomy():
    while True:
        # Module Name
        print("Module 27: Annotations Taxonomy Extract")
        print("")
        # Parsing Annotations Taxonomy
        print("Step 1: Parsing Annotations Taxonomy, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_taxonomy_list_with_none = []
        for record in handle_parse_db:
            annotations_taxonomy_list_with_none.append(record.annotations["taxonomy"])
        annotations_taxonomy_list_with_none_list = [['None'] if annotations_taxonomy_value is None else annotations_taxonomy_value for annotations_taxonomy_value in annotations_taxonomy_list_with_none]
        annotations_taxonomy_list = []
        for annotations_taxonomy_none in annotations_taxonomy_list_with_none_list:
            annotations_taxonomy_list.append(annotations_taxonomy_none)
        print(annotations_taxonomy_list)
        print("")
        # Writing Annotations Taxonomy
        print("Step 2: Writing Annotations Taxonomy, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_taxonomy_string = '\n'.join(str(annotations_taxonomy_value) for annotations_taxonomy_value in annotations_taxonomy_list)
        print(annotations_taxonomy_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_taxonomy.txt"), "w")
        handle_outputfile.write(annotations_taxonomy_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 27: Annotations Taxonomy Extract Complete")
        print("")
        break

# Annotations Title Extract ##################################################
def annotationstitle():
    while True:
        # Module Name
        print("Module 28: Annotations Title Extract")
        print("")
        # Parsing Annotations Title
        print("Step 1: Parsing Annotations Title, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_title_list_with_none = []
        for record in handle_parse_db:
            annotations_title_list_with_none.append(record.annotations["references"][0].title)
        annotations_title_list_with_none_list = [['None'] if annotations_title_value is None else annotations_title_value for annotations_title_value in annotations_title_list_with_none]
        annotations_title_list = []
        for annotations_title_none in annotations_title_list_with_none_list:
            annotations_title_list.append(annotations_title_none)
        print(annotations_title_list)
        print("")
        # Writing Annotations Title
        print("Step 2: Writing Annotations Title, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_title_string = '\n'.join(str(annotations_title_value) for annotations_title_value in annotations_title_list)
        print(annotations_title_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_title.txt"), "w")
        handle_outputfile.write(annotations_title_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 28: Annotations Title Extract Complete")
        print("")
        break

# Annotations Topology Extract ##################################################
def annotationstopology():
    while True:
        # Module Name
        print("Module 29: Annotations Topology Extract")
        print("")
        # Parsing Annotations Topology
        print("Step 1: Parsing Annotations Topology, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_topology_list_with_none = []
        for record in handle_parse_db:
            annotations_topology_list_with_none.append(record.annotations["topology"])
        annotations_topology_list_with_none_list = [['None'] if annotations_topology_value is None else annotations_topology_value for annotations_topology_value in annotations_topology_list_with_none]
        annotations_topology_list = []
        for annotations_topology_none in annotations_topology_list_with_none_list:
            annotations_topology_list.append(annotations_topology_none)
        print(annotations_topology_list)
        print("")
        # Writing Annotations Topology
        print("Step 2: Writing Annotations Topology, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_topology_string = '\n'.join(str(annotations_topology_value) for annotations_topology_value in annotations_topology_list)
        print(annotations_topology_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_topology.txt"), "w")
        handle_outputfile.write(annotations_topology_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 29: Annotations Topology Extract Complete")
        print("")
        break

# Annotations Keys Extract ##################################################
def annotationskeys():
    while True:
        # Module Name
        print("Module 30: Annotations Keys Extract")
        print("")
        # Parsing Annotations Keys
        print("Step 1: Parsing Annotations Keys, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_keys_list_with_none = []
        for record in handle_parse_db:
            annotations_keys_list_with_none.append(record.annotations.keys)
        annotations_keys_list_with_none_list = [['None'] if annotations_keys_value is None else annotations_keys_value for annotations_keys_value in annotations_keys_list_with_none]
        annotations_keys_list = []
        for annotations_keys_none in annotations_keys_list_with_none_list:
            annotations_keys_list.append(annotations_keys_none)
        print(annotations_keys_list)
        print("")
        # Writing Annotations Keys
        print("Step 2: Writing Annotations Keys, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_keys_string = '\n'.join(str(annotations_keys_value) for annotations_keys_value in annotations_keys_list)
        print(annotations_keys_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_keys.txt"), "w")
        handle_outputfile.write(annotations_keys_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 30: Annotations Keys Extract Complete")
        print("")
        break

# Annotations Values Extract ##################################################
def annotationsvalues():
    while True:
        # Module Name
        print("Module 31: Annotations Values Extract")
        print("")
        # Parsing Annotations Values
        print("Step 1: Parsing Annotations Values, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_values_list_with_none = []
        for record in handle_parse_db:
            annotations_values_list_with_none.append(record.annotations.values)
        annotations_values_list_with_none_list = [['None'] if annotations_values_value is None else annotations_values_value for annotations_values_value in annotations_values_list_with_none]
        annotations_values_list = []
        for annotations_values_none in annotations_values_list_with_none_list:
            annotations_values_list.append(annotations_values_none)
        print(annotations_values_list)
        print("")
        # Writing Annotations Values
        print("Step 2: Writing Annotations Values, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        annotations_values_string = '\n'.join(str(annotations_values_value) for annotations_values_value in annotations_values_list)
        print(annotations_values_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_annotations_values.txt"), "w")
        handle_outputfile.write(annotations_values_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 31: Annotations Values Extract Complete")
        print("")
        break

# DB XREFS Extract ##################################################
def dbxrefs():
    while True:
        # Module Name
        print("Module 32: DB XREFS Extract")
        print("")
        # Parsing DB XREFS
        print("Step 1: Parsing DB XREFS, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        db_xrefs_list_with_none = []
        for record in handle_parse_db:
            db_xrefs_list_with_none.append(record.dbxrefs)
        db_xrefs_list_with_none_list = [['None'] if db_xrefs_value is None else db_xrefs_value for db_xrefs_value in db_xrefs_list_with_none]
        db_xrefs_list = []
        for db_xrefs_none in db_xrefs_list_with_none_list:
            db_xrefs_list.append(db_xrefs_none)
        print(db_xrefs_list)
        print("")
        # Writing DB XREFS
        print("Step 2: Writing DB XREFS, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        db_xrefs_string = '\n'.join(str(db_xrefs_value) for db_xrefs_value in db_xrefs_list)
        print(db_xrefs_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_db_xrefs.txt"), "w")
        handle_outputfile.write(db_xrefs_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 32: DB XREFS Extract Complete")
        print("")
        break

# Definition Extract ##################################################
def definition():
    while True:
        # Module Name
        print("Module 33: Definition Extract")
        print("")
        # Parsing Definition
        print("Step 1: Parsing Definition, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        definition_list_with_none = []
        for record in handle_parse_db:
            definition_list_with_none.append(record.description)
        definition_list_with_none_list = [['None'] if definition_value is None else definition_value for definition_value in definition_list_with_none]
        definition_list = []
        for definition_none in definition_list_with_none_list:
            definition_list.append(definition_none)
        print(definition_list)
        print("")
        # Writing Definition
        print("Step 2: Writing Definition, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        definition_string = '\n'.join(str(definition_value) for definition_value in definition_list)
        print(definition_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_definition.txt"), "w")
        handle_outputfile.write(definition_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 33: Definition Extract Complete")
        print("")
        break

# Features Extract ##################################################
def features():
    while True:
        # Module Name
        print("Module 34: Features Extract")
        print("")
        # Parsing Features
        print("Step 1: Parsing Features, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        features_list_with_none = []
        for record in handle_parse_db:
            features_list_with_none.append(record.features)
        features_list_with_none_list = [['None'] if features_value is None else features_value for features_value in features_list_with_none]
        features_list = []
        for features_none in features_list_with_none_list:
            features_list.append(features_none)
        print(features_list)
        print("")
        # Writing Features
        print("Step 2: Writing Features, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        features_string = '\n'.join(str(features_value) for features_value in features_list)
        print(features_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_features.txt"), "w")
        handle_outputfile.write(features_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 34: Features Extract Complete")
        print("")
        break

# Features Length Extract ##################################################
def featureslength():
    while True:
        # Module Name
        print("Module 35: Features Length Extract")
        print("")
        # Parsing Features
        print("Step 1: Parsing Features Length, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        features_length_list_with_none = []
        for record in handle_parse_db:
            features_length_list_with_none.append(len(record.features))
        features_length_list_with_none_list = [['None'] if features_length_value is None else features_length_value for features_length_value in features_length_list_with_none]
        features_length_list = []
        for features_length_none in features_length_list_with_none_list:
            features_length_list.append(features_length_none)
        print(features_length_list)
        print("")
        # Writing Features Length
        print("Step 2: Writing Features Length, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        features_length_string = '\n'.join(str(features_length_value) for features_length_value in features_length_list)
        print(features_length_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_features_length.txt"), "w")
        handle_outputfile.write(features_length_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 35: Features Length Extract Complete")
        print("")
        break

# Format Extract ##################################################
def formatdescription():
    while True:
        # Module Name
        print("Module 36: Format Extract")
        print("")
        # Parsing Format
        print("Step 1: Parsing Format, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        format_list_with_none = []
        for record in handle_parse_db:
            format_list_with_none.append(record.format)
        format_list_with_none_list = [['None'] if format_value is None else format_value for format_value in format_list_with_none]
        format_list = []
        for format_none in format_list_with_none_list:
            format_list.append(format_none)
        print(format_list)
        print("")
        # Writing Format
        print("Step 2: Writing Format, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        format_string = '\n'.join(str(format_value) for format_value in format_list)
        print(format_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_format.txt"), "w")
        handle_outputfile.write(format_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 36: Format Extract Complete")
        print("")
        break

# Length Extract ##################################################
def length():
    while True:
        # Module Name
        print("Module 37: Length Extract")
        print("")
        # Parsing Length
        print("Step 1: Parsing Length, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        length_list_with_none = []
        for record in handle_parse_db:
            length_list_with_none.append(len(record))
        length_list_with_none_list = [['None'] if length_value is None else length_value for length_value in length_list_with_none]
        length_list = []
        for length_none in length_list_with_none_list:
            length_list.append(length_none)
        print(length_list)
        print("")
        # Writing Length
        print("Step 2: Writing Length, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        length_string = '\n'.join(str(length_value) for length_value in length_list)
        print(length_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_length.txt"), "w")
        handle_outputfile.write(length_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 37: Length Extract Complete")
        print("")
        break

# Length Annotations Extract ##################################################
def lengthannotations():
    while True:
        # Module Name
        print("Module 38: Length Annotations Extract")
        print("")
        # Parsing Length Annotations
        print("Step 1: Parsing Length Annotations, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        length_annotations_list_with_none = []
        for record in handle_parse_db:
            length_annotations_list_with_none.append(len(record.annotations))
        length_annotations_list_with_none_list = [['None'] if length_annotations_value is None else length_annotations_value for length_annotations_value in length_annotations_list_with_none]
        length_annotations_list = []
        for length_annotations_none in length_annotations_list_with_none_list:
            length_annotations_list.append(length_annotations_none)
        print(length_annotations_list)
        print("")
        # Writing Length Annotations
        print("Step 2: Writing Length Annotations, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        length_annotations_string = '\n'.join(str(length_annotations_value) for length_annotations_value in length_annotations_list)
        print(length_annotations_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_length_annotations.txt"), "w")
        handle_outputfile.write(length_annotations_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 38: Length Annotations Extract Complete")
        print("")
        break

# Letter Annotations Extract ##################################################
def letterannotations():
    while True:
        # Module Name
        print("Module 39: Letter Annotations Extract")
        print("")
        # Parsing Letter Annotations
        print("Step 1: Parsing Letter Annotations, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        letter_annotations_list_with_none = []
        for record in handle_parse_db:
            letter_annotations_list_with_none.append(record.letter_annotations)
        letter_annotations_list_with_none_list = [['None'] if letter_annotations_value is None else letter_annotations_value for letter_annotations_value in letter_annotations_list_with_none]
        letter_annotations_list = []
        for letter_annotations_none in letter_annotations_list_with_none_list:
            letter_annotations_list.append(letter_annotations_none)
        print(letter_annotations_list)
        print("")
        # Writing Letter Annotations
        print("Step 2: Writing Letter Annotations, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        letter_annotations_string = '\n'.join(str(letter_annotations_value) for letter_annotations_value in letter_annotations_list)
        print(letter_annotations_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_letter_annotations.txt"), "w")
        handle_outputfile.write(letter_annotations_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 39: Letter Annotations Extract Complete")
        print("")
        break

# Letter Annotations Phred Quality Extract ##################################################
def letterannotationsphredquality():
    while True:
        # Module Name
        print("Module 40: Letter Annotations Phred Quality Extract")
        print("")
        # Parsing Letter Annotations Phred Quality
        print("Step 1: Parsing Letter Annotations Phred Quality, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        letter_annotations_phred_quality_list_with_none = []
        for record in handle_parse_db:
            letter_annotations_quality_list_with_none.append(record.letter_annotations["phred_quality"])
        letter_annotations_quality_list_with_none_list = [['None'] if letter_annotations_quality_value is None else letter_annotations_quality_value for letter_annotations_quality_value in letter_annotations_quality_list_with_none]
        letter_annotations_quality_list = []
        for letter_annotations_quality_none in letter_annotations_quality_list_with_none_list:
            letter_annotations_quality_list.append(letter_annotations_quality_none)
        print(letter_annotations_quality_list)
        print("")
        # Writing Letter Annotations Phred Quality
        print("Step 2: Writing Letter Annotations Phred Quality, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        letter_annotations_quality_string = '\n'.join(str(letter_annotations_quality_value) for letter_annotations_quality_value in letter_annotations_quality_list)
        print(letter_annotations_quality_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_letter_annotations_quality.txt"), "w")
        handle_outputfile.write(letter_annotations_quality_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 40: Letter Annotations Phred Quality Extract Complete")
        print("")
        break

# Locus Extract ##################################################
def locus():
    while True:
        # Module Name
        print("Module 41: Locus Extract")
        print("")
        # Parsing Locus
        print("Step 1: Parsing Locus, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        locus_list_with_none = []
        for record in handle_parse_db:
            locus_list_with_none.append(record.name)
        locus_list_with_none_list = [['None'] if locus_value is None else locus_value for locus_value in locus_list_with_none]
        locus_list = []
        for locus_none in locus_list_with_none_list:
            locus_list.append(locus_none)
        print(locus_list)
        print("")
        # Writing Locus
        print("Step 2: Writing Locus, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        locus_string = '\n'.join(str(locus_value) for locus_value in locus_list)
        print(locus_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_locus.txt"), "w")
        handle_outputfile.write(locus_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 41: Locus Extract Complete")
        print("")
        break

# Origin Extract ##################################################
def origin():
    while True:
        # Module Name
        print("Module 42: Origin Extract")
        print("")
        # Parsing Origin
        print("Step 1: Parsing Origin, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        origin_list_with_none = []
        for record in handle_parse_db:
            origin_list_with_none.append(repr(record.seq))
        origin_list_with_none_list = [['None'] if origin_value is None else origin_value for origin_value in origin_list_with_none]
        origin_list = []
        for origin_none in origin_list_with_none_list:
            origin_list.append(origin_none)
        print(origin_list)
        print("")
        # Writing Origin
        print("Step 2: Writing Origin, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        origin_string = '\n'.join(str(origin_value) for origin_value in origin_list)
        print(origin_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_origin.txt"), "w")
        handle_outputfile.write(origin_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 42: Origin Extract Complete")
        print("")
        break

# Origin All Extract ##################################################
def originall():
    while True:
        # Module Name
        print("Module 43: Origin All Extract")
        print("")
        # Parsing Origin All
        print("Step 1: Parsing Origin All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        origin_list_with_none = []
        for record in handle_parse_db:
            origin_list_with_none.append(record.seq)
        origin_list_with_none_list = [['None'] if origin_value is None else origin_value for origin_value in origin_list_with_none]
        origin_list = []
        for origin_none in origin_list_with_none_list:
            origin_list.append(origin_none)
        print(origin_list)
        print("")
        # Writing Origin All
        print("Step 2: Writing Origin All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        origin_string = '\n'.join(str(origin_value) for origin_value in origin_list)
        print(origin_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_origin_all.txt"), "w")
        handle_outputfile.write(origin_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 43: Origin All Extract Complete")
        print("")
        break

# ORF5 Filter ##################################################
def orf5filter():
    while True:
        # Module Name
        print("Module 44: ORF5 Filter")
        print("")
        # Parsing ORF5
        print("Step 1: Parsing ORF5, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list = []
        #organism_list = []
        #gene_list = []
        #length_list = []
        #pubmedid_list = []
        #date_list = []
        #note_list = []
        country_list = []
        #origin_list = []
        for record in handle_parse_db:
            if record.annotations["organism"] == 'Porcine reproductive and respiratory syndrome virus':
                if len(record) == 600 or len(record) == 603 or len(record) == 606:
                    if record.features:
                        for feature in record.features:
                            if feature.type == "source":
                                accession_list.append(record.id)
                                #organism_list.append(record.annotations["organism"])
                                #gene_list.append('ORF5')
                                #length_list.append(str(len(record)))
                                #pubmedid_list.append(record.annotations["references"][0].pubmed_id)
                                #date_list.append(record.annotations["date"])
                                #note_list.append(feature.qualifiers.get("note"))
                                country_list.append(feature.qualifiers.get("country"))
                                #origin_list.append(str(record.seq))
        # pubmedid
        #organism_list_with_none = [['None'] if organism_value is None else organism_value for organism_value in organism_list]
        #organism_none_list = []
        #for organism_none in organism_list_with_none:
        #    organism_none_list.append(organism_none)
        #date_list_with_none = [['None'] if date_value is None else date_value for date_value in date_list]
        #date_none_list = []
        #for date_none in date_list_with_none:
        #    date_none_list.append(date_none)
        #note_list_with_none = [['None'] if note_value is None else note_value for note_value in note_list]
        #note_none_list = []
        #for note_none in note_list_with_none:
        #    note_none_list.append('\n'.join(note_none))
        country_list_with_none = [['None'] if country_value is None else country_value for country_value in country_list]
        country_none_list = []
        for country_none in country_list_with_none:
            country_none_list.append('\n'.join(country_none))
        #origin_list_with_none = [['None'] if origin_value is None else origin_value for origin_value in origin_list]
        #origin_none_list = []
        #for origin_none in origin_list_with_none:
        #    origin_none_list.append(origin_none)
        orf5_filter_zip = zip(accession_list, country_none_list)
        print(orf5_filter_zip)
    
        orf5_list = []
        for orf5_filter in orf5_filter_zip:
            orf5_list.append("\t".join(orf5_filter))
        # print(ORF5_list)
        print("")
        # Writing ORF5
        print("Step 2: Writing ORF5, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        orf5_string = '\n'.join(str(orf5_value) for orf5_value in orf5_list)
        print(orf5_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_orf5_filter.txt"), "w")
        handle_outputfile.write(orf5_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 44: ORF5 Filter Complete")
        print("")
        break

# ORF5 EU Extract ##################################################
def orf5euextract():
    while True:
        # Module Name
        print("Module 45: ORF5 EU Extract")
        print("")
        # Parsing ORF5 EU
        print("Step 1: Parsing ORF5 EU, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        # 1 #################################################
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        country_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        country = feature.qualifiers.get('country')
                        country_list_with_none.append(country)
        country_list_with_none_list = [['None'] if country_value is None else country_value for country_value in country_list_with_none]
        country_extract_list = []
        for country_none in country_list_with_none_list:
            country_extract_list.append('\n'.join(country_none))
        print(country_extract_list)
        country_string = "\n".join(country_extract_list)
        print(country_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_country.txt"), "w")
        handle_outputfile.write(country_string)
        handle_outputfile.close()
        # 2 ##################################################
        handle_inputfile = open(os.path.expanduser("~/extractor/local/export/source_country.txt"))
        country_define_list = []
        latitude_define_list = []
        longitude_define_list = []
        continent_define_list = []
        for row in handle_inputfile:
            row = row.rstrip()
            if row.startswith('Albania'):
                country_define_list.append('Albania')
                latitude_define_list.append('41.3275000')
                longitude_define_list.append('19.8188896')
                continent_define_list.append('Europe')
            if row.startswith('Australia'):
                country_define_list.append('Australia')
                latitude_define_list.append('-33.8678500')
                longitude_define_list.append('133.7751360')
                continent_define_list.append('Australia/Oceania')
            if row.startswith('Austria'):
                country_define_list.append('Austria')
                latitude_define_list.append('48.2084878')
                longitude_define_list.append('16.3720760')
                continent_define_list.append('Europe')
            if row.startswith('Belarus'):
                country_define_list.append('Belarus')
                latitude_define_list.append('53.9000000')
                longitude_define_list.append('27.5666676')
                continent_define_list.append('Europe')
            if row.startswith('Bulgaria'):
                country_define_list.append('Bulgaria')
                latitude_define_list.append('42.6975135')
                longitude_define_list.append('23.3241463')
                continent_define_list.append('Europe')
            if row.startswith('Cambodia'):
                country_define_list.append('Cambodia')
                latitude_define_list.append('11.5500000')
                longitude_define_list.append('104.9166641')
                continent_define_list.append('Asia')
            if row == 'Canada':
                country_define_list.append('Canada')
                latitude_define_list.append('43.7001138')
                longitude_define_list.append('-79.4163055')
                continent_define_list.append('North American')
            if row == 'Canada: ON':
                country_define_list.append('Canada: Ontario')
                latitude_define_list.append('51.253775')
                longitude_define_list.append('-85.323214')
                continent_define_list.append('North American')
            if row == 'Canada: QC':
                country_define_list.append('Canada: Quebec')
                latitude_define_list.append('46.8122791')
                longitude_define_list.append('-71.2145386')
                continent_define_list.append('North American')
            if row == 'Canada: SK':
                country_define_list.append('Canada: Saskatchewan')
                latitude_define_list.append('52.939916')
                longitude_define_list.append('-106.450864')
                continent_define_list.append('North American')
            if row == 'Chile':
                country_define_list.append('Chile')
                latitude_define_list.append('-33.4500000')
                longitude_define_list.append('-70.6666641')
                continent_define_list.append('South American')
            if row == 'Chile: Maipo':
                country_define_list.append('Chile: Maipo River')
                latitude_define_list.append('-33.8683360')
                longitude_define_list.append('-70.8111840')
                continent_define_list.append('South American')
            if row == 'China':
                country_define_list.append('China')
                latitude_define_list.append('31.2222222')
                longitude_define_list.append('121.4580536')
                continent_define_list.append('Asia')
            if row == 'China: Anhui':
                country_define_list.append('China: Anhui')
                latitude_define_list.append('30.6006770')
                longitude_define_list.append('117.9249000')
                continent_define_list.append('Asia')
            if row == 'China: Boxing City':
                country_define_list.append('China: Boxing City')
                latitude_define_list.append('37.1502260')
                longitude_define_list.append('118.1318150')
                continent_define_list.append('Asia')
            if row == 'China: Chongqing':
                country_define_list.append('China: Chongqing')
                latitude_define_list.append('29.5627778')
                longitude_define_list.append('106.5527802')
                continent_define_list.append('Asia')
            if row == 'China: Dezhou City':
                country_define_list.append('China: Dezhou City')
                latitude_define_list.append('37.4340920')
                longitude_define_list.append('116.3574640')
                continent_define_list.append('Asia')
            if row == 'China: Fujian, Fuzhou':
                country_define_list.append('China: Fujian, Fuzhou')
                latitude_define_list.append('26.0613889')
                longitude_define_list.append('119.3061142')
                continent_define_list.append('Asia')
            if row == 'China: Fujian, Nanping':
                country_define_list.append('China: Fujian, Nanping')
                latitude_define_list.append('26.6417740')
                longitude_define_list.append('118.1777100')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong':
                country_define_list.append('China: Guangdong')
                latitude_define_list.append('23.3790330')
                longitude_define_list.append('113.7632830')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Foshan':
                country_define_list.append('China: Guangdong, Foshan')
                latitude_define_list.append('23.0333333')
                longitude_define_list.append('113.1166687')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Guangzhou':
                country_define_list.append('China: Guangdong, Guangzhou')
                latitude_define_list.append('23.1166667')
                longitude_define_list.append('113.2500000')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Guangzhou, Huadu':
                country_define_list.append('China: Guangdong, Guangzhou, Huadu')
                latitude_define_list.append('23.4041650')
                longitude_define_list.append('113.2202180')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Heyuan':
                country_define_list.append('China: Guangdong, Heyuan')
                latitude_define_list.append('23.7436850')
                longitude_define_list.append('114.7009610')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Huizhou':
                country_define_list.append('China: Guangdong, Huizhou')
                latitude_define_list.append('23.1122570')
                longitude_define_list.append('114.4158010')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Jiangmen':
                country_define_list.append('China: Guangdong, Jiangmen')
                latitude_define_list.append('22.5791170')
                longitude_define_list.append('113.0815080')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Jiangmeng':
                country_define_list.append('China: Guangdong, Jiangmeng')
                latitude_define_list.append('22.5791170')
                longitude_define_list.append('113.0815080')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Jieyang':
                country_define_list.append('China: Guangdong, Jieyang')
                latitude_define_list.append('23.5297222')
                longitude_define_list.append('116.3655548')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Kaiping':
                country_define_list.append('China: Guangdong, Kaiping')
                latitude_define_list.append('22.3763950')
                longitude_define_list.append('112.6985450')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Maoming':
                country_define_list.append('China: Guangdong, Maoming')
                latitude_define_list.append('21.6629910')
                longitude_define_list.append('110.9254390')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Meizhou':
                country_define_list.append('China: Guangdong, Meizhou')
                latitude_define_list.append('24.2885780')
                longitude_define_list.append('116.1225230')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Qingyuan':
                country_define_list.append('China: Guangdong, Qingyuan')
                latitude_define_list.append('23.6817740')
                longitude_define_list.append('113.0560420')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Sanshui':
                country_define_list.append('China: Guangdong, Sanshui')
                latitude_define_list.append('23.1560450')
                longitude_define_list.append('112.8966060')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shantou':
                country_define_list.append('China: Guangdong, Shantou')
                latitude_define_list.append('23.3600000')
                longitude_define_list.append('116.6783371')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shanwei':
                country_define_list.append('China: Guangdong, Shanwei')
                latitude_define_list.append('22.7861860')
                longitude_define_list.append('115.3751580')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shaoguan':
                country_define_list.append('China: Guangdong, Shaoguan')
                latitude_define_list.append('24.8000000')
                longitude_define_list.append('113.5833359')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shenzhen':
                country_define_list.append('China: Guangdong, Shenzhen')
                latitude_define_list.append('22.5455377')
                longitude_define_list.append('114.0682983')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shiling':
                country_define_list.append('China: Guangdong, Shiling')
                latitude_define_list.append('23.4607280')
                longitude_define_list.append('113.1515960')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Sihui':
                country_define_list.append('China: Guangdong, Sihui')
                latitude_define_list.append('23.3270010')
                longitude_define_list.append('112.7341030')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Yangchun':
                country_define_list.append('China: Guangdong, Yangchun')
                latitude_define_list.append('22.1704370')
                longitude_define_list.append('111.7915390')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Yangjiang':
                country_define_list.append('China: Guangdong, Yangjiang')
                latitude_define_list.append('21.8579580')
                longitude_define_list.append('111.9822320')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Yunfu':
                country_define_list.append('China: Guangdong, Yunfu')
                latitude_define_list.append('22.9333333')
                longitude_define_list.append('112.0333328')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Zhanjiang':
                country_define_list.append('China: Guangdong, Zhanjiang')
                latitude_define_list.append('21.2000000')
                longitude_define_list.append('110.3833313')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Zhaoqing':
                country_define_list.append('China: Guangdong, Zhaoqing')
                latitude_define_list.append('23.0471910')
                longitude_define_list.append('112.4650910')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Zhongshan':
                country_define_list.append('China: Guangdong, Zhongshan')
                latitude_define_list.append('22.5175850')
                longitude_define_list.append('113.3927700')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong,Guangzhou':
                country_define_list.append('China: Guangdong, Guangzhou')
                latitude_define_list.append('23.1166667')
                longitude_define_list.append('113.2500000')
                continent_define_list.append('Asia')
            if row == 'China: Guangxi':
                country_define_list.append('China: Guangxi')
                latitude_define_list.append('23.7247600')
                longitude_define_list.append('108.8076190')
                continent_define_list.append('Asia')
            if row == 'China: Guodian City':
                country_define_list.append('China: Guodian City')
                latitude_define_list.append('36.0839480')
                longitude_define_list.append('103.6665500')
                continent_define_list.append('Asia')
            if row == 'China: Hainan':
                country_define_list.append('China: Hainan')
                latitude_define_list.append('19.5663950')
                longitude_define_list.append('109.9496860')
                continent_define_list.append('Asia')
            if row == 'China: Hebei Province':
                country_define_list.append('China: Hebei Province')
                latitude_define_list.append('37.8956590')
                longitude_define_list.append('114.9042210')
                continent_define_list.append('Asia')
            if row == 'China: Heilongjiang':
                country_define_list.append('China: Heilongjiang')
                latitude_define_list.append('47.1216470')
                longitude_define_list.append('128.7382310')
                continent_define_list.append('Asia')
            if row == 'China: Henan':
                country_define_list.append('China: Henan')
                latitude_define_list.append('34.2904300')
                longitude_define_list.append('113.3823550')
                continent_define_list.append('Asia')
            if row == 'China: Hennan,Yuoyang':
                country_define_list.append('China: Hennan,Yuoyang')
                latitude_define_list.append('34.6202020')
                longitude_define_list.append('112.4539260')
                continent_define_list.append('Asia')
            if row == 'China: Hnan,Luoyang':
                country_define_list.append('China: Hnan,Luoyang')
                latitude_define_list.append('34.6836111')
                longitude_define_list.append('112.4536133')
                continent_define_list.append('Asia')
            if row == 'China: Hunan':
                country_define_list.append('China: Hunan')
                latitude_define_list.append('27.6253000')
                longitude_define_list.append('111.8568590')
                continent_define_list.append('Asia')
            if row == 'China: Jiangxi':
                country_define_list.append('China: Jiangxi')
                latitude_define_list.append('27.0874560')
                longitude_define_list.append('114.9042210')
                continent_define_list.append('Asia')
            if row == 'China: Jilin, Changchun':
                country_define_list.append('China: Jilin, Changchun')
                latitude_define_list.append('43.8800000')
                longitude_define_list.append('125.3227768')
                continent_define_list.append('Asia')
            if row == 'China: Jinan City':
                country_define_list.append('China: Jinan City')
                latitude_define_list.append('36.6683333')
                longitude_define_list.append('116.9972229')
                continent_define_list.append('Asia')
            if row == 'China: Liaocheng City':
                country_define_list.append('China: Liaocheng City')
                latitude_define_list.append('36.4570300')
                longitude_define_list.append('115.9854600')
                continent_define_list.append('Asia')
            if row == 'China: SD Province':
                country_define_list.append('China: SD Province')
                latitude_define_list.append('35.8939570')
                longitude_define_list.append('117.9249000')
                continent_define_list.append('Asia')
            if row == 'China: Shandong':
                country_define_list.append('China: Shandong')
                latitude_define_list.append('35.8939570')
                longitude_define_list.append('117.9249000')
                continent_define_list.append('Asia')
            if row == 'China: Sichuan':
                country_define_list.append('China: Sichuan')
                latitude_define_list.append('30.6516520')
                longitude_define_list.append('104.0759310')
                continent_define_list.append('Asia')
            if row == 'China: Wenzu City':
                country_define_list.append('China: Wenzu City')
                latitude_define_list.append('27.9938280')
                longitude_define_list.append('120.6993610')
                continent_define_list.append('Asia')
            if row == 'China: Wuzhou':
                country_define_list.append('China: Wuzhou')
                latitude_define_list.append('23.4769620')
                longitude_define_list.append('111.2791150')
                continent_define_list.append('Asia')
            if row == 'China: Xinjiang':
                country_define_list.append('China: Xinjiang')
                latitude_define_list.append('42.5246360')
                longitude_define_list.append('87.5395850')
                continent_define_list.append('Asia')
            if row == 'China: Yun Nan':
                country_define_list.append('China: Yun Nan')
                latitude_define_list.append('24.4752850')
                longitude_define_list.append('101.3431060')
                continent_define_list.append('Asia')
            if row == 'China: Yunnan':
                country_define_list.append('China: Yunnan')
                latitude_define_list.append('24.4752850')
                longitude_define_list.append('101.3431060')
                continent_define_list.append('Asia')
            if row == 'China: Yunnan, Qujing':
                country_define_list.append('China: Yunnan, Qujing')
                latitude_define_list.append('25.4899990')
                longitude_define_list.append('103.7961670')
                continent_define_list.append('Asia')
            if row == 'China: Zhangqiu City':
                country_define_list.append('China: Zhangqiu City')
                latitude_define_list.append('36.6812590')
                longitude_define_list.append('117.5262280')
                continent_define_list.append('Asia')
            if row.startswith('Croatia'):
                country_define_list.append('Croatia')
                latitude_define_list.append('45.8131276')
                longitude_define_list.append('15.9775333')
                continent_define_list.append('Europe')
            if row.startswith('Czech Republic'):
                country_define_list.append('Czech Republic')
                latitude_define_list.append('50.0878368')
                longitude_define_list.append('14.4241323')
                continent_define_list.append('Europe')
            if row.startswith('Denmark'):
                country_define_list.append('Denmark')
                latitude_define_list.append('55.6776812')
                longitude_define_list.append('12.5709343')
                continent_define_list.append('Europe')
            if row.startswith('France'):
                country_define_list.append('France')
                latitude_define_list.append('48.8534100')
                longitude_define_list.append('2.3487999')
                continent_define_list.append('Europe')
            if row.startswith('Germany'):
                country_define_list.append('Germany')
                latitude_define_list.append('52.5166667')
                longitude_define_list.append('13.3999996')
                continent_define_list.append('Europe')
            if row.startswith('Hungary'):
                country_define_list.append('Hungary')
                latitude_define_list.append('47.5000000')
                longitude_define_list.append('19.0833340')
                continent_define_list.append('Europe')
            if row.startswith('India'):
                country_define_list.append('India')
                latitude_define_list.append('19.0144100')
                longitude_define_list.append('72.8479385')
                continent_define_list.append('Asia')
            if row.startswith('Italy'):
                country_define_list.append('Italy')
                latitude_define_list.append('41.9000000')
                longitude_define_list.append('12.4833336')
                continent_define_list.append('Europe')
            if row == 'Japan':
                country_define_list.append('Japan')
                latitude_define_list.append('35.6895266')
                longitude_define_list.append('139.6916809')
                continent_define_list.append('Asia')
            if row == 'Japan:Aichi':
                country_define_list.append('Japan: Aichi')
                latitude_define_list.append('34.9666667')
                longitude_define_list.append('136.6166687')
                continent_define_list.append('Asia')
            if row == 'Japan:Aomori':
                country_define_list.append('Japan: Aomori')
                latitude_define_list.append('40.8211111')
                longitude_define_list.append('140.7511139')
                continent_define_list.append('Asia')
            if row == 'Japan:Chiba':
                country_define_list.append('Japan: Chiba')
                latitude_define_list.append('35.6000000')
                longitude_define_list.append('140.1166687')
                continent_define_list.append('Asia')
            if row == 'Japan:Ehime':
                country_define_list.append('Japan: Ehime')
                latitude_define_list.append('33.8416240')
                longitude_define_list.append('132.7656810')
                continent_define_list.append('Asia')
            if row == 'Japan:Fukushima':
                country_define_list.append('Japan: Fukushima')
                latitude_define_list.append('37.7500000')
                longitude_define_list.append('140.4666595')
                continent_define_list.append('Asia')
            if row == 'Japan:Gunma':
                country_define_list.append('Japan: Gunma')
                latitude_define_list.append('36.3906670')
                longitude_define_list.append('139.0604060')
                continent_define_list.append('Asia')
            if row == 'Japan:Hokkaido':
                country_define_list.append('Japan: Hokkaido')
                latitude_define_list.append('43.2203270')
                longitude_define_list.append('142.8634740')
                continent_define_list.append('Asia')
            if row == 'Japan:Ibaraki':
                country_define_list.append('Japan: Ibaraki')
                latitude_define_list.append('34.8164106')
                longitude_define_list.append('135.5682831')
                continent_define_list.append('Asia')
            if row == 'Japan:Ishikawa':
                country_define_list.append('Japan: Ishikawa')
                latitude_define_list.append('36.5946820')
                longitude_define_list.append('136.6255730')
                continent_define_list.append('Asia')
            if row == 'Japan:Iwate':
                country_define_list.append('Japan: Iwate')
                latitude_define_list.append('39.7036190')
                longitude_define_list.append('141.1526840')
                continent_define_list.append('Asia')
            if row == 'Japan:Kagoshima':
                country_define_list.append('Japan: Kagoshima')
                latitude_define_list.append('31.6000000')
                longitude_define_list.append('130.5500031')
                continent_define_list.append('Asia')
            if row == 'Japan:Kyoto':
                country_define_list.append('Japan: Kyoto')
                latitude_define_list.append('35.0210700')
                longitude_define_list.append('135.7538452')
                continent_define_list.append('Asia')
            if row == 'Japan:Nagasaki':
                country_define_list.append('Japan: Nagasaki')
                latitude_define_list.append('32.7550000')
                longitude_define_list.append('129.8683319')
                continent_define_list.append('Asia')
            if row == 'Japan:Niigata':
                country_define_list.append('Japan: Niigata')
                latitude_define_list.append('37.9166667')
                longitude_define_list.append('139.0500031')
                continent_define_list.append('Asia')
            if row == 'Japan:Osaka':
                country_define_list.append('Japan: Osaka')
                latitude_define_list.append('34.6937398')
                longitude_define_list.append('135.5021820')
                continent_define_list.append('Asia')
            if row == 'Japan:Saga':
                country_define_list.append('Japan: Saga')
                latitude_define_list.append('33.2634820')
                longitude_define_list.append('130.3008580')
                continent_define_list.append('Asia')
            if row == 'Japan:Saitama':
                country_define_list.append('Japan: Saitama')
                latitude_define_list.append('35.8617290')
                longitude_define_list.append('139.6454820')
                continent_define_list.append('Asia')
            if row == 'Japan:Shiga':
                country_define_list.append('Japan: Shiga')
                latitude_define_list.append('35.0045310')
                longitude_define_list.append('135.8685900')
                continent_define_list.append('Asia')
            if row == 'Japan:Shizuoka':
                country_define_list.append('Japan: Shizuoka')
                latitude_define_list.append('34.9666667')
                longitude_define_list.append('138.3833313')
                continent_define_list.append('Asia')
            if row == 'Japan:Tochigi':
                country_define_list.append('Japan: Tochigi')
                latitude_define_list.append('36.5657250')
                longitude_define_list.append('139.8835650')
                continent_define_list.append('Asia')
            if row == 'Japan:Yamagata':
                country_define_list.append('Japan: Yamagata')
                latitude_define_list.append('38.2527778')
                longitude_define_list.append('140.3374939')
                continent_define_list.append('Asia')
            if row == 'Japan:Yamaguchi':
                country_define_list.append('Japan: Yamaguchi')
                latitude_define_list.append('34.1859560')
                longitude_define_list.append('131.4706490')
                continent_define_list.append('Asia')
            if row == 'Japan:Yamanashi':
                country_define_list.append('Japan: Yamanashi')
                latitude_define_list.append('35.6641580')
                longitude_define_list.append('138.5684490')
                continent_define_list.append('Asia')
            if row.startswith('Laos'):
                country_define_list.append('Laos')
                latitude_define_list.append('19.856270')
                longitude_define_list.append('102.495496')
                continent_define_list.append('Asia')
            if row.startswith('Lithuania'):
                country_define_list.append('Lithuania')
                latitude_define_list.append('54.6833333')
                longitude_define_list.append('25.3166676')
                continent_define_list.append('Europe')
            if row == 'Malaysia':
                country_define_list.append('Malaysia')
                latitude_define_list.append('3.1666667')
                longitude_define_list.append('101.6999969')
                continent_define_list.append('Asia')
            if row == 'Malaysia: Sarawak':
                country_define_list.append('Malaysia: Sarawak')
                latitude_define_list.append('1.5532780')
                longitude_define_list.append('110.3592130')
                continent_define_list.append('Asia')
            if row == 'Malaysia: Selangor':
                country_define_list.append('Malaysia: Selangor')
                latitude_define_list.append('3.3500000')
                longitude_define_list.append('101.2500000')
                continent_define_list.append('Asia')
            if row == 'Mexico':
                country_define_list.append('Mexico')
                latitude_define_list.append('19.4341667')
                longitude_define_list.append('-99.1386108')
                continent_define_list.append('North American')
            if row == 'Mexico: Morelos':
                country_define_list.append('Mexico: Morelos')
                latitude_define_list.append('18.8000000')
                longitude_define_list.append('-98.9499969')
                continent_define_list.append('North American')
            if row.startswith('Myanmar'):
                country_define_list.append('Myanmar')
                latitude_define_list.append('21.9162210')
                longitude_define_list.append('95.9559740')
                continent_define_list.append('Asia')
            if row.startswith('Poland'):
                country_define_list.append('Poland')
                latitude_define_list.append('52.2500000')
                longitude_define_list.append('21.0000000')
                continent_define_list.append('Europe')
            if row.startswith('Romania'):
                country_define_list.append('Romania')
                latitude_define_list.append('44.4333333')
                longitude_define_list.append('26.1000004')
                continent_define_list.append('Europe')
            if row.startswith('Russia'):
                country_define_list.append('Russia')
                latitude_define_list.append('55.7522222')
                longitude_define_list.append('37.6155548')
                continent_define_list.append('Europe')
            if row.startswith('Serbia'):
                country_define_list.append('Serbia')
                latitude_define_list.append('44.7865680')
                longitude_define_list.append('20.4489220')
                continent_define_list.append('Europe')
            if row.startswith('Singapore'):
                country_define_list.append('Singapore')
                latitude_define_list.append('1.3520830')
                longitude_define_list.append('103.8198360')
                continent_define_list.append('Asia')
            if row.startswith('Slovakia'):
                country_define_list.append('Slovakia')
                latitude_define_list.append('48.1485960')
                longitude_define_list.append('17.1077480')
                continent_define_list.append('Europe')
            if row.startswith('Slovenia'):
                country_define_list.append('Slovenia')
                latitude_define_list.append('46.0552778')
                longitude_define_list.append('14.5144444')
                continent_define_list.append('Europe')
            if row.startswith('South Korea'):
                country_define_list.append('South Korea')
                latitude_define_list.append('37.5663889')
                longitude_define_list.append('126.9997253')
                continent_define_list.append('Asia')
            if row.startswith('Spain'):
                country_define_list.append('Spain')
                latitude_define_list.append('40.4165021')
                longitude_define_list.append('-3.7025642')
                continent_define_list.append('Europe')
            if row.startswith('Taiwan'):
                country_define_list.append('Taiwan')
                latitude_define_list.append('23.6978100')
                longitude_define_list.append('120.9605150')
                continent_define_list.append('Asia')
            if row.startswith('Thailand'):
                country_define_list.append('Thailand')
                latitude_define_list.append('13.7500000')
                longitude_define_list.append('100.5166702')
                continent_define_list.append('Asia')
            if row.startswith('United Kingdom'):
                country_define_list.append('United Kingdom')
                latitude_define_list.append('51.5084153')
                longitude_define_list.append('-0.1255327')
                continent_define_list.append('Europe')
            if row == 'USA':
                country_define_list.append('USA')
                latitude_define_list.append('40.7142691')
                longitude_define_list.append('-74.0059738')
                continent_define_list.append('North American')
            if row == 'USA: AZ':
                country_define_list.append('USA: Arizona')
                latitude_define_list.append('15.6333333')
                longitude_define_list.append('-87.3166656')
                continent_define_list.append('North American')
            if row == 'USA: Can':
                country_define_list.append('USA: Can Lane')
                latitude_define_list.append('29.1949620')
                longitude_define_list.append('-81.0073860')
                continent_define_list.append('North American')
            if row == 'USA: CAN':
                country_define_list.append('USA: Can Lane')
                latitude_define_list.append('29.1949620')
                longitude_define_list.append('-81.0073860')
                continent_define_list.append('North American')
            if row == 'USA: CO':
                country_define_list.append('USA: Colorado')
                latitude_define_list.append('38.8338816')
                longitude_define_list.append('-104.8213654')
                continent_define_list.append('North American')
            if row == 'USA: IA':
                country_define_list.append('USA: Iowa')
                latitude_define_list.append('41.8780030')
                longitude_define_list.append('-93.0977020')
                continent_define_list.append('North American')
            if row == 'USA: IL':
                country_define_list.append('USA: Illinois')
                latitude_define_list.append('40.6331250')
                longitude_define_list.append('-89.3985280')
                continent_define_list.append('North American')
            if row == 'USA: Illinois':
                country_define_list.append('USA: Illinois')
                latitude_define_list.append('40.6331250')
                longitude_define_list.append('-89.3985280')
                continent_define_list.append('North American')
            if row == 'USA: IN':
                country_define_list.append('USA: Indiana')
                latitude_define_list.append('39.7683765')
                longitude_define_list.append('-86.1580429')
                continent_define_list.append('North American')
            if row == 'USA: KS':
                country_define_list.append('USA: Kansas')
                latitude_define_list.append('39.0997266')
                longitude_define_list.append('-94.5785675')
                continent_define_list.append('North American')
            if row == 'USA: KY':
                country_define_list.append('USA: Kentucky')
                latitude_define_list.append('37.8393330')
                longitude_define_list.append('-84.2700180')
                continent_define_list.append('North American')
            if row == 'USA: Minnesota':
                country_define_list.append('USA: Minnesota')
                latitude_define_list.append('46.7295530')
                longitude_define_list.append('-94.6859000')
                continent_define_list.append('North American')
            if row == 'USA: MN':
                country_define_list.append('USA: Minnesota')
                latitude_define_list.append('46.7295530')
                longitude_define_list.append('-94.6859000')
                continent_define_list.append('North American')
            if row == 'USA: MO':
                country_define_list.append('USA: Missouri')
                latitude_define_list.append('37.9642530')
                longitude_define_list.append('-91.8318330')
                continent_define_list.append('North American')
            if row == 'USA: MS':
                country_define_list.append('USA: Mississippi')
                latitude_define_list.append('32.3546680')
                longitude_define_list.append('-89.3985280')
                continent_define_list.append('North American')
            if row == 'USA: NC':
                country_define_list.append('USA: North Carolina')
                latitude_define_list.append('35.7595730')
                longitude_define_list.append('-79.0193000')
                continent_define_list.append('North American')
            if row == 'USA: ND':
                country_define_list.append('USA: North Dakota')
                latitude_define_list.append('47.5514930')
                longitude_define_list.append('-101.0020120')
                continent_define_list.append('North American')
            if row == 'USA: NE':
                country_define_list.append('USA: Nebraska')
                latitude_define_list.append('41.4925370')
                longitude_define_list.append('-99.9018130')
                continent_define_list.append('North American')
            if row == 'USA: OH':
                country_define_list.append('USA: Ohio')
                latitude_define_list.append('40.4172870')
                longitude_define_list.append('-82.90712O0')
                continent_define_list.append('North American')
            if row == 'USA: OK':
                country_define_list.append('USA: Oklahoma')
                latitude_define_list.append('35.4675602')
                longitude_define_list.append('-97.5164261')
                continent_define_list.append('North American')
            if row == 'USA: OK and TX':
                country_define_list.append('USA: Oklahoma, Texas County')
                latitude_define_list.append('36.8301320')
                longitude_define_list.append('-101.4339150')
                continent_define_list.append('North American')
            if row == 'USA: PA':
                country_define_list.append('USA: Pennsylvania')
                latitude_define_list.append('41.2033220')
                longitude_define_list.append('-77.1945250')
                continent_define_list.append('North American')
            if row == 'USA: SC':
                country_define_list.append('USA: South Carolina')
                latitude_define_list.append('33.8360810')
                longitude_define_list.append('-81.1637250')
                continent_define_list.append('North American')
            if row == 'USA: SD':
                country_define_list.append('USA: South Dakota')
                latitude_define_list.append('43.9695150')
                longitude_define_list.append('-99.9018130')
                continent_define_list.append('North American')
            if row == 'USA: TN':
                country_define_list.append('USA: Tennessee')
                latitude_define_list.append('35.5174910')
                longitude_define_list.append('-86.5804470')
                continent_define_list.append('North American')
            if row == 'USA: TX':
                country_define_list.append('USA: Texas')
                latitude_define_list.append('31.9685990')
                longitude_define_list.append('-99.901813')
                continent_define_list.append('North American0')
            if row == 'USA: UT':
                country_define_list.append('USA: Utah')
                latitude_define_list.append('39.3209800')
                longitude_define_list.append('-111.0937310')
                continent_define_list.append('North American')
            if row == 'USA: VA':
                country_define_list.append('USA: Virginia')
                latitude_define_list.append('36.8529263')
                longitude_define_list.append('-78.656894')
                continent_define_list.append('North American')
            if row == 'USA: WI':
                country_define_list.append('USA: Wisconsin')
                latitude_define_list.append('43.784440')
                longitude_define_list.append('-75.9779816')
                continent_define_list.append('North American')
            if row == 'USA: WY':
                country_define_list.append('USA: Wyoming, USA')
                latitude_define_list.append('43.0759680')
                longitude_define_list.append('-107.2902840')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Adams County, (Latitude, Longitude = 39.932888, 91.165397)':
                country_define_list.append('USA: Illinois, Adams County')
                latitude_define_list.append('39.9328880')
                longitude_define_list.append('91.1653970')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Bureau County, (Latitude, Longitude = 41.193264, 89.608245)':
                country_define_list.append('USA: Illinois, Bureau County')
                latitude_define_list.append('41.1932640')
                longitude_define_list.append('89.6082450')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Bureau County, (Latitude, Longitude = 41.297188, 89.750043)':
                country_define_list.append('USA: Illinois, Bureau County')
                latitude_define_list.append('41.2971880')
                longitude_define_list.append('89.7500430')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Calhoun County, (Latitude, Longitude = 39.336947, 90.7058)':
                country_define_list.append('USA: Illinois, Calhoun County')
                latitude_define_list.append('39.3369470')
                longitude_define_list.append('90.7058000')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Cass County, (Latitude, Longitude = 40.015423, 90.445059)':
                country_define_list.append('USA: Illinois, Cass County')
                latitude_define_list.append('40.0154230')
                longitude_define_list.append('90.4450590')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, DeKalb County, (Latitude, Longitude = 42.018187, 88.911136)':
                country_define_list.append('USA: Illinois, DeKalb County')
                latitude_define_list.append('42.0181870')
                longitude_define_list.append('88.9111360')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, eastern Iowa (Latitude, Longitude = 40.961538, 91.272947)':
                country_define_list.append('USA: Illinois, eastern Iowa')
                latitude_define_list.append('40.9615380')
                longitude_define_list.append('91.2729470')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, eastern Iowa (Latitude, Longitude = 41.686447, 90.76462)':
                country_define_list.append('USA: Illinois, eastern Iowa')
                latitude_define_list.append('41.6864470')
                longitude_define_list.append('90.764620')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Edwards County, (Latitude, Longitude = 38.292699, 88.080045)':
                country_define_list.append('USA: Illinois, Edwards County')
                latitude_define_list.append('38.2926990')
                longitude_define_list.append('88.0800450')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Effingham County, (Latitude, Longitude = 39.00111, 88.494312)':
                country_define_list.append('USA: Illinois, Effingham County')
                latitude_define_list.append('39.0011100')
                longitude_define_list.append('88.4943120')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Effingham County, (Latitude, Longitude = 39.122349, 88.434839)':
                country_define_list.append('USA: Illinois, Effingham County')
                latitude_define_list.append('39.1223490')
                longitude_define_list.append('88.4348390')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Fayette County, (Latitude, Longitude = 39.056735, 88.956398)':
                country_define_list.append('USA: Illinois, Fayette County')
                latitude_define_list.append('39.0567350')
                longitude_define_list.append('88.9563980')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Hancock County, (Latitude, Longitude = 40.355539, 91.020382)':
                country_define_list.append('USA: Illinois, Hancock County')
                latitude_define_list.append('40.3555390')
                longitude_define_list.append('91.0203820')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.152444, 90.080107)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.1524440')
                longitude_define_list.append('90.0801070')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.375614, 90.095164)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.3756140')
                longitude_define_list.append('90.0951640')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.387142, 90.32882)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.3871420')
                longitude_define_list.append('90.3288200')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.415749, 90.392462)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.4157490')
                longitude_define_list.append('90.3924620')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.49018, 90.01181)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.4901800')
                longitude_define_list.append('90.0118100')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.017453, 90.174874)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0174530')
                longitude_define_list.append('90.1748740')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.034799, 90.174874)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0347990')
                longitude_define_list.append('90.1748740')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.063306, 89.988752)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0633060')
                longitude_define_list.append('89.9887520')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.068925, 90.11095)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0689250')
                longitude_define_list.append('90.1109500')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Livingston County, (Latitude, Longitude = 40.87108, 88.900607)':
                country_define_list.append('USA: Illinois, Livingston County')
                latitude_define_list.append('40.87108000')
                longitude_define_list.append('88.90060700')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Livingston County, (Latitude, Longitude = 40.933942, 88.490876)':
                country_define_list.append('USA: Illinois, Livingston County')
                latitude_define_list.append('40.9339420')
                longitude_define_list.append('88.4908760')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Logan County, (Latitude, Longitude = 40.063501, 89.285535)':
                country_define_list.append('USA: Illinois, Logan County')
                latitude_define_list.append('40.0635010')
                longitude_define_list.append('89.2855350')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Logan County, (Latitude, Longitude = 40.138741, 89.557122)':
                country_define_list.append('USA: Illin0ois, Logan County')
                latitude_define_list.append('40.138741')
                longitude_define_list.append('89.5571220')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mason County, (Latitude, Longitude = 40.301321, 89.970103)':
                country_define_list.append('USA: Illinois, Mason County')
                latitude_define_list.append('40.301321')
                longitude_define_list.append('89.9701030')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, McHenry County, (Latitude, Longitude = 42.291314, 88.538917)':
                country_define_list.append('USA: Illinois, McHenry County')
                latitude_define_list.append('42.2913140')
                longitude_define_list.append('88.5389170')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, McLean County, (Latitude, Longitude = 40.676635, 88.595756)':
                country_define_list.append('USA: Illinois, McLean County')
                latitude_define_list.append('40.676635')
                longitude_define_list.append('88.5957560')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, McLean County, (Latitude, Longitude = 40.716803, 88.882006)':
                country_define_list.append('USA: Illinois, McLean County')
                latitude_define_list.append('40.7168030')
                longitude_define_list.append('88.8820060')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Menard County, (Latitude, Longitude = 40.08722, 89.710619)':
                country_define_list.append('USA: Illinois, Menard County')
                latitude_define_list.append('40.0872200')
                longitude_define_list.append('89.710619')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Menard County, (Latitude, Longitude = 40.093363, 89.760698)':
                country_define_list.append('USA: Illinois, Menard County')
                latitude_define_list.append('40.0933630')
                longitude_define_list.append('89.7606980')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mercer County, (Latitude, Longitude = 41.06919, 90.743887)':
                country_define_list.append('USA: Illinois, Mercer County')
                latitude_define_list.append('41.0691900')
                longitude_define_list.append('90.7438870')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mercer County, (Latitude, Longitude = 41.085009, 90.496348)':
                country_define_list.append('USA: Illinois, Mercer County')
                latitude_define_list.append('41.0850090')
                longitude_define_list.append('90.4963480')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mercer County, (Latitude, Longitude = 41.106902, 90.730675)':
                country_define_list.append('USA: Illinois, Mercer County')
                latitude_define_list.append('41.1069020')
                longitude_define_list.append('90.7306750')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Ogle County, (Latitude, Longitude = 41.983211, 89.245471)':
                country_define_list.append('USA: Illinois, Ogle County')
                latitude_define_list.append('41.9832110')
                longitude_define_list.append('89.2454710')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Schuyler County, (Latitude, Longitude = 40.244924, 90.668438)':
                country_define_list.append('USA: Illinois, Schuyler County')
                latitude_define_list.append('40.2449240')
                longitude_define_list.append('90.6684380')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Scott County, (Latitude, Longitude = 39.557258, 90.5532)':
                country_define_list.append('USA: Illinois, Scott County')
                latitude_define_list.append('39.5572580')
                longitude_define_list.append('90.5532000')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Shelby County, (Latitude, Longitude = 39.274053, 88.655896)':
                country_define_list.append('USA: Illinois, Shelby County')
                latitude_define_list.append('39.2740530')
                longitude_define_list.append('88.6558960')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Shelby County, (Latitude, Longitude = 41.034799, 90.174874)':
                country_define_list.append('USA: Illinois, Shelby County')
                latitude_define_list.append('41.0347990')
                longitude_define_list.append('90.1748740')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, St. Clair County, (Latiude, Longitude = 38.456916, 89.928058)':
                country_define_list.append('USA: Illinois, St. Clair County')
                latitude_define_list.append('38.4569160')
                longitude_define_list.append('89.9280580')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Warren County, (Latitude, Longitude = 40.747517, 90.570944)':
                country_define_list.append('USA: Illinois, Warren County')
                latitude_define_list.append('40.7475170')
                longitude_define_list.append('90.5709440')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Warren County, (Latitude, Longitude = 40.79257, 90.505294)':
                country_define_list.append('USA: Illinois, Warren County')
                latitude_define_list.append('40.7925700')
                longitude_define_list.append('90.5052940')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Warren County, (Latitude, Longitude = 40.882615, 90.465536)':
                country_define_list.append('USA: Illinois, Warren County')
                latitude_define_list.append('40.8826150')
                longitude_define_list.append('90.4655360')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Whiteside County, (Latitude, Longitude = 41.664499, 89.912534)':
                country_define_list.append('USA: Illinois, Whiteside County')
                latitude_define_list.append('41.6644990')
                longitude_define_list.append('89.9125340')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Whiteside County, (Latitude, Longitude = 41.715868, 89.966289)':
                country_define_list.append('USA: Illinois, Whiteside County')
                latitude_define_list.append('41.7158680')
                longitude_define_list.append('89.966289')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Woodford County, (Latitude, Longitude = 40.765446, 89.352716)':
                country_define_list.append('USA: Illinois, Woodford County')
                latitude_define_list.append('40.7654460')
                longitude_define_list.append('89.3527160')
                continent_define_list.append('North American')
            if row == 'Viet Nam':
                country_define_list.append('Viet Nam')
                latitude_define_list.append('10.7500000')
                longitude_define_list.append('106.6666641')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: An Giang':
                country_define_list.append('Viet Nam: An Giang')
                latitude_define_list.append('10.5215840')
                longitude_define_list.append('105.1258960')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Bac Lieu province':
                country_define_list.append('Viet Nam: Bac Lieu province')
                latitude_define_list.append('9.2850000')
                longitude_define_list.append('105.7244415')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Bac Ninh':
                country_define_list.append('Viet Nam: Bac Ninh')
                latitude_define_list.append('21.1833333')
                longitude_define_list.append('106.0500031')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Binh Duong':
                country_define_list.append('Viet Nam: Binh Duong')
                latitude_define_list.append('21.3833333')
                longitude_define_list.append('103.0166702')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Dien Bien':
                country_define_list.append('Viet Nam: Dien Bien')
                latitude_define_list.append('21.3833333')
                longitude_define_list.append('103.0166702')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Dong Nai':
                country_define_list.append('Viet Nam: Dong Nai')
                latitude_define_list.append('11.0686300')
                longitude_define_list.append('107.1675980')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hai Duong':
                country_define_list.append('Viet Nam: Hai Duong')
                latitude_define_list.append('20.9333333')
                longitude_define_list.append('106.3166656')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hai Phong':
                country_define_list.append('Viet Nam: Hai Phong')
                latitude_define_list.append('20.8449110')
                longitude_define_list.append('106.6880840')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hanoi':
                country_define_list.append('Viet Nam: Hanoi')
                latitude_define_list.append('21.0277640')
                longitude_define_list.append('105.8341600')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Ho Chi Minh':
                country_define_list.append('Viet Nam: Ho Chi Minh')
                latitude_define_list.append('10.75')
                longitude_define_list.append('106.6666641')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hung Yen':
                country_define_list.append('Viet Nam: Hung Yen')
                latitude_define_list.append('20.8525710')
                longitude_define_list.append('106.0169970')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hung Yen province':
                country_define_list.append('Viet Nam: Hung Yen province')
                latitude_define_list.append('20.8525710')
                longitude_define_list.append('106.0169970')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Lao Cai':
                country_define_list.append('Viet Nam: Lao Cai')
                latitude_define_list.append('22.48333330')
                longitude_define_list.append('103.9499969')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Nghe An':
                country_define_list.append('Viet Nam: Nghe An')
                latitude_define_list.append('19.2342490')
                longitude_define_list.append('104.9200360')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Nghe An province':
                country_define_list.append('Viet Nam: Nghe An province')
                latitude_define_list.append('19.2342490')
                longitude_define_list.append('104.9200360')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Quang Ninh province':
                country_define_list.append('Viet Nam: Quang Ninh province')
                latitude_define_list.append('21.0063820')
                longitude_define_list.append('107.2925140')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Son La':
                country_define_list.append('Viet Nam: Son La')
                latitude_define_list.append('21.3166667')
                longitude_define_list.append('103.728917')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Tay Ninh':
                country_define_list.append('Viet Nam: Tay Ninh')
                latitude_define_list.append('11.335155')
                longitude_define_list.append('103.9000015')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Thai Binh':
                country_define_list.append('Viet Nam: Thai Binh')
                latitude_define_list.append('20.5386940')
                longitude_define_list.append('20.5386940')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Tien Giang province':
                country_define_list.append('Viet Nam: Tien Giang province')
                latitude_define_list.append('10.4493320')
                longitude_define_list.append('106.3420500')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Yen Bai':
                country_define_list.append('Viet Nam: Yen Bai')
                latitude_define_list.append('21.683792')
                longitude_define_list.append('21.7000000')
                continent_define_list.append('Asia')
            if row == 'Viet Nam:Ho Chi Minh':
                country_define_list.append('Viet Nam:Ho Chi Minh')
                latitude_define_list.append('10.7500000')
                longitude_define_list.append('106.6666641')
                continent_define_list.append('Asia')
            if row == 'Viet Nam:Nam Dinh':
                country_define_list.append('Viet Nam:Nam Dinh')
                latitude_define_list.append('20.4388230')
                longitude_define_list.append('106.1621050')
                continent_define_list.append('Asia')
            if row == 'Viet Nam:Thai Binh':
                country_define_list.append('Viet Nam: Thai Binh')
                latitude_define_list.append('20.4500000')
                longitude_define_list.append('106.393478')
                continent_define_list.append('Asia')
        #extract_zip = zip(country_list,latitude_list,longitude_list,continent_list)
        #print(extract_zip)
        #extract_list = []
        #for extract in extract_zip:
        #    extract_list.append("\t".join(extract))
        #country_string = '\n'.join(str(country_value) for country_value in extract_list)
        #print(country_string)
        #handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_country_define.txt", "w")
        #handle_outputfile.write(str(country_string))
        #handle_outputfile.close()
        # 3 ##################################################
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk", "genbank"))
        # pummedid
        accession_list = []
        organism_list = []
        gene_list = []
        length_list = []
        pubmedid_list = []
        publishdate_list = []
        note_list = []
        genotype_list = []
        country_list = []
        collectiondate_list = []
        origin_list = []
        for record in handle_parse_db:
            if record.annotations["organism"] == 'Porcine reproductive and respiratory syndrome virus':
                if len(record) == 600 or len(record) == 603 or len(record) == 606:
                    if record.features:
                        for feature in record.features:
                            if feature.type == "source":
                                accession_list.append(record.id)
                                organism_list.append(record.annotations["organism"])
                                gene_list.append('ORF5')
                                length_list.append(str(len(record)))
                                pubmedid_list.append(record.annotations["references"][0].pubmed_id)
                                publishdate_list.append(record.annotations["date"])
                                note_list.append(feature.qualifiers.get("note"))
                                genotype_list.append('European')
                                country_list.append(feature.qualifiers.get("country"))
                                collectiondate_list.append(feature.qualifiers.get("collection_date"))
                                origin_list.append(str(record.seq))
        # pubmedid
        organism_list_with_none = [['None'] if organism_value is None else organism_value for organism_value in organism_list]
        organism_none_list = []
        for organism_none in organism_list_with_none:
            organism_none_list.append(organism_none)
        publishdate_list_with_none = [['None'] if publishdate_value is None else publishdate_value for publishdate_value in publishdate_list]
        publishdate_none_list = []
        for publishdate_none in publishdate_list_with_none:
            publishdate_none_list.append(publishdate_none)
        note_list_with_none = [['None'] if note_value is None else note_value for note_value in note_list]
        note_none_list = []
        for note_none in note_list_with_none:
            note_none_list.append('\n'.join(note_none))
        country_list_with_none = [['None'] if country_value is None else country_value for country_value in country_list]
        country_none_list = []
        for country_none in country_list_with_none:
            country_none_list.append('\n'.join(country_none))
        collectiondate_list_with_none = [['None'] if collectiondate_value is None else collectiondate_value for collectiondate_value in collectiondate_list]
        collectiondate_none_list = []
        for collectiondate_none in collectiondate_list_with_none:
            collectiondate_none_list.append('\n'.join(collectiondate_none))
        origin_list_with_none = [['None'] if origin_value is None else origin_value for origin_value in origin_list]
        origin_none_list = []
        for origin_none in origin_list_with_none:
            origin_none_list.append(origin_none)
        pubmedid_zip = zip(accession_list,organism_none_list,gene_list,length_list,collectiondate_none_list,pubmedid_list,publishdate_none_list,note_none_list,genotype_list,country_none_list,country_define_list,latitude_define_list,longitude_define_list,continent_define_list,origin_none_list)
        print(pubmedid_zip)
        orf5_pubmedid_list = []
        for orf5_pubmedid in pubmedid_zip:
            orf5_pubmedid_list.append("\t".join(orf5_pubmedid))
        # print(ORF5_list)
        print("")
        # Writing ORF5
        print("Step 2: Writing ORF5 EU, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        orf5_pubmedid_string = '\n'.join(str(orf5_value) for orf5_value in orf5_pubmedid_list)
        print(orf5_pubmedid_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_orf5_eu_extract.txt"), "w")
        handle_outputfile.write(orf5_pubmedid_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 45: ORF5 EU Extract Complete")
        print("")
        break

# ORF5 NA Extract ##################################################
def orf5naextract():
    while True:
        # Module Name
        print("Module 46: ORF5 NA Extract")
        print("")
        # Parsing ORF5 NA
        print("Step 1: Parsing ORF5 NA, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        # 1 #################################################
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        country_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        country = feature.qualifiers.get('country')
                        country_list_with_none.append(country)
        country_list_with_none_list = [['None'] if country_value is None else country_value for country_value in country_list_with_none]
        country_extract_list = []
        for country_none in country_list_with_none_list:
            country_extract_list.append('\n'.join(country_none))
        print(country_extract_list)
        country_string = "\n".join(country_extract_list)
        print(country_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_country.txt"), "w")
        handle_outputfile.write(country_string)
        handle_outputfile.close()
        # 2 ##################################################
        handle_inputfile = open(os.path.expanduser("~/extractor/local/export/source_country.txt"))
        country_define_list = []
        latitude_define_list = []
        longitude_define_list = []
        continent_define_list = []
        for row in handle_inputfile:
            row = row.rstrip()
            if row.startswith('Albania'):
                country_define_list.append('Albania')
                latitude_define_list.append('41.3275000')
                longitude_define_list.append('19.8188896')
                continent_define_list.append('Europe')
            if row.startswith('Australia'):
                country_define_list.append('Australia')
                latitude_define_list.append('-33.8678500')
                longitude_define_list.append('133.7751360')
                continent_define_list.append('Australia/Oceania')
            if row.startswith('Austria'):
                country_define_list.append('Austria')
                latitude_define_list.append('48.2084878')
                longitude_define_list.append('16.3720760')
                continent_define_list.append('Europe')
            if row.startswith('Belarus'):
                country_define_list.append('Belarus')
                latitude_define_list.append('53.9000000')
                longitude_define_list.append('27.5666676')
                continent_define_list.append('Europe')
            if row.startswith('Bulgaria'):
                country_define_list.append('Bulgaria')
                latitude_define_list.append('42.6975135')
                longitude_define_list.append('23.3241463')
                continent_define_list.append('Europe')
            if row.startswith('Cambodia'):
                country_define_list.append('Cambodia')
                latitude_define_list.append('11.5500000')
                longitude_define_list.append('104.9166641')
                continent_define_list.append('Asia')
            if row == 'Canada':
                country_define_list.append('Canada')
                latitude_define_list.append('43.7001138')
                longitude_define_list.append('-79.4163055')
                continent_define_list.append('North American')
            if row == 'Canada: ON':
                country_define_list.append('Canada: Ontario')
                latitude_define_list.append('51.253775')
                longitude_define_list.append('-85.323214')
                continent_define_list.append('North American')
            if row == 'Canada: QC':
                country_define_list.append('Canada: Quebec')
                latitude_define_list.append('46.8122791')
                longitude_define_list.append('-71.2145386')
                continent_define_list.append('North American')
            if row == 'Canada: SK':
                country_define_list.append('Canada: Saskatchewan')
                latitude_define_list.append('52.939916')
                longitude_define_list.append('-106.450864')
                continent_define_list.append('North American')
            if row == 'Chile':
                country_define_list.append('Chile')
                latitude_define_list.append('-33.4500000')
                longitude_define_list.append('-70.6666641')
                continent_define_list.append('South American')
            if row == 'Chile: Maipo':
                country_define_list.append('Chile: Maipo River')
                latitude_define_list.append('-33.8683360')
                longitude_define_list.append('-70.8111840')
                continent_define_list.append('South American')
            if row == 'China':
                country_define_list.append('China')
                latitude_define_list.append('31.2222222')
                longitude_define_list.append('121.4580536')
                continent_define_list.append('Asia')
            if row == 'China: Anhui':
                country_define_list.append('China: Anhui')
                latitude_define_list.append('30.6006770')
                longitude_define_list.append('117.9249000')
                continent_define_list.append('Asia')
            if row == 'China: Boxing City':
                country_define_list.append('China: Boxing City')
                latitude_define_list.append('37.1502260')
                longitude_define_list.append('118.1318150')
                continent_define_list.append('Asia')
            if row == 'China: Chongqing':
                country_define_list.append('China: Chongqing')
                latitude_define_list.append('29.5627778')
                longitude_define_list.append('106.5527802')
                continent_define_list.append('Asia')
            if row == 'China: Dezhou City':
                country_define_list.append('China: Dezhou City')
                latitude_define_list.append('37.4340920')
                longitude_define_list.append('116.3574640')
                continent_define_list.append('Asia')
            if row == 'China: Fujian, Fuzhou':
                country_define_list.append('China: Fujian, Fuzhou')
                latitude_define_list.append('26.0613889')
                longitude_define_list.append('119.3061142')
                continent_define_list.append('Asia')
            if row == 'China: Fujian, Nanping':
                country_define_list.append('China: Fujian, Nanping')
                latitude_define_list.append('26.6417740')
                longitude_define_list.append('118.1777100')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong':
                country_define_list.append('China: Guangdong')
                latitude_define_list.append('23.3790330')
                longitude_define_list.append('113.7632830')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Foshan':
                country_define_list.append('China: Guangdong, Foshan')
                latitude_define_list.append('23.0333333')
                longitude_define_list.append('113.1166687')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Guangzhou':
                country_define_list.append('China: Guangdong, Guangzhou')
                latitude_define_list.append('23.1166667')
                longitude_define_list.append('113.2500000')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Guangzhou, Huadu':
                country_define_list.append('China: Guangdong, Guangzhou, Huadu')
                latitude_define_list.append('23.4041650')
                longitude_define_list.append('113.2202180')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Heyuan':
                country_define_list.append('China: Guangdong, Heyuan')
                latitude_define_list.append('23.7436850')
                longitude_define_list.append('114.7009610')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Huizhou':
                country_define_list.append('China: Guangdong, Huizhou')
                latitude_define_list.append('23.1122570')
                longitude_define_list.append('114.4158010')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Jiangmen':
                country_define_list.append('China: Guangdong, Jiangmen')
                latitude_define_list.append('22.5791170')
                longitude_define_list.append('113.0815080')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Jiangmeng':
                country_define_list.append('China: Guangdong, Jiangmeng')
                latitude_define_list.append('22.5791170')
                longitude_define_list.append('113.0815080')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Jieyang':
                country_define_list.append('China: Guangdong, Jieyang')
                latitude_define_list.append('23.5297222')
                longitude_define_list.append('116.3655548')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Kaiping':
                country_define_list.append('China: Guangdong, Kaiping')
                latitude_define_list.append('22.3763950')
                longitude_define_list.append('112.6985450')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Maoming':
                country_define_list.append('China: Guangdong, Maoming')
                latitude_define_list.append('21.6629910')
                longitude_define_list.append('110.9254390')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Meizhou':
                country_define_list.append('China: Guangdong, Meizhou')
                latitude_define_list.append('24.2885780')
                longitude_define_list.append('116.1225230')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Qingyuan':
                country_define_list.append('China: Guangdong, Qingyuan')
                latitude_define_list.append('23.6817740')
                longitude_define_list.append('113.0560420')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Sanshui':
                country_define_list.append('China: Guangdong, Sanshui')
                latitude_define_list.append('23.1560450')
                longitude_define_list.append('112.8966060')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shantou':
                country_define_list.append('China: Guangdong, Shantou')
                latitude_define_list.append('23.3600000')
                longitude_define_list.append('116.6783371')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shanwei':
                country_define_list.append('China: Guangdong, Shanwei')
                latitude_define_list.append('22.7861860')
                longitude_define_list.append('115.3751580')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shaoguan':
                country_define_list.append('China: Guangdong, Shaoguan')
                latitude_define_list.append('24.8000000')
                longitude_define_list.append('113.5833359')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shenzhen':
                country_define_list.append('China: Guangdong, Shenzhen')
                latitude_define_list.append('22.5455377')
                longitude_define_list.append('114.0682983')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Shiling':
                country_define_list.append('China: Guangdong, Shiling')
                latitude_define_list.append('23.4607280')
                longitude_define_list.append('113.1515960')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Sihui':
                country_define_list.append('China: Guangdong, Sihui')
                latitude_define_list.append('23.3270010')
                longitude_define_list.append('112.7341030')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Yangchun':
                country_define_list.append('China: Guangdong, Yangchun')
                latitude_define_list.append('22.1704370')
                longitude_define_list.append('111.7915390')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Yangjiang':
                country_define_list.append('China: Guangdong, Yangjiang')
                latitude_define_list.append('21.8579580')
                longitude_define_list.append('111.9822320')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Yunfu':
                country_define_list.append('China: Guangdong, Yunfu')
                latitude_define_list.append('22.9333333')
                longitude_define_list.append('112.0333328')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Zhanjiang':
                country_define_list.append('China: Guangdong, Zhanjiang')
                latitude_define_list.append('21.2000000')
                longitude_define_list.append('110.3833313')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Zhaoqing':
                country_define_list.append('China: Guangdong, Zhaoqing')
                latitude_define_list.append('23.0471910')
                longitude_define_list.append('112.4650910')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong, Zhongshan':
                country_define_list.append('China: Guangdong, Zhongshan')
                latitude_define_list.append('22.5175850')
                longitude_define_list.append('113.3927700')
                continent_define_list.append('Asia')
            if row == 'China: Guangdong,Guangzhou':
                country_define_list.append('China: Guangdong, Guangzhou')
                latitude_define_list.append('23.1166667')
                longitude_define_list.append('113.2500000')
                continent_define_list.append('Asia')
            if row == 'China: Guangxi':
                country_define_list.append('China: Guangxi')
                latitude_define_list.append('23.7247600')
                longitude_define_list.append('108.8076190')
                continent_define_list.append('Asia')
            if row == 'China: Guodian City':
                country_define_list.append('China: Guodian City')
                latitude_define_list.append('36.0839480')
                longitude_define_list.append('103.6665500')
                continent_define_list.append('Asia')
            if row == 'China: Hainan':
                country_define_list.append('China: Hainan')
                latitude_define_list.append('19.5663950')
                longitude_define_list.append('109.9496860')
                continent_define_list.append('Asia')
            if row == 'China: Hebei Province':
                country_define_list.append('China: Hebei Province')
                latitude_define_list.append('37.8956590')
                longitude_define_list.append('114.9042210')
                continent_define_list.append('Asia')
            if row == 'China: Heilongjiang':
                country_define_list.append('China: Heilongjiang')
                latitude_define_list.append('47.1216470')
                longitude_define_list.append('128.7382310')
                continent_define_list.append('Asia')
            if row == 'China: Henan':
                country_define_list.append('China: Henan')
                latitude_define_list.append('34.2904300')
                longitude_define_list.append('113.3823550')
                continent_define_list.append('Asia')
            if row == 'China: Hennan,Yuoyang':
                country_define_list.append('China: Hennan,Yuoyang')
                latitude_define_list.append('34.6202020')
                longitude_define_list.append('112.4539260')
                continent_define_list.append('Asia')
            if row == 'China: Hnan,Luoyang':
                country_define_list.append('China: Hnan,Luoyang')
                latitude_define_list.append('34.6836111')
                longitude_define_list.append('112.4536133')
                continent_define_list.append('Asia')
            if row == 'China: Hunan':
                country_define_list.append('China: Hunan')
                latitude_define_list.append('27.6253000')
                longitude_define_list.append('111.8568590')
                continent_define_list.append('Asia')
            if row == 'China: Jiangxi':
                country_define_list.append('China: Jiangxi')
                latitude_define_list.append('27.0874560')
                longitude_define_list.append('114.9042210')
                continent_define_list.append('Asia')
            if row == 'China: Jilin, Changchun':
                country_define_list.append('China: Jilin, Changchun')
                latitude_define_list.append('43.8800000')
                longitude_define_list.append('125.3227768')
                continent_define_list.append('Asia')
            if row == 'China: Jinan City':
                country_define_list.append('China: Jinan City')
                latitude_define_list.append('36.6683333')
                longitude_define_list.append('116.9972229')
                continent_define_list.append('Asia')
            if row == 'China: Liaocheng City':
                country_define_list.append('China: Liaocheng City')
                latitude_define_list.append('36.4570300')
                longitude_define_list.append('115.9854600')
                continent_define_list.append('Asia')
            if row == 'China: SD Province':
                country_define_list.append('China: SD Province')
                latitude_define_list.append('35.8939570')
                longitude_define_list.append('117.9249000')
                continent_define_list.append('Asia')
            if row == 'China: Shandong':
                country_define_list.append('China: Shandong')
                latitude_define_list.append('35.8939570')
                longitude_define_list.append('117.9249000')
                continent_define_list.append('Asia')
            if row == 'China: Sichuan':
                country_define_list.append('China: Sichuan')
                latitude_define_list.append('30.6516520')
                longitude_define_list.append('104.0759310')
                continent_define_list.append('Asia')
            if row == 'China: Wenzu City':
                country_define_list.append('China: Wenzu City')
                latitude_define_list.append('27.9938280')
                longitude_define_list.append('120.6993610')
                continent_define_list.append('Asia')
            if row == 'China: Wuzhou':
                country_define_list.append('China: Wuzhou')
                latitude_define_list.append('23.4769620')
                longitude_define_list.append('111.2791150')
                continent_define_list.append('Asia')
            if row == 'China: Xinjiang':
                country_define_list.append('China: Xinjiang')
                latitude_define_list.append('42.5246360')
                longitude_define_list.append('87.5395850')
                continent_define_list.append('Asia')
            if row == 'China: Yun Nan':
                country_define_list.append('China: Yun Nan')
                latitude_define_list.append('24.4752850')
                longitude_define_list.append('101.3431060')
                continent_define_list.append('Asia')
            if row == 'China: Yunnan':
                country_define_list.append('China: Yunnan')
                latitude_define_list.append('24.4752850')
                longitude_define_list.append('101.3431060')
                continent_define_list.append('Asia')
            if row == 'China: Yunnan, Qujing':
                country_define_list.append('China: Yunnan, Qujing')
                latitude_define_list.append('25.4899990')
                longitude_define_list.append('103.7961670')
                continent_define_list.append('Asia')
            if row == 'China: Zhangqiu City':
                country_define_list.append('China: Zhangqiu City')
                latitude_define_list.append('36.6812590')
                longitude_define_list.append('117.5262280')
                continent_define_list.append('Asia')
            if row.startswith('Croatia'):
                country_define_list.append('Croatia')
                latitude_define_list.append('45.8131276')
                longitude_define_list.append('15.9775333')
                continent_define_list.append('Europe')
            if row.startswith('Czech Republic'):
                country_define_list.append('Czech Republic')
                latitude_define_list.append('50.0878368')
                longitude_define_list.append('14.4241323')
                continent_define_list.append('Europe')
            if row.startswith('Denmark'):
                country_define_list.append('Denmark')
                latitude_define_list.append('55.6776812')
                longitude_define_list.append('12.5709343')
                continent_define_list.append('Europe')
            if row.startswith('France'):
                country_define_list.append('France')
                latitude_define_list.append('48.8534100')
                longitude_define_list.append('2.3487999')
                continent_define_list.append('Europe')
            if row.startswith('Germany'):
                country_define_list.append('Germany')
                latitude_define_list.append('52.5166667')
                longitude_define_list.append('13.3999996')
                continent_define_list.append('Europe')
            if row.startswith('Hungary'):
                country_define_list.append('Hungary')
                latitude_define_list.append('47.5000000')
                longitude_define_list.append('19.0833340')
                continent_define_list.append('Europe')
            if row.startswith('India'):
                country_define_list.append('India')
                latitude_define_list.append('19.0144100')
                longitude_define_list.append('72.8479385')
                continent_define_list.append('Asia')
            if row.startswith('Italy'):
                country_define_list.append('Italy')
                latitude_define_list.append('41.9000000')
                longitude_define_list.append('12.4833336')
                continent_define_list.append('Europe')
            if row == 'Japan':
                country_define_list.append('Japan')
                latitude_define_list.append('35.6895266')
                longitude_define_list.append('139.6916809')
                continent_define_list.append('Asia')
            if row == 'Japan:Aichi':
                country_define_list.append('Japan: Aichi')
                latitude_define_list.append('34.9666667')
                longitude_define_list.append('136.6166687')
                continent_define_list.append('Asia')
            if row == 'Japan:Aomori':
                country_define_list.append('Japan: Aomori')
                latitude_define_list.append('40.8211111')
                longitude_define_list.append('140.7511139')
                continent_define_list.append('Asia')
            if row == 'Japan:Chiba':
                country_define_list.append('Japan: Chiba')
                latitude_define_list.append('35.6000000')
                longitude_define_list.append('140.1166687')
                continent_define_list.append('Asia')
            if row == 'Japan:Ehime':
                country_define_list.append('Japan: Ehime')
                latitude_define_list.append('33.8416240')
                longitude_define_list.append('132.7656810')
                continent_define_list.append('Asia')
            if row == 'Japan:Fukushima':
                country_define_list.append('Japan: Fukushima')
                latitude_define_list.append('37.7500000')
                longitude_define_list.append('140.4666595')
                continent_define_list.append('Asia')
            if row == 'Japan:Gunma':
                country_define_list.append('Japan: Gunma')
                latitude_define_list.append('36.3906670')
                longitude_define_list.append('139.0604060')
                continent_define_list.append('Asia')
            if row == 'Japan:Hokkaido':
                country_define_list.append('Japan: Hokkaido')
                latitude_define_list.append('43.2203270')
                longitude_define_list.append('142.8634740')
                continent_define_list.append('Asia')
            if row == 'Japan:Ibaraki':
                country_define_list.append('Japan: Ibaraki')
                latitude_define_list.append('34.8164106')
                longitude_define_list.append('135.5682831')
                continent_define_list.append('Asia')
            if row == 'Japan:Ishikawa':
                country_define_list.append('Japan: Ishikawa')
                latitude_define_list.append('36.5946820')
                longitude_define_list.append('136.6255730')
                continent_define_list.append('Asia')
            if row == 'Japan:Iwate':
                country_define_list.append('Japan: Iwate')
                latitude_define_list.append('39.7036190')
                longitude_define_list.append('141.1526840')
                continent_define_list.append('Asia')
            if row == 'Japan:Kagoshima':
                country_define_list.append('Japan: Kagoshima')
                latitude_define_list.append('31.6000000')
                longitude_define_list.append('130.5500031')
                continent_define_list.append('Asia')
            if row == 'Japan:Kyoto':
                country_define_list.append('Japan: Kyoto')
                latitude_define_list.append('35.0210700')
                longitude_define_list.append('135.7538452')
                continent_define_list.append('Asia')
            if row == 'Japan:Nagasaki':
                country_define_list.append('Japan: Nagasaki')
                latitude_define_list.append('32.7550000')
                longitude_define_list.append('129.8683319')
                continent_define_list.append('Asia')
            if row == 'Japan:Niigata':
                country_define_list.append('Japan: Niigata')
                latitude_define_list.append('37.9166667')
                longitude_define_list.append('139.0500031')
                continent_define_list.append('Asia')
            if row == 'Japan:Osaka':
                country_define_list.append('Japan: Osaka')
                latitude_define_list.append('34.6937398')
                longitude_define_list.append('135.5021820')
                continent_define_list.append('Asia')
            if row == 'Japan:Saga':
                country_define_list.append('Japan: Saga')
                latitude_define_list.append('33.2634820')
                longitude_define_list.append('130.3008580')
                continent_define_list.append('Asia')
            if row == 'Japan:Saitama':
                country_define_list.append('Japan: Saitama')
                latitude_define_list.append('35.8617290')
                longitude_define_list.append('139.6454820')
                continent_define_list.append('Asia')
            if row == 'Japan:Shiga':
                country_define_list.append('Japan: Shiga')
                latitude_define_list.append('35.0045310')
                longitude_define_list.append('135.8685900')
                continent_define_list.append('Asia')
            if row == 'Japan:Shizuoka':
                country_define_list.append('Japan: Shizuoka')
                latitude_define_list.append('34.9666667')
                longitude_define_list.append('138.3833313')
                continent_define_list.append('Asia')
            if row == 'Japan:Tochigi':
                country_define_list.append('Japan: Tochigi')
                latitude_define_list.append('36.5657250')
                longitude_define_list.append('139.8835650')
                continent_define_list.append('Asia')
            if row == 'Japan:Yamagata':
                country_define_list.append('Japan: Yamagata')
                latitude_define_list.append('38.2527778')
                longitude_define_list.append('140.3374939')
                continent_define_list.append('Asia')
            if row == 'Japan:Yamaguchi':
                country_define_list.append('Japan: Yamaguchi')
                latitude_define_list.append('34.1859560')
                longitude_define_list.append('131.4706490')
                continent_define_list.append('Asia')
            if row == 'Japan:Yamanashi':
                country_define_list.append('Japan: Yamanashi')
                latitude_define_list.append('35.6641580')
                longitude_define_list.append('138.5684490')
                continent_define_list.append('Asia')
            if row.startswith('Laos'):
                country_define_list.append('Laos')
                latitude_define_list.append('19.856270')
                longitude_define_list.append('102.495496')
                continent_define_list.append('Asia')
            if row.startswith('Lithuania'):
                country_define_list.append('Lithuania')
                latitude_define_list.append('54.6833333')
                longitude_define_list.append('25.3166676')
                continent_define_list.append('Europe')
            if row == 'Malaysia':
                country_define_list.append('Malaysia')
                latitude_define_list.append('3.1666667')
                longitude_define_list.append('101.6999969')
                continent_define_list.append('Asia')
            if row == 'Malaysia: Sarawak':
                country_define_list.append('Malaysia: Sarawak')
                latitude_define_list.append('1.5532780')
                longitude_define_list.append('110.3592130')
                continent_define_list.append('Asia')
            if row == 'Malaysia: Selangor':
                country_define_list.append('Malaysia: Selangor')
                latitude_define_list.append('3.3500000')
                longitude_define_list.append('101.2500000')
                continent_define_list.append('Asia')
            if row == 'Mexico':
                country_define_list.append('Mexico')
                latitude_define_list.append('19.4341667')
                longitude_define_list.append('-99.1386108')
                continent_define_list.append('North American')
            if row == 'Mexico: Morelos':
                country_define_list.append('Mexico: Morelos')
                latitude_define_list.append('18.8000000')
                longitude_define_list.append('-98.9499969')
                continent_define_list.append('North American')
            if row.startswith('Myanmar'):
                country_define_list.append('Myanmar')
                latitude_define_list.append('21.9162210')
                longitude_define_list.append('95.9559740')
                continent_define_list.append('Asia')
            if row.startswith('Poland'):
                country_define_list.append('Poland')
                latitude_define_list.append('52.2500000')
                longitude_define_list.append('21.0000000')
                continent_define_list.append('Europe')
            if row.startswith('Romania'):
                country_define_list.append('Romania')
                latitude_define_list.append('44.4333333')
                longitude_define_list.append('26.1000004')
                continent_define_list.append('Europe')
            if row.startswith('Russia'):
                country_define_list.append('Russia')
                latitude_define_list.append('55.7522222')
                longitude_define_list.append('37.6155548')
                continent_define_list.append('Europe')
            if row.startswith('Serbia'):
                country_define_list.append('Serbia')
                latitude_define_list.append('44.7865680')
                longitude_define_list.append('20.4489220')
                continent_define_list.append('Europe')
            if row.startswith('Singapore'):
                country_define_list.append('Singapore')
                latitude_define_list.append('1.3520830')
                longitude_define_list.append('103.8198360')
                continent_define_list.append('Asia')
            if row.startswith('Slovakia'):
                country_define_list.append('Slovakia')
                latitude_define_list.append('48.1485960')
                longitude_define_list.append('17.1077480')
                continent_define_list.append('Europe')
            if row.startswith('Slovenia'):
                country_define_list.append('Slovenia')
                latitude_define_list.append('46.0552778')
                longitude_define_list.append('14.5144444')
                continent_define_list.append('Europe')
            if row.startswith('South Korea'):
                country_define_list.append('South Korea')
                latitude_define_list.append('37.5663889')
                longitude_define_list.append('126.9997253')
                continent_define_list.append('Asia')
            if row.startswith('Spain'):
                country_define_list.append('Spain')
                latitude_define_list.append('40.4165021')
                longitude_define_list.append('-3.7025642')
                continent_define_list.append('Europe')
            if row.startswith('Taiwan'):
                country_define_list.append('Taiwan')
                latitude_define_list.append('23.6978100')
                longitude_define_list.append('120.9605150')
                continent_define_list.append('Asia')
            if row.startswith('Thailand'):
                country_define_list.append('Thailand')
                latitude_define_list.append('13.7500000')
                longitude_define_list.append('100.5166702')
                continent_define_list.append('Asia')
            if row.startswith('United Kingdom'):
                country_define_list.append('United Kingdom')
                latitude_define_list.append('51.5084153')
                longitude_define_list.append('-0.1255327')
                continent_define_list.append('Europe')
            if row == 'USA':
                country_define_list.append('USA')
                latitude_define_list.append('40.7142691')
                longitude_define_list.append('-74.0059738')
                continent_define_list.append('North American')
            if row == 'USA: AZ':
                country_define_list.append('USA: Arizona')
                latitude_define_list.append('15.6333333')
                longitude_define_list.append('-87.3166656')
                continent_define_list.append('North American')
            if row == 'USA: Can':
                country_define_list.append('USA: Can Lane')
                latitude_define_list.append('29.1949620')
                longitude_define_list.append('-81.0073860')
                continent_define_list.append('North American')
            if row == 'USA: CAN':
                country_define_list.append('USA: Can Lane')
                latitude_define_list.append('29.1949620')
                longitude_define_list.append('-81.0073860')
                continent_define_list.append('North American')
            if row == 'USA: CO':
                country_define_list.append('USA: Colorado')
                latitude_define_list.append('38.8338816')
                longitude_define_list.append('-104.8213654')
                continent_define_list.append('North American')
            if row == 'USA: IA':
                country_define_list.append('USA: Iowa')
                latitude_define_list.append('41.8780030')
                longitude_define_list.append('-93.0977020')
                continent_define_list.append('North American')
            if row == 'USA: IL':
                country_define_list.append('USA: Illinois')
                latitude_define_list.append('40.6331250')
                longitude_define_list.append('-89.3985280')
                continent_define_list.append('North American')
            if row == 'USA: Illinois':
                country_define_list.append('USA: Illinois')
                latitude_define_list.append('40.6331250')
                longitude_define_list.append('-89.3985280')
                continent_define_list.append('North American')
            if row == 'USA: IN':
                country_define_list.append('USA: Indiana')
                latitude_define_list.append('39.7683765')
                longitude_define_list.append('-86.1580429')
                continent_define_list.append('North American')
            if row == 'USA: KS':
                country_define_list.append('USA: Kansas')
                latitude_define_list.append('39.0997266')
                longitude_define_list.append('-94.5785675')
                continent_define_list.append('North American')
            if row == 'USA: KY':
                country_define_list.append('USA: Kentucky')
                latitude_define_list.append('37.8393330')
                longitude_define_list.append('-84.2700180')
                continent_define_list.append('North American')
            if row == 'USA: Minnesota':
                country_define_list.append('USA: Minnesota')
                latitude_define_list.append('46.7295530')
                longitude_define_list.append('-94.6859000')
                continent_define_list.append('North American')
            if row == 'USA: MN':
                country_define_list.append('USA: Minnesota')
                latitude_define_list.append('46.7295530')
                longitude_define_list.append('-94.6859000')
                continent_define_list.append('North American')
            if row == 'USA: MO':
                country_define_list.append('USA: Missouri')
                latitude_define_list.append('37.9642530')
                longitude_define_list.append('-91.8318330')
                continent_define_list.append('North American')
            if row == 'USA: MS':
                country_define_list.append('USA: Mississippi')
                latitude_define_list.append('32.3546680')
                longitude_define_list.append('-89.3985280')
                continent_define_list.append('North American')
            if row == 'USA: NC':
                country_define_list.append('USA: North Carolina')
                latitude_define_list.append('35.7595730')
                longitude_define_list.append('-79.0193000')
                continent_define_list.append('North American')
            if row == 'USA: ND':
                country_define_list.append('USA: North Dakota')
                latitude_define_list.append('47.5514930')
                longitude_define_list.append('-101.0020120')
                continent_define_list.append('North American')
            if row == 'USA: NE':
                country_define_list.append('USA: Nebraska')
                latitude_define_list.append('41.4925370')
                longitude_define_list.append('-99.9018130')
                continent_define_list.append('North American')
            if row == 'USA: OH':
                country_define_list.append('USA: Ohio')
                latitude_define_list.append('40.4172870')
                longitude_define_list.append('-82.90712O0')
                continent_define_list.append('North American')
            if row == 'USA: OK':
                country_define_list.append('USA: Oklahoma')
                latitude_define_list.append('35.4675602')
                longitude_define_list.append('-97.5164261')
                continent_define_list.append('North American')
            if row == 'USA: OK and TX':
                country_define_list.append('USA: Oklahoma, Texas County')
                latitude_define_list.append('36.8301320')
                longitude_define_list.append('-101.4339150')
                continent_define_list.append('North American')
            if row == 'USA: PA':
                country_define_list.append('USA: Pennsylvania')
                latitude_define_list.append('41.2033220')
                longitude_define_list.append('-77.1945250')
                continent_define_list.append('North American')
            if row == 'USA: SC':
                country_define_list.append('USA: South Carolina')
                latitude_define_list.append('33.8360810')
                longitude_define_list.append('-81.1637250')
                continent_define_list.append('North American')
            if row == 'USA: SD':
                country_define_list.append('USA: South Dakota')
                latitude_define_list.append('43.9695150')
                longitude_define_list.append('-99.9018130')
                continent_define_list.append('North American')
            if row == 'USA: TN':
                country_define_list.append('USA: Tennessee')
                latitude_define_list.append('35.5174910')
                longitude_define_list.append('-86.5804470')
                continent_define_list.append('North American')
            if row == 'USA: TX':
                country_define_list.append('USA: Texas')
                latitude_define_list.append('31.9685990')
                longitude_define_list.append('-99.901813')
                continent_define_list.append('North American0')
            if row == 'USA: UT':
                country_define_list.append('USA: Utah')
                latitude_define_list.append('39.3209800')
                longitude_define_list.append('-111.0937310')
                continent_define_list.append('North American')
            if row == 'USA: VA':
                country_define_list.append('USA: Virginia')
                latitude_define_list.append('36.8529263')
                longitude_define_list.append('-78.656894')
                continent_define_list.append('North American')
            if row == 'USA: WI':
                country_define_list.append('USA: Wisconsin')
                latitude_define_list.append('43.784440')
                longitude_define_list.append('-75.9779816')
                continent_define_list.append('North American')
            if row == 'USA: WY':
                country_define_list.append('USA: Wyoming, USA')
                latitude_define_list.append('43.0759680')
                longitude_define_list.append('-107.2902840')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Adams County, (Latitude, Longitude = 39.932888, 91.165397)':
                country_define_list.append('USA: Illinois, Adams County')
                latitude_define_list.append('39.9328880')
                longitude_define_list.append('91.1653970')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Bureau County, (Latitude, Longitude = 41.193264, 89.608245)':
                country_define_list.append('USA: Illinois, Bureau County')
                latitude_define_list.append('41.1932640')
                longitude_define_list.append('89.6082450')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Bureau County, (Latitude, Longitude = 41.297188, 89.750043)':
                country_define_list.append('USA: Illinois, Bureau County')
                latitude_define_list.append('41.2971880')
                longitude_define_list.append('89.7500430')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Calhoun County, (Latitude, Longitude = 39.336947, 90.7058)':
                country_define_list.append('USA: Illinois, Calhoun County')
                latitude_define_list.append('39.3369470')
                longitude_define_list.append('90.7058000')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Cass County, (Latitude, Longitude = 40.015423, 90.445059)':
                country_define_list.append('USA: Illinois, Cass County')
                latitude_define_list.append('40.0154230')
                longitude_define_list.append('90.4450590')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, DeKalb County, (Latitude, Longitude = 42.018187, 88.911136)':
                country_define_list.append('USA: Illinois, DeKalb County')
                latitude_define_list.append('42.0181870')
                longitude_define_list.append('88.9111360')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, eastern Iowa (Latitude, Longitude = 40.961538, 91.272947)':
                country_define_list.append('USA: Illinois, eastern Iowa')
                latitude_define_list.append('40.9615380')
                longitude_define_list.append('91.2729470')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, eastern Iowa (Latitude, Longitude = 41.686447, 90.76462)':
                country_define_list.append('USA: Illinois, eastern Iowa')
                latitude_define_list.append('41.6864470')
                longitude_define_list.append('90.764620')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Edwards County, (Latitude, Longitude = 38.292699, 88.080045)':
                country_define_list.append('USA: Illinois, Edwards County')
                latitude_define_list.append('38.2926990')
                longitude_define_list.append('88.0800450')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Effingham County, (Latitude, Longitude = 39.00111, 88.494312)':
                country_define_list.append('USA: Illinois, Effingham County')
                latitude_define_list.append('39.0011100')
                longitude_define_list.append('88.4943120')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Effingham County, (Latitude, Longitude = 39.122349, 88.434839)':
                country_define_list.append('USA: Illinois, Effingham County')
                latitude_define_list.append('39.1223490')
                longitude_define_list.append('88.4348390')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Fayette County, (Latitude, Longitude = 39.056735, 88.956398)':
                country_define_list.append('USA: Illinois, Fayette County')
                latitude_define_list.append('39.0567350')
                longitude_define_list.append('88.9563980')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Hancock County, (Latitude, Longitude = 40.355539, 91.020382)':
                country_define_list.append('USA: Illinois, Hancock County')
                latitude_define_list.append('40.3555390')
                longitude_define_list.append('91.0203820')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.152444, 90.080107)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.1524440')
                longitude_define_list.append('90.0801070')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.375614, 90.095164)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.3756140')
                longitude_define_list.append('90.0951640')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.387142, 90.32882)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.3871420')
                longitude_define_list.append('90.3288200')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.415749, 90.392462)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.4157490')
                longitude_define_list.append('90.3924620')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Henry County, (Latitude, Longitude = 41.49018, 90.01181)':
                country_define_list.append('USA: Illinois, Henry County')
                latitude_define_list.append('41.4901800')
                longitude_define_list.append('90.0118100')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.017453, 90.174874)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0174530')
                longitude_define_list.append('90.1748740')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.034799, 90.174874)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0347990')
                longitude_define_list.append('90.1748740')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.063306, 89.988752)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0633060')
                longitude_define_list.append('89.9887520')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Knox County, (Latitude, Longitude = 41.068925, 90.11095)':
                country_define_list.append('USA: Illinois, Knox County')
                latitude_define_list.append('41.0689250')
                longitude_define_list.append('90.1109500')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Livingston County, (Latitude, Longitude = 40.87108, 88.900607)':
                country_define_list.append('USA: Illinois, Livingston County')
                latitude_define_list.append('40.87108000')
                longitude_define_list.append('88.90060700')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Livingston County, (Latitude, Longitude = 40.933942, 88.490876)':
                country_define_list.append('USA: Illinois, Livingston County')
                latitude_define_list.append('40.9339420')
                longitude_define_list.append('88.4908760')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Logan County, (Latitude, Longitude = 40.063501, 89.285535)':
                country_define_list.append('USA: Illinois, Logan County')
                latitude_define_list.append('40.0635010')
                longitude_define_list.append('89.2855350')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Logan County, (Latitude, Longitude = 40.138741, 89.557122)':
                country_define_list.append('USA: Illin0ois, Logan County')
                latitude_define_list.append('40.138741')
                longitude_define_list.append('89.5571220')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mason County, (Latitude, Longitude = 40.301321, 89.970103)':
                country_define_list.append('USA: Illinois, Mason County')
                latitude_define_list.append('40.301321')
                longitude_define_list.append('89.9701030')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, McHenry County, (Latitude, Longitude = 42.291314, 88.538917)':
                country_define_list.append('USA: Illinois, McHenry County')
                latitude_define_list.append('42.2913140')
                longitude_define_list.append('88.5389170')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, McLean County, (Latitude, Longitude = 40.676635, 88.595756)':
                country_define_list.append('USA: Illinois, McLean County')
                latitude_define_list.append('40.676635')
                longitude_define_list.append('88.5957560')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, McLean County, (Latitude, Longitude = 40.716803, 88.882006)':
                country_define_list.append('USA: Illinois, McLean County')
                latitude_define_list.append('40.7168030')
                longitude_define_list.append('88.8820060')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Menard County, (Latitude, Longitude = 40.08722, 89.710619)':
                country_define_list.append('USA: Illinois, Menard County')
                latitude_define_list.append('40.0872200')
                longitude_define_list.append('89.710619')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Menard County, (Latitude, Longitude = 40.093363, 89.760698)':
                country_define_list.append('USA: Illinois, Menard County')
                latitude_define_list.append('40.0933630')
                longitude_define_list.append('89.7606980')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mercer County, (Latitude, Longitude = 41.06919, 90.743887)':
                country_define_list.append('USA: Illinois, Mercer County')
                latitude_define_list.append('41.0691900')
                longitude_define_list.append('90.7438870')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mercer County, (Latitude, Longitude = 41.085009, 90.496348)':
                country_define_list.append('USA: Illinois, Mercer County')
                latitude_define_list.append('41.0850090')
                longitude_define_list.append('90.4963480')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Mercer County, (Latitude, Longitude = 41.106902, 90.730675)':
                country_define_list.append('USA: Illinois, Mercer County')
                latitude_define_list.append('41.1069020')
                longitude_define_list.append('90.7306750')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Ogle County, (Latitude, Longitude = 41.983211, 89.245471)':
                country_define_list.append('USA: Illinois, Ogle County')
                latitude_define_list.append('41.9832110')
                longitude_define_list.append('89.2454710')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Schuyler County, (Latitude, Longitude = 40.244924, 90.668438)':
                country_define_list.append('USA: Illinois, Schuyler County')
                latitude_define_list.append('40.2449240')
                longitude_define_list.append('90.6684380')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Scott County, (Latitude, Longitude = 39.557258, 90.5532)':
                country_define_list.append('USA: Illinois, Scott County')
                latitude_define_list.append('39.5572580')
                longitude_define_list.append('90.5532000')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Shelby County, (Latitude, Longitude = 39.274053, 88.655896)':
                country_define_list.append('USA: Illinois, Shelby County')
                latitude_define_list.append('39.2740530')
                longitude_define_list.append('88.6558960')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Shelby County, (Latitude, Longitude = 41.034799, 90.174874)':
                country_define_list.append('USA: Illinois, Shelby County')
                latitude_define_list.append('41.0347990')
                longitude_define_list.append('90.1748740')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, St. Clair County, (Latiude, Longitude = 38.456916, 89.928058)':
                country_define_list.append('USA: Illinois, St. Clair County')
                latitude_define_list.append('38.4569160')
                longitude_define_list.append('89.9280580')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Warren County, (Latitude, Longitude = 40.747517, 90.570944)':
                country_define_list.append('USA: Illinois, Warren County')
                latitude_define_list.append('40.7475170')
                longitude_define_list.append('90.5709440')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Warren County, (Latitude, Longitude = 40.79257, 90.505294)':
                country_define_list.append('USA: Illinois, Warren County')
                latitude_define_list.append('40.7925700')
                longitude_define_list.append('90.5052940')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Warren County, (Latitude, Longitude = 40.882615, 90.465536)':
                country_define_list.append('USA: Illinois, Warren County')
                latitude_define_list.append('40.8826150')
                longitude_define_list.append('90.4655360')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Whiteside County, (Latitude, Longitude = 41.664499, 89.912534)':
                country_define_list.append('USA: Illinois, Whiteside County')
                latitude_define_list.append('41.6644990')
                longitude_define_list.append('89.9125340')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Whiteside County, (Latitude, Longitude = 41.715868, 89.966289)':
                country_define_list.append('USA: Illinois, Whiteside County')
                latitude_define_list.append('41.7158680')
                longitude_define_list.append('89.966289')
                continent_define_list.append('North American')
            if row == 'USA:Illinois, Woodford County, (Latitude, Longitude = 40.765446, 89.352716)':
                country_define_list.append('USA: Illinois, Woodford County')
                latitude_define_list.append('40.7654460')
                longitude_define_list.append('89.3527160')
                continent_define_list.append('North American')
            if row == 'Viet Nam':
                country_define_list.append('Viet Nam')
                latitude_define_list.append('10.7500000')
                longitude_define_list.append('106.6666641')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: An Giang':
                country_define_list.append('Viet Nam: An Giang')
                latitude_define_list.append('10.5215840')
                longitude_define_list.append('105.1258960')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Bac Lieu province':
                country_define_list.append('Viet Nam: Bac Lieu province')
                latitude_define_list.append('9.2850000')
                longitude_define_list.append('105.7244415')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Bac Ninh':
                country_define_list.append('Viet Nam: Bac Ninh')
                latitude_define_list.append('21.1833333')
                longitude_define_list.append('106.0500031')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Binh Duong':
                country_define_list.append('Viet Nam: Binh Duong')
                latitude_define_list.append('21.3833333')
                longitude_define_list.append('103.0166702')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Dien Bien':
                country_define_list.append('Viet Nam: Dien Bien')
                latitude_define_list.append('21.3833333')
                longitude_define_list.append('103.0166702')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Dong Nai':
                country_define_list.append('Viet Nam: Dong Nai')
                latitude_define_list.append('11.0686300')
                longitude_define_list.append('107.1675980')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hai Duong':
                country_define_list.append('Viet Nam: Hai Duong')
                latitude_define_list.append('20.9333333')
                longitude_define_list.append('106.3166656')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hai Phong':
                country_define_list.append('Viet Nam: Hai Phong')
                latitude_define_list.append('20.8449110')
                longitude_define_list.append('106.6880840')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hanoi':
                country_define_list.append('Viet Nam: Hanoi')
                latitude_define_list.append('21.0277640')
                longitude_define_list.append('105.8341600')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Ho Chi Minh':
                country_define_list.append('Viet Nam: Ho Chi Minh')
                latitude_define_list.append('10.75')
                longitude_define_list.append('106.6666641')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hung Yen':
                country_define_list.append('Viet Nam: Hung Yen')
                latitude_define_list.append('20.8525710')
                longitude_define_list.append('106.0169970')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Hung Yen province':
                country_define_list.append('Viet Nam: Hung Yen province')
                latitude_define_list.append('20.8525710')
                longitude_define_list.append('106.0169970')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Lao Cai':
                country_define_list.append('Viet Nam: Lao Cai')
                latitude_define_list.append('22.48333330')
                longitude_define_list.append('103.9499969')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Nghe An':
                country_define_list.append('Viet Nam: Nghe An')
                latitude_define_list.append('19.2342490')
                longitude_define_list.append('104.9200360')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Nghe An province':
                country_define_list.append('Viet Nam: Nghe An province')
                latitude_define_list.append('19.2342490')
                longitude_define_list.append('104.9200360')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Quang Ninh province':
                country_define_list.append('Viet Nam: Quang Ninh province')
                latitude_define_list.append('21.0063820')
                longitude_define_list.append('107.2925140')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Son La':
                country_define_list.append('Viet Nam: Son La')
                latitude_define_list.append('21.3166667')
                longitude_define_list.append('103.728917')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Tay Ninh':
                country_define_list.append('Viet Nam: Tay Ninh')
                latitude_define_list.append('11.335155')
                longitude_define_list.append('103.9000015')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Thai Binh':
                country_define_list.append('Viet Nam: Thai Binh')
                latitude_define_list.append('20.5386940')
                longitude_define_list.append('20.5386940')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Tien Giang province':
                country_define_list.append('Viet Nam: Tien Giang province')
                latitude_define_list.append('10.4493320')
                longitude_define_list.append('106.3420500')
                continent_define_list.append('Asia')
            if row == 'Viet Nam: Yen Bai':
                country_define_list.append('Viet Nam: Yen Bai')
                latitude_define_list.append('21.683792')
                longitude_define_list.append('21.7000000')
                continent_define_list.append('104.8666687')
            if row == 'Viet Nam:Ho Chi Minh':
                country_define_list.append('Viet Nam:Ho Chi Minh')
                latitude_define_list.append('10.7500000')
                longitude_define_list.append('106.6666641')
                continent_define_list.append('Asia')
            if row == 'Viet Nam:Nam Dinh':
                country_define_list.append('Viet Nam:Nam Dinh')
                latitude_define_list.append('20.4388230')
                longitude_define_list.append('106.1621050')
                continent_define_list.append('Asia')
            if row == 'Viet Nam:Thai Binh':
                country_define_list.append('Viet Nam: Thai Binh')
                latitude_define_list.append('20.4500000')
                longitude_define_list.append('106.393478')
                continent_define_list.append('Asia')
        #extract_zip = zip(country_list,latitude_list,longitude_list,continent_list)
        #print(extract_zip)
        #extract_list = []
        #for extract in extract_zip:
        #    extract_list.append("\t".join(extract))
        #country_string = '\n'.join(str(country_value) for country_value in extract_list)
        #print(country_string)
        #handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_country_define.txt", "w")
        #handle_outputfile.write(str(country_string))
        #handle_outputfile.close()
        # 3 ##################################################
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk", "genbank"))
        # pummedid
        accession_list = []
        organism_list = []
        gene_list = []
        length_list = []
        pubmedid_list = []
        publishdate_list = []
        note_list = []
        genotype_list = []
        country_list = []
        collectiondate_list = []
        origin_list = []
        for record in handle_parse_db:
            if record.annotations["organism"] == 'Porcine reproductive and respiratory syndrome virus':
                if len(record) == 600 or len(record) == 603 or len(record) == 606:
                    if record.features:
                        for feature in record.features:
                            if feature.type == "source":
                                accession_list.append(record.id)
                                organism_list.append(record.annotations["organism"])
                                gene_list.append('ORF5')
                                length_list.append(str(len(record)))
                                pubmedid_list.append(record.annotations["references"][0].pubmed_id)
                                publishdate_list.append(record.annotations["date"])
                                note_list.append(feature.qualifiers.get("note"))
                                genotype_list.append('North American')
                                country_list.append(feature.qualifiers.get("country"))
                                collectiondate_list.append(feature.qualifiers.get("collection_date"))
                                origin_list.append(str(record.seq))
        # pubmedid
        organism_list_with_none = [['None'] if organism_value is None else organism_value for organism_value in organism_list]
        organism_none_list = []
        for organism_none in organism_list_with_none:
            organism_none_list.append(organism_none)
        publishdate_list_with_none = [['None'] if publishdate_value is None else publishdate_value for publishdate_value in publishdate_list]
        publishdate_none_list = []
        for publishdate_none in publishdate_list_with_none:
            publishdate_none_list.append(publishdate_none)
        note_list_with_none = [['None'] if note_value is None else note_value for note_value in note_list]
        note_none_list = []
        for note_none in note_list_with_none:
            note_none_list.append('\n'.join(note_none))
        country_list_with_none = [['None'] if country_value is None else country_value for country_value in country_list]
        country_none_list = []
        for country_none in country_list_with_none:
            country_none_list.append('\n'.join(country_none))
        collectiondate_list_with_none = [['None'] if collectiondate_value is None else collectiondate_value for collectiondate_value in collectiondate_list]
        collectiondate_none_list = []
        for collectiondate_none in collectiondate_list_with_none:
            collectiondate_none_list.append('\n'.join(collectiondate_none))
        origin_list_with_none = [['None'] if origin_value is None else origin_value for origin_value in origin_list]
        origin_none_list = []
        for origin_none in origin_list_with_none:
            origin_none_list.append(origin_none)
        pubmedid_zip = zip(accession_list,organism_none_list,gene_list,length_list,collectiondate_none_list,pubmedid_list,publishdate_none_list,note_none_list,genotype_list,country_none_list,country_define_list,latitude_define_list,longitude_define_list,continent_define_list,origin_none_list)
        print(pubmedid_zip)
        orf5_pubmedid_list = []
        for orf5_pubmedid in pubmedid_zip:
            orf5_pubmedid_list.append("\t".join(orf5_pubmedid))
        # print(ORF5_NA_list)
        print("")
        # Writing ORF5 NA
        print("Step 2: Writing ORF5 NA, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        orf5_pubmedid_string = '\n'.join(str(orf5_value) for orf5_value in orf5_pubmedid_list)
        print(orf5_pubmedid_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_orf5_na_extract.txt"), "w")
        handle_outputfile.write(orf5_pubmedid_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 46: ORF5 NA Extract Complete")
        print("")
        break

# All Extract ##################################################
def allextract():
    while True:
        # Module Name
        print("Module 47: All Extract")
        print("")
        # Parsing All
        print("Step 1: Parsing All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        # Accession Number
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_number_list_with_none = []
        for record in handle_parse_db:
            accession_number_list_with_none.append(record.id)
        accession_number_list_with_none_list = [['None'] if accession_number_value is None else accession_number_value for accession_number_value in accession_number_list_with_none]
        accession_number_list = []
        for accession_number_none in accession_number_list_with_none_list:
            accession_number_list.append(accession_number_none)
        # Annotations Accessions
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_accessions_list_with_none = []
        for record in handle_parse_db:
            annotations_accessions_list_with_none.append(record.annotations["accessions"])
        annotations_accessions_list_with_none_list = [['None'] if annotations_accessions_value is None else annotations_accessions_value for annotations_accessions_value in annotations_accessions_list_with_none]
        annotations_accessions_list = []
        for annotations_accessions_none in annotations_accessions_list_with_none_list:
            annotations_accessions_list.append(annotations_accessions_none)
        # Annotations Comment
        #handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk", "genbank")
        #annotations_comment_list_with_none = []
        #for record in handle_parse_db:
        #    annotations_comment_list_with_none.append(record.annotations["comment"])
        #annotations_comment_list_with_none_list = [['None'] if annotations_comment_value is None else annotations_comment_value for annotations_comment_value in annotations_comment_list_with_none]
        #annotations_comment_list = []
        #for annotations_comment_none in annotations_comment_list_with_none_list:
        #    annotations_comment_list.append(annotations_comment_none)
        # Annotations Date
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_date_list_with_none = []
        for record in handle_parse_db:
            annotations_date_list_with_none.append(record.annotations["date"])
        annotations_date_list_with_none_list = [['None'] if annotations_date_value is None else annotations_date_value for annotations_date_value in annotations_date_list_with_none]
        annotations_date_list = []
        for annotations_date_none in annotations_date_list_with_none_list:
            annotations_date_list.append(annotations_date_none)
        # Annotations Data FIle Division
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_data_file_division_list_with_none = []
        for record in handle_parse_db:
            annotations_data_file_division_list_with_none.append(record.annotations["data_file_division"])
        annotations_data_file_division_list_with_none_list = [['None'] if annotations_data_file_division_value is None else annotations_data_file_division_value for annotations_data_file_division_value in annotations_data_file_division_list_with_none]
        annotations_data_file_division_list = []
        for annotations_data_file_division_none in annotations_data_file_division_list_with_none_list:
            annotations_data_file_division_list.append(annotations_data_file_division_none)
        # Annotations Evidence
        #handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk", "genbank")
        #annotations_evidence_list_with_none = []
        #for record in handle_parse_db:
        #    annotations_evidence_list_with_none.append(record.annotations["evidence"])
        #annotations_evidence_list_with_none_list = [['None'] if annotations_evidence_value is None else annotations_evidence_value for annotations_evidence_value in annotations_evidence_list_with_none]
        #annotations_evidence_list = []
        #for annotations_evidence_none in annotations_evidence_list_with_none_list:
        #    annotations_evidence_list.append(annotations_evidence_none)
        # Annotations GI
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_gi_list_with_none = []
        for record in handle_parse_db:
            annotations_gi_list_with_none.append(record.annotations["gi"])
        annotations_gi_list_with_none_list = [['None'] if annotations_gi_value is None else annotations_gi_value for annotations_gi_value in annotations_gi_list_with_none]
        annotations_gi_list = []
        for annotations_gi_none in annotations_gi_list_with_none_list:
            annotations_gi_list.append(annotations_gi_none)
        # Annotations Keywords
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_keywords_list_with_none = []
        for record in handle_parse_db:
            annotations_keywords_list_with_none.append(record.annotations["keywords"])
        annotations_keywords_list_with_none_list = [['None'] if annotations_keywords_value is None else annotations_keywords_value for annotations_keywords_value in annotations_keywords_list_with_none]
        annotations_keywords_list = []
        for annotations_keywords_none in annotations_keywords_list_with_none_list:
            annotations_keywords_list.append(annotations_keywords_none)
        # Annotations Molecule Type
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_molecule_type_list_with_none = []
        for record in handle_parse_db:
            annotations_molecule_type_list_with_none.append(record.annotations["molecule_type"])
        annotations_molecule_type_list_with_none_list = [['None'] if annotations_molecule_type_value is None else annotations_molecule_type_value for annotations_molecule_type_value in annotations_molecule_type_list_with_none]
        annotations_molecule_type_list = []
        for annotations_molecule_type_none in annotations_molecule_type_list_with_none_list:
            annotations_molecule_type_list.append(annotations_molecule_type_none)
        # Annotations Organism
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_organism_list_with_none = []
        for record in handle_parse_db:
            annotations_organism_list_with_none.append(record.annotations["organism"])
        annotations_organism_list_with_none_list = [['None'] if annotations_organism_value is None else annotations_organism_value for annotations_organism_value in annotations_organism_list_with_none]
        annotations_organism_list = []
        for annotations_organism_none in annotations_organism_list_with_none_list:
            annotations_organism_list.append(annotations_organism_none)
        # Annotations References
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_references_list_with_none = []
        for record in handle_parse_db:
            annotations_references_list_with_none.append(record.annotations["references"])
        annotations_references_list_with_none_list = [['None'] if annotations_references_value is None else annotations_references_value for annotations_references_value in annotations_references_list_with_none]
        annotations_references_list = []
        for annotations_references_none in annotations_references_list_with_none_list:
            annotations_references_list.append(annotations_references_none)
        # Annotations Sequence Version
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_sequence_version_list_with_none = []
        for record in handle_parse_db:
            annotations_sequence_version_list_with_none.append(record.annotations["sequence_version"])
        annotations_sequence_version_list_with_none_list = [['None'] if annotations_sequence_version_value is None else annotations_sequence_version_value for annotations_sequence_version_value in annotations_sequence_version_list_with_none]
        annotations_sequence_version_list = []
        for annotations_sequence_version_none in annotations_sequence_version_list_with_none_list:
            annotations_sequence_version_list.append(annotations_sequence_version_none)
        # Annotation Source
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_source_list_with_none = []
        for record in handle_parse_db:
            annotations_source_list_with_none.append(record.annotations["source"])
        annotations_source_list_with_none_list = [['None'] if annotations_source_value is None else annotations_source_value for annotations_source_value in annotations_source_list_with_none]
        annotations_source_list = []
        for annotations_source_none in annotations_source_list_with_none_list:
            annotations_source_list.append(annotations_source_none)
        # Annotations Taxonomy
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_taxonomy_list_with_none = []
        for record in handle_parse_db:
            annotations_taxonomy_list_with_none.append(record.annotations["taxonomy"])
        annotations_taxonomy_list_with_none_list = [['None'] if annotations_taxonomy_value is None else annotations_taxonomy_value for annotations_taxonomy_value in annotations_taxonomy_list_with_none]
        annotations_taxonomy_list = []
        for annotations_taxonomy_none in annotations_taxonomy_list_with_none_list:
            annotations_taxonomy_list.append(annotations_taxonomy_none)
        # Annotations Topology
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_topology_list_with_none = []
        for record in handle_parse_db:
            annotations_topology_list_with_none.append(record.annotations["topology"])
        annotations_topology_list_with_none_list = [['None'] if annotations_topology_value is None else annotations_topology_value for annotations_topology_value in annotations_topology_list_with_none]
        annotations_topology_list = []
        for annotations_topology_none in annotations_topology_list_with_none_list:
            annotations_topology_list.append(annotations_topology_none)
        # Annotations Keys
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_keys_list_with_none = []
        for record in handle_parse_db:
            annotations_keys_list_with_none.append(record.annotations.keys)
        annotations_keys_list_with_none_list = [['None'] if annotations_keys_value is None else annotations_keys_value for annotations_keys_value in annotations_keys_list_with_none]
        annotations_keys_list = []
        for annotations_keys_none in annotations_keys_list_with_none_list:
            annotations_keys_list.append(annotations_keys_none)
        # Annotations Values
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        annotations_values_list_with_none = []
        for record in handle_parse_db:
            annotations_values_list_with_none.append(record.annotations.values)
        annotations_values_list_with_none_list = [['None'] if annotations_values_value is None else annotations_values_value for annotations_values_value in annotations_values_list_with_none]
        annotations_values_list = []
        for annotations_values_none in annotations_values_list_with_none_list:
            annotations_values_list.append(annotations_values_none)
        # DB XREFS
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        db_xrefs_list_with_none = []
        for record in handle_parse_db:
            db_xrefs_list_with_none.append(record.dbxrefs)
        db_xrefs_list_with_none_list = [['None'] if db_xrefs_value is None else db_xrefs_value for db_xrefs_value in db_xrefs_list_with_none]
        db_xrefs_list = []
        for db_xrefs_none in db_xrefs_list_with_none_list:
            db_xrefs_list.append(db_xrefs_none)
        # Definition
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        definition_list_with_none = []
        for record in handle_parse_db:
            definition_list_with_none.append(record.description)
        definition_list_with_none_list = [['None'] if definition_value is None else definition_value for definition_value in definition_list_with_none]
        definition_list = []
        for definition_none in definition_list_with_none_list:
            definition_list.append(definition_none)
        # Features
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        features_list_with_none = []
        for record in handle_parse_db:
            features_list_with_none.append(record.features)
        features_list_with_none_list = [['None'] if features_value is None else features_value for features_value in features_list_with_none]
        features_list = []
        for features_none in features_list_with_none_list:
            features_list.append(features_none)
        # Features Length
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        features_length_list_with_none = []
        for record in handle_parse_db:
            features_length_list_with_none.append(len(record.features))
        features_length_list_with_none_list = [['None'] if features_length_value is None else features_length_value for features_length_value in features_length_list_with_none]
        features_length_list = []
        for features_length_none in features_length_list_with_none_list:
            features_length_list.append(features_length_none)
        # Format
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        format_list_with_none = []
        for record in handle_parse_db:
            format_list_with_none.append(record.format)
        format_list_with_none_list = [['None'] if format_value is None else format_value for format_value in format_list_with_none]
        format_list = []
        for format_none in format_list_with_none_list:
            format_list.append(format_none)
        # Length
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        length_list_with_none = []
        for record in handle_parse_db:
            length_list_with_none.append(len(record))
        length_list_with_none_list = [['None'] if length_value is None else length_value for length_value in length_list_with_none]
        length_list = []
        for length_none in length_list_with_none_list:
            length_list.append(length_none)
        # Lenght Annotations
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        length_annotations_list_with_none = []
        for record in handle_parse_db:
            length_annotations_list_with_none.append(len(record.annotations))
        length_annotations_list_with_none_list = [['None'] if length_annotations_value is None else length_annotations_value for length_annotations_value in length_annotations_list_with_none]
        length_annotations_list = []
        for length_annotations_none in length_annotations_list_with_none_list:
            length_annotations_list.append(length_annotations_none)
        # Letter Annotations
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        letter_annotations_list_with_none = []
        for record in handle_parse_db:
            letter_annotations_list_with_none.append(record.letter_annotations)
        letter_annotations_list_with_none_list = [['None'] if letter_annotations_value is None else letter_annotations_value for letter_annotations_value in letter_annotations_list_with_none]
        letter_annotations_list = []
        for letter_annotations_none in letter_annotations_list_with_none_list:
            letter_annotations_list.append(letter_annotations_none)
        # Letter Annotations Phred Quality
        #handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk", "genbank")
        #letter_annotations_phred_quality_list_with_none = []
        #for record in handle_parse_db:
        #    letter_annotations_quality_list_with_none.append(record.letter_annotations["phred_quality"])
        #letter_annotations_quality_list_with_none_list = [['None'] if letter_annotations_quality_value is None else letter_annotations_quality_value for letter_annotations_quality_value in letter_annotations_quality_list_with_none]
        #letter_annotations_quality_list = []
        #for letter_annotations_quality_none in letter_annotations_quality_list_with_none_list:
        #    letter_annotations_quality_list.append(letter_annotations_quality_none)
        # Locus
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        locus_list_with_none = []
        for record in handle_parse_db:
            locus_list_with_none.append(record.name)
        locus_list_with_none_list = [['None'] if locus_value is None else locus_value for locus_value in locus_list_with_none]
        locus_list = []
        for locus_none in locus_list_with_none_list:
            locus_list.append(locus_none)
        # Origin
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        origin_list_with_none = []
        for record in handle_parse_db:
            origin_list_with_none.append(repr(record.seq))
        origin_list_with_none_list = [['None'] if origin_value is None else origin_value for origin_value in origin_list_with_none]
        origin_list = []
        for origin_none in origin_list_with_none_list:
            origin_list.append(origin_none)
        # Zipping All Extract
        #all_extract_zip = zip(accession_number_list,annotations_accessions_list,annotations_comment_list,annotations_date_list,annotations_data_file_division_list,annotations_evidence_list,annotations_gi_list,annotations_keywords_list,annotations_molecule_type_list,annotations_organism_list,annotations_references_list,annotations_sequence_version_list,annotations_source_list,annotations_taxonomy_list,annotations_topology_list,annotations_keys_list,annotations_values_list,db_xrefs_list,definition_list,features_list,features_length_list,format_list,length_list,length_annotations_list,letter_annotations_list,letter_annotations_phred_quality_list,locus_list,origin_list)
        all_extract_zip = zip(accession_number_list,annotations_accessions_list,annotations_date_list,annotations_data_file_division_list,annotations_gi_list,annotations_keywords_list,annotations_molecule_type_list,annotations_organism_list,annotations_references_list,annotations_sequence_version_list,annotations_source_list,annotations_taxonomy_list,annotations_topology_list,annotations_keys_list,annotations_values_list,db_xrefs_list,definition_list,features_list,features_length_list,format_list,length_list,length_annotations_list,letter_annotations_list,locus_list,origin_list)
        print(all_extract_zip)
        all_extract_list = []
        for all_feature in all_extract_zip:
            all_extract_list.append(all_feature)
        print("")
        # Writing All
        print("Step 2: Writing All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        all_string = '\t'.join(str(all_value) for all_value in all_extract_list)
        print(all_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/description_all_extract.txt"), "w")
        handle_outputfile.write(str(all_string))
        handle_outputfile.close()
        print("")
        # All Complete
        print("Module 47: All Extract Complete")
        print("")
        break

# Module Main ##################################################
def maindescription():
    while True:
        print("")
        print("\tMenu: Extract Description\t\t\t# Update later.")
        print("")
        print("\t\tModule 1.  Accession Number \t\t: Button[1]")
        print("\t\tModule 2.  Annotations (Accessions) \t: Button[2]")
        print("\t\tModule 3.  Annotations (Author) \t: Button[3]")
        print("\t\tModule 4.  Annotations (Comment) \t: Button[4]")
        print("\t\tModule 5.  Annotations (Consrtm) \t: Button[5]")
        print("\t\tModule 6.  Annotations (Date) \t\t: Button[6]")
        print("\t\tModule 7.  Annotations (Data File Div.) : Button[7]")
        print("\t\tModule 8.  Annotations (Evidence) \t: Button[8] #")       
        print("\t\tModule 9.  Annotations (GI) \t\t: Button[9]")
        print("\t\tModule 10. Annotations (Journal[0]) \t: Button[10]")
        print("\t\tModule 11. Annotations (Keywords) \t: Button[11]")
        print("\t\tModule 12. Annotations (Location) \t: Button[12]")
        print("\t\tModule 13. Annotations (Medline ID) \t: Button[13]")
        print("\t\tModule 14. Annotations (Molecule Type) \t: Button[14]")    
        print("\t\tModule 15. Annotations (Organism) \t: Button[15]")
        print("\t\tModule 16. Annotations (Reference Length) : Button[16]")
        print("\t\tModule 17. Annotations (Ref[0].Authors) : Button[17]")
        print("\t\tModule 18. Annotations (Ref[0].Title) \t: Button[18]")
        print("\t\tModule 19. Annotations (Ref[0].Journal) : Button[19]")
        print("\t\tModule 20. Annotations (Ref[0].PubmedID): Button[20]")
        print("\t\tModule 21. Annotations (Ref[All].Authors) : Button[21]")
        print("\t\tModule 22. Annotations (Ref[All].Title)   : Button[22]")
        print("\t\tModule 23. Annotations (Ref[All].Journal) : Button[23]")
        print("\t\tModule 24. Annotations (Ref[All].PubmedID): Button[24]")
        print("\t\tModule 25. Annotations (Sequence Ver.) \t: Button[25]")     
        print("\t\tModule 26. Annotations (Source) \t: Button[26]")
        print("\t\tModule 27. Annotations (Taxonomy) \t: Button[27]")
        print("\t\tModule 28. Annotations (Title[0]) \t: Button[28]")
        print("\t\tModule 29. Annotations (Topology) \t: Button[29]")     
        print("\t\tModule 30. Annotations Keys \t\t: Button[30]")
        print("\t\tModule 31. Annotations Values \t\t: Button[31]")       
        print("\t\tModule 32. DB XREFS \t\t\t: Button[32]")
        print("\t\tModule 33. Definition \t\t\t: Button[33]")
        print("\t\tModule 34. Features \t\t\t: Button[34]")
        print("\t\tModule 35. Features Length \t\t: Button[35]")
        print("\t\tModule 36. Format \t\t\t: Button[36]")
        print("\t\tModule 37. Length \t\t\t: Button[37]")
        print("\t\tModule 38. Length Annotations\t\t: Button[38]")
        print("\t\tModule 39. Letter Annotations\t\t: Button[39]")
        print("\t\tModule 40. Letter Anno. (Phred Quality)\t: Button[40] #")       
        print("\t\tModule 41. Locus \t\t\t: Button[41]")
        print("\t\tModule 42. Origin \t\t\t: Button[42]")
        print("\t\tModule 43. Origin[All] \t\t\t: Button[43]")
        print("\t\tModule 44. ORF5 Filter \t\t\t: Button[44]")
        print("\t\tModule 45. ORF5 EU extract \t\t: Button[45]")
        print("\t\tModule 46. ORF5 NA extract \t\t: Button[46]")
        print("\t\tModule 47. All Extract \t\t\t: Button[47] #")
        print("\t\tModule 48. Main Menu \t\t\t: Button[48]")
        print("\t\tModule 49. Exit. \t\t\t: Button[49]")
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
                accessionnumber()
            elif button == 2:
                print("")
                annotationsaccessions()               
            elif button == 3:
                print("")
                annotationsauthors()
            elif button == 4:
                print("")
                annotationscomment()
            elif button == 5:
                print("")
                annotationsconsrtm()               
            elif button == 6:
                print("")
                annotationsdate()
            elif button == 7:
                print("")
                annotationsdatafiledivision()
            #elif button == 8:
            #    print("")
            #    annotationsevidence()               
            elif button == 9:
                print("")
                annotationsgi()
            elif button == 10:
                print("")
                annotationsjournal()
            elif button == 11:
                print("")
                annotationskeywords()
            elif button == 12:
                print("")
                annotationslocation()
            elif button == 13:
                print("")
                annotationsmedlineid()
            elif button == 14:
                print("")
                annotationsmoleculetype()
            elif button == 15:
                print("")
                annotationsorganism()
            elif button == 16:
                print("")
                annotationsreferenceslength()
            elif button == 17:
                print("")
                annotationsref0authors()
            elif button == 18:
                print("")
                annotationsref0title()
            elif button == 19:
                print("")
                annotationsref0journal()
            elif button == 20:
                print("")
                annotationsref0pubmedid()
            elif button == 21:
                print("")
                annotationsrefallauthors()
            elif button == 22:
                print("")
                annotationsrefalltitle()
            elif button == 23:
                print("")
                annotationsrefalljournal()
            elif button == 24:
                print("")
                annotationsrefallpubmedid()
            elif button == 25:
                print("")
                annotationssequenceversion()
            elif button == 26:
                print("")
                annotationssource()
            elif button == 27:
                print("")
                annotationstaxonomy()
            elif button == 28:
                print("")
                annotationstitle()               
            elif button == 29:
                print("")
                annotationstopology()               
            elif button == 30:
                print("")
                annotationskeys()
            elif button == 31:
                print("")
                annotationsvalues()
            elif button == 32:
                print("")
                dbxrefs()
            elif button == 33:
                print("")
                definition()
            elif button == 34:
                print("")
                features()
            elif button == 35:
                print("")
                featureslength()
            elif button == 36:
                print("")
                formatdescription()
            elif button == 37:
                print("")
                length()
            elif button == 38:
                print("")
                lengthannotations()
            elif button == 39:
                print("")
                letterannotations()
            #elif button == 40:
            #    print("")
            #    letterannotationsphredquality()               
            elif button == 41:
                print("")
                locus()
            elif button == 42:
                print("")
                origin()
            elif button == 43:
                print("")
                originall()
            elif button == 44:
                print("")
                orf5filter()
            elif button == 45:
                print("")
                orf5euextract()
            elif button == 46:
                print("")
                orf5naextract()        
            #elif button == 47:
            #    print("")
            #    allextract()
            elif button == 48:
                print("")
                import main
                main.main()
            elif button == 49:
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

maindescription()
