#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from geopy.geocoders import Nominatim
from Bio.Alphabet import generic_dna, generic_protein

# Module: Source Features Extract
button = ""

# Acronym Extract ##################################################
def acronym():
    while True:
        # Module Name
        print("Module 1: Acronym Extract")
        print("")
        # Parsing Acronym
        print("Step 1: Parsing Acronym, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        acronym_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        acronym = feature.qualifiers.get('acronym')
                        acronym_list_with_none.append(acronym)
        acronym_list_with_none_list = [['None'] if acronym_value is None else acronym_value for acronym_value in acronym_list_with_none]
        acronym_list = []
        for acronym_none in acronym_list_with_none_list:
            acronym_list.append('\n'.join(acronym_none))
        print(acronym_list)
        print("")
        # Writing Acronym
        print("Step 2: Writing Acronym, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        acronym_string = "\n".join(acronym_list)
        print(acronym_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_acronym.txt"), "w")
        handle_outputfile.write(acronym_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 1: Acronym Extract Complete")
        print("")
        break

# Altitude Extract ##################################################
def altitude():
    while True:
        # Module Name
        print("Module 2: Altitude Extract")
        print("")
        # Parsing Altitude
        print("Step 1: Parsing Altitude, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        altitude_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        altitude = feature.qualifiers.get('altitude')
                        altitude_list_with_none.append(altitude)
        altitude_list_with_none_list = [['None'] if altitude_value is None else altitude_value for altitude_value in altitude_list_with_none]
        altitude_list = []
        for altitude_none in altitude_list_with_none_list:
            altitude_list.append('\n'.join(altitude_none))
        print(altitude_list)
        print("")
        # Writing Altitude
        print("Step 2: Writing Altitude, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        altitude_string = "\n".join(altitude_list)
        print(altitude_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_altitude.txt"), "w")
        handle_outputfile.write(altitude_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 2: Altitude Extract Complete")
        print("")
        break

# Anamorph Extract ##################################################
def anamorph():
    while True:
        # Module Name
        print("Module 3: Anamorph Extract")
        print("")
        # Parsing Anamorph
        print("Step 1: Parsing Anamorph, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        anamorph_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        anamorph = feature.qualifiers.get('anamorph')
                        anamorph_list_with_none.append(anamorph)
        anamorph_list_with_none_list = [['None'] if anamorph_value is None else anamorph_value for anamorph_value in anamorph_list_with_none]
        anamorph_list = []
        for anamorph_none in anamorph_list_with_none_list:
            anamorph_list.append('\n'.join(anamorph_none))
        print(anamorph_list)
        print("")
        # Writing Anamorph
        print("Step 2: Writing Anamorph, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        anamorph_string = "\n".join(anamorph_list)
        print(anamorph_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_anamorph.txt"), "w")
        handle_outputfile.write(anamorph_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 3: Anamorph Extract Complete")
        print("")
        break

# Authority Extract ##################################################
def authority():
    while True:
        # Module Name
        print("Module 4: Authority Extract")
        print("")
        # Parsing Authority
        print("Step 1: Parsing Authority, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        authority_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        authority = feature.qualifiers.get('authority')
                        authority_list_with_none.append(authority)
        authority_list_with_none_list = [['None'] if authority_value is None else authority_value for authority_value in authority_list_with_none]
        authority_list = []
        for authority_none in authority_list_with_none_list:
            authority_list.append('\n'.join(authority_none))
        print(authority_list)
        print("")
        # Writing Authority
        print("Step 2: Writing Authority, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        authority_string = "\n".join(authority_list)
        print(authority_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_authority.txt"), "w")
        handle_outputfile.write(authority_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 4: Authority Extract Complete")
        print("")
        break

# Bio Material Extract ##################################################
def biomaterial():
    while True:
        # Module Name
        print("Module 5: Bio Material Extract")
        print("")
        # Parsing Bio Material
        print("Step 1: Parsing Bio Material, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        bio_material_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        bio_material = feature.qualifiers.get('bio_material')
                        bio_material_list_with_none.append(bio_material)
        bio_material_list_with_none_list = [['None'] if bio_material_value is None else bio_material_value for bio_material_value in bio_material_list_with_none]
        bio_material_list = []
        for bio_material_none in bio_material_list_with_none_list:
            bio_material_list.append('\n'.join(bio_material_none))
        print(bio_material_list)
        print("")
        # Writing Bio Material
        print("Step 2: Writing Bio Material, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        bio_material_string = "\n".join(bio_material_list)
        print(bio_material_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_bio_material.txt"), "w")
        handle_outputfile.write(bio_material_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 5: Bio Material Extract Complete")
        print("")
        break

# Biotype Extract ##################################################
def biotype():
    while True:
        # Module Name
        print("Module 6: Biotype Extract")
        print("")
        # Parsing Biotype
        print("Step 1: Parsing Biotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        biotype_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        biotype = feature.qualifiers.get('biotype')
                        biotype_list_with_none.append(biotype)
        biotype_list_with_none_list = [['None'] if biotype_value is None else biotype_value for biotype_value in biotype_list_with_none]
        biotype_list = []
        for biotype_none in biotype_list_with_none_list:
            biotype_list.append('\n'.join(biotype_none))
        print(biotype_list)
        print("")
        # Writing Biotype
        print("Step 2 : Writing Biotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        biotype_string = "\n".join(biotype_list)
        print(biotype_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_biotype.txt"), "w")
        handle_outputfile.write(biotype_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 6 : Biotype Extract Complete")
        print("")
        break

# Biovar Extract ##################################################
def biovar():
    while True:
        # Module Name
        print("Module 7: Biovar Extract")
        print("")
        # Parsing Biotype
        print("Step 1: Parsing Biovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        biovar_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        biovar = feature.qualifiers.get('biovar')
                        biovar_list_with_none.append(biovar)
        biovar_list_with_none_list = [['None'] if biovar_value is None else biovar_value for biovar_value in biovar_list_with_none]
        biovar_list = []
        for biovar_none in biovar_list_with_none_list:
            biovar_list.append('\n'.join(biovar_none))
        print(biovar_list)
        print("")
        # Writing Biovar
        print("Step 2: Writing Biovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        biovar_string = "\n".join(biovar_list)
        print(biovar_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_biovar.txt"), "w")
        handle_outputfile.write(biovar_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 7: Biovar Extract Complete")
        print("")
        break

# Breed Extract ##################################################
def breed():
    while True:
        # Module Name
        print("Module 8: Breed Extract")
        print("")
        # Parsing Breed
        print("Step 1: Parsing Breed, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        breed_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        breed = feature.qualifiers.get('breed')
                        breed_list_with_none.append(breed)
        breed_list_with_none_list = [['None'] if breed_value is None else breed_value for breed_value in breed_list_with_none]
        breed_list = []
        for breed_none in breed_list_with_none_list:
            breed_list.append('\n'.join(breed_none))
        print(breed_list)
        print("")
        # Writing Breed
        print("Step 2: Writing Breed, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        breed_string = "\n".join(breed_list)
        print(breed_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_breed.txt"), "w")
        handle_outputfile.write(breed_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 8: Breed Extract Complete")
        print("")
        break

# Cell Line Extract ##################################################
def cellline():
    while True:
        # Module Name
        print("Module 9: Cell Line Extract")
        print("")
        # Parsing Cell Line
        print("Step 1: Parsing Cell Line, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        cell_line_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        cell_line = feature.qualifiers.get('cell_line')
                        cell_line_list_with_none.append(cell_line)
        cell_line_list_with_none_list = [['None'] if cell_line_value is None else cell_line_value for cell_line_value in cell_line_list_with_none]
        cell_line_list = []
        for cell_line_none in cell_line_list_with_none_list:
            cell_line_list.append('\n'.join(cell_line_none))
        print(cell_line_list)
        print("")
        # Writing Cell Line
        print("Step 2: Writing Cell Line, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        cell_line_string = "\n".join(cell_line_list)
        print(cell_line_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_cell_line.txt"), "w")
        handle_outputfile.write(cell_line_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 9: Cell Line Extract Complete")
        print("")
        break

# Cell Type Extract ##################################################
def celltype():
    while True:
        # Module Name
        print("Module 10: Cell Type Extract")
        print("")
        # Parsing Cell Type
        print("Step 1: Parsing Cell Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        cell_type_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        cell_type = feature.qualifiers.get('cell_type')
                        cell_type_list_with_none.append(cell_type)
        cell_type_list_with_none_list = [['None'] if cell_type_value is None else cell_type_value for cell_type_value in cell_type_list_with_none]
        cell_type_list = []
        for cell_type_none in cell_type_list_with_none_list:
            cell_type_list.append('\n'.join(cell_type_none))
        print(cell_type_list)
        print("")
        # Writing Cell Type
        print("Step 2: Writing Cell Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        cell_type_string = "\n".join(cell_type_list)
        print(cell_type_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_cell_type.txt"), "w")
        handle_outputfile.write(cell_type_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 10: Cell Type Extract Complete")
        print("")
        break

# Chemovar Extract ##################################################
def chemovar():
    while True:
        # Module Name
        print("Module 11: Chemovar Extract")
        print("")
        # Parsing Chemovar
        print("Step 1: Parsing Chemovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        chemovar_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        chemovar = feature.qualifiers.get('chemovar')
                        chemovar_list_with_none.append(chemovar)
        chemovar_list_with_none_list = [['None'] if chemovar_value is None else chemovar_value for chemovar_value in chemovar_list_with_none]
        chemovar_list = []
        for chemovar_none in chemovar_list_with_none_list:
            chemovar_list.append('\n'.join(chemovar_none))
        print(chemovar_list)
        print("")
        # Writing Chemovar
        print("Step 2: Writing Chemovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        chemovar_string = "\n".join(chemovar_list)
        print(chemovar_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_chemovar.txt"), "w")
        handle_outputfile.write(chemovar_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 11: Chemovar Extract Complete")
        print("")
        break

# Chromosome Extract ##################################################
def chromosome():
    while True:
        # Module Name
        print("Module 12: Chromosome Extract")
        print("")
        # Parsing Chromosome
        print("Step 1: Parsing Chromosome, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        chromosome_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        chromosome = feature.qualifiers.get('chromosome')
                        chromosome_list_with_none.append(chromosome)
        chromosome_list_with_none_list = [['None'] if chromosome_value is None else chromosome_value for chromosome_value in chromosome_list_with_none]
        chromosome_list = []
        for chromosome_none in chromosome_list_with_none_list:
            chromosome_list.append('\n'.join(chromosome_none))
        print(chromosome_list)
        print("")
        # Writing Chromosome
        print("Step 2: Writing Chromosome, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        chromosome_string = "\n".join(chromosome_list)
        print(chromosome_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_chromosome.txt"), "w")
        handle_outputfile.write(chromosome_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 12: Chromosome Extract Complete")
        print("")
        break

# Clone Extract ##################################################
def clone():
    while True:
        # Module Name
        print("Module 13: Clone Extract")
        print("")
        # Parsing Clone
        print("Step 1: Parsing Clone, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        clone_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        clone = feature.qualifiers.get('clone')
                        clone_list_with_none.append(clone)
        clone_list_with_none_list = [['None'] if clone_value is None else clone_value for clone_value in clone_list_with_none]
        clone_list = []
        for clone_none in clone_list_with_none_list:
            clone_list.append('\n'.join(clone_none))
        print(clone_list)
        print("")
        # Writing Clone
        print("Step 2: Writing Clone, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        clone_string = "\n".join(clone_list)
        print(clone_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_clone.txt"), "w")
        handle_outputfile.write(clone_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 13: Clone Extract Complete")
        print("")
        break

# Clone Lib Extract ##################################################
def clonelib():
    while True:
        # Module Name
        print("Module 14: Clone Lib Extract")
        print("")
        # Parsing Clone Lib
        print("Step 1: Parsing Clone Lib, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        clone_lib_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        clone_lib = feature.qualifiers.get('clone_lib')
                        clone_lib_list_with_none.append(clone_lib)
        clone_lib_list_with_none_list = [['None'] if clone_lib_value is None else clone_lib_value for clone_lib_value in clone_lib_list_with_none]
        clone_lib_list = []
        for clone_lib_none in clone_lib_list_with_none_list:
            clone_lib_list.append('\n'.join(clone_lib_none))
        print(clone_lib_list)
        print("")
        # Writing Clone Lib
        print("Step 2: Writing Clone Lib, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        clone_lib_string = "\n".join(clone_lib_list)
        print(clone_lib_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_clone_lib.txt"), "w")
        handle_outputfile.write(clone_lib_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 14: Clone Lib Extract Complete")
        print("")
        break

# Collected By Extract ##################################################
def collectedby():
    while True:
        # Module Name
        print("Module 15: Collected By Extract")
        print("")
        # Parsing Collected By
        print("Step 1: Parsing Collected By, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        collected_by_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        collected_by = feature.qualifiers.get('collected_by')
                        collected_by_list_with_none.append(collected_by)
        collected_by_list_with_none_list = [['None'] if collected_by_value is None else collected_by_value for collected_by_value in collected_by_list_with_none]
        collected_by_list = []
        for collected_by_none in collected_by_list_with_none_list:
            collected_by_list.append('\n'.join(collected_by_none))
        print(collected_by_list)
        print("")
        # Writing Collected By
        print("Step 2: Writing Collected By, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        collected_by_string = "\n".join(collected_by_list)
        print(collected_by_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_collected_by.txt"), "w")
        handle_outputfile.write(collected_by_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 15: Collected By Extract Complete")
        print("")
        break

# Collection Date Extract ##################################################
def collectiondate():
    while True:
        # Module Name
        print("Module 16: Collection Date Extract")
        print("")
        # Parsing Collection Date
        print("Step 1: Parsing Collection Date, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        collection_date_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        collection_date = feature.qualifiers.get('collection_date')
                        collection_date_list_with_none.append(collection_date)
        collection_date_list_with_none_list = [['None'] if collection_date_value is None else collection_date_value for collection_date_value in collection_date_list_with_none]
        collection_date_list = []
        for collection_date_none in collection_date_list_with_none_list:
            collection_date_list.append('\n'.join(collection_date_none))
        print(collection_date_list)
        print("")
        # Writing Collection Date
        print("Step 2: Writing Collection Date, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        collection_date_string = "\n".join(collection_date_list)
        print(collection_date_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_collection_date.txt"), "w")
        handle_outputfile.write(collection_date_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 16: Collection Date Extract Complete")
        print("")
        break

# Common Extract ##################################################
def common():
    while True:
        # Module Name
        print("Module 17: Common Extract")
        print("")
        # Parsing Common
        print("Step 1: Parsing Common, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        common_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        common = feature.qualifiers.get('common')
                        common_list_with_none.append(common)
        common_list_with_none_list = [['None'] if common_value is None else common_value for common_value in common_list_with_none]
        common_list = []
        for common_none in common_list_with_none_list:
            common_list.append('\n'.join(common_none))
        print(common_list)
        print("")
        # Writing Common
        print("Step 2: Writing Common, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        common_string = "\n".join(common_list)
        print(common_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_common.txt"), "w")
        handle_outputfile.write(common_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 17: Common Extract Complete")
        print("")
        break

# Country Extract ##################################################
def country():
    while True:
        # Module Name
        print("Module 18: Country Extract")
        print("")
        # Parsing Country
        print("Step 1: Parsing Country, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        country_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        country = feature.qualifiers.get('country')
                        country_list_with_none.append(country)
        country_list_with_none_list = [['None'] if country_value is None else country_value for country_value in country_list_with_none]
        country_list = []
        for country_none in country_list_with_none_list:
            country_list.append('\n'.join(country_none))
        print(country_list)
        print("")
        # Writing Country
        print("Step 2: Writing Country, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        country_string = "\n".join(country_list)
        print(country_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_country.txt"), "w")
        handle_outputfile.write(country_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 18: Country Extract Complete")
        print("")
        break

# Country ACCN Extract ##################################################
def countryaccn():
    while True:
        # Module Name
        print("Module 18.1: Country ACCN Extract")
        print("")
        # Parsing Country ACCN
        print("Step 1: Parsing Country ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        country_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        country = feature.qualifiers.get('country')
                        accession_list_with_none.append(record.id)
                        country_list_with_none.append(country)
        country_list_with_none_list = [['None'] if country_value is None else country_value for country_value in country_list_with_none]
        #country_list = []
        country_list_with_none_list_none = []
        for country_none in country_list_with_none_list:
            country_list_with_none_list_none.append('\n'.join(country_none))
        country_zip = zip(accession_list_with_none,country_list_with_none_list_none)
        print(country_zip)
        country_list = []
        for accession_country in country_zip:
            country_list.append("\t".join(accession_country))
        #print(country_list)
        print("")
        # Writing Country ACCN
        print("Step 2: Writing Country ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #country_string = "\n".join(country_list)
        country_string = "\n".join(str(accesion_country_value) for accesion_country_value in country_list)
        print(country_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_country_accn.txt"), "w")
        handle_outputfile.write(str(country_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 18.1: Country ACCN Extract Complete")
        print("")
        break

# Cultivar Extract ##################################################
def cultivar():
    while True:
        # Module Name
        print("Module 19: Cultivar Extract")
        print("")
        # Parsing Cultivar
        print("Step 1: Parsing Cultivar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        cultivar_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        cultivar = feature.qualifiers.get('cultivar')
                        cultivar_list_with_none.append(cultivar)
        cultivar_list_with_none_list = [['None'] if cultivar_value is None else cultivar_value for cultivar_value in cultivar_list_with_none]
        cultivar_list = []
        for cultivar_none in cultivar_list_with_none_list:
            cultivar_list.append('\n'.join(cultivar_none))
        print(cultivar_list)
        print("")
        # Writing Cultivar
        print("Step 2: Writing Cultivar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        cultivar_string = "\n".join(cultivar_list)
        print(cultivar_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_cultivar.txt"), "w")
        handle_outputfile.write(cultivar_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 19: Cultivar Extract Complete")
        print("")
        break

# Culture Collection Extract ##################################################
def culturecollection():
    while True:
        # Module Name
        print("Module 20: Culture Collection Extract")
        print("")
        # Parsing Culture Collection
        print("Step 1: Parsing Culture Collection, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        culture_collection_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        culture_collection = feature.qualifiers.get('culture_collection')
                        culture_collection_list_with_none.append(culture_collection)
        culture_collection_list_with_none_list = [['None'] if culture_collection_value is None else culture_collection_value for culture_collection_value in culture_collection_list_with_none]
        culture_collection_list = []
        for culture_collection_none in culture_collection_list_with_none_list:
            culture_collection_list.append('\n'.join(culture_collection_none))
        print(culture_collection_list)
        print("")
        # Writing Culture Collection
        print("Step 2: Writing Culture Collection, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        culture_collection_string = "\n".join(culture_collection_list)
        print(culture_collection_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_culture_collection.txt"), "w")
        handle_outputfile.write(culture_collection_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 20: Culture Collection Extract Complete")
        print("")
        break

# Dev Stage Extract ##################################################
def devstage():
    while True:
        # Module Name
        print("Module 21: Dev Stage Extract")
        print("")
        # Parsing Dev Stage
        print("Step 1: Parsing Dev Stage, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        dev_stage_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        dev_stage = feature.qualifiers.get('dev_stage')
                        dev_stage_list_with_none.append(dev_stage)
        dev_stage_list_with_none_list = [['None'] if dev_stage_value is None else dev_stage_value for dev_stage_value in dev_stage_list_with_none]
        dev_stage_list = []
        for dev_stage_none in dev_stage_list_with_none_list:
            dev_stage_list.append('\n'.join(dev_stage_none))
        print(dev_stage_list)
        print("")
        # Writing Dev Stage
        print("Step 2: Writing Dev Stage, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        dev_stage_string = "\n".join(dev_stage_list)
        print(dev_stage_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_dev_stage.txt"), "w")
        handle_outputfile.write(dev_stage_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 21: Dev Stage Extract Complete")
        print("")
        break

# Ecotype Extract ##################################################
def ecotype():
    while True:
        # Module Name
        print("Module 22: Ecotype Extract")
        print("")
        # Parsing Ecotype
        print("Step 1: Parsing Ecotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        ecotype_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        ecotype = feature.qualifiers.get('ecotype')
                        ecotype_list_with_none.append(ecotype)
        ecotype_list_with_none_list = [['None'] if ecotype_value is None else ecotype_value for ecotype_value in ecotype_list_with_none]
        ecotype_list = []
        for ecotype_none in ecotype_list_with_none_list:
            ecotype_list.append('\n'.join(ecotype_none))
        print(ecotype_list)
        print("")
        # Writing Ecotype
        print("Step 2: Writing Ecotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        ecotype_string = "\n".join(ecotype_list)
        print(ecotype_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_ecotype.txt"), "w")
        handle_outputfile.write(ecotype_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 22: Ecotype Extract Complete")
        print("")
        break

# Endogenous Virus Name Extract ##################################################
def endogenousvirusname():
    while True:
        # Module Name
        print("Module 23: Endogenous Virus Name Extract")
        print("")
        # Parsing Endogenous Virus Name
        print("Step 1: Parsing Endogenous Virus Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        endogenous_virus_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        endogenous_virus_name = feature.qualifiers.get('endogenous_virus_name')
                        endogenous_virus_name_list_with_none.append(endogenous_virus_name)
        endogenous_virus_name_list_with_none_list = [['None'] if endogenous_virus_name_value is None else endogenous_virus_name_value for endogenous_virus_name_value in endogenous_virus_name_list_with_none]
        endogenous_virus_name_list = []
        for endogenous_virus_name_none in endogenous_virus_name_list_with_none_list:
            endogenous_virus_name_list.append('\n'.join(endogenous_virus_name_none))
        print(endogenous_virus_name_list)
        print("")
        # Writing Endogenous Virus Name
        print("Step 2: Writing Endogenous Virus Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        endogenous_virus_name_string = "\n".join(endogenous_virus_name_list)
        print(endogenous_virus_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_endogenous_virus_name.txt"), "w")
        handle_outputfile.write(endogenous_virus_name_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 2 : Ecotype Endogenous Virus Name Complete")
        print("")
        break

# Environmental Sample Extract ##################################################
def environmentalsample():
    while True:
        # Module Name
        print("Module 24: Environmental Sample Extract")
        print("")
        # Parsing Environmental Sample
        print("Step 1: Parsing Environmental Sample, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        environmental_sample_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        environmental_sample = feature.qualifiers.get('environmental_sample')
                        environmental_sample_list_with_none.append(environmental_sample)
        environmental_sample_list_with_none_list = [['None'] if environmental_sample_value is None else environmental_sample_value for environmental_sample_value in environmental_sample_list_with_none]
        environmental_sample_list = []
        for environmental_sample_none in environmental_sample_list_with_none_list:
            environmental_sample_list.append('\n'.join(environmental_sample_none))
        print(environmental_sample_list)
        print("")
        # Writing Environmental Sample
        print("Step 2: Writing Environmental Sample, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        environmental_sample_string = "\n".join(environmental_sample_list)
        print(environmental_sample_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_environmental_sample.txt"), "w")
        handle_outputfile.write(environmental_sample_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 24: Environmental Sample Complete")
        print("")
        break

# Focus Extract ##################################################
def focus():
    while True:
        # Module Name
        print("Module 25: Focus Extract")
        print("")
        # Parsing Focus
        print("Step 1: Parsing Focus, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        focus_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        focus = feature.qualifiers.get('focus')
                        focus_list_with_none.append(focus)
        focus_list_with_none_list = [['None'] if focus_value is None else focus_value for focus_value in focus_list_with_none]
        focus_list = []
        for focus_none in focus_list_with_none_list:
            focus_list.append('\n'.join(focus_none))
        print(focus_list)
        print("")
        # Writing Focus
        print("Step 2: Writing Focus, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        focus_string = "\n".join(focus_list)
        print(focus_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_focus.txt"), "w")
        handle_outputfile.write(focus_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 25: Focus Complete")
        print("")
        break

# Forma Extract ##################################################
def forma():
    while True:
        # Module Name
        print("Module 26: Forma Extract")
        print("")
        # Parsing Forma
        print("Step 1: Parsing Forma, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        forma_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        forma = feature.qualifiers.get('forma')
                        forma_list_with_none.append(forma)
        forma_list_with_none_list = [['None'] if forma_value is None else forma_value for forma_value in forma_list_with_none]
        forma_list = []
        for forma_none in forma_list_with_none_list:
            forma_list.append('\n'.join(forma_none))
        print(forma_list)
        print("")
        # Writing Forma
        print("Step 2: Writing Forma, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        forma_string = "\n".join(forma_list)
        print(forma_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_forma.txt"), "w")
        handle_outputfile.write(forma_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 26: Forma Complete")
        print("")
        break

# Forma Specialis Extract ##################################################
def formaspecialis():
    while True:
        # Module Name
        print("Module 27: Forma Specialis Extract")
        print("")
        # Parsing Forma Specialis
        print("Step 1: Parsing Forma Specialis, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        forma_specialis_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        forma_specialis = feature.qualifiers.get('forma_specialis')
                        forma_specialis_list_with_none.append(forma_specialis)
        forma_specialis_list_with_none_list = [['None'] if forma_specialis_value is None else forma_specialis_value for forma_specialis_value in forma_specialis_list_with_none]
        forma_specialis_list = []
        for forma_specialis_none in forma_specialis_list_with_none_list:
            forma_specialis_list.append('\n'.join(forma_specialis_none))
        print(forma_specialis_list)
        print("")
        # Writing Forma Specialis
        print("Step 2: Writing Forma Specialis, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        forma_specialis_string = "\n".join(forma_specialis_list)
        print(forma_specialis_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_forma_specialis.txt"), "w")
        handle_outputfile.write(forma_specialis_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 27: Forma Specialis Complete")
        print("")
        break

# Fwd PCR Primer Name Extract ##################################################
def fwdpcrprimername():
    while True:
        # Module Name
        print("Module 28: Fwd PCR Primer Name Extract")
        print("")
        # Parsing Fwd PCR Primer Name
        print("Step 1: Parsing Fwd PCR Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        fwd_pcr_primer_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        fwd_pcr_primer_name = feature.qualifiers.get('fwd_pcr_primer_name')
                        fwd_pcr_primer_name_list_with_none.append(fwd_pcr_primer_name)
        fwd_pcr_primer_name_list_with_none_list = [['None'] if fwd_pcr_primer_name_value is None else fwd_pcr_primer_name_value for fwd_pcr_primer_name_value in fwd_pcr_primer_name_list_with_none]
        fwd_pcr_primer_name_list = []
        for fwd_pcr_primer_name_none in fwd_pcr_primer_name_list_with_none_list:
            fwd_pcr_primer_name_list.append('\n'.join(fwd_pcr_primer_name_none))
        print(fwd_pcr_primer_name_list)
        print("")
        # Writing Fwd PCR Primer Name
        print("Step 2: Writing Fwd PCR Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        fwd_pcr_primer_name_string = "\n".join(fwd_pcr_primer_name_list)
        print(fwd_pcr_primer_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_fwd_pcr_primer_name.txt"), "w")
        handle_outputfile.write(fwd_pcr_primer_name_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 28: Fwd PCR Primer Name Complete")
        print("")
        break

# Fwd PCR Primer Seq Extract ##################################################
def fwdpcrprimerseq():
    while True:
        # Module Name
        print("Module 29: Fwd PCR Primer Seq Extract")
        print("")
        # Parsing Fwd PCR Primer Seq
        print("Step 1: Parsing Fwd PCR Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        fwd_pcr_primer_seq_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        fwd_pcr_primer_seq = feature.qualifiers.get('fwd_pcr_primer_seq')
                        fwd_pcr_primer_seq_list_with_none.append(fwd_pcr_primer_seq)
        fwd_pcr_primer_seq_list_with_none_list = [['None'] if fwd_pcr_primer_seq_value is None else fwd_pcr_primer_seq_value for fwd_pcr_primer_seq_value in fwd_pcr_primer_seq_list_with_none]
        fwd_pcr_primer_seq_list = []
        for fwd_pcr_primer_seq_none in fwd_pcr_primer_seq_list_with_none_list:
            fwd_pcr_primer_seq_list.append('\n'.join(fwd_pcr_primer_seq_none))
        print(fwd_pcr_primer_seq_list)
        print("")
        # Writing Fwd PCR Primer Seq
        print("Step 2: Writing Fwd PCR Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        fwd_pcr_primer_seq_string = "\n".join(fwd_pcr_primer_seq_list)
        print(fwd_pcr_primer_seq_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_fwd_pcr_primer_seq.txt"), "w")
        handle_outputfile.write(fwd_pcr_primer_seq_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 29: Fwd PCR Primer Seq Complete")
        print("")
        break

# Fwd Primer Name Extract ##################################################
def fwdprimername():
    while True:
        # Module Name
        print("Module 30: Fwd Primer Name Extract")
        print("")
        # Parsing Fwd Primer Name
        print("Step 1: Parsing Fwd Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        fwd_primer_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        fwd_primer_name = feature.qualifiers.get('fwd_primer_name')
                        fwd_primer_name_list_with_none.append(fwd_primer_name)
        fwd_primer_name_list_with_none_list = [['None'] if fwd_primer_name_value is None else fwd_primer_name_value for fwd_primer_name_value in fwd_primer_name_list_with_none]
        fwd_primer_name_list = []
        for fwd_primer_name_none in fwd_primer_name_list_with_none_list:
            fwd_primer_name_list.append('\n'.join(fwd_primer_name_none))
        print(fwd_primer_name_list)
        print("")
        # Writing Fwd Primer Name
        print("Step 2: Writing Fwd Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        fwd_primer_name_string = "\n".join(fwd_primer_name_list)
        print(fwd_primer_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_fwd_primer_name.txt"), "w")
        handle_outputfile.write(fwd_primer_name_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 30: Fwd Primer Name Complete")
        print("")
        break

# Fwd Primer Seq Extract ##################################################
def fwdprimerseq():
    while True:
        # Module Name
        print("Module 31: Fwd Primer Seq Extract")
        print("")
        # Parsing Fwd Primer Seq
        print("Step 1: Parsing Fwd Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        fwd_primer_seq_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        fwd_primer_seq = feature.qualifiers.get('fwd_primer_seq')
                        fwd_primer_seq_list_with_none.append(fwd_primer_seq)
        fwd_primer_seq_list_with_none_list = [['None'] if fwd_primer_seq_value is None else fwd_primer_seq_value for fwd_primer_seq_value in fwd_primer_seq_list_with_none]
        fwd_primer_seq_list = []
        for fwd_primer_seq_none in fwd_primer_seq_list_with_none_list:
            fwd_primer_seq_list.append('\n'.join(fwd_primer_seq_none))
        print(fwd_primer_seq_list)
        print("")
        # Writing Fwd Primer Seq
        print("Step 2: Writing Fwd Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        fwd_primer_seq_string = "\n".join(fwd_primer_seq_list)
        print(fwd_primer_seq_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_fwd_primer_seq.txt"), "w")
        handle_outputfile.write(fwd_primer_seq_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 31: Fwd Primer Seq Complete")
        print("")
        break

# Gene Extract ##################################################
def gene():
    while True:
        # Module Name
        print("Module 32: Gene Extract")
        print("")
        # Parsing Gene
        print("Step 1: Parsing Gene, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        gene_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        gene = feature.qualifiers.get('gene')
                        gene_list_with_none.append(gene)
        gene_list_with_none_list = [['None'] if gene_value is None else gene_value for gene_value in gene_list_with_none]
        gene_list = []
        for gene_none in gene_list_with_none_list:
            gene_list.append('\n'.join(gene_none))
        print(gene_list)
        print("")
        # Writing Gene
        print("Step 2: Writing Gene, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        gene_string = "\n".join(gene_list)
        print(gene_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_gene.txt"), "w")
        handle_outputfile.write(gene_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 32: Gene Complete")
        print("")
        break

# Genotype Extract ##################################################
def genotype():
    while True:
        # Module Name
        print("Module 33: Genotype Extract")
        print("")
        # Parsing Gene
        print("Step 1: Parsing Genotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        genotype_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        genotype = feature.qualifiers.get('genotype')
                        genotype_list_with_none.append(genotype)
        genotype_list_with_none_list = [['None'] if genotype_value is None else genotype_value for genotype_value in genotype_list_with_none]
        genotype_list = []
        for genotype_none in genotype_list_with_none_list:
            genotype_list.append('\n'.join(genotype_none))
        print(genotype_list)
        print("")
        # Writing Genotype
        print("Step 2: Writing Genotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        genotype_string = "\n".join(genotype_list)
        print(genotype_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_genotype.txt"), "w")
        handle_outputfile.write(genotype_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 33: Genotype Complete")
        print("")
        break

# Germline Extract ##################################################
def germline():
    while True:
        # Module Name
        print("Module 34: Germline Extract")
        print("")
        # Parsing Germline
        print("Step 1: Parsing Germline, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        germline_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        germline = feature.qualifiers.get('germline')
                        germline_list_with_none.append(germline)
        germline_list_with_none_list = [['None'] if germline_value is None else germline_value for germline_value in germline_list_with_none]
        germline_list = []
        for germline_none in germline_list_with_none_list:
            germline_list.append('\n'.join(germline_none))
        print(germline_list)
        print("")
        # Writing Germline
        print("Step 2: Writing Germline, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        germline_string = "\n".join(germline_list)
        print(germline_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_germline.txt"), "w")
        handle_outputfile.write(germline_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 34: Germline Complete")
        print("")
        break

# Group Extract ##################################################
def group():
    while True:
        # Module Name
        print("Module 35: Group Extract")
        print("")
        # Parsing Group
        print("Step 1: Parsing Group, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        group_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        group = feature.qualifiers.get('group')
                        group_list_with_none.append(group)
        group_list_with_none_list = [['None'] if group_value is None else group_value for group_value in group_list_with_none]
        group_list = []
        for group_none in group_list_with_none_list:
            group_list.append('\n'.join(group_none))
        print(group_list)
        print("")
        # Writing Group
        print("Step 2: Writing Group, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        group_string = "\n".join(group_list)
        print(group_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_group.txt"), "w")
        handle_outputfile.write(group_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 35: Group Complete")
        print("")
        break

# Haplogroup Extract ##################################################
def haplogroup():
    while True:
        # Module Name
        print("Module 36: Haplogroup Extract")
        print("")
        # Parsing Haplogroup
        print("Step 1: Parsing Haplogroup, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        haplogroup_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        haplogroup = feature.qualifiers.get('haplogroup')
                        haplogroup_list_with_none.append(haplogroup)
        haplogroup_list_with_none_list = [['None'] if haplogroup_value is None else haplogroup_value for haplogroup_value in haplogroup_list_with_none]
        haplogroup_list = []
        for haplogroup_none in haplogroup_list_with_none_list:
            haplogroup_list.append('\n'.join(haplogroup_none))
        print(haplogroup_list)
        print("")
        # Writing Haplogroup
        print("Step 2: Writing Haplogroup, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        haplogroup_string = "\n".join(haplogroup_list)
        print(haplogroup_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_haplogroup.txt"), "w")
        handle_outputfile.write(haplogroup_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 36: Haplogroup Complete")
        print("")
        break

# Haplotype Extract ##################################################
def haplotype():
    while True:
        # Module Name
        print("Module 37: Haplotype Extract")
        print("")
        # Parsing Haplotype
        print("Step 1: Parsing Haplotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        haplotype_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        haplotype = feature.qualifiers.get('haplotype')
                        haplotype_list_with_none.append(haplotype)
        haplotype_list_with_none_list = [['None'] if haplotype_value is None else haplotype_value for haplotype_value in haplotype_list_with_none]
        haplotype_list = []
        for haplotype_none in haplotype_list_with_none_list:
            haplotype_list.append('\n'.join(haplotype_none))
        print(haplotype_list)
        print("")
        # Writing Haplotype
        print("Step 2: Writing Haplotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        haplotype_string = "\n".join(haplotype_list)
        print(haplotype_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_haplotype.txt"), "w")
        handle_outputfile.write(haplotype_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 37: Haplotype Complete")
        print("")
        break

# Host Extract ##################################################
def host():
    while True:
        # Module Name
        print("Module 38: Host Extract")
        print("")
        # Parsing Host
        print("Step 1: Parsing Host, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        host_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        host = feature.qualifiers.get('host')
                        host_list_with_none.append(host)
        host_list_with_none_list = [['None'] if host_value is None else host_value for host_value in host_list_with_none]
        host_list = []
        for host_none in host_list_with_none_list:
            host_list.append('\n'.join(host_none))
        print(host_list)
        print("")
        # Writing Host
        print("Step 2: Writing Host, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        host_string = "\n".join(host_list)
        print(host_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_host.txt"), "w")
        handle_outputfile.write(host_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 38: Host Complete")
        print("")
        break

# Host ACCN Extract ##################################################
def hostaccn():
    while True:
        # Module Name
        print("Module 38.1: Host ACCN Extract")
        print("")
        # Parsing Host ACCN
        print("Step 1: Parsing Host ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        host_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        host = feature.qualifiers.get('host')
                        accession_list_with_none.append(record.id)
                        host_list_with_none.append(host)
        host_list_with_none_list = [['None'] if host_value is None else host_value for host_value in host_list_with_none]
        #host_list = []
        host_list_with_none_list_none = []
        for host_none in host_list_with_none_list:
            host_list_with_none_list_none.append('\n'.join(host_none))
        host_zip = zip(accession_list_with_none,host_list_with_none_list_none)
        print(host_zip)
        host_list = []
        for accession_host in host_zip:
            host_list.append("\t".join(accession_host))
        #print(host_list)
        print("")
        # Writing Host
        print("Step 2: Writing Host ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #host_string = "\n".join(host_list)
        host_string = "\n".join(str(accesion_host_value) for accesion_host_value in host_list)
        print(host_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_host_accn.txt"), "w")
        handle_outputfile.write(str(host_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 38.1: Host ACCN Complete")
        print("")
        break

# Identified By Extract ##################################################
def identifiedby():
    while True:
        # Module Name
        print("Module 39: Identified By Extract")
        print("")
        # Parsing Identified By
        print("Step 1: Parsing Identified By, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        identified_by_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        identified_by = feature.qualifiers.get('identified_by')
                        identified_by_list_with_none.append(identified_by)
        identified_by_list_with_none_list = [['None'] if identified_by_value is None else identified_by_value for identified_by_value in identified_by_list_with_none]
        identified_by_list = []
        for identified_by_none in identified_by_list_with_none_list:
            identified_by_list.append('\n'.join(identified_by_none))
        print(identified_by_list)
        print("")
        # Writing Identified By
        print("Step 2: Writing Identified By, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        identified_by_string = "\n".join(identified_by_list)
        print(identified_by_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_identified_by.txt"), "w")
        handle_outputfile.write(identified_by_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 39: Identified By Complete")
        print("")
        break

# Isolate Extract ##################################################
def isolate():
    while True:
        # Module Name
        print("Module 40: Isolate Extract")
        print("")
        # Parsing Isolate
        print("Step 1: Parsing Isolate, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        isolate_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        isolate = feature.qualifiers.get('isolate')
                        isolate_list_with_none.append(isolate)
        isolate_list_with_none_list = [['None'] if isolate_value is None else isolate_value for isolate_value in isolate_list_with_none]
        isolate_list = []
        for isolate_none in isolate_list_with_none_list:
            isolate_list.append('\n'.join(isolate_none))
        print(isolate_list)
        print("")
        # Writing Isolate
        print("Step 2: Writing Isolate, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        isolate_string = "\n".join(isolate_list)
        print(isolate_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_isolate.txt"), "w")
        handle_outputfile.write(isolate_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 40: Isolate Complete")
        print("")
        break

# Isolate ACCN Extract ##################################################
def isolateaccn():
    while True:
        # Module Name
        print("Module 40.1: Isolate ACCN Extract")
        print("")
        # Parsing Isolate ACCN
        print("Step 1: Parsing Isolate ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        isolate_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        isolate = feature.qualifiers.get('isolate')
                        accession_list_with_none.append(record.id)
                        isolate_list_with_none.append(isolate)
        isolate_list_with_none_list = [['None'] if isolate_value is None else isolate_value for isolate_value in isolate_list_with_none]
        #isolate_list = []
        isolate_list_with_none_list_none = []
        for isolate_none in isolate_list_with_none_list:
            isolate_list_with_none_list_none.append('\n'.join(isolate_none))
        isolate_zip = zip(accession_list_with_none,isolate_list_with_none_list_none)
        print(isolate_zip)
        isolate_list = []
        for accession_isolate in isolate_zip:
            isolate_list.append("\t".join(accession_isolate))
        #print(isolate_list)
        print("")
        # Writing Isolate ACCN
        print("Step 2: Writing Isolate ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #isolate_string = "\n".join(isolate_list)
        isolate_string = "\n".join(str(accesion_isolate_value) for accesion_isolate_value in isolate_list)
        print(isolate_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_isolate_accn.txt"), "w")
        handle_outputfile.write(str(isolate_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 40.1: Isolate ACCN Complete")
        print("")
        break

# Isolation Source Extract ##################################################
def isolationsource():
    while True:
        # Module Name
        print("Module 41: Isolation Source Extract")
        print("")
        # Parsing Isolation Source
        print("Step 1: Parsing Isolation Source, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        isolation_source_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        isolation_source = feature.qualifiers.get('isolation_source')
                        isolation_source_list_with_none.append(isolation_source)
        isolation_source_list_with_none_list = [['None'] if isolation_source_value is None else isolation_source_value for isolation_source_value in isolation_source_list_with_none]
        isolation_source_list = []
        for isolation_source_none in isolation_source_list_with_none_list:
            isolation_source_list.append('\n'.join(isolation_source_none))
        print(isolation_source_list)
        print("")
        # Writing Isolation Source
        print("Step 2: Writing Isolation Source, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        isolation_source_string = "\n".join(isolation_source_list)
        print(isolation_source_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_isolation_source.txt"), "w")
        handle_outputfile.write(isolation_source_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 41: Isolation Source Complete")
        print("")
        break

# Isolation Source ACCN Extract ##################################################
def isolationsourceaccn():
    while True:
        # Module Name
        print("Module 41.1: Isolation Source ACCN Extract")
        print("")
        # Parsing Isolation Source ACCN
        print("Step 1: Parsing Isolation Source ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        isolation_source_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        isolation_source = feature.qualifiers.get('isolation_source')
                        accession_list_with_none.append(record.id)
                        isolation_source_list_with_none.append(isolation_source)
        isolation_source_list_with_none_list = [['None'] if isolation_source_value is None else isolation_source_value for isolation_source_value in isolation_source_list_with_none]
        #isolation_source_list = []
        isolation_source_list_with_none_list_none = []
        for isolation_source_none in isolation_source_list_with_none_list:
            isolation_source_list_with_none_list_none.append('\n'.join(isolation_source_none))
        isolation_source_zip = zip(accession_list_with_none,isolation_source_list_with_none_list_none)
        print(isolation_source_zip)
        isolation_source_list = []
        for accession_isolation_source in isolation_source_zip:
            isolation_source_list.append("\t".join(accession_isolation_source))
        #print(isolation_source_list)
        print("")
        # Writing Isolation Source ACCN
        print("Step 2: Writing Isolation Source ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #isolation_source_string = "\n".join(isolation_source_list)
        isolation_source_string = "\n".join(str(accesion_isolation_source_value) for accesion_isolation_source_value in isolation_source_list)
        print(isolation_source_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_isolation_source_accn.txt"), "w")
        handle_outputfile.write(str(isolation_source_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 41.1: Isolation Source ACCN Complete")
        print("")
        break

# Lab Host Extract ##################################################
def labhost():
    while True:
        # Module Name
        print("Module 42: Lab Host Extract")
        print("")
        # Parsing Lab Host
        print("Step 1: Parsing Lab Host, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        lab_host_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        lab_host = feature.qualifiers.get('lab_host')
                        lab_host_list_with_none.append(lab_host)
        lab_host_list_with_none_list = [['None'] if lab_host_value is None else lab_host_value for lab_host_value in lab_host_list_with_none]
        lab_host_list = []
        for lab_host_none in lab_host_list_with_none_list:
            lab_host_list.append('\n'.join(lab_host_none))
        print(lab_host_list)
        print("")
        # Writing Lab Host
        print("Step 2: Writing Lab Host, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        lab_host_string = "\n".join(lab_host_list)
        print(lab_host_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_lab_host.txt"), "w")
        handle_outputfile.write(lab_host_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 42: Lab Host Complete")
        print("")
        break

# Lat Lon Extract ##################################################
def latlon():
    while True:
        # Module Name
        print("Module 43: Lat Lon Extract")
        print("")
        # Parsing Lat Lon
        print("Step 1: Parsing Lat Lon, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        lat_lon_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        lat_lon = feature.qualifiers.get('lat_lon')
                        lat_lon_list_with_none.append(lat_lon)
        lat_lon_list_with_none_list = [['None'] if lat_lon_value is None else lat_lon_value for lat_lon_value in lat_lon_list_with_none]
        lat_lon_list = []
        for lat_lon_none in lat_lon_list_with_none_list:
            lat_lon_list.append('\n'.join(lat_lon_none))
        print(lat_lon_list)
        print("")
        # Writing Lat Lon
        print("Step 2: Writing Lat Lon, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        lat_lon_string = "\n".join(lat_lon_list)
        print(lat_lon_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_lat_lon.txt"), "w")
        handle_outputfile.write(lat_lon_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 43: Lat Lon Complete")
        print("")
        break

# Latitude\Longitude Extract ##################################################
def latitudelongitude():
    while True:
        # Module Name
        print("Module 44: Latitude\Longitude (GeoPy) Extract")
        print("")
        # Parsing Latitude\Longitude (GeoPy) 
        print("Step 1: Parsing Latitude\Longitude (GeoPy) , please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_inputfile = open(os.path.expanduser("~/extractor/local/export/source_country.txt"))
        latitude_list = []
        longitude_list = []
        for row in handle_inputfile:
            row = row.rstrip()
            geolocator = Nominatim()
            latitude_longitude = geolocator.geocode(row)
            latitude_list.append(str(latitude_longitude.latitude))
            longitude_list.append(str(latitude_longitude.longitude))
        # Zipping Latitude/Longitude (GeoPy) 
        latitude_longitude_zip = zip(latitude_list,longitude_list)
        print(latitude_longitude_zip)
        latitude_longitude_list = []
        for latitude_longitude in latitude_longitude_zip:
            latitude_longitude_list.append("\t".join(latitude_longitude))
        print("")
        # Writing Latitude/Longitude (GeoPy) 
        print("Step 2: Writing Latitude/Longitude (GeoPy) , please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        latitude_longitude_string = '\n'.join(str(latitude_longitude_value) for latitude_longitude_value in latitude_longitude_list)
        print(latitude_longitude_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_latitude_longitude_geopy.txt"), "w")
        handle_outputfile.write(str(latitude_longitude_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 44: Latitud/Longitude (GeoPy) Extract Complete")
        print("")
        break

# Map Extract ##################################################
def mapsource():
    while True:
        # Module Name
        print("Module 45: Map Extract")
        print("")
        # Parsing Map
        print("Step 1: Parsing Map, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        map_source_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        map_source = feature.qualifiers.get('map')
                        map_source_list_with_none.append(map_source)
        map_source_list_with_none_list = [['None'] if map_source_value is None else map_source_value for map_source_value in map_source_list_with_none]
        map_source_list = []
        for map_source_none in map_source_list_with_none_list:
            map_source_list.append('\n'.join(map_source_none))
        print(map_source_list)
        print("")
        # Writing Map
        print("Step 2: Writing Map, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        map_source_string = "\n".join(map_source_list)
        print(map_source_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_map.txt"), "w")
        handle_outputfile.write(map_source_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 45: Map Complete")
        print("")
        break

# Metagenome Source Extract ##################################################
def metagenomesource():
    while True:
        # Module Name
        print("Module 46: Metagenome Source Extract")
        print("")
        # Parsing Metagenome Source
        print("Step 1: Parsing Metagenome Source, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        metagenome_source_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        metagenome_source = feature.qualifiers.get('metagenome_source')
                        metagenome_source_list_with_none.append(metagenome_source)
        metagenome_source_list_with_none_list = [['None'] if metagenome_source_value is None else metagenome_source_value for metagenome_source_value in metagenome_source_list_with_none]
        metagenome_source_list = []
        for metagenome_source_none in metagenome_source_list_with_none_list:
            metagenome_source_list.append('\n'.join(metagenome_source_none))
        print(metagenome_source_list)
        print("")
        # Writing Metagenome Source
        print("Step 2: Writing Metagenome Source, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        metagenome_source_string = "\n".join(metagenome_source_list)
        print(metagenome_source_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_metagenome_source.txt"), "w")
        handle_outputfile.write(metagenome_source_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 46: Metagenome Source Complete")
        print("")
        break

# Metagenomic Extract ##################################################
def metagenomic():
    while True:
        # Module Name
        print("Module 47: Metagenomic Extract")
        print("")
        # Parsing Metagenomic
        print("Step 1: Parsing Metagenomic, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        metagenomic_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        metagenomic = feature.qualifiers.get('metagenomic')
                        metagenomic_list_with_none.append(metagenomic)
        metagenomic_list_with_none_list = [['None'] if metagenomic_value is None else metagenomic_value for metagenomic_value in metagenomic_list_with_none]
        metagenomic_list = []
        for metagenomic_none in metagenomic_list_with_none_list:
            metagenomic_list.append('\n'.join(metagenomic_none))
        print(metagenomic_list)
        print("")
        # Writing Metagenomic
        print("Step 2: Writing Metagenomic, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        metagenomic_string = "\n".join(metagenomic_list)
        print(metagenomic_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_metagenomic.txt"), "w")
        handle_outputfile.write(metagenomic_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 47: Metagenomic Complete")
        print("")
        break

# Note Extract ##################################################
def note():
    while True:
        # Module Name
        print("Module 48: Note Extract")
        print("")
        # Parsing Note
        print("Step 1: Parsing Note, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        note_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        note = feature.qualifiers.get('note')
                        note_list_with_none.append(note)
        note_list_with_none_list = [['None'] if note_value is None else note_value for note_value in note_list_with_none]
        note_list = []
        for note_none in note_list_with_none_list:
            note_list.append('\n'.join(note_none))
        print(note_list)
        print("")
        # Writing Note
        print("Step 2: Writing Note, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        note_string = "\n".join(note_list)
        print(note_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_note.txt"), "w")
        handle_outputfile.write(note_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 48: Note Extract Complete")
        print("")
        break

# Organism Extract ##################################################
def organism():
    while True:
        # Module Name
        print("Module 49: Organism Extract")
        print("")
        # Parsing Organism
        print("Step 1: Parsing Organism, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        organism_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        organism = feature.qualifiers.get('organism')
                        organism_list_with_none.append(organism)
        organism_list_with_none_list = [['None'] if organism_value is None else organism_value for organism_value in organism_list_with_none]
        organism_list = []
        for organism_none in organism_list_with_none_list:
            organism_list.append('\n'.join(organism_none))
        print(organism_list)
        print("")
        # Writing Organism
        print("Step 2: Writing Organism, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        organism_string = "\n".join(organism_list)
        print(organism_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_organism.txt"), "w")
        handle_outputfile.write(organism_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 49: Organism Complete")
        print("")
        break

# Pathovar Extract ##################################################
def pathovar():
    while True:
        # Module Name
        print("Module 50: Pathovar Extract")
        print("")
        # Parsing Pathovar
        print("Step 1: Parsing Pathovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        pathovar_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        pathovar = feature.qualifiers.get('pathovar')
                        pathovar_list_with_none.append(pathovar)
        pathovar_list_with_none_list = [['None'] if pathovar_value is None else pathovar_value for pathovar_value in pathovar_list_with_none]
        pathovar_list = []
        for pathovar_none in pathovar_list_with_none_list:
            pathovar_list.append('\n'.join(pathovar_none))
        print(pathovar_list)
        print("")
        # Writing Pathovar
        print("Step 2: Writing Pathovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)

        print("")
        pathovar_string = "\n".join(pathovar_list)
        print(pathovar_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_pathovar.txt"), "w")
        handle_outputfile.write(pathovar_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 50: Pathovar Complete")
        print("")
        break

# Plasmid Name Extract ##################################################
def plasmidname():
    while True:
        # Module Name
        print("Module 51: Plasmid Name Extract")
        print("")
        # Parsing Plasmid Name
        print("Step 1: Parsing Plasmid Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        plasmid_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        plasmid_name = feature.qualifiers.get('plasmid_name')
                        plasmid_name_list_with_none.append(plasmid_name)
        plasmid_name_list_with_none_list = [['None'] if plasmid_name_value is None else plasmid_name_value for plasmid_name_value in plasmid_name_list_with_none]
        plasmid_name_list = []
        for plasmid_name_none in plasmid_name_list_with_none_list:
            plasmid_name_list.append('\n'.join(plasmid_name_none))
        print(plasmid_name_list)
        print("")
        # Writing Plasmid Name
        print("Step 2: Writing Plasmid Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        plasmid_name_string = "\n".join(plasmid_name_list)
        print(plasmid_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_plasmid_name.txt"), "w")
        handle_outputfile.write(plasmid_name_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 51: Plasmid Name Complete")
        print("")
        break

# Plastid Name Extract ##################################################
def plastidname():
    while True:
        # Module Name
        print("Module 52: Plastid Name Extract")
        print("")
        # Parsing Plastid Name
        print("Step 1: Parsing Plastid Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        plastid_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        plastid_name = feature.qualifiers.get('plastid_name')
                        plastid_name_list_with_none.append(plastid_name)
        plastid_name_list_with_none_list = [['None'] if plastid_name_value is None else plastid_name_value for plastid_name_value in plastid_name_list_with_none]
        plastid_name_list = []
        for plastid_name_none in plastid_name_list_with_none_list:
            plastid_name_list.append('\n'.join(plastid_name_none))
        print(plastid_name_list)
        print("")
        # Writing Plastid Name
        print("Step 2: Writing Plastid Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        plastid_name_string = "\n".join(plastid_name_list)
        print(plastid_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_plastid_name.txt"), "w")
        handle_outputfile.write(plastid_name_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 52: Plastid Name Complete")
        print("")
        break

# Pop Variant Extract ##################################################
def popvariant():
    while True:
        # Module Name
        print("Module 53: Pop Variant Extract")
        print("")
        # Parsing Pop Variant
        print("Step 1: Parsing Pop Variant, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        pop_variant_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        pop_variant = feature.qualifiers.get('pop_variant')
                        pop_variant_list_with_none.append(pop_variant)
        pop_variant_list_with_none_list = [['None'] if pop_variant_value is None else pop_variant_value for pop_variant_value in pop_variant_list_with_none]
        pop_variant_list = []
        for pop_variant_none in pop_variant_list_with_none_list:
            pop_variant_list.append('\n'.join(pop_variant_none))
        print(pop_variant_list)
        print("")
        # Writing Pop Variant
        print("Step 2: Writing Pop Variant, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        pop_variant_string = "\n".join(pop_variant_list)
        print(pop_variant_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_pop_variant.txt"), "w")
        handle_outputfile.write(pop_variant_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 53: Pop Variant Complete")
        print("")
        break

# Protein Extract ##################################################
def protein():
    while True:
        # Module Name
        print("Module 54: Protein Extract")
        print("")
        # Parsing Protein
        print("Step 1: Parsing Protein, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        protein_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        protein = feature.qualifiers.get('protein')
                        protein_list_with_none.append(protein)
        protein_list_with_none_list = [['None'] if protein_value is None else protein_value for protein_value in protein_list_with_none]
        protein_list = []
        for protein_none in protein_list_with_none_list:
            protein_list.append('\n'.join(protein_none))
        print(protein_list)
        print("")
        # Writing Protein
        print("Step 2: Writing Protein, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        protein_string = "\n".join(protein_list)
        print(protein_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_protein.txt"), "w")
        handle_outputfile.write(protein_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 54: Protein Complete")
        print("")
        break

# Prot Desc Extract ##################################################
def protdesc():
    while True:
        # Module Name
        print("Module 55: Prot Desc Extract")
        print("")
        # Parsing Prot Desc
        print("Step 1: Parsing Prot Desc, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        prot_desc_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        prot_desc = feature.qualifiers.get('prot_desc')
                        prot_desc_list_with_none.append(prot_desc)
        prot_desc_list_with_none_list = [['None'] if prot_desc_value is None else prot_desc_value for prot_desc_value in prot_desc_list_with_none]
        prot_desc_list = []
        for prot_desc_none in prot_desc_list_with_none_list:
            prot_desc_list.append('\n'.join(prot_desc_none))
        print(prot_desc_list)
        print("")
        # Writing Prot Desc
        print("Step 2: Writing Prot Desc, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        prot_desc_string = "\n".join(prot_desc_list)
        print(prot_desc_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_prot_desc.txt"), "w")
        handle_outputfile.write(prot_desc_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 55: Prot Desc Complete")
        print("")
        break

# Rearranged Extract ##################################################
def rearranged():
    while True:
        # Module Name
        print("Module 56: Rearranged Extract")
        print("")
        # Parsing Rearranged
        print("Step 1: Parsing Rearranged, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        rearranged_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        rearranged = feature.qualifiers.get('rearranged')
                        rearranged_list_with_none.append(rearranged)
        rearranged_list_with_none_list = [['None'] if rearranged_value is None else rearranged_value for rearranged_value in rearranged_list_with_none]
        rearranged_list = []
        for rearranged_none in rearranged_list_with_none_list:
            rearranged_list.append('\n'.join(rearranged_none))
        print(rearranged_list)
        print("")
        # Writing Rearranged
        print("Step 2: Writing Rearranged, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        rearranged_string = "\n".join(rearranged_list)
        print(rearranged_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_rearranged.txt"), "w")
        handle_outputfile.write(rearranged_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 56: Rearranged Complete")
        print("")
        break

# Rev PCR Primer Name Extract ##################################################
def revpcrprimername():
    while True:
        # Module Name
        print("Module 57: Rev PCR Primer Name Extract")
        print("")
        # Parsing Rev PCR Primer Name
        print("Step 1: Parsing Rev PCR Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        rev_pcr_primer_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        rev_pcr_primer_name = feature.qualifiers.get('rev_pcr_primer_name')
                        rev_pcr_primer_name_list_with_none.append(rev_pcr_primer_name)
        rev_pcr_primer_name_list_with_none_list = [['None'] if rev_pcr_primer_name_value is None else rev_pcr_primer_name_value for rev_pcr_primer_name_value in rev_pcr_primer_name_list_with_none]
        rev_pcr_primer_name_list = []
        for rev_pcr_primer_name_none in rev_pcr_primer_name_list_with_none_list:
            rev_pcr_primer_name_list.append('\n'.join(rev_pcr_primer_name_none))
        print(rev_pcr_primer_name_list)
        print("")
        # Writing Rev PCR Primer Name
        print("Step 2: Writing Rev PCR Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        rev_pcr_primer_name_string = "\n".join(rev_pcr_primer_name_list)
        print(rev_pcr_primer_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_rev_pcr_primer_name.txt"), "w")
        handle_outputfile.write(rev_pcr_primer_name_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 57: Rev PCR Primer Name Complete")
        print("")
        break

# Rev PCR Primer Seq Extract ##################################################
def revpcrprimerseq():
    while True:
        # Module Name
        print("Module 58: Rev PCR Primer Seq Extract")
        print("")
        # Parsing Rev PCR Primer Seq
        print("Step 1: Parsing Rev PCR Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        rev_pcr_primer_seq_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        rev_pcr_primer_seq = feature.qualifiers.get('rev_pcr_primer_seq')
                        rev_pcr_primer_seq_list_with_none.append(rev_pcr_primer_seq)
        rev_pcr_primer_seq_list_with_none_list = [['None'] if rev_pcr_primer_seq_value is None else rev_pcr_primer_seq_value for rev_pcr_primer_seq_value in rev_pcr_primer_seq_list_with_none]
        rev_pcr_primer_seq_list = []
        for rev_pcr_primer_seq_none in rev_pcr_primer_seq_list_with_none_list:
            rev_pcr_primer_seq_list.append('\n'.join(rev_pcr_primer_seq_none))
        print(rev_pcr_primer_seq_list)
        print("")
        # Writing Rev PCR Primer Seq
        print("Step 2: Writing Rev PCR Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        rev_pcr_primer_seq_string = "\n".join(rev_pcr_primer_seq_list)
        print(rev_pcr_primer_seq_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_rev_pcr_primer_seq.txt"), "w")
        handle_outputfile.write(rev_pcr_primer_seq_string)
        handle_outputfile.close()
        print("")

        # Module Complete
        print("Module 58: Rev PCR Primer Seq Complete")
        print("")
        break

# Rev Primer Name Extract ##################################################
def revprimername():
    while True:
        # Module Name
        print("Module 59: Rev Primer Name Extract")
        print("")
        # Parsing Rev Primer Name
        print("Step 1: Parsing Rev Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        rev_primer_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        rev_primer_name = feature.qualifiers.get('rev_primer_name')
                        rev_primer_name_list_with_none.append(rev_primer_name)
        rev_primer_name_list_with_none_list = [['None'] if rev_primer_name_value is None else rev_primer_name_value for rev_primer_name_value in rev_primer_name_list_with_none]
        rev_primer_name_list = []
        for rev_primer_name_none in rev_primer_name_list_with_none_list:
            rev_primer_name_list.append('\n'.join(rev_primer_name_none))
        print(rev_primer_name_list)
        print("")
        # Writing Rev Primer Name
        print("Step 2: Writing Rev Primer Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        rev_primer_name_string = "\n".join(rev_primer_name_list)
        print(rev_primer_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_rev_primer_name.txt"), "w")
        handle_outputfile.write(rev_primer_name_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 59: Rev Primer Name Complete")
        print("")
        break

# Rev Primer Seq Extract ##################################################
def revprimerseq():
    while True:
        # Module Name
        print("Module 60: Rev Primer Seq Extract")
        print("")
        # Parsing Rev Primer Seq
        print("Step 1: Parsing Rev Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        rev_primer_seq_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        rev_primer_seq = feature.qualifiers.get('rev_primer_seq')
                        rev_primer_seq_list_with_none.append(rev_primer_seq)
        rev_primer_seq_list_with_none_list = [['None'] if rev_primer_seq_value is None else rev_primer_seq_value for rev_primer_seq_value in rev_primer_seq_list_with_none]
        rev_primer_seq_list = []
        for rev_primer_seq_none in rev_primer_seq_list_with_none_list:
            rev_primer_seq_list.append('\n'.join(rev_primer_seq_none))
        print(rev_primer_seq_list)
        print("")
        # Writing Rev Primer Seq
        print("Step 2: Writing Rev Primer Seq, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        rev_primer_seq_string = "\n".join(rev_primer_seq_list)
        print(rev_primer_seq_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_rev_primer_seq.txt"), "w")
        handle_outputfile.write(rev_primer_seq_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 60: Rev Primer Seq Complete")
        print("")
        break

# Segment Extract ##################################################
def segment():
    while True:
        # Module Name
        print("Module 61: Segment Extract")
        print("")
        # Parsing Segment
        print("Step 1: Parsing Segment, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        segment_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        segment = feature.qualifiers.get('segment')
                        segment_list_with_none.append(segment)
        segment_list_with_none_list = [['None'] if segment_value is None else segment_value for segment_value in segment_list_with_none]
        segment_list = []
        for segment_none in segment_list_with_none_list:
            segment_list.append('\n'.join(segment_none))
        print(segment_list)
        print("")
        # Writing Segment
        print("Step 2: Writing Segment, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        segment_string = "\n".join(segment_list)
        print(segment_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_segment.txt"), "w")
        handle_outputfile.write(segment_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 61: Segment Complete")
        print("")
        break

# Serogroup Extract ##################################################
def serogroup():
    while True:
        # Module Name
        print("Module 62: Serogroup Extract")
        print("")
        # Parsing Serogroup
        print("Step 1 : Parsing Serogroup, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        serogroup_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        serogroup = feature.qualifiers.get('serogroup')
                        serogroup_list_with_none.append(serogroup)
        serogroup_list_with_none_list = [['None'] if serogroup_value is None else serogroup_value for serogroup_value in serogroup_list_with_none]
        serogroup_list = []
        for serogroup_none in serogroup_list_with_none_list:
            serogroup_list.append('\n'.join(serogroup_none))
        print(serogroup_list)
        print("")
        # Writing Serogroup
        print("Step 2: Writing Serogroup, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        serogroup_string = "\n".join(serogroup_list)
        print(serogroup_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_serogroup.txt"), "w")
        handle_outputfile.write(serogroup_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 62: Serogroup Complete")
        print("")
        break

# Serotype Extract ##################################################
def serotype():
    while True:
        # Module Name
        print("Module 63: Serotype Extract")
        print("")
        # Parsing Serotype
        print("Step 1: Parsing Serotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        serotype_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        serotype = feature.qualifiers.get('serotype')
                        serotype_list_with_none.append(serotype)
        serotype_list_with_none_list = [['None'] if serotype_value is None else serotype_value for serotype_value in serotype_list_with_none]
        serotype_list = []
        for serotype_none in serotype_list_with_none_list:
            serotype_list.append('\n'.join(serotype_none))
        print(serotype_list)
        print("")
        # Writing Serotype
        print("Step 2: Writing Serotype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        serotype_string = "\n".join(serotype_list)
        print(serotype_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_serotype.txt"), "w")
        handle_outputfile.write(serotype_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 63: Serotype Complete")
        print("")
        break

# Serotype ACCN Extract ##################################################
def serotypeaccn():
    while True:
        # Module Name
        print("Module 63.1: Serotype ACCN Extract")
        print("")
        # Parsing Serotype ACCN
        print("Step 1: Parsing Serotype ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        serotype_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        serotype = feature.qualifiers.get('serotype')
                        accession_list_with_none.append(record.id)
                        serotype_list_with_none.append(serotype)
        serotype_list_with_none_list = [['None'] if serotype_value is None else serotype_value for serotype_value in serotype_list_with_none]
        #serotype_list = []
        serotype_list_with_none_list_none = []
        for serotype_none in serotype_list_with_none_list:
            serotype_list_with_none_list_none.append('\n'.join(serotype_none))
        serotype_zip = zip(accession_list_with_none,serotype_list_with_none_list_none)
        print(serotype_zip)
        serotype_list = []
        for accession_serotype in serotype_zip:
            serotype_list.append("\t".join(accession_serotype))
        #print(serotype_list)
        print("")
        # Writing Serotype ACCN
        print("Step 2: Writing Serotype ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #serotype_string = "\n".join(serotype_list)
        serotype_string = "\n".join(str(accesion_serotype_value) for accesion_serotype_value in serotype_list)
        print(serotype_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_serotype_accn.txt"), "w")
        handle_outputfile.write(str(serotype_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 63.1: Serotype ACCN Complete")
        print("")
        break

# Serovar Extract ##################################################
def serovar():
    while True:
        # Module Name
        print("Module 64: Serovar Extract")
        print("")
        # Parsing Serovar
        print("Step 1: Parsing Serovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        serovar_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        serovar = feature.qualifiers.get('serovar')
                        serovar_list_with_none.append(serovar)
        serovar_list_with_none_list = [['None'] if serovar_value is None else serovar_value for serovar_value in serovar_list_with_none]
        serovar_list = []
        for serovar_none in serovar_list_with_none_list:
            serovar_list.append('\n'.join(serovar_none))
        print(serovar_list)
        print("")
        # Writing Serovar
        print("Step 2: Writing Serovar, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        serovar_string = "\n".join(serovar_list)
        print(serovar_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_serovar.txt"), "w")
        handle_outputfile.write(serovar_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 64: Serovar Complete")
        print("")
        break

# Sex Extract ##################################################
def sex():
    while True:
        # Module Name
        print("Module 65: Sex Extract")
        print("")
        # Parsing Sex
        print("Step 1: Parsing Sex, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        sex_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        sex = feature.qualifiers.get('sex')
                        sex_list_with_none.append(sex)
        sex_list_with_none_list = [['None'] if sex_value is None else sex_value for sex_value in sex_list_with_none]
        sex_list = []
        for sex_none in sex_list_with_none_list:
            sex_list.append('\n'.join(sex_none))
        print(sex_list)
        print("")
        # Writing Sex
        print("Step 2: Writing Sex, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        sex_string = "\n".join(sex_list)
        print(sex_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_sex.txt"), "w")
        handle_outputfile.write(sex_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 65: Sex Complete")
        print("")
        break

# Specific Host Extract ##################################################
def specifichost():
    while True:
        # Module Name
        print("Module 66: Specific Host Extract")
        print("")
        # Parsing Specific Host
        print("Step 1: Parsing Specific Host, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        specific_host_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        specific_host = feature.qualifiers.get('specific_host')
                        specific_host_list_with_none.append(specific_host)
        specific_host_list_with_none_list = [['None'] if specific_host_value is None else specific_host_value for specific_host_value in specific_host_list_with_none]
        specific_host_list = []
        for specific_host_none in specific_host_list_with_none_list:
            specific_host_list.append('\n'.join(specific_host_none))
        print(specific_host_list)
        print("")
        # Writing Specific Host
        print("Step 2: Writing Specific Host, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        specific_host_string = "\n".join(specific_host_list)
        print(specific_host_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_specific_host.txt"), "w")
        handle_outputfile.write(specific_host_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 66: Specific Host Complete")
        print("")
        break

# Specimen Voucher Extract ##################################################
def specimenvoucher():
    while True:
        # Module Name
        print("Module 67: Specimen Voucher Extract")
        print("")
        # Parsing Specimen Voucher
        print("Step 1: Parsing Specimen Voucher, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        specimen_voucher_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        specimen_voucher = feature.qualifiers.get('specimen_voucher')
                        specimen_voucher_list_with_none.append(specimen_voucher)
        specimen_voucher_list_with_none_list = [['None'] if specimen_voucher_value is None else specimen_voucher_value for specimen_voucher_value in specimen_voucher_list_with_none]
        specimen_voucher_list = []
        for specimen_voucher_none in specimen_voucher_list_with_none_list:
            specimen_voucher_list.append('\n'.join(specimen_voucher_none))
        print(specimen_voucher_list)
        print("")
        # Writing Specimen Voucher
        print("Step 2: Writing Specimen Voucher, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        specimen_voucher_string = "\n".join(specimen_voucher_list)
        print(specimen_voucher_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_specimen_voucher.txt"), "w")
        handle_outputfile.write(specimen_voucher_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 67: Specimen Voucher Extract Complete")
        print("")
        break

# Strain Extract ##################################################
def strain():
    while True:
        # Module Name
        print("Module 68: Strain Extract")
        print("")
        # Parsing Strain
        print("Step 1: Parsing Strain, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        strain_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        strain = feature.qualifiers.get('strain')
                        strain_list_with_none.append(strain)
        strain_list_with_none_list = [['None'] if strain_value is None else strain_value for strain_value in strain_list_with_none]
        strain_list = []
        for strain_none in strain_list_with_none_list:
            strain_list.append('\n'.join(strain_none))
        print(strain_list)
        print("")
        # Writing strain
        print("Step 2: Writing Strain, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        strain_string = "\n".join(strain_list)
        print(strain_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_strain.txt"), "w")
        handle_outputfile.write(strain_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 68: strain Extract Complete")
        print("")
        break

# Strain ACCN Extract ##################################################
def strainaccn():
    while True:
        # Module Name
        print("Module 68.1: Strain ACCN Extract")
        print("")
        # Parsing Strain ACCN
        print("Step 1: Parsing Strain ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        strain_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        strain = feature.qualifiers.get('strain')
                        accession_list_with_none.append(record.id)
                        strain_list_with_none.append(strain)
        strain_list_with_none_list = [['None'] if strain_value is None else strain_value for strain_value in strain_list_with_none]
        #strain_list = []
        strain_list_with_none_list_none = []
        for strain_none in strain_list_with_none_list:
            strain_list_with_none_list_none.append('\n'.join(strain_none))
        strain_zip = zip(accession_list_with_none,strain_list_with_none_list_none)
        print(strain_zip)
        strain_list = []
        for accession_strain in strain_zip:
            strain_list.append("\t".join(accession_strain))
        #print(strain_list)
        print("")
        # Writing strain ACCN
        print("Step 2: Writing Strain ACCN, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        #strain_string = "\n".join(strain_list)
        strain_string = "\n".join(str(accesion_strain_value) for accesion_strain_value in strain_list)
        print(strain_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_strain_accn.txt"), "w")
        handle_outputfile.write(str(strain_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 68.1: strain ACCN Extract Complete")
        print("")
        break

# Sub Species Extract ##################################################
def subspecies():
    while True:
        # Module Name
        print("Module 69: Sub Species Extract")
        print("")
        # Parsing Sub Species
        print("Step 1: Parsing Sub Species, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        sub_species_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        sub_species = feature.qualifiers.get('sub_species')
                        sub_species_list_with_none.append(sub_species)
        sub_species_list_with_none_list = [['None'] if sub_species_value is None else sub_species_value for sub_species_value in sub_species_list_with_none]
        sub_species_list = []
        for sub_species_none in sub_species_list_with_none_list:
            sub_species_list.append('\n'.join(sub_species_none))
        print(sub_species_list)
        print("")
        # Writing Sub Species
        print("Step 2: Writing Sub Species, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        sub_species_string = "\n".join(sub_species_list)
        print(sub_species_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_sub_species.txt"), "w")
        handle_outputfile.write(sub_species_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 69: Sub Species Complete")
        print("")
        break

# Subclone Extract ##################################################
def subclone():
    while True:
        # Module Name
        print("Module 70: Subclone Extract")
        print("")
        # Parsing Subclone
        print("Step 1: Parsing Subclone, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        subclone_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        subclone = feature.qualifiers.get('subclone')
                        subclone_list_with_none.append(subclone)
        subclone_list_with_none_list = [['None'] if subclone_value is None else subclone_value for subclone_value in subclone_list_with_none]
        subclone_list = []
        for subclone_none in subclone_list_with_none_list:
            subclone_list.append('\n'.join(subclone_none))
        print(subclone_list)
        print("")
        # Writing Subclone
        print("Step 2: Writing Subclone, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        subclone_string = "\n".join(subclone_list)
        print(subclone_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_subclone.txt"), "w")
        handle_outputfile.write(subclone_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 70: Subclone Extract Complete")
        print("")
        break

# Subgroup Extract ##################################################
def subgroup():
    while True:
        # Module Name
        print("Module 71: Subgroup Extract")
        print("")
        # Parsing Subgroup
        print("Step 1: Parsing Subgroup, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        subgroup_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        subgroup = feature.qualifiers.get('subgroup')
                        subgroup_list_with_none.append(subgroup)
        subgroup_list_with_none_list = [['None'] if subgroup_value is None else subgroup_value for subgroup_value in subgroup_list_with_none]
        subgroup_list = []
        for subgroup_none in subgroup_list_with_none_list:
            subgroup_list.append('\n'.join(subgroup_none))
        print(subgroup_list)
        print("")
        # Writing Subgroup
        print("Step 2: Writing Subgroup, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        subgroup_string = "\n".join(subgroup_list)
        print(subgroup_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_subgroup.txt"), "w")
        handle_outputfile.write(subgroup_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 71: Subgroup Extract Complete")
        print("")
        break

# Substrain Extract ##################################################
def substrain():
    while True:
        # Module Name
        print("Module 72: Substrain Extract")
        print("")
        # Parsing Substrain
        print("Step 1: Parsing Substrain, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        substrain_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        substrain = feature.qualifiers.get('substrain')
                        substrain_list_with_none.append(substrain)
        substrain_list_with_none_list = [['None'] if substrain_value is None else substrain_value for substrain_value in substrain_list_with_none]
        substrain_list = []
        for substrain_none in substrain_list_with_none_list:
            substrain_list.append('\n'.join(substrain_none))
        print(substrain_list)
        print("")
        # Writing Substrain
        print("Step 2: Writing Substrain, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        substrain_string = "\n".join(substrain_list)
        print(substrain_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_substrain.txt"), "w")
        handle_outputfile.write(substrain_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 72: Substrain Extract Complete")
        print("")
        break

# Subtype Extract ##################################################
def subtype():
    while True:
        # Module Name
        print("Module 73: Subtype Extract")
        print("")
        # Parsing Subtype
        print("Step 1: Parsing Subtype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        subtype_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        subtype = feature.qualifiers.get('subtype')
                        subtype_list_with_none.append(subtype)
        subtype_list_with_none_list = [['None'] if subtype_value is None else subtype_value for subtype_value in subtype_list_with_none]
        subtype_list = []
        for subtype_none in subtype_list_with_none_list:
            subtype_list.append('\n'.join(subtype_none))
        print(subtype_list)
        print("")
        # Writing Subtype
        print("Step 2: Writing Subtype, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        subtype_string = "\n".join(subtype_list)
        print(subtype_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_subtype.txt"), "w")
        handle_outputfile.write(subtype_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 73: Subtype Extract Complete")
        print("")
        break

# Synonym Extract ##################################################
def synonym():
    while True:
        # Module Name
        print("Module 74: Synonym Extract")
        print("")
        # Parsing Synonym
        print("Step 1: Parsing Synonym, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        synonym_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        synonym = feature.qualifiers.get('synonym')
                        synonym_list_with_none.append(synonym)
        synonym_list_with_none_list = [['None'] if synonym_value is None else synonym_value for synonym_value in synonym_list_with_none]
        synonym_list = []
        for synonym_none in synonym_list_with_none_list:
            synonym_list.append('\n'.join(synonym_none))
        print(synonym_list)
        print("")
        # Writing Synonym
        print("Step 2: Writing Synonym, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        synonym_string = "\n".join(synonym_list)
        print(synonym_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_synonym.txt"), "w")
        handle_outputfile.write(synonym_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 74: Synonym Extract Complete")
        print("")
        break

# Taxon Extract ##################################################
def taxon():
    while True:
        # Module Name
        print("Module 75: Taxon Extract")
        print("")
        # Parsing Taxon
        print("Step 1: Parsing Taxon, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        taxon_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        #teleomorph = feature.qualifiers.get('teleomorph')
                        taxon = feature.qualifiers["db_xref"][0].strip('taxon:')
                        taxon_list_with_none.append(taxon)
        taxon_list_with_none_list = [['None'] if taxon_value is None else taxon_value for taxon_value in taxon_list_with_none]
        taxon_list = []
        for taxon_none in taxon_list_with_none_list:
            taxon_list.append('\n'.join(taxon_none))
        print(taxon_list)
        print("")
        # Writing Taxon
        print("Step 2: Writing Taxon, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        taxon_string = "\n".join(taxon_list)
        print(taxon_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_taxon.txt"), "w")
        handle_outputfile.write(taxon_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 75: Taxon Extract Complete")
        print("")
        break

# Teleomorph Extract ##################################################
def teleomorph():
    while True:
        # Module Name
        print("Module 76: Teleomorph Extract")
        print("")
        # Parsing Teleomorph
        print("Step 1: Parsing Teleomorph, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        teleomorph_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        teleomorph = feature.qualifiers.get('teleomorph')
                        teleomorph_list_with_none.append(teleomorph)
        teleomorph_list_with_none_list = [['None'] if teleomorph_value is None else teleomorph_value for teleomorph_value in teleomorph_list_with_none]
        teleomorph_list = []
        for teleomorph_none in teleomorph_list_with_none_list:
            teleomorph_list.append('\n'.join(teleomorph_none))
        print(teleomorph_list)
        print("")
        # Writing Teleomorph
        print("Step 2: Writing Teleomorph, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        teleomorph_string = "\n".join(teleomorph_list)
        print(teleomorph_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_teleomorph.txt"), "w")
        handle_outputfile.write(teleomorph_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 76: Teleomorph Extract Complete")
        print("")
        break

# Tissue Lib Extract ##################################################
def tissuelib():
    while True:
        # Module Name
        print("Module 77: Tissue Lib Extract")
        print("")
        # Parsing Tissue Lib
        print("Step 1: Parsing Tissue Lib, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        tissue_lib_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        tissue_lib = feature.qualifiers.get('tissue_lib')
                        tissue_lib_list_with_none.append(tissue_lib)
        tissue_lib_list_with_none_list = [['None'] if tissue_lib_value is None else tissue_lib_value for tissue_lib_value in tissue_lib_list_with_none]
        tissue_lib_list = []
        for tissue_lib_none in tissue_lib_list_with_none_list:
            tissue_lib_list.append('\n'.join(tissue_lib_none))
        print(tissue_lib_list)
        print("")
        # Writing Lab Host
        print("Step 2: Writing Tissue Lib, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        tissue_lib_string = "\n".join(tissue_lib_list)
        print(tissue_lib_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_tissue_lib.txt"), "w")
        handle_outputfile.write(tissue_lib_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 77: Tissue Lib Extract Complete")
        print("")
        break

# Tissue Type Extract ##################################################
def tissuetype():
    while True:
        # Module Name
        print("Module 78: Tissue Type Extract")
        print("")
        # Parsing Tissue Type
        print("Step 1: Parsing Tissue Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        tissue_type_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        tissue_type = feature.qualifiers.get('tissue_type')
                        tissue_type_list_with_none.append(tissue_type)
        tissue_type_list_with_none_list = [['None'] if tissue_type_value is None else tissue_type_value for tissue_type_value in tissue_type_list_with_none]
        tissue_type_list = []
        for tissue_type_none in tissue_type_list_with_none_list:
            tissue_type_list.append('\n'.join(tissue_type_none))
        print(tissue_type_list)
        print("")
        # Writing Lab Host
        print("Step 2: Writing Tissue Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        tissue_type_string = "\n".join(tissue_type_list)
        print(tissue_type_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_tissue_type.txt"), "w")
        handle_outputfile.write(tissue_type_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 78: Tissue Type Extract Complete")
        print("")
        break

# Transgenic Extract ##################################################
def transgenic():
    while True:
        # Module Name
        print("Module 79: Transgenic Extract")
        print("")
        # Parsing Transgenic
        print("Step 1: Parsing Transgenic, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        transgenic_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        transgenic = feature.qualifiers.get('transgenic')
                        transgenic_list_with_none.append(transgenic)
        transgenic_list_with_none_list = [['None'] if transgenic_value is None else transgenic_value for transgenic_value in transgenic_list_with_none]
        transgenic_list = []
        for transgenic_none in transgenic_list_with_none_list:
            transgenic_list.append('\n'.join(transgenic_none))
        print(transgenic_list)
        print("")
        # Writing Transgenic
        print("Step 2: Writing Transgenic, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        transgenic_string = "\n".join(transgenic_list)
        print(transgenic_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_transgenic.txt"), "w")
        handle_outputfile.write(transgenic_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 79: Transgenic Extract Complete")
        print("")
        break

# Type Extract ##################################################
def sourcetype():
    while True:
        # Module Name
        print("Module 80: Type Extract")
        print("")
        # Parsing Type
        print("Step 1: Parsing Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        type_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        source_type = feature.qualifiers.get('type')
                        type_list_with_none.append(source_type)
        type_list_with_none_list = [['None'] if type_value is None else type_value for type_value in type_list_with_none]
        type_list = []
        for type_none in type_list_with_none_list:
            type_list.append('\n'.join(type_none))
        print(type_list)
        print("")
        # Writing Type
        print("Step 2: Writing Type, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        type_string = "\n".join(type_list)
        print(type_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_type.txt"), "w")
        handle_outputfile.write(type_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 80: Type Extract Complete")
        print("")
        break

# Variety Extract ##################################################
def variety():
    while True:
        # Module Name
        print("Module 81: Variety Extract")
        print("")
        # Parsing Variety
        print("Step 1: Parsing Variety, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        variety_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "source":
                        variety = feature.qualifiers.get('variety')
                        variety_list_with_none.append(variety)
        variety_list_with_none_list = [['None'] if variety_value is None else variety_value for variety_value in variety_list_with_none]
        variety_list = []
        for variety_none in variety_list_with_none_list:
            variety_list.append('\n'.join(variety_none))
        print(variety_list)
        print("")
        # Writing Variety
        print("Step 2: Writing Variety, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        variety_string = "\n".join(variety_list)
        print(variety_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/source_variety.txt"), "w")
        handle_outputfile.write(variety_string)
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 81: Variety Extract Complete")
        print("")
        break

# Module Main ##################################################
def mainsource():
    while True:
        print("")
        print("\tMenu: Extract Source")
        print("")
        print("\t\tModule 1.  Acronym \t\t\t: Button[1]")
        print("\t\tModule 2.  Altitude \t\t\t: Button[2]")
        print("\t\tModule 3.  Anamorph \t\t\t: Button[3]")
        print("\t\tModule 4.  Authority \t\t\t: Button[4]")
        print("\t\tModule 5.  Bio Material \t\t: Button[5]")
        print("\t\tModule 6.  Biotype \t\t\t: Button[6]")
        print("\t\tModule 7.  Biovar \t\t\t: Button[7]")
        print("\t\tModule 8.  Breed \t\t\t: Button[8]")
        print("\t\tModule 9.  Cell Line \t\t\t: Button[9]")
        print("\t\tModule 10. Cell Type \t\t\t: Button[10]")
        print("\t\tModule 11. Chemovar \t\t\t: Button[11]")
        print("\t\tModule 12. Chromosome \t\t\t: Button[12]")        
        print("\t\tModule 13. Clone \t\t\t: Button[13]")
        print("\t\tModule 14. Clone Lib \t\t\t: Button[14]")
        print("\t\tModule 15. Collected By \t\t: Button[15]")
        print("\t\tModule 16. Collection Date \t\t: Button[16]")
        print("\t\tModule 17. Common \t\t\t: Button[17]")
        print("\t\tModule 18. Country \t\t\t: Button[18]")
        print("\t\tModule 18.1 Country ACCN \t\t: Button[181]")
        print("\t\tModule 19. Cultivar \t\t\t: Button[19]")
        print("\t\tModule 20. Culture Collection \t\t: Button[20]")
        print("\t\tModule 21. Dev Stage \t\t\t: Button[21]")
        print("\t\tModule 22. Ecotype \t\t\t: Button[22]")
        print("\t\tModule 23. Endogenous Virus Name \t: Button[23]")
        print("\t\tModule 24. Environmental Sample \t: Button[24]")
        print("\t\tModule 25. Focus \t\t\t: Button[25]")
        print("\t\tModule 26. Forma \t\t\t: Button[26]")
        print("\t\tModule 27. Forma Specialis \t\t: Button[27]")
        print("\t\tModule 28. Fwd PCR Primer Name \t\t: Button[28]")
        print("\t\tModule 29. Fwd PCR Primer Seq \t\t: Button[29]")
        print("\t\tModule 30. Fwd Primer Name \t\t: Button[30]")
        print("\t\tModule 31. Fwd Primer Seq \t\t: Button[32]")
        print("\t\tModule 32. Gene \t\t\t: Button[32]")
        print("\t\tModule 33. Genotype \t\t\t: Button[33]")       
        print("\t\tModule 34. Germline \t\t\t: Button[34]")
        print("\t\tModule 35. Group \t\t\t: Button[35]")
        print("\t\tModule 36. Haplogroup \t\t\t: Button[36]")
        print("\t\tModule 37. Haplotype \t\t\t: Button[37]")
        print("\t\tModule 38. Host \t\t\t: Button[38]")
        print("\t\tModule 38.1 Host ACCN \t\t\t: Button[381]")
        print("\t\tModule 39. Identified By \t\t: Button[39]")
        print("\t\tModule 40. Isolate \t\t\t: Button[40]")
        print("\t\tModule 40.1 Isolate ACCN \t\t: Button[401]")
        print("\t\tModule 41. Isolation Source \t\t: Button[41]")
        print("\t\tModule 41.1 Isolation Source ACCN \t: Button[411]")        
        print("\t\tModule 42. Lab Host \t\t\t: Button[42]")
        print("\t\tModule 43. Lat Lon \t\t\t: Button[43]")
        print("\t\tModule 44. Latitude\Longitude (GeoPy) \t: Button[44]")
        print("\t\tModule 45. Map \t\t\t\t: Button[45]")
        print("\t\tModule 46. Metagenome Source \t\t: Button[46]")
        print("\t\tModule 47. Metagenomic \t\t\t: Button[47]")
        print("\t\tModule 48. Note \t\t\t: Button[48]")
        print("\t\tModule 49. Organism \t\t\t: Button[49]")
        print("\t\tModule 50. Pathovar \t\t\t: Button[50]")
        print("\t\tModule 51. Plasmid Name \t\t: Button[51]")
        print("\t\tModule 52. Plastid Name \t\t: Button[52]")
        print("\t\tModule 53. Pop Variant \t\t\t: Button[53]")
        print("\t\tModule 54. Protein \t\t\t: Button[54]")
        print("\t\tModule 55. Prot Desc \t\t\t: Button[55]")        
        print("\t\tModule 56. Rearranged \t\t\t: Button[56]")
        print("\t\tModule 57. Rev PCR Primer Name \t\t: Button[57]")
        print("\t\tModule 58. Rev PCR Primer Seq \t\t: Button[58]")        
        print("\t\tModule 59. Rev Primer Name \t\t: Button[59]")
        print("\t\tModule 60. Rev Primer Seq \t\t: Button[60]")
        print("\t\tModule 61. Segment \t\t\t: Button[61]")
        print("\t\tModule 62. Serogroup \t\t\t: Button[62]")
        print("\t\tModule 63. Serotype \t\t\t: Button[63]")
        print("\t\tModule 63.1 Serotype ACCN \t\t: Button[631]")        
        print("\t\tModule 64. Serovar \t\t\t: Button[64]")
        print("\t\tModule 65. Sex \t\t\t\t: Button[65]")
        print("\t\tModule 66. Specific Host \t\t: Button[66]")
        print("\t\tModule 67. Specimen Voucher \t\t: Button[67]")
        print("\t\tModule 68. Strain \t\t\t: Button[68]")
        print("\t\tModule 68.1 Strain ACCN \t\t: Button[681]")        
        print("\t\tModule 69. Sub Species \t\t\t: Button[69]")
        print("\t\tModule 70. Subclone \t\t\t: Button[70]")
        print("\t\tModule 71. Subgroup \t\t\t: Button[71]")
        print("\t\tModule 72. Substrain \t\t\t: Button[72]")
        print("\t\tModule 73. Subtype \t\t\t: Button[73]")
        print("\t\tModule 74. Synonym \t\t\t: Button[74]")
        print("\t\tModule 75. Taxon \t\t\t: Button[75]")
        print("\t\tModule 76. Teleomorph \t\t\t: Button[76]")
        print("\t\tModule 77. Tissue Lib \t\t\t: Button[77]")
        print("\t\tModule 78. Tissue Type \t\t\t: Button[78]")
        print("\t\tModule 79. Transgenic \t\t\t: Button[79]")        
        print("\t\tModule 80. Type \t\t\t: Button[80]")
        print("\t\tModule 81. Variety \t\t\t: Button[81]")
        print("\t\tModule 82. Main Menu \t\t\t: Button[82]")
        print("\t\tModule 83. Exit. \t\t\t: Button[83]")
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
                acronym()
            elif button == 2:
                print("")
                altitude()
            elif button == 3:
                print("")
                anamorph()
            elif button == 4:
                print("")
                authority()
            elif button == 5:
                print("")
                biomaterial()
            elif button == 6:
                print("")
                biotype()
            elif button == 7:
                print("")
                biovar()
            elif button == 8:
                print("")
                breed()               
            elif button == 9:
                print("")
                cellline()
            elif button == 10:
                print("")
                celltype()               
            elif button == 11:
                print("")
                chemovar()
            elif button == 12:
                print("")
                chromosome()
            elif button == 13:
                print("")
                clone()
            elif button == 14:
                print("")
                clonelib()
            elif button == 15:
                print("")
                collectedby()               
            elif button == 16:
                print("")
                collectiondate()
            elif button == 17:
                print("")
                common()
            elif button == 18:
                print("")
                country()
            elif button == 181:
                print("")
                countryaccn()
            elif button == 19:
                print("")
                cultivar()
            elif button == 20:
                print("")
                culturecollection()
            elif button == 21:
                print("")
                devstage()
            elif button == 22:
                print("")
                ecotype()
            elif button == 23:
                print("")
                endogenousvirusname()
            elif button == 24:
                print("")
                environmentalsample()
            elif button == 25:
                print("")
                focus()
            elif button == 26:
                print("")
                forma()
            elif button == 27:
                print("")
                formaspecialis()
            elif button == 28:
                print("")
                fwdpcrprimername()
            elif button == 29:
                print("")
                fwdpcrprimerseq()
            elif button == 30:
                print("")
                fwdprimername()
            elif button == 31:
                print("")
                fwdprimerseq()               
            elif button == 32:
                print("")
                gene()
            elif button == 33:
                print("")
                genotype()               
            elif button == 34:
                print("")
                germline()
            elif button == 35:
                print("")
                group()
            elif button == 36:
                print("")
                haplogroup()
            elif button == 37:
                print("")
                haplotype()
            elif button == 38:
                print("")
                host()
            elif button == 381:
                print("")
                hostaccn()
            elif button == 39:
                print("")
                identifiedby()
            elif button == 40:
                print("")
                isolate()
            elif button == 401:
                print("")
                isolateaccn()
            elif button == 41:
                print("")
                isolationsource()
            elif button == 411:
                print("")
                isolationsourceaccn()
            elif button == 42:
                print("")
                labhost()
            elif button == 43:
                print("")
                latlon()
            elif button == 44:
                print("")
                latitudelongitude()
            elif button == 45:
                print("")
                mapsource()
            elif button == 46:
                print("")
                metagenomesource()
            elif button == 47:
                print("")
                metagenomic()               
            elif button == 48:
                print("")
                note()
            elif button == 49:
                print("")
                organism()
            elif button == 50:
                print("")
                pathovar()
            elif button == 51:
                print("")
                plasmidname()               
            elif button == 52:
                print("")
                plastidname()
            elif button == 53:
                print("")
                popvariant()
            elif button == 54:
                print("")
                protein()
            elif button == 55:
                print("")
                protdesc()
            elif button == 56:
                print("")
                rearranged()               
            elif button == 57:
                print("")
                revpcrprimername()
            elif button == 58:
                print("")
                revpcrprimerseq()
            elif button == 59:
                print("")
                revprimername()
            elif button == 60:
                print("")
                revprimerseq()
            elif button == 61:
                print("")
                segment()
            elif button == 62:
                print("")
                serogroup()
            elif button == 63:
                print("")
                serotype()
            elif button == 631:
                print("")
                serotypeaccn()
            elif button == 64:
                print("")
                serovar()
            elif button == 65:
                print("")
                sex()
            elif button == 66:
                print("")
                specifichost()
            elif button == 67:
                print("")
                specimenvoucher()
            elif button == 68:
                print("")
                strain()
            elif button == 681:
                print("")
                strainaccn()
            elif button == 69:
                print("")
                subspecies()
            elif button == 70:
                print("")
                subclone()
            elif button == 71:
                print("")
                subgroup()
            elif button == 72:
                print("")
                substrain()
            elif button == 73:
                print("")
                subtype()
            elif button == 74:
                print("")
                synonym()
            elif button == 75:
                print("")
                taxon()
            elif button == 76:
                print("")
                teleomorph()
            elif button == 77:
                print("")
                tissuelib()
            elif button == 78:
                print("")
                tissuetype()
            elif button == 79:
                print("")
                transgenic()
            elif button == 80:
                print("")
                sourcetype()
            elif button == 81:
                print("")
                variety()
            elif button == 82:
                print("")
                import main
                main.main()
            elif button == 83:
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

mainsource()
