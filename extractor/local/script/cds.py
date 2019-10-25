#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from geopy.geocoders import Nominatim
from Bio.Alphabet import generic_dna, generic_protein

# Module: CDS Features Extract
button = ""

# Codon Start Extract ##################################################
def codonstart():
    while True:
        # Module Name
        print("Module 1: Codon Start Extract")
        print("")
        # Parsing Codon Start
        print("Step 1: Parsing Codon Start, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        codon_start_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        codon_start = feature.qualifiers.get('codon_start')
                        accession_list_with_none.append(record.id)
                        codon_start_list_with_none.append(codon_start)
        codon_start_list_with_none_list = [['None'] if codon_start_value is None else codon_start_value for codon_start_value in codon_start_list_with_none]
        codon_start_list_with_none_list_none = []
        for codon_start_none in codon_start_list_with_none_list:
            codon_start_list_with_none_list_none.append('\n'.join(codon_start_none))
        codon_start_zip = zip(accession_list_with_none,codon_start_list_with_none_list_none)
        print(codon_start_zip)
        codon_start_list = []
        for accession_codon_start in codon_start_zip:
            codon_start_list.append("\t".join(accession_codon_start))
        print("")
        # Writing Codon Start
        print("Step 2: Writing Codon Start, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        codon_start_string = "\n".join(str(accesion_codon_start_value) for accesion_codon_start_value in codon_start_list)
        print(codon_start_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_codon_start.txt"), "w")
        handle_outputfile.write(str(codon_start_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 1 : Codon Start Extract Complete")
        print("")
        break

# DB XREF Extract ##################################################
def dbxref():
    while True:
        # Module Name
        print("Module 2: DB XREF Extract")
        print("")
        # Parsing DB XREF
        print("Step 1: Parsing DB XREF, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        db_xref_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        db_xref = feature.qualifiers.get('db_xref')
                        accession_list_with_none.append(record.id)
                        db_xref_list_with_none.append(db_xref)
        db_xref_list_with_none_list = [['None'] if db_xref_value is None else db_xref_value for db_xref_value in db_xref_list_with_none]
        db_xref_list_with_none_list_none = []
        for db_xref_none in db_xref_list_with_none_list:
            db_xref_list_with_none_list_none.append('\n'.join(db_xref_none))
        db_xref_zip = zip(accession_list_with_none,db_xref_list_with_none_list_none)
        print(db_xref_zip)
        db_xref_list = []
        for accession_db_xref in db_xref_zip:
            db_xref_list.append("\t".join(accession_db_xref))
        print("")
        # Writing DB XREF
        print("Step 2: Writing DB XREF, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        db_xref_string = "\n".join(str(accesion_db_xref_value) for accesion_db_xref_value in db_xref_list)
        print(db_xref_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_db_xref.txt"), "w")
        handle_outputfile.write(str(db_xref_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 2: DB XREF Extract Complete")
        print("")
        break

# EC Number Extract ##################################################
def ecnumber():
    while True:
        # Module Name
        print("Module 3: EC Number Extract")
        print("")
        # Parsing EC Number
        print("Step 1: Parsing EC Number, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        ec_number_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        ec_number = feature.qualifiers.get('ec_number')
                        accession_list_with_none.append(record.id)
                        ec_number_list_with_none.append(ec_number)
        ec_number_list_with_none_list = [['None'] if ec_number_value is None else ec_number_value for ec_number_value in ec_number_list_with_none]
        ec_number_list_with_none_list_none = []
        for ec_number_none in ec_number_list_with_none_list:
            ec_number_list_with_none_list_none.append('\n'.join(ec_number_none))
        ec_number_zip = zip(accession_list_with_none,ec_number_list_with_none_list_none)
        print(ec_number_zip)
        ec_number_list = []
        for accession_ec_number in ec_number_zip:
            ec_number_list.append("\t".join(accession_ec_number))
        print("")
        # Writing EC Number
        print("Step 2: Writing EC Number, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        ec_number_string = "\n".join(str(accession_ec_number_value) for accession_ec_number_value in ec_number_list)
        print(ec_number_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_ec_number.txt"), "w")
        handle_outputfile.write(str(ec_number_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 3: EC Number Extract Complete")
        print("")
        break

# Evidence Extract ##################################################
def evidence():
    while True:
        # Module Name
        print("Module 4: Evidence Extract")
        print("")
        # Parsing Evidence
        print("Step 1: Parsing Evidence, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        evidence_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        evidence = feature.qualifiers.get('evidence')
                        accession_list_with_none.append(record.id)
                        evidence_list_with_none.append(evidence)
        evidence_list_with_none_list = [['None'] if evidence_value is None else evidence_value for evidence_value in evidence_list_with_none]
        evidence_list_with_none_list_none = []
        for evidence_none in evidence_list_with_none_list:
            evidence_list_with_none_list_none.append('\n'.join(evidence_none))
        evidence_zip = zip(accession_list_with_none,evidence_list_with_none_list_none)
        print(evidence_zip)
        evidence_list = []
        for accession_evidence in evidence_zip:
            evidence_list.append("\t".join(accession_evidence))
        print("")
        # Writing Evidence
        print("Step 2: Writing Evidence, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        evidence_string = "\n".join(str(accession_evidence_value) for accession_evidence_value in evidence_list)
        print(evidence_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_evidence.txt"), "w")
        handle_outputfile.write(str(evidence_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 4: Evidence Extract Complete")
        print("")
        break

# Exception Extract ##################################################
def exception():
    while True:
        # Module Name
        print("Module 5: Exception Extract")
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
                    if feature.type == "CDS":
                        exception = feature.qualifiers.get('exception')
                        accession_list_with_none.append(record.id)
                        exception_list_with_none.append(exception)
        exception_list_with_none_list = [['None'] if exception_value is None else exception_value for exception_value in exception_list_with_none]
        exception_list_with_none_list_none = []
        for exception_none in exception_list_with_none_list:
            exception_list_with_none_list_none.append('\n'.join(exception_none))
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
        exception_string = "\n".join(str(accession_exception_value) for accession_exception_value in exception_list)
        print(exception_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_exception.txt"), "w")
        handle_outputfile.write(str(exception_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 5: Exception Extract Complete")
        print("")
        break

# Experiment Extract ##################################################
def experiment():
    while True:
        # Module Name
        print("Module 6: Experiment Extract")
        print("")
        # Parsing Experiment
        print("Step 1: Parsing Experiment, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        experiment_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        experiment = feature.qualifiers.get('experiment')
                        accession_list_with_none.append(record.id)
                        experiment_list_with_none.append(experiment)
        experiment_list_with_none_list = [['None'] if experiment_value is None else experiment_value for experiment_value in experiment_list_with_none]
        experiment_list_with_none_list_none = []
        for experiment_none in experiment_list_with_none_list:
            experiment_list_with_none_list_none.append('\n'.join(experiment_none))
        experiment_zip = zip(accession_list_with_none,experiment_list_with_none_list_none)
        print(experiment_zip)
        experiment_list = []
        for accession_experiment in experiment_zip:
            experiment_list.append("\t".join(accession_experiment))
        print("")
        # Writing Experiment
        print("Step 2: Writing Experiment, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        experiment_string = "\n".join(str(accession_experiment_value) for accession_experiment_value in experiment_list)
        print(experiment_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_experiment.txt"), "w")
        handle_outputfile.write(str(experiment_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 6: Experiment Extract Complete")
        print("")
        break

# Function Extract ##################################################
def function():
    while True:
        # Module Name
        print("Module 7: Function Extract")
        print("")
        # Parsing Function
        print("Step 1: Parsing Function, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        function_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        function = feature.qualifiers.get('function')
                        accession_list_with_none.append(record.id)
                        function_list_with_none.append(function)
        function_list_with_none_list = [['None'] if function_value is None else function_value for function_value in function_list_with_none]
        function_list_with_none_list_none = []
        for function_none in function_list_with_none_list:
            function_list_with_none_list_none.append('\n'.join(function_none))
        function_zip = zip(accession_list_with_none,function_list_with_none_list_none)
        print(function_zip)
        function_list = []
        for accession_function in function_zip:
            function_list.append("\t".join(accession_function))
        print("")
        # Writing Function
        print("Step 2: Writing Function, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        function_string = "\n".join(str(accession_function_value) for accession_function_value in function_list)
        print(function_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_function.txt"), "w")
        handle_outputfile.write(str(function_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 7: Function Extract Complete")
        print("")
        break

# Gene Extract ##################################################
def gene():
    while True:
        # Module Name
        print("Module 8: Gene Extract")
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
                    if feature.type == "CDS":
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
        gene_string = '\n'.join(str(accession_gene_value) for accession_gene_value in gene_list)
        print(gene_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_gene.txt"), "w")
        handle_outputfile.write(str(gene_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 8: Gene Extract Complete")
        print("")
        break

# Gene Synonym Extract ##################################################
def genesynonym():
    while True:
        # Module Name
        print("Module 9: Gene Synonym Extract")
        print("")
        # Parsing Gene Synonym
        print("Step 1: Parsing Gene Synonym, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        gene_synonym_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        gene_synonym = feature.qualifiers.get('gene_synonym')
                        accession_list_with_none.append(record.id)
                        gene_synonym_list_with_none.append(gene_synonym)
        gene_synonym_list_with_none_list = [['None'] if gene_synonym_value is None else gene_synonym_value for gene_synonym_value in gene_synonym_list_with_none]
        gene_synonym_list_with_none_list_none = []
        for gene_synonym_none in gene_synonym_list_with_none_list:
            gene_synonym_list_with_none_list_none.append('\n'.join(gene_synonym_none))
        # Zipping Gene Synonym
        gene_synonym_zip = zip(accession_list_with_none,gene_synonym_list_with_none_list_none)
        print(gene_synonym_zip)
        gene_synonym_list = []
        for accession_gene_synonym in gene_synonym_zip:
            gene_synonym_list.append("\t".join(accession_gene_synonym))
        print("")
        # Writing Gene Synonym
        print("Step 2: Writing Gene Synonym, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        gene_synonym_string = '\n'.join(str(accession_gene_synonym_value) for accession_gene_synonym_value in gene_synonym_list)
        print(gene_synonym_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_gene_synonym.txt"), "w")
        handle_outputfile.write(str(gene_synonym_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 9: Gene Synonym Extract Complete")
        print("")
        break

# Go Component Extract ##################################################
def gocomponent():
    while True:
        # Module Name
        print("Module 10: Go Component Extract")
        print("")
        # Parsing Go Component
        print("Step 1: Parsing Go Component, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        go_component_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        go_component = feature.qualifiers.get('go_component')
                        accession_list_with_none.append(record.id)
                        go_component_list_with_none.append(go_component)
        go_component_list_with_none_list = [['None'] if go_component_value is None else go_component_value for go_component_value in go_component_list_with_none]
        go_component_list_with_none_list_none = []
        for go_component_none in go_component_list_with_none_list:
            go_component_list_with_none_list_none.append('\n'.join(go_component_none))
        # Zipping Go Component
        go_component_zip = zip(accession_list_with_none,go_component_list_with_none_list_none)
        print(go_component_zip)
        go_component_list = []
        for accession_go_component in go_component_zip:
            go_component_list.append("\t".join(accession_go_component))
        print("")
        # Writing Go Component
        print("Step 2: Writing Go Component, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        go_component_string = '\n'.join(str(accession_go_component_value) for accession_go_component_value in go_component_list)
        print(go_component_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_go_component.txt"), "w")
        handle_outputfile.write(str(go_component_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 10: Go Component Extract Complete")
        print("")
        break

# Go Function Extract ##################################################
def gofunction():
    while True:
        # Module Name
        print("Module 11: Go Function Extract")
        print("")
        # Parsing Go Function
        print("Step 1: Parsing Go Function, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        go_function_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        go_function = feature.qualifiers.get('go_function')
                        accession_list_with_none.append(record.id)
                        go_function_list_with_none.append(go_function)
        go_function_list_with_none_list = [['None'] if go_function_value is None else go_function_value for go_function_value in go_function_list_with_none]
        go_function_list_with_none_list_none = []
        for go_function_none in go_function_list_with_none_list:
            go_function_list_with_none_list_none.append('\n'.join(go_function_none))
        # Zipping Go Function
        go_function_zip = zip(accession_list_with_none,go_function_list_with_none_list_none)
        print(go_function_zip)
        go_function_list = []
        for accession_go_function in go_function_zip:
            go_function_list.append("\t".join(accession_go_function))
        print("")
        # Writing Go Function
        print("Step 2: Writing Go Component, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        go_function_string = '\n'.join(str(accession_go_function_value) for accession_go_function_value in go_function_list)
        print(go_function_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_go_function.txt"), "w")
        handle_outputfile.write(str(go_function_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 11: Go Function Extract Complete")
        print("")
        break

# Go Process Extract ##################################################
def goprocess():
    while True:
        # Module Name
        print("Module 12: Go Process Extract")
        print("")
        # Parsing Go Process
        print("Step 1: Parsing Go Process, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        go_process_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        go_process = feature.qualifiers.get('go_process')
                        accession_list_with_none.append(record.id)
                        go_process_list_with_none.append(go_process)
        go_process_list_with_none_list = [['None'] if go_process_value is None else go_process_value for go_process_value in go_process_list_with_none]
        go_process_list_with_none_list_none = []
        for go_process_none in go_process_list_with_none_list:
            go_process_list_with_none_list_none.append('\n'.join(go_process_none))
        # Zipping Go Process
        go_process_zip = zip(accession_list_with_none,go_process_list_with_none_list_none)
        print(go_process_zip)
        go_process_list = []
        for accession_go_process in go_process_zip:
            go_process_list.append("\t".join(accession_go_process))
        print("")
        # Writing Go Process
        print("Step 2: Writing Go Process, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        go_process_string = '\n'.join(str(accession_go_process_value) for accession_go_process_value in go_process_list)
        print(go_process_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_go_process.txt"), "w")
        handle_outputfile.write(str(go_process_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 12: Go Process Extract Complete")
        print("")
        break

# Inference Extract ##################################################
def inference():
    while True:
        # Module Name
        print("Module 13: Inference Extract")
        print("")
        # Parsing Inference
        print("Step 1: Parsing Inference, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        inference_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        inference = feature.qualifiers.get('inference')
                        accession_list_with_none.append(record.id)
                        inference_list_with_none.append(inference)
        inference_list_with_none_list = [['None'] if inference_value is None else inference_value for inference_value in inference_list_with_none]
        inference_list_with_none_list_none = []
        for inference_none in inference_list_with_none_list:
            inference_list_with_none_list_none.append('\n'.join(inference_none))
        # Zipping Inference
        inference_zip = zip(accession_list_with_none,inference_list_with_none_list_none)
        print(inference_zip)
        inference_list = []
        for accession_inference in inference_zip:
            inference_list.append("\t".join(accession_inference))
        print("")
        # Writing Inference
        print("Step 2: Writing Label, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        inference_string = '\n'.join(str(accession_inference_value) for accession_inference_value in inference_list)
        print(inference_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_inference.txt"), "w")
        handle_outputfile.write(str(inference_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 13: Inference Extract Complete")
        print("")
        break

# Label Extract ##################################################
def label():
    while True:
        # Module Name
        print("Module 14: Label Extract")
        print("")
        # Parsing Label
        print("Step 1: Parsing Label, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        label_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        label = feature.qualifiers.get('label')
                        accession_list_with_none.append(record.id)
                        label_list_with_none.append(label)
        label_list_with_none_list = [['None'] if label_value is None else label_value for label_value in label_list_with_none]
        label_list_with_none_list_none = []
        for label_none in label_list_with_none_list:
            label_list_with_none_list_none.append('\n'.join(label_none))
        # Zipping Label
        label_zip = zip(accession_list_with_none,label_list_with_none_list_none)
        print(label_zip)
        label_list = []
        for accession_label in label_zip:
            label_list.append("\t".join(accession_label))
        print("")
        # Writing Label
        print("Step 2: Writing Label, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        label_string = '\n'.join(str(accession_label_value) for accession_label_value in label_list)
        print(label_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_label.txt"), "w")
        handle_outputfile.write(str(label_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 14: Label Extract Complete")
        print("")
        break

# Location Extract ##################################################
def location():
    while True:
        # Module Name
        print("Module 15: Location Extract")
        print("")
        # Parsing Location
        print("Step 1: Parsing Location, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        location_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        #location = feature.qualifiers.get('location')
                        location = feature.location
                        accession_list_with_none.append(record.id)
                        location_list_with_none.append(feature.location)
        location_list_with_none_list = [['None'] if location_value is None else location_value for location_value in location_list_with_none]
        location_list_with_none_list_none = []
        for location_none in location_list_with_none_list:
            #location_list_with_none_list_none.append('\n'.join(location_none))
            location_list_with_none_list_none.append(location_none)
        # Zipping Location
        location_zip = zip(accession_list_with_none,location_list_with_none_list_none)
        print(location_zip)
        location_list = []
        for accession_location in location_zip:
            #location_list.append("\t".join(accession_location))
            location_list.append(accession_location)
        print("")
        # Writing Location
        print("Step 2: Writing Label, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        location_string = '\n'.join(str(accession_location_value) for accession_location_value in location_list)
        print(location_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_location.txt"), "w")
        handle_outputfile.write(str(location_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 15: Location Extract Complete")
        print("")
        break

# Locus Tag Extract ##################################################
def locustag():
    while True:
        # Module Name
        print("Module 16: Locus Tag Extract")
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
                    if feature.type == "CDS":
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
        locus_tag_string = '\n'.join(str(accession_locus_tag_value) for accession_locus_tag_value in locus_tag_list)
        print(locus_tag_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_locus_tag.txt"), "w")
        handle_outputfile.write(str(locus_tag_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 16: Locus Tag Extract Complete")
        print("")
        break

# Note Extract ##################################################
def note():
    while True:
        # Module Name
        print("Module 17: Note Extract")
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
                    if feature.type == "CDS":
                        note = feature.qualifiers.get('note')
                        accession_list_with_none.append(record.id)
                        note_list_with_none.append(note)
        note_list_with_none_list = [['None'] if note_value is None else note_value for note_value in note_list_with_none]
        note_list_with_none_list_none = []
        for note_none in note_list_with_none_list:
            note_list_with_none_list_none.append('\n'.join(note_none))
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
        note_string = "\n".join(str(accesion_note_value) for accesion_note_value in note_list)
        print(note_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_note.txt"), "w")
        handle_outputfile.write(str(note_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 17: Note Extract Complete")
        print("")
        break

# Product Extract ##################################################
def product():
    while True:
        # Module Name
        print("Module 18: Product Extract")
        print("")
        # Parsing Product
        print("Step 1: Parsing Product, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        product_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        product = feature.qualifiers.get('product')
                        accession_list_with_none.append(record.id)
                        product_list_with_none.append(product)
        product_list_with_none_list = [['None'] if product_value is None else product_value for product_value in product_list_with_none]
        product_list_with_none_list_none = []
        for product_none in product_list_with_none_list:
            product_list_with_none_list_none.append('\n'.join(product_none))
        product_zip = zip(accession_list_with_none,product_list_with_none_list_none)
        print(product_zip)
        product_list = []
        for accession_product in product_zip:
            product_list.append("\t".join(accession_product))
        print("")
        # Writing Product
        print("Step 2: Writing Product, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        product_string = "\n".join(str(accession_product_value) for accession_product_value in product_list)
        print(product_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_product.txt"), "w")
        handle_outputfile.write(str(product_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 18: Product Extract Complete")
        print("")
        break

# Prot DESC Extract ##################################################
def protdesc():
    while True:
        # Module Name
        print("Module 19: Prot DESC Extract")
        print("")
        # Parsing Prot DESC
        print("Step 1: Parsing Prot DESC, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        prot_desc_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        prot_desc = feature.qualifiers.get('prot_desc')
                        accession_list_with_none.append(record.id)
                        prot_desc_list_with_none.append(prot_desc)
        prot_desc_list_with_none_list = [['None'] if prot_desc_value is None else prot_desc_value for prot_desc_value in prot_desc_list_with_none]
        prot_desc_list_with_none_list_none = []
        for prot_desc_none in prot_desc_list_with_none_list:
            prot_desc_list_with_none_list_none.append('\n'.join(prot_desc_none))
        prot_desc_zip = zip(accession_list_with_none,prot_desc_list_with_none_list_none)
        print(prot_desc_zip)
        prot_desc_list = []
        for accession_prot_desc in prot_desc_zip:
            prot_desc_list.append("\t".join(accession_prot_desc))
        print("")
        # Writing Prot DESC
        print("Step 2: Writing Prot DESC, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        prot_desc_string = "\n".join(str(accesion_prot_desc_value) for accesion_prot_desc_value in prot_desc_list)
        print(prot_desc_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_prot_desc.txt"), "w")
        handle_outputfile.write(str(prot_desc_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 19: Prot DESC Extract Complete")
        print("")
        break

# Protein ID Extract ##################################################
def proteinid():
    while True:
        # Module Name
        print("Module 20: Protein ID Extract")
        print("")
        # Parsing Protein ID
        print("Step 1: Parsing Protein ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        protein_id_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        protein_id = feature.qualifiers.get('protein_id')
                        accession_list_with_none.append(record.id)
                        protein_id_list_with_none.append(protein_id)
        protein_id_list_with_none_list = [['None'] if protein_id_value is None else protein_id_value for protein_id_value in protein_id_list_with_none]
        protein_id_list_with_none_list_none = []
        for protein_id_none in protein_id_list_with_none_list:
            protein_id_list_with_none_list_none.append('\n'.join(protein_id_none))
        protein_id_zip = zip(accession_list_with_none,protein_id_list_with_none_list_none)
        print(protein_id_zip)
        protein_id_list = []
        for accession_protein_id in protein_id_zip:
            protein_id_list.append("\t".join(accession_protein_id))
        print("")
        # Writing Protein ID
        print("Step 2: Writing Protein ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        protein_id_string = "\n".join(str(accesion_protein_id_value) for accesion_protein_id_value in protein_id_list)
        print(protein_id_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_protein_id.txt"), "w")
        handle_outputfile.write(str(protein_id_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 20: Protein ID Extract Complete")
        print("")
        break

# Pseudo Extract ##################################################
def pseudo():
    while True:
        # Module Name
        print("Module 21: Pseudo Extract")
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
                    if feature.type == "CDS":
                        pseudo = feature.qualifiers.get('pseudo')
                        accession_list_with_none.append(record.id)
                        pseudo_list_with_none.append(pseudo)
        pseudo_list_with_none_list = [['None'] if pseudo_value is None else pseudo_value for pseudo_value in pseudo_list_with_none]
        pseudo_list_with_none_list_none = []
        for pseudo_none in pseudo_list_with_none_list:
            pseudo_list_with_none_list_none.append('\n'.join(pseudo_none))
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
        pseudo_string = "\n".join(str(accesion_pseudo_value) for accesion_pseudo_value in pseudo_list)
        print(pseudo_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_pseudo.txt"), "w")
        handle_outputfile.write(str(pseudo_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 21: Pseudo Extract Complete")
        print("")
        break

# Ribosomal Slippage Extract ##################################################
def ribosomalslippage():
    while True:
        # Module Name
        print("Module 22: Ribosomal Slippage Extract")
        print("")
        # Parsing Ribosomal Slippage
        print("Step 1: Parsing Ribosomal Slippage, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        ribosomal_slippage_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        ribosomal_slippage = feature.qualifiers.get('ribosomal_slippage')
                        accession_list_with_none.append(record.id)
                        ribosomal_slippage_list_with_none.append(ribosomal_slippage)
        ribosomal_slippage_list_with_none_list = [['None'] if ribosomal_slippage_value is None else ribosomal_slippage_value for ribosomal_slippage_value in ribosomal_slippage_list_with_none]
        ribosomal_slippage_list_with_none_list_none = []
        for ribosomal_slippage_none in ribosomal_slippage_list_with_none_list:
            ribosomal_slippage_list_with_none_list_none.append('\n'.join(ribosomal_slippage_none))
        # Zipping Ribosomal Slippage
        ribosomal_slippage_zip = zip(accession_list_with_none,ribosomal_slippage_list_with_none_list_none)
        print(ribosomal_slippage_zip)
        ribosomal_slippage_list = []
        for accession_ribosomal_slippage in ribosomal_slippage_zip:
            ribosomal_slippage_list.append("\t".join(accession_ribosomal_slippage))
        print("")
        # Writing Ribosomal Slippage
        print("Step 2: Writing Ribosomal Slippage, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        ribosomal_slippage_string = "\n".join(str(accesion_ribosomal_slippage_value) for accesion_ribosomal_slippage_value in ribosomal_slippage_list)
        print(ribosomal_slippage_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_ribosomal_slippage.txt"), "w")
        handle_outputfile.write(str(ribosomal_slippage_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 22: Ribosomal Slippage Extract Complete")
        print("")
        break

# Standard Name Extract ##################################################
def standardname():
    while True:
        # Module Name
        print("Module 23: Standard Name Extract")
        print("")
        # Parsing Standard Name
        print("Step 1: Parsing Pseudo, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        standard_name_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        standard_name = feature.qualifiers.get('standard_name')
                        accession_list_with_none.append(record.id)
                        standard_name_list_with_none.append(standard_name)
        standard_name_list_with_none_list = [['None'] if standard_name_value is None else standard_name_value for standard_name_value in standard_name_list_with_none]
        standard_name_list_with_none_list_none = []
        for standard_name_none in standard_name_list_with_none_list:
            standard_name_list_with_none_list_none.append('\n'.join(standard_name_none))
        standard_name_zip = zip(accession_list_with_none,standard_name_list_with_none_list_none)
        print(standard_name_zip)
        standard_name_list = []
        for accession_standard_name in standard_name_zip:
            standard_name_list.append("\t".join(accession_standard_name))
        print("")
        # Writing Standard Name
        print("Step 2: Writing Standard Name, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        standard_name_string = "\n".join(str(accesion_standard_name_value) for accesion_standard_name_value in standard_name_list)
        print(standard_name_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_standard_name.txt"), "w")
        handle_outputfile.write(str(standard_name_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 23: Standard Name Extract Complete")
        print("")
        break

# Transcript ID Extract ##################################################
def transcriptid():
    while True:
        # Module Name
        print("Module 24: Transcript ID Extract")
        print("")
        # Parsing Transcript ID
        print("Step 1: Parsing Transcript ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        transcript_id_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        transcript_id = feature.qualifiers.get('transcript_id')
                        accession_list_with_none.append(record.id)
                        transcript_id_list_with_none.append(transcript_id)
        transcript_id_list_with_none_list = [['None'] if transcript_id_value is None else transcript_id_value for transcript_id_value in transcript_id_list_with_none]
        transcript_id_list_with_none_list_none = []
        for transcript_id_none in transcript_id_list_with_none_list:
            transcript_id_list_with_none_list_none.append('\n'.join(transcript_id_none))
        transcript_id_zip = zip(accession_list_with_none,transcript_id_list_with_none_list_none)
        print(transcript_id_zip)
        transcript_id_list = []
        for accession_transcript_id in transcript_id_zip:
            transcript_id_list.append("\t".join(accession_transcript_id))
        print("")
        # Writing Transcript ID
        print("Step 2: Writing Transcript ID, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        transcript_id_string = "\n".join(str(accesion_transcript_id_value) for accesion_transcript_id_value in transcript_id_list)
        print(transcript_id_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_transcript_id.txt"), "w")
        handle_outputfile.write(str(transcript_id_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 24: Transcript ID Extract Complete")
        print("")
        break

# Transl Except Extract ##################################################
def translexcept():
    while True:
        # Module Name
        print("Module 25: Transl Except Extract")
        print("")
        # Parsing Transl Except
        print("Step 1: Parsing Transl Except, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        transl_except_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        transl_except = feature.qualifiers.get('transl_except')
                        accession_list_with_none.append(record.id)
                        transl_except_list_with_none.append(transl_except)
        transl_except_list_with_none_list = [['None'] if transl_except_value is None else transl_except_value for transl_except_value in transl_except_list_with_none]
        transl_except_list_with_none_list_none = []
        for transl_except_none in transl_except_list_with_none_list:
            transl_except_list_with_none_list_none.append('\n'.join(transl_except_none))
        transl_except_zip = zip(accession_list_with_none,transl_except_list_with_none_list_none)
        print(transl_except_zip)
        transl_except_list = []
        for accession_transl_except in transl_except_zip:
            transl_except_list.append("\t".join(accession_transl_except))
        print("")
        # Writing Transl Except
        print("Step 2: Writing Transl Except, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        transl_except_string = "\n".join(str(accesion_transl_except_value) for accesion_transl_except_value in transl_except_list)
        print(transl_except_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_transl_except.txt"), "w")
        handle_outputfile.write(str(transl_except_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 25: Transl Except Extract Complete")
        print("")
        break

# Transl Table Extract ##################################################
def transltable():
    while True:
        # Module Name
        print("Module 26: Transl Table Extract")
        print("")
        # Parsing Transl Table
        print("Step 1: Parsing Transl Table, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        transl_table_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        transl_table = feature.qualifiers.get('transl_table')
                        accession_list_with_none.append(record.id)
                        transl_table_list_with_none.append(transl_table)
        transl_table_list_with_none_list = [['None'] if transl_table_value is None else transl_table_value for transl_table_value in transl_table_list_with_none]
        transl_table_list_with_none_list_none = []
        for transl_table_none in transl_table_list_with_none_list:
            transl_table_list_with_none_list_none.append('\n'.join(transl_table_none))
        transl_table_zip = zip(accession_list_with_none,transl_table_list_with_none_list_none)
        print(transl_table_zip)
        transl_table_list = []
        for accession_transl_table in transl_table_zip:
            transl_table_list.append("\t".join(accession_transl_table))
        print("")
        # Writing Transl Table
        print("Step 2: Writing Transl Table, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        transl_table_string = "\n".join(str(accesion_transl_table_value) for accesion_transl_table_value in transl_table_list)
        print(transl_table_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_transl_table.txt"), "w")
        handle_outputfile.write(str(transl_table_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 26: Transl Table Extract Complete")
        print("")
        break

# Translation Extract ##################################################
def translation():
    while True:
        # Module Name
        print("Module 27: Translation Extract")
        print("")
        # Parsing Translation
        print("Step 1: Parsing Translation, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        translation_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        translation = feature.qualifiers.get('translation')
                        accession_list_with_none.append(record.id)
                        translation_list_with_none.append(translation)
        translation_list_with_none_list = [['None'] if translation_value is None else translation_value for translation_value in translation_list_with_none]
        translation_list_with_none_list_none = []
        for translation_none in translation_list_with_none_list:
            translation_list_with_none_list_none.append('\n'.join(translation_none))
        translation_zip = zip(accession_list_with_none,translation_list_with_none_list_none)
        print(translation_zip)
        translation_list = []
        for accession_translation in translation_zip:
            translation_list.append("\t".join(accession_translation))
        print("")
        # Writing Translation
        print("Step 2: Writing Translation, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        translation_string = "\n".join(str(accesion_translation_value) for accesion_translation_value in translation_list)
        print(translation_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_translation.txt"), "w")
        handle_outputfile.write(str(translation_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 27: Translation Extract Complete")
        print("")
        break

# Usedin Extract ##################################################
def usedin():
    while True:
        # Module Name
        print("Module 28: Usedin Extract")
        print("")
        # Parsing Usedin
        print("Step 1: Parsing Usedin, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        usedin_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        usedin = feature.qualifiers.get('usedin')
                        accession_list_with_none.append(record.id)
                        usedin_list_with_none.append(usedin)
        usedin_list_with_none_list = [['None'] if usedin_value is None else usedin_value for usedin_value in usedin_list_with_none]
        usedin_list_with_none_list_none = []
        for usedin_none in usedin_list_with_none_list:
            usedin_list_with_none_list_none.append('\n'.join(usedin_none))
        usedin_zip = zip(accession_list_with_none,usedin_list_with_none_list_none)
        print(usedin_zip)
        usedin_list = []
        for accession_usedin in usedin_zip:
            usedin_list.append("\t".join(accession_usedin))
        print("")
        # Writing Usedin
        print("Step 2: Writing Usedin, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        usedin_string = "\n".join(str(accesion_usedin_value) for accesion_usedin_value in usedin_list)
        print(usedin_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_usedin.txt"), "w")
        handle_outputfile.write(str(usedin_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module 28: Usedin Extract Complete")
        print("")
        break

# All Extract ##################################################
def allextract():
    while True:
        # Module Name
        print("Module -: All Extract")
        print("")
        # Parsing All
        print("Step 1: Parsing All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        # Codon Start
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        codon_start_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        codon_start = feature.qualifiers.get('codon_start')
                        accession_list_with_none.append(record.id)
                        codon_start_list_with_none.append(codon_start)
        codon_start_list_with_none_list = [['None'] if codon_start_value is None else codon_start_value for codon_start_value in codon_start_list_with_none]
        codon_start_list = []
        for codon_start_none in codon_start_list_with_none_list:
            codon_start_list.append('\n'.join(codon_start_none))
        # DB XREF
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        db_xref_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        db_xref = feature.qualifiers.get('db_xref')
                        accession_list_with_none.append(record.id)
                        db_xref_list_with_none.append(db_xref)
        db_xref_list_with_none_list = [['None'] if db_xref_value is None else db_xref_value for db_xref_value in db_xref_list_with_none]
        db_xref_list = []
        for db_xref_none in db_xref_list_with_none_list:
            db_xref_list.append('\n'.join(db_xref_none))
        # EC Number
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        ec_number_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        ec_number = feature.qualifiers.get('ec_number')
                        accession_list_with_none.append(record.id)
                        ec_number_list_with_none.append(ec_number)
        ec_number_list_with_none_list = [['None'] if ec_number_value is None else ec_number_value for ec_number_value in ec_number_list_with_none]
        ec_number_list = []
        for ec_number_none in ec_number_list_with_none_list:
            ec_number_list.append('\n'.join(ec_number_none))
        # Evidence
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        evidence_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        evidence = feature.qualifiers.get('evidence')
                        accession_list_with_none.append(record.id)
                        evidence_list_with_none.append(evidence)
        evidence_list_with_none_list = [['None'] if evidence_value is None else evidence_value for evidence_value in evidence_list_with_none]
        evidence_list = []
        for evidence_none in evidence_list_with_none_list:
            evidence_list.append('\n'.join(evidence_none))
        # Exception
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        exception_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        exception = feature.qualifiers.get('exception')
                        accession_list_with_none.append(record.id)
                        exception_list_with_none.append(exception)
        exception_list_with_none_list = [['None'] if exception_value is None else exception_value for exception_value in exception_list_with_none]
        exception_list = []
        for exception_none in exception_list_with_none_list:
            exception_list.append('\n'.join(exception_none))
        # Experiment
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        experiment_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        experiment = feature.qualifiers.get('experiment')
                        accession_list_with_none.append(record.id)
                        experiment_list_with_none.append(experiment)
        experiment_list_with_none_list = [['None'] if experiment_value is None else experiment_value for experiment_value in experiment_list_with_none]
        experiment_list = []
        for experiment_none in experiment_list_with_none_list:
            experiment_list.append('\n'.join(experiment_none))
        # Function
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        function_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        function = feature.qualifiers.get('function')
                        accession_list_with_none.append(record.id)
                        function_list_with_none.append(function)
        function_list_with_none_list = [['None'] if function_value is None else function_value for function_value in function_list_with_none]
        function_list = []
        for function_none in function_list_with_none_list:
            function_list.append('\n'.join(function_none))
        # Gene
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        gene_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        gene = feature.qualifiers.get('gene')
                        accession_list_with_none.append(record.id)
                        gene_list_with_none.append(gene)
        gene_list_with_none_list = [['None'] if gene_value is None else gene_value for gene_value in gene_list_with_none]
        gene_list = []
        for gene_none in gene_list_with_none_list:
            gene_list.append('\n'.join(gene_none))
        # Gene Synonym
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        gene_synonym_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        gene_synonym = feature.qualifiers.get('gene_synonym')
                        accession_list_with_none.append(record.id)
                        gene_synonym_list_with_none.append(gene_synonym)
        gene_synonym_list_with_none_list = [['None'] if gene_synonym_value is None else gene_synonym_value for gene_synonym_value in gene_synonym_list_with_none]
        gene_synonym_list = []
        for gene_synonym_none in gene_synonym_list_with_none_list:
            gene_synonym_list.append('\n'.join(gene_synonym_none))
        # Go Component
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        go_component_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        go_component = feature.qualifiers.get('go_component')
                        accession_list_with_none.append(record.id)
                        go_component_list_with_none.append(go_component)
        go_component_list_with_none_list = [['None'] if go_component_value is None else go_component_value for go_component_value in go_component_list_with_none]
        go_component_list = []
        for go_component_none in go_component_list_with_none_list:
            go_component_list.append('\n'.join(go_component_none))
        # Go Function
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        go_function_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        go_function = feature.qualifiers.get('go_function')
                        accession_list_with_none.append(record.id)
                        go_function_list_with_none.append(go_function)
        go_function_list_with_none_list = [['None'] if go_function_value is None else go_function_value for go_function_value in go_function_list_with_none]
        go_function_list = []
        for go_function_none in go_function_list_with_none_list:
            go_function_list.append('\n'.join(go_function_none))
        # Go Process
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        go_process_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        go_process = feature.qualifiers.get('go_process')
                        accession_list_with_none.append(record.id)
                        go_process_list_with_none.append(go_process)
        go_process_list_with_none_list = [['None'] if go_process_value is None else go_process_value for go_process_value in go_process_list_with_none]
        go_process_list = []
        for go_process_none in go_process_list_with_none_list:
            go_process_list.append('\n'.join(go_process_none))
        # Inference
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        inference_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        inference = feature.qualifiers.get('inference')
                        accession_list_with_none.append(record.id)
                        inference_list_with_none.append(inference)
        inference_list_with_none_list = [['None'] if inference_value is None else inference_value for inference_value in inference_list_with_none]
        inference_list = []
        for inference_none in inference_list_with_none_list:
            inference_list.append('\n'.join(inference_none))
        # Label
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        label_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        label = feature.qualifiers.get('label')
                        accession_list_with_none.append(record.id)
                        label_list_with_none.append(label)
        label_list_with_none_list = [['None'] if label_value is None else label_value for label_value in label_list_with_none]
        label_list = []
        for label_none in label_list_with_none_list:
            label_list.append('\n'.join(label_none))
        # Location
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        location_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        #location = feature.qualifiers.get('location')
                        location = feature.location
                        accession_list_with_none.append(record.id)
                        location_list_with_none.append(feature.location)
        location_list_with_none_list = [['None'] if location_value is None else location_value for location_value in location_list_with_none]
        location_list = []
        for location_none in location_list_with_none_list:
            #location_list_with_none_list_none.append('\n'.join(location_none))
            location_list.append(location_none)
        # Locus Tag
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        locus_tag_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
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
                    if feature.type == "CDS":
                        note = feature.qualifiers.get('note')
                        note_list_with_none.append(note)
        note_list_with_none_list = [['None'] if note_value is None else note_value for note_value in note_list_with_none]
        note_list = []
        for note_none in note_list_with_none_list:
            note_list.append('\n'.join(note_none))
        # Product
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        product_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        product = feature.qualifiers.get('product')
                        product_list_with_none.append(product)
        product_list_with_none_list = [['None'] if product_value is None else product_value for product_value in product_list_with_none]
        product_list = []
        for product_none in product_list_with_none_list:
            product_list.append('\n'.join(product_none))
        # Prot DESC
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        prot_desc_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        prot_desc = feature.qualifiers.get('prot_desc')
                        accession_list_with_none.append(record.id)
                        prot_desc_list_with_none.append(prot_desc)
        prot_desc_list_with_none_list = [['None'] if prot_desc_value is None else prot_desc_value for prot_desc_value in prot_desc_list_with_none]
        prot_desc_list = []
        for prot_desc_none in prot_desc_list_with_none_list:
            prot_desc_list.append('\n'.join(prot_desc_none))
        # Protein ID
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        protein_id_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        protein_id = feature.qualifiers.get('protein_id')
                        protein_id_list_with_none.append(protein_id)
        protein_id_list_with_none_list = [['None'] if protein_id_value is None else protein_id_value for protein_id_value in protein_id_list_with_none]
        protein_id_list = []
        for protein_id_none in protein_id_list_with_none_list:
            protein_id_list.append('\n'.join(protein_id_none))
        # Pseudo
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        pseudo_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        pseudo = feature.qualifiers.get('pseudo')
                        accession_list_with_none.append(record.id)
                        pseudo_list_with_none.append(pseudo)
        pseudo_list_with_none_list = [['None'] if pseudo_value is None else pseudo_value for pseudo_value in pseudo_list_with_none]
        pseudo_list = []
        for pseudo_none in pseudo_list_with_none_list:
            pseudo_list.append('\n'.join(pseudo_none))
        # Ribosomal Slippage
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        ribosomal_slippage_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        ribosomal_slippage = feature.qualifiers.get('ribosomal_slippage')
                        ribosomal_slippage_list_with_none.append(ribosomal_slippage)
        ribosomal_slippage_list_with_none_list = [['None'] if ribosomal_slippage_value is None else ribosomal_slippage_value for ribosomal_slippage_value in ribosomal_slippage_list_with_none]
        ribosomal_slippage_list = []
        for ribosomal_slippage_none in ribosomal_slippage_list_with_none_list:
            ribosomal_slippage_list.append('\n'.join(ribosomal_slippage_none))
        # Transcript ID
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        transcript_id_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        transcript_id = feature.qualifiers.get('transcript_id')
                        accession_list_with_none.append(record.id)
                        transcript_id_list_with_none.append(transcript_id)
        transcript_id_list_with_none_list = [['None'] if transcript_id_value is None else transcript_id_value for transcript_id_value in transcript_id_list_with_none]
        transcript_id_list = []
        for transcript_id_none in transcript_id_list_with_none_list:
            transcript_id_list.append('\n'.join(transcript_id_none))
        # Transl Except
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        transl_except_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        transl_except = feature.qualifiers.get('transl_except')
                        accession_list_with_none.append(record.id)
                        transl_except_list_with_none.append(transl_except)
        transl_except_list_with_none_list = [['None'] if transl_except_value is None else transl_except_value for transl_except_value in transl_except_list_with_none]
        transl_except_list = []
        for transl_except_none in transl_except_list_with_none_list:
            transl_except_list.append('\n'.join(transl_except_none))
        # Transl Table
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        transl_table_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        transl_table = feature.qualifiers.get('transl_table')
                        accession_list_with_none.append(record.id)
                        transl_table_list_with_none.append(transl_table)
        transl_table_list_with_none_list = [['None'] if transl_table_value is None else transl_table_value for transl_table_value in transl_table_list_with_none]
        transl_table_list = []
        for transl_table_none in transl_table_list_with_none_list:
            transl_table_list.append('\n'.join(transl_table_none))
         # Translation
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        translation_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        translation = feature.qualifiers.get('translation')
                        translation_list_with_none.append(translation)
        translation_list_with_none_list = [['None'] if translation_value is None else translation_value for translation_value in translation_list_with_none]
        translation_list = []
        for translation_none in translation_list_with_none_list:
            translation_list.append('\n'.join(translation_none))
        # Usedin
        handle_parse_db = SeqIO.parse(os.path.expanduser("~/extractor/local/import/sequence.gbk"), "genbank")
        accession_list_with_none = []
        usedin_list_with_none = []
        for record in handle_parse_db:
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        usedin = feature.qualifiers.get('usedin')
                        accession_list_with_none.append(record.id)
                        usedin_list_with_none.append(usedin)
        usedin_list_with_none_list = [['None'] if usedin_value is None else usedin_value for usedin_value in usedin_list_with_none]
        usedin_list = []
        for usedin_none in usedin_list_with_none_list:
            usedin_list.append('\n'.join(usedin_none))
        # Zipping All Extract
        all_extract_zip = zip(accession_list_with_none,codon_start_list,db_xref_list,ec_number_list,evidence_list,exception_list,experiment_list,function_list,gene_list,gene_synonym_list,go_component_list,go_function_list,go_process_list,inference_list,label_list,location_list,locus_tag_list,note_list,product_list,prot_desc_list,protein_id_list,pseudo_list,ribosomal_slippage_list,transcript_id_list,transl_except_list,transl_table_list,translation_list,usedin_list)
        print(all_extract_zip)
        all_extract_list = []
        for all_feature in all_extract_zip:
            all_extract_list.append("\t".join(all_feature))
        print("")
        # Writing All
        print("Step 2 : Writing All, please wait...")
        with progressbar.ProgressBar(max_value=100) as bar:
            for i in range(100):
                time.sleep(0.1)
                bar.update(i)
        print("")
        all_string = '\n'.join(str(all_value) for all_value in all_extract_list)
        print(all_string)
        handle_outputfile = open(os.path.expanduser("~/extractor/local/export/cds_all_extract.txt"), "w")
        handle_outputfile.write(str(all_string))
        handle_outputfile.close()
        print("")
        # Module Complete
        print("Module -: All Extract Complete")
        print("")
        break

# Module Main ##################################################
def maincds():
    while True:
        print("")
        print("\tMenu: Extract CDS")
        print("")
        print("\t\tModule 1.  Codon Start \t\t\t: Button[1]")
        print("\t\tModule 2.  DB XREF \t\t\t: Button[2]")
        print("\t\tModule 3.  EC Number \t\t\t: Button[3]")
        print("\t\tModule 4.  Evidence \t\t\t: Button[4]")
        print("\t\tModule 5.  Exception \t\t\t: Button[5]")
        print("\t\tModule 6.  Experiment \t\t\t: Button[6]")
        print("\t\tModule 7.  Function \t\t\t: Button[7]")
        print("\t\tModule 8.  Gene \t\t\t: Button[8]")
        print("\t\tModule 9.  Gene Synonym \t\t: Button[9]")
        print("\t\tModule 10. Go Component \t\t: Button[10]")
        print("\t\tModule 11. Go Function \t\t\t: Button[11]")
        print("\t\tModule 12. Go Process \t\t\t: Button[12]")
        print("\t\tModule 13. Inference \t\t\t: Button[13]")
        print("\t\tModule 14. Label \t\t\t: Button[14]")
        print("\t\tModule 15. Location \t\t\t: Button[15]")
        print("\t\tModule 16. Locus Tag \t\t\t: Button[16]")
        print("\t\tModule 17. Note \t\t\t: Button[17]")
        print("\t\tModule 18. Product \t\t\t: Button[18]")
        print("\t\tModule 19. Prot DESC \t\t\t: Button[19]")
        print("\t\tModule 20. Protein ID \t\t\t: Button[20]")
        print("\t\tModule 21. Pseudo \t\t\t: Button[21]")
        print("\t\tModule 22. Ribosomal Slippage \t\t: Button[22]")
        print("\t\tModule 23. Standard Name \t\t: Button[23]")
    
        print("\t\tModule 24. Transcript ID \t\t: Button[24]")
        print("\t\tModule 25. Transl Except \t\t: Button[25]")
        print("\t\tModule 26. Transl Table \t\t: Button[26]")
        print("\t\tModule 27. Translation \t\t\t: Button[27]")
        print("\t\tModule 28. Usedin \t\t\t: Button[28]")
        #print("\t\tModule -. All Extract \t\t\t: Button[-]")
        print("\t\tModule 29. Main Menu \t\t\t: Button[29]")
        print("\t\tModule 30. Exit. \t\t\t: Button[30]")
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
                codonstart()
            elif button == 2:
                print("")
                dbxref()           
            elif button == 3:
                print("")
                ecnumber()
            elif button == 4:
                print("")
                evidence()
            elif button == 5:
                print("")
                exception()               
            elif button == 6:
                print("")
                experiment()
            elif button == 7:
                print("")
                function()
            elif button == 8:
                print("")
                gene()
            elif button == 9:
                print("")
                genesynonym()
            elif button == 10:
                print("")
                gocomponent()
            elif button == 11:
                print("")
                gofunction()
            elif button == 12:
                print("")
                goprocess()
            elif button == 13:
                print("")
                inference()               
            elif button == 14:
                print("")
                label()
            elif button == 15:
                print("")
                location()
            elif button == 16:
                print("")
                locustag()
            elif button == 17:
                print("")
                note()
            elif button == 18:
                print("")
                product()
            elif button == 19:
                print("")
                protdesc()          
            elif button == 20:
                print("")
                proteinid()
            elif button == 21:
                print("")
                pseudo()               
            elif button == 22:
                print("")
                ribosomalslippage()
            elif button == 23:
                print("")
                standardname()
            elif button == 24:
                print("")
                transcriptid()
            elif button == 25:
                print("")
                translexcept()               
            elif button == 26:
                print("")
                transltable()
            elif button == 27:
                print("")
                translation()
            elif button == 28:
                print("")
                usedin()
            #elif button == 27:
                #print("")
                #allextract()
            elif button == 29:
                print("")
                import main
                main.main()
            elif button == 30:
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

maincds()
