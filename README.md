# NCBI-Extractor

Using for simple NCBI sequences filters.

Support:
  1. Python3

Requires:
  1. Biopython
  2. GeoPy
  3. Progressbar2
  4. Requests

Setup:
  1. Putting "NCBI-Extractor" to home directory.
  2. Putting "sequence.gbk" or "sequence.gb.txt" or "sequence.fasta" to import folder.

Compiler:
  1. cd NCBI-Extractor/extractor/local/script/
  2. python3 main.py

Documents:
  1. Main Menu: Feature Extractor 6 modules includes Convert File, Extract Description, Extract Source, Extract Gene, Extract CDS and Extract MISC Features.
  2. Module 1 using file "sequence.gb.txt" convert to .fasta, .gbk, using file "sequence.fasta" convert to .gbk.
  3. Module 2-6 using file "sequence.gbk" to extract info of Description, Source, Gene, CDS and MISC Features.
  *** Field of NCBI Sequence see in references ***

References:
  1. https://www.ncbi.nlm.nih.gov/books/NBK53701/
  2. https://www.ncbi.nlm.nih.gov/genbank/samplerecord/
  3. https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html

Bug report:
chalermpong.int@biotec.or.th
