#!/usr/bin/env python3

# This script is a modification of the script found in Peter Cock's site (http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta/).
# Script borrowed from fhsantanna (https://github.com/fhsantanna/bioinfo_scripts/blob/master/gbk2faa.py) edited by valzip (https://github.com/valzip) to work with python3 and added more qualifiers to output
# Script is able to process multiple gbk entries in single file, e.g.: output from PROKKA_*.gbk. Tested on Python 3.6.9.
# Script requires biopython in python3 (see: https://biopython.org/wiki/Download)
# Usage: python3 py3_gbk2faa.py <input> <output>



import sys
from Bio import GenBank
from Bio import SeqIO

input_handle  = open(sys.argv[1], "r")
output_handle = open(sys.argv[2], "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    #print("Dealing with GenBank record %s" % seq_record.id) #uncomment start of this line if you want to see which record has been processed
    for seq_feature in seq_record.features :
        try: # Without "try", it crashes when it finds a CDS without translation (pseudogene).
            if seq_feature.type=="CDS" :
                assert len(seq_feature.qualifiers['translation'])==1
                seq_qualifiers = seq_feature.qualifiers.keys()                
                if "gene" in seq_qualifiers:
                    gene_nm = seq_feature.qualifiers['gene'][0]
                else:
                    gene_nm = "no_name"
                if "EC_number" in seq_qualifiers:
                    EC_no = seq_feature.qualifiers['EC_number'][0]
                else:
                    EC_no = "no_EC"
                if "db_xref" in seq_qualifiers:
                    ref_db = seq_feature.qualifiers['db_xref'][0]
                else:
                    ref_db = "no_xref"

            
                output_handle.write(">%s,%s,%s,%s,%s,%s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_feature.qualifiers['product'][0],
                    gene_nm,
                    EC_no,
                    ref_db,
                    seq_record.description,
                    seq_feature.qualifiers['translation'][0]))
                pass
        except:
            continue

output_handle.close()
input_handle.close()
print("gbk2faa done")
