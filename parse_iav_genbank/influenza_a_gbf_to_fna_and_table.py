#!/usr/bin/env python

## meant for getting sequences and info tables from influenza a genome assemblies
## searches subdirectories in input directory to look for genbank files
## such as those downloaded with ncbi datasets
## makes info table for each segment in input directory
## makes a .fna of each assembly in individual assembly directories

from Bio import GenBank as GenBank
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import sys
from pathlib import Path
import pandas as pd

input_dir = sys.argv[1]

gbf_list = []
for path, subdirs, files in os.walk(input_dir):
    for name in files:
        gbff = os.path.join(path, name)
        if gbff.endswith('.gbff') and os.path.isfile(gbff) and os.path.getsize(gbff) > 0:
                gbf_list.append(gbff)


seq_data_list = []

for gbf_file in gbf_list:

    file_path = str(Path(gbf_file).parents[0])

    assembly = file_path.split('/')[-1:]

    seq_output_file = os.path.join(file_path, f'{assembly[0]}.fna')

    if os.path.isfile(seq_output_file):
        os.remove(seq_output_file)

    for seq_record in SeqIO.parse(gbf_file, "genbank"):

        if seq_record.seq and seq_record.name and seq_record.features:
            serotype = []
            segment = []
            organism = []
            host = []
            db_xref = []

            for feat in seq_record.features:
                if feat.qualifiers.get('serotype'):
                    serotype.append(feat.qualifiers.get('serotype')[0])
                if feat.qualifiers.get('segment'):
                    segment.append(feat.qualifiers.get('segment')[0])
                if feat.qualifiers.get('organism'):
                    organism.append(feat.qualifiers.get('organism')[0])
                if feat.qualifiers.get('host'):
                    host.append(feat.qualifiers.get('host')[0])
                if feat.qualifiers.get('db_xref'):
                    db_xref.append(feat.qualifiers.get('db_xref')[0])
            
            try:
                seq_data_list.append([assembly[0], seq_record.name, serotype[0], 
                                    segment[0], host[0], db_xref[0]])
                
                print(f">{seq_record.name}\n{seq_record.seq}", file = open(seq_output_file, "a"))
            except Exception as e:
                print(seq_record.name, seq_output_file, e) 

seq_data_df = pd.DataFrame(seq_data_list, 
                           columns=["assembly", "accession", "serotype", 
                                    "segment", "host", "db_xref"])


if not seq_data_df.empty:
    seq_data_output_file = os.path.join(input_dir, "Influenza_A_segment_info1.tsv")

    seq_data_df.to_csv(seq_data_output_file,
                            sep = "\t", index = False)