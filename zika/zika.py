from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.append('')  # need to import from base
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.process import process, get_parser
from base.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime

attribute_nesting = {'geographic location':['region', 'country', 'city']}

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()

    lineage = 'zika'
    input_data_path = '../fauna/data/'+lineage
    store_data_path = 'store/'+lineage + '_'
    build_data_path = 'build/'+lineage + '_'

    zika = process(input_data_path = input_data_path, store_data_path = store_data_path, build_data_path = build_data_path,
                   reference='zika/metadata/zika_outgroup.gb',
                   proteins=['CA', 'PRO', 'MP', 'ENV', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', 'NS4B'],
                   method='SLSQP')

    fasta_fields = {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country', 8:'db', 10:'authors'}
    zika.load_sequences(fields=fasta_fields)
    zika.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                               s.attributes['date']< datetime(2017,1,1).date())
    dropped_strains = ["THA/PLCal_ZV/2013", "PLCal_ZV"]
    zika.seqs.filter(lambda s: s.id not in dropped_strains)
    zika.seqs.subsample(category = lambda x:(x.attributes['region'],
                                             x.attributes['date'].year,
                                             x.attributes['date'].month), threshold=100)

    zika.align()
    zika.build_tree()
    zika.clock_filter(n_iqd=3, plot=True)
    zika.annotate_tree(Tc=0.005, timetree=True, reroot='best')
    zika.tree.geo_inference('region')
    zika.tree.geo_inference('country')
    zika.export(controls = attribute_nesting)
