# build a database using cellmarker for... something?
from collections import defaultdict
import sqlite3
import pandas as pd

data = pd.read_table('all_cell_markers.txt')

conn = sqlite3.connect('cellmarker/data/cell_marker.db')
c = conn.cursor()
try:
    data.to_sql('all_cell_markers', conn)
except ValueError as e:
    pass
try:
    c.execute('CREATE INDEX gene_symbol_index ON all_cell_markers("geneSymbol")')
except:
    pass
try:
    c.execute('CREATE INDEX entry_index ON all_cell_markers("index_labels")')
except:
    pass

# build a mapping of gene markers to cell types, and a map of cell types to gene markers.
genes_to_cells = defaultdict(lambda: [])
cells_to_genes = defaultdict(lambda: [])
genes_to_indices = defaultdict(lambda: [])
for i, row in data.iterrows():
    # one of 'Human' or 'Mouse'
    species = row['speciesType']
    tissue = row['tissueType']
    # one of 'Normal cell' or 'Cancer cell'
    cell_type = row['cellType']
    cell_name = row['cellName']
    gene_symbol = row['geneSymbol']
    print(gene_symbol)
    if not isinstance(gene_symbol, str):
        continue
    gene_symbols = [gene_symbol]
    if ',' in gene_symbol:
        gene_symbols = [x.strip(' []') for x in gene_symbol.split(',')]
    uberon_id = row['UberonOntologyID']
    cell_ontology_id = row['CellOntologyID']
    for gene_symbol in gene_symbols:
        genes_to_indices[gene_symbol].append(i)
        cells_to_genes[cell_name].append(gene_symbol)
        genes_to_cells[gene_symbol].append(cell_name)

# create a table representing a gene-index mapping...
try:
    c.execute('CREATE TABLE gene_indices (gene text, row_id integer)')
    for gene, indices in genes_to_indices.items():
        for index in indices:
            c.execute('INSERT INTO gene_indices VALUES (?, ?)', (gene, index))
except:
    pass
try:
    c.execute('CREATE INDEX gene_indices_index ON gene_indices(gene)')
except:
    pass

# TODO: do some sort of hypergeometric test???
# given a list of genes, we want to find both tissue types and cell names that are overrepresented.
# how do we do that???
# how do we calculate p-values for this thing?
