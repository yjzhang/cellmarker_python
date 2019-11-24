# build a database using cellmarker for... something?
from collections import defaultdict
import sqlite3
import pandas as pd

data_table = pd.read_csv('data/PanglaoDB_markers_11_Nov_2019.tsv', sep='\t')

panglao_conn = sqlite3.connect('data/panglao.db')
panglao_c = panglao_conn.cursor()

panglao_c.execute('SELECT * FROM cell_dict')
results = panglao_c.fetchall()


conn = sqlite3.connect('cellmarker/data/panglao.db')
c = conn.cursor()


# build a mapping of gene markers to cell types, and a map of cell types to gene markers.
genes_to_cells = defaultdict(lambda: set())
genes_to_tissues = defaultdict(lambda: set())
cells_to_genes = defaultdict(lambda: set())
tissues_to_genes = defaultdict(lambda: set())
genes_to_indices = defaultdict(lambda: [])
cells_genes_to_pmids = defaultdict(lambda: set())
cells_to_cellonto = dict()
species_map = {'Hs': 'Human', 'Mm': 'Mouse'}
for i, row in data_table.iterrows():
    species = row['species']
    species = species.split()
    gene_symbol = row['official gene symbol']
    cell_type = row['cell type']
    tissue_type = row['organ']
    print(species, gene_symbol, cell_type, tissue_type)
    cells_to_cellonto[cell_type] = 0
    genes_to_indices[gene_symbol].append(i)
    for s in species:
        if s not in species_map:
            continue
        s = species_map[s]
        genes_to_cells[gene_symbol, s].add((cell_type))
        genes_to_tissues[gene_symbol, s].add((tissue_type))
        cells_to_genes[cell_type, s].add((gene_symbol))
        tissues_to_genes[tissue_type, s].add((gene_symbol))
        cells_genes_to_pmids[(cell_type, gene_symbol, s)].add(1)
        cells_genes_to_pmids[(tissue_type, gene_symbol, s)].add(1)

print('Calculating tfidf')

# TODO: tfidf transform
# load corpus
import numpy as np
mouse_gene_symbols = [x[1] for x in cells_genes_to_pmids.keys() if x[2] == 'Mouse']
mouse_cells = [x[0] for x in cells_genes_to_pmids.keys() if x[2] == 'Mouse']
human_gene_symbols = [x[1] for x in cells_genes_to_pmids.keys() if x[2] == 'Human']
human_cells = [x[0] for x in cells_genes_to_pmids.keys() if x[2] == 'Human']
mouse_gene_symbols = np.array(list(set(mouse_gene_symbols)))
human_gene_symbols = np.array(list(set(human_gene_symbols)))
mouse_cells = np.array(list(set(mouse_cells)))
human_cells = np.array(list(set(human_cells)))
gene_to_id_human = {x:i for i, x in enumerate(human_gene_symbols)}
gene_to_id_mouse = {x:i for i, x in enumerate(mouse_gene_symbols)}
cell_to_id_human = {x:i for i, x in enumerate(human_cells)}
cell_to_id_mouse = {x:i for i, x in enumerate(mouse_cells)}

human_data = np.zeros((len(human_cells), len(human_gene_symbols)))
mouse_data = np.zeros((len(mouse_cells), len(mouse_gene_symbols)))
for cell_gene, pmids in cells_genes_to_pmids.items():
    cell, gene, species = cell_gene
    print(cell, gene, species)
    if not isinstance(cell, str):
        continue
    if species == 'Human':
        human_data[cell_to_id_human[cell], gene_to_id_human[gene]] = len(pmids)
    if species == 'Mouse':
        mouse_data[cell_to_id_mouse[cell], gene_to_id_mouse[gene]] = len(pmids)


from scipy import sparse
corpus = sparse.csr_matrix(human_data)
corpus_data = corpus.toarray()
cells, genes = corpus.shape
# TODO: create a tf-idf matrix from corpus.mm
# each cell is like a document, each gene is like a word
from sklearn.feature_extraction.text import TfidfTransformer
tfidf = TfidfTransformer()
data_human_tfidf = tfidf.fit_transform(corpus).toarray()
corpus_mouse = sparse.csr_matrix(mouse_data)
data_mouse_tfidf = tfidf.fit_transform(corpus_mouse).toarray()


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

# TODO: separate human and mouse datasets...
# create a table representing a cell type - gene mapping.
try:
    c.execute('CREATE TABLE cell_gene(cellName text, species text, gene text)')
    for cell, genes in cells_to_genes.items():
        cell, species = cell
        for gene in genes:
            c.execute('INSERT INTO cell_gene VALUES (?, ?, ?)', (cell, species, gene))
except:
    pass
try:
    c.execute('CREATE INDEX cell_names_index ON cell_gene(cellName, species)')
except:
    pass

# create a table representing a tissue - gene mapping.
try:
    c.execute('CREATE TABLE tissue_gene (tissueType text, species text, gene text)')
    for tissue, genes in tissues_to_genes.items():
        tissue, species = tissue
        for gene in genes:
            c.execute('INSERT INTO tissue_gene VALUES (?, ?, ?)', (tissue, species, gene))
except:
    pass
try:
    c.execute('CREATE INDEX tissue_type_index ON tissue_gene(tissueType, species)')
except:
    pass

# create a table represent a cell name - cell ontology ID mapping
try:
    c.execute('CREATE TABLE cell_cl (cellName text, CellOntologyID text)')
    for cell_name, onto_id in cells_to_cellonto.items():
        c.execute('INSERT INTO cell_cl VALUES (?, ?)', (cell_name, onto_id))
except:
    pass
try:
    c.execute('CREATE INDEX cell_cl_index ON cell_cl(cellName)')
except:
    pass

# create a table representing a cellName-gene : PMID mapping.
try:
    c.execute('CREATE TABLE cell_gene_pmid (cellName text, gene text, species text, pmid text)')
    for cell_gene, pmids in cells_genes_to_pmids.items():
        cell, gene, species = cell_gene
        for pmid in pmids:
            c.execute('INSERT INTO cell_gene_pmid VALUES (?, ?, ?, ?)', (cell, gene, species, pmid))
except:
    pass
try:
    c.execute('CREATE INDEX cell_gene_pmid_index ON cell_gene_pmid(cellName, gene, species)')
except:
    pass

try:
    c.execute('CREATE TABLE cell_gene_tfidf (cellName text, gene text, species text, tfidf real)')
    for cell_gene, pmids in cells_genes_to_pmids.items():
        cell, gene, species = cell_gene
        if species == 'Human':
            tfidf_val = data_human_tfidf[cell_to_id_human[cell], gene_to_id_human[gene]]
        if species == 'Mouse':
            tfidf_val = data_mouse_tfidf[cell_to_id_mouse[cell], gene_to_id_mouse[gene]]
        c.execute('INSERT INTO cell_gene_tfidf VALUES (?, ?, ?, ?)', (cell, gene, species, tfidf_val))
except:
    pass
try:
    c.execute('CREATE INDEX cell_gene_tfidf_index ON cell_gene_tfidf(cellName, gene, species)')
except:
    pass


conn.commit()
conn.close()
