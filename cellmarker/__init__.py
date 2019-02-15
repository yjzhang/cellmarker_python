import os
import sqlite3

PATH = os.path.dirname(__file__)
DB_DIR = os.path.join(PATH, 'data', 'cell_marker.db')

def get_all_genes():
    """
    Returns a list of all unique gene symbols.
    """
    conn = sqlite3.connect(DB_DIR)
    c = conn.cursor()
    c.execute('SELECT DISTINCT gene FROM gene_indices')
    results = c.fetchall()
    return [x[0] for x in results]

def get_genes_indices(genes):
    """
    Given a list of genes, this returns a dict that maps gene symbols
    to all the indices that correspond to the given gene.
    """
    conn = sqlite3.connect(DB_DIR)
    c = conn.cursor()
    gene_indices = {}
    for gene in genes:
        c.execute('SELECT row_id FROM gene_indices WHERE gene=?', (gene,))
        results = c.fetchall()
        gene_indices[gene] = [x[0] for x in results]
    return gene_indices

