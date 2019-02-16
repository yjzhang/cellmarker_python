import os
import sqlite3

PATH = os.path.dirname(__file__)
DB_DIR = os.path.join(PATH, 'data', 'cell_marker.db')

def get_all_genes():
    """
    Returns a list of all unique gene symbols.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT DISTINCT gene FROM gene_indices')
    results = C.fetchall()
    conn.close()
    return [x[0] for x in results]

def get_all_cells():
    """
    Returns a list of all unique cell names.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT DISTINCT cellName FROM cell_gene')
    results = C.fetchall()
    conn.close()
    return [x[0] for x in results]

def get_genes_indices(genes):
    """
    Given a list of genes, this returns a dict that maps gene symbols
    to all the indices that correspond to the given gene.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    gene_indices = {}
    for gene in genes:
        C.execute('SELECT row_id FROM gene_indices WHERE gene=?', (gene,))
        results = C.fetchall()
        gene_indices[gene] = [x[0] for x in results]
    conn.close()
    return gene_indices

def get_cell_genes(cell):
    """
    Given the name of a cell, this returns a list of all genes associated with that cell.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT gene FROM cell_gene WHERE cellName=?', (cell,))
    results = C.fetchall()
    results = [x[0] for x in results]
    conn.close()
    return results

def get_papers_cell_genes(cell, genes):
    """
    """

def hypergeometric_test(genes, cells_or_tissues='cells'):
    """
    Uses a hypergeometric test to identify the most relevant cell types.
    """
    # TODO: cells_or_tissues is unused right now.
    from scipy import stats
    all_cells = get_all_cells()
    all_genes = get_all_genes()
    cell_p_vals = {}
    genes = set(genes)
    for cell in all_cells:
        cell_genes = set(get_cell_genes(cell))
        k = len(genes.intersection(cell_genes))
        pv = stats.hypergeom.cdf(k, len(all_genes), len(cell_genes), len(genes))
        cell_p_vals[cell] = 1 - pv
    cell_p_vals = list(cell_p_vals.items())
    cell_p_vals.sort(key=lambda x: x[1])
    return cell_p_vals

# TODO: what is a more sophisticated test that accounts for the same genes being present in many different cell types?
