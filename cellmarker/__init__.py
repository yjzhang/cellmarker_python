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

def get_all_cells(cells_or_tissues='cells'):
    """
    Returns a list of all unique cell names, or tissue names.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    if cells_or_tissues == 'cells':
        C.execute('SELECT DISTINCT cellName FROM cell_gene')
    elif cells_or_tissues == 'tissues':
        C.execute('SELECT DISTINCT tissueType FROM tissue_gene')
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

def get_cell_genes(cell, cells_or_tissues='cells'):
    """
    Given the name of a cell, this returns a list of all genes associated with that cell.
    Alternatively, if cells_or_tissues is 'tissues', this returns a list of
    all genes associated with that tissue.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    if cells_or_tissues == 'cells':
        C.execute('SELECT gene FROM cell_gene WHERE cellName=?', (cell,))
    elif cells_or_tissues == 'tissues':
        C.execute('SELECT gene FROM tissue_gene WHERE tissueType=?', (cell,))
    results = C.fetchall()
    results = [x[0] for x in results]
    conn.close()
    return results

def get_papers_cell_gene(cell, gene):
    """
    Returns all PMIDs associated with the cell type - gene combination.
    """
    # pubmed link format: https://www.ncbi.nlm.nih.gov/pubmed/<pmid>
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT pmid FROM cell_gene_pmid WHERE cellName=? AND gene=?', (cell, gene))
    results = C.fetchall()
    conn.close()
    if len(results) > 0:
        return [x[0] for x in results]
    else:
        return []

def hypergeometric_test(genes, cells_or_tissues='cells'):
    """
    Uses a hypergeometric test to identify the most relevant cell types.
    """
    from scipy import stats
    all_cells = get_all_cells(cells_or_tissues)
    all_genes = get_all_genes()
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for cell in all_cells:
        cell_genes = set(get_cell_genes(cell, cells_or_tissues))
        overlapping_genes = list(genes.intersection(cell_genes))
        if len(overlapping_genes) == 0:
            continue
        pmids = {}
        for gene in overlapping_genes:
            pmids[gene] = get_papers_cell_gene(cell, gene)
        k = len(overlapping_genes)
        pv = stats.hypergeom.cdf(k - 1, len(all_genes), len(cell_genes), len(genes))
        cell_p_vals[cell] = (1 - pv, overlapping_genes, pmids)
    cell_p_vals = list(cell_p_vals.items())
    cell_p_vals.sort(key=lambda x: x[1][0])
    # merge items
    cell_p_vals = [(x[0],) + x[1] for x in cell_p_vals]
    return cell_p_vals

# TODO: what is a more sophisticated test that accounts for the same genes being present in many different cell types?
