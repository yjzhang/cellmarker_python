import os
import sqlite3

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

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

def get_all_cell_cls(cells_or_tissues='cells'):
    """
    Returns a list of (cell, cellontology id) tuples.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT DISTINCT cellName, CellOntologyID FROM cell_cl')
    results = C.fetchall()
    conn.close()
    return results

def get_all_species():
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT DISTINCT species FROM cell_gene')
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

@lru_cache(maxsize=None)
def get_cell_genes(cell, cells_or_tissues='cells', species='all'):
    """
    Given the name of a cell, this returns a list of all genes associated with that cell.
    Alternatively, if cells_or_tissues is 'tissues', this returns a list of
    all genes associated with that tissue.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    if cells_or_tissues == 'cells':
        if species != 'all':
            C.execute('SELECT gene FROM cell_gene WHERE cellName=? AND species=?', (cell, species))
        else:
            C.execute('SELECT gene FROM cell_gene WHERE cellName=?', (cell,))
    elif cells_or_tissues == 'tissues':
        if species != 'all':
            C.execute('SELECT gene FROM tissue_gene WHERE tissueType=? AND species=?', (cell, species))
        else:
            C.execute('SELECT gene FROM tissue_gene WHERE tissueType=?', (cell,))
    results = C.fetchall()
    results = [x[0] for x in results]
    conn.close()
    return results

@lru_cache(maxsize=None)
def get_papers_cell_gene(cell, gene, species='all'):
    """
    Returns all PMIDs associated with the cell type - gene combination.
    """
    # pubmed link format: https://www.ncbi.nlm.nih.gov/pubmed/<pmid>
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    if species == 'all':
        C.execute('SELECT pmid FROM cell_gene_pmid WHERE cellName=? AND gene=?', (cell, gene))
    else:
        C.execute('SELECT pmid FROM cell_gene_pmid WHERE cellName=? AND gene=? AND species=?', (cell, gene, species))
    results = C.fetchall()
    conn.close()
    if len(results) > 0:
        return [x[0] for x in results]
    else:
        return []

def hypergeometric_test(genes, cells_or_tissues='cells', species='all', return_header=False, return_cl=False):
    """
    Uses a hypergeometric test to identify the most relevant cell types.

    Returns:
        list of 4-tuples: cell name, p-value, overlapping genes, pmids
        in order of ascending p-value.
    """
    from scipy import stats
    genes = [x.upper() for x in genes]
    if return_cl:
        all_cell_cls = get_all_cell_cls(cells_or_tissues)
        all_cells = [x[0] for x in all_cell_cls]
        cells_to_cls = {c[0]: c[1] for c in all_cell_cls}
    else:
        all_cells = get_all_cells(cells_or_tissues)
    all_genes = get_all_genes()
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for cell in all_cells:
        cell_genes = set(get_cell_genes(cell, cells_or_tissues, species=species))
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
    if return_header:
        header = ['Cell', 'P-value', 'Overlapping Genes', 'PMIDs']
        cell_p_vals = [header] + cell_p_vals
    if return_cl:
        start_id = 0
        if return_header:
            cell_p_vals[0] = ['Cell Ontology ID'] + cell_p_vals[0]
            start_id = 1
        for i in range(start_id, len(cell_p_vals)):
            cell_p_vals[i] = (cells_to_cls[cell_p_vals[i][0]],) + cell_p_vals[i]
    return cell_p_vals

@lru_cache(maxsize=None)
def get_all_child_cells():
    """
    Args:
        cell_cls (dict): dict of cell type name : ontology thing

    Returns a dict of cell type : list of child cell types.
    """
    from owlready2 import get_ontology, default_world
    cl_db_dir = os.path.join(PATH, 'data', 'cl.db')
    default_world.set_backend(filename=cl_db_dir)
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT * FROM cell_cl')
    results = C.fetchall()
    cell_cl = {r[0]: r[1] for r in results}
    conn.close()
    # open cell ontology
    onto = get_ontology('cl.owl').load()
    # http://owlready.8326.n8.nabble.com/Accessing-class-by-its-name-in-owlready2-td457.html 
    namespace = onto.get_namespace("http://purl.obolibrary.org/obo/")
    cell_cls = dict()
    for cell, cl in cell_cl.items():
        cell_cls[cell] = namespace[cl]
    cl_cells = {cl.get_name(cl) : cell for cell,  cl in cell_cls.items()}

    cells_children = dict()
    for cell_name, cl in cell_cls.items():
        children = list(cl.descendants())
        cells_children[cell_name] = []
        for child in children:
            child_name = child.get_name(child)
            if child_name in cl_cells:
                cell_type = cl_cells[child.get_name(child)]
                cells_children[cell_name].append(cell_type)
    return cells_children


def hierarchical_hypergeom_test(genes, cells_or_tissues='cells', species='all', return_header=False):
    """
    Returns a list and a dict: {cell: (pval, overlapping genes, pmids, child cell types), ...]
    """
    from scipy import stats
    all_cells = get_all_cells(cells_or_tissues)
    all_genes = get_all_genes()
    cells_children = get_all_child_cells()
    # TODO: add all 'child' genes?
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    # cells_genes_mapping contains
    cells_genes_mapping = {}
    for cell in all_cells:
        cells_genes_mapping[cell] = {cell: get_cell_genes(cell, cells_or_tissues, species)}
    for cell in all_cells:
        if cell in cells_children:
            for child in cells_children[cell]:
                cells_genes_mapping[cell][child] = get_cell_genes(child, cells_or_tissues, species)
    for cell, cell_genes in cells_genes_mapping.items():
        # all_cell_genes contains all genes for this and child cell types.
        all_cell_genes = set()
        for gene_list in cell_genes.values():
            all_cell_genes.update(gene_list)
        overlapping_genes = list(genes.intersection(all_cell_genes))
        if len(overlapping_genes) == 0:
            continue
        # TODO: deal with child cell types?
        pmids = {}
        for child, gene_list in cell_genes.items():
            sub_overlapping_genes = genes.intersection(gene_list)
            for gene in sub_overlapping_genes:
                pmids[gene] = get_papers_cell_gene(cell, gene, species)
        k = len(overlapping_genes)
        pv = stats.hypergeom.cdf(k - 1, len(all_genes), len(all_cell_genes), len(genes))
        cell_p_vals[cell] = (1 - pv, overlapping_genes, pmids)
    cell_p_vals = list(cell_p_vals.items())
    cell_p_vals.sort(key=lambda x: x[1][0])
    # merge items
    cell_p_vals = [(x[0],) + x[1] for x in cell_p_vals]
    if return_header:
        header = ['Cell', 'P-value', 'Overlapping Genes', 'PMIDs', 'Child Cells']
        cell_p_vals = [header] + cell_p_vals
    return cell_p_vals


def heuristic_test(genes, cells_or_tissues='cells', return_header=False):
    # TODO
    pass

# TODO: what is a more sophisticated test that accounts for the same genes being present in many different cell types?
