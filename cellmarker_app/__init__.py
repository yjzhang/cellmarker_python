import json

from flask import Flask, render_template, request
from flask_bootstrap import Bootstrap
from flask_caching import Cache

import cellmarker

cache = Cache()
app = Flask(__name__)
Bootstrap(app)

# maximum file length is 1000MB
app.config['MAX_CONTENT_LENGTH'] = 1000 * 1024 * 1024
# default args to pass to uncurl.run_state_estimation

def pmid_to_link(pmid):
    return '<a href="https://www.ncbi.nlm.nih.gov/pubmed/{0}">{0}</a>'.format(pmid)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/input', methods=['POST'])
def input():
    genes = request.form['genes']
    print(genes)
    method = request.form['test-type']
    print(method)
    cells_or_tissues = request.form['cells-or-tissues']
    print(cells_or_tissues)
    genes = [x.strip().upper() for x in genes.split()]
    print(genes)
    if method == 'hypergeom':
        result = cellmarker.hypergeometric_test(genes, cells_or_tissues)
        cell_types = []
        for i in range(10):
            ri = result[i]
            genes = ri[2]
            gene_pmids = []
            for g in genes:
                gene_pmids.append('{0}: {1}'.format(g, ','.join(pmid_to_link(x) for x in ri[3][g])))
            cell_types.append((ri[0], ri[1], ', '.join(ri[2]), ', '.join(gene_pmids)))
        return render_template('response.html', header=['Cell', 'P-value', 'Overlapping Genes', 'PMIDs'],
                cell_types=cell_types)


@app.route('/api_post', methods=['POST'])
def api_post():
    """
    Returns the json directly...
    """
    genes = request.form['genes']
    method = request.form['test-type']
    cells_or_tissues = request.form['cells-or-tissues']
    genes = [x.strip().upper() for x in genes.split()]
    if method == 'hypergeom':
        result = cellmarker.hypergeometric_test(genes, cells_or_tissues)
        return result[:10]

if __name__ == '__main__':
    app.run(debug=True)
