import numpy as np
import scipy.io
import scipy.sparse
import cellmarker

path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/10x_pure_pooled/'

data = scipy.io.mmread(path + 'data_400_cells.mtx')
data = scipy.sparse.csc_matrix(data)
genes = np.loadtxt(path + 'gene_names_400.tsv', dtype=str)
labels = np.loadtxt(path + 'labs_400_cells.txt').astype(int)
true_labels = {
        0:   'CD19+ b cells',
        1:   'CD14+ monocytes',
        2:   'CD34+',
        3:   'CD4+ t helper',
        4:   'CD56+ nk',
        5:   'CD8+ cytotoxic t',
        6:   'CD4+/CD45RO+ memory t',
        7:   'CD8+/CD45RA+ naive cytotoxic',
        8:   'CD4+/CD45RA+/CD25- naive t',
        9:   'CD4+/CD25 regulatory t'
}

# potentially correct labels
true_labels_cell_ontology = {
        0:   ['B cell'],
        1:   ['Monocyte', 'Classical monocyte', 'Non-classical monocyte'],
        2:   ['Hematopoietic precursor cell', 'Hematopoietic progenitor cell', 'Hematopoietic cell'],
        3:   ['T cell', 'T helper cell', 'CD4+ T cell'],
        4:   ['Natural killer cell', 'Natural killer T (NKT) cell'],
        5:   ['T cell', 'CD8+ cytotoxic T cell'],
        6:   ['T cell', 'CD4+ memory T cell', 'Memory T cell'],
        7:   ['T cell', 'Naive T cell'],
        8:   ['T cell', 'Immature T cell', 'Naive T cell', 'Naive CD4+ T cell'],
        9:   ['T cell', 'Regulatory T (Treg) cell']
}
# TODO: do a diffexp based on labels
from uncurl_analysis import gene_extraction
scores_u, pvals_u = gene_extraction.one_vs_rest_t(data, labels, eps=0.1, test='u')
scores_t, pvals_t = gene_extraction.one_vs_rest_t(data, labels, eps=0.1, test='t')

# TODO: run a systematic experiment
import matplotlib.pyplot as plt
import random
plt.figure(figsize=(12, 7))
# TODO: random gene vs random cell
methods = ['t', 'u', 'ratio', 'random_genes', 'random_cells']
n_genes = [20, 100]
ticks = ['x', 'o', '*', '+', 's']
ranks = 100
all_cell_types = cellmarker.get_all_cells()
for num_genes in n_genes:
    for method, tick in zip(methods, ticks):
        cell_label_results = {}
        for lab in sorted(list(set(labels))):
            if method == 't':
                top_genes = [genes[x[0]] for x in pvals_t[lab][:num_genes]]
            elif method == 'u':
                top_genes = [genes[x[0]] for x in pvals_u[lab][:num_genes]]
            elif method == 'ratio':
                top_genes = [genes[x[0]] for x in scores_t[lab][:num_genes]]
            elif method == 'random_genes':
                top_genes = random.sample(list(genes), num_genes)
            if method != 'random_cells':
                results = cellmarker.hypergeometric_test(top_genes)
                cell_label_results[lab] = results
            elif method == 'random_cells':
                results = [(cell_type, 1.0) for cell_type in random.sample(list(all_cell_types), ranks)]
                cell_label_results[lab] = results
        # get accuracy at all ranks
        result_ranks = []
        prev_rank = -1
        for labi, results in cell_label_results.items():
            true_results = true_labels_cell_ontology[labi]
            rank = ranks
            for i, xr in enumerate(results):
                if xr[0] in true_results:
                    rank = i
                    break
            result_ranks.append(rank)
        result_ranks = np.array(result_ranks)
        accuracies = []
        ranks_list = []
        prev_accuracy = -1
        prev_rank = 0
        # plot accuracies only when change happens?
        for i in range(ranks):
            ac = float(sum(result_ranks <= i))/len(result_ranks)
            if ac != prev_accuracy:
                #if prev_accuracy != -1:
                #    accuracies.append(prev_accuracy)
                #    ranks_list.append(i-1)
                accuracies.append(ac)
                ranks_list.append(i)
                prev_accuracy = ac
                prev_rank = i
        plt.plot(ranks_list, accuracies, '--' + tick, label=method + ' ' + str(num_genes) + ' genes')

plt.grid()
plt.title('Accuracy vs rank on 10x_400')
plt.xlabel('rank')
plt.ylabel('accuracy')
plt.xticks(range(0, ranks, int(ranks/10)))
plt.legend()
plt.savefig('query_accuracy_10x_400.png', dpi=200)


# try a different dataset... Zeisel dataset?
path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/zeisel/'
data = scipy.io.loadmat(path + 'Zeisel.mat')
