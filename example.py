
import cellmarker


all_cells = cellmarker.get_all_cells()

cell0 = all_cells[0]
print(cell0)
cell0_genes = cellmarker.get_cell_genes(cell0)
print(cell0_genes)
pvals = cellmarker.hypergeometric_test(cell0_genes)
print(pvals)

# top genes from cluster 0 of http://uncurl-app.yjzhang.com:8888/user/test_10x_400_new/view
genes_ex = ['KCNB1', 'LINC00856', 'REST', 'RP11-35G22.1', 'LCN6', 'LYRM7', 'CEP97', 'FGF18', 'CTD-2514K5.2']
pvals_ex = cellmarker.hypergeometric_test(genes_ex)
print(pvals_ex[:10])

print('Panglao DB test')
print(cellmarker.hypergeometric_test(genes_ex, db_dir=cellmarker.PANGLAO_DB_DIR))
# cluster 10 of http://uncurl-app.yjzhang.com:8888/user/d909e044-63cc-4b2a-942b-0611b0be4837-cerebellum-k10-all-cells/view
genes_brain = """
HDAC9
LMX1A
GRIA2
CACNA2D1
DGKB
PHACTR2
FAM19A2
MID1
PTPRO
CELF4
""".strip().split('\n')
pvals_brain = cellmarker.hypergeometric_test(genes_brain)
print('basic')
print(pvals_brain[:10])


