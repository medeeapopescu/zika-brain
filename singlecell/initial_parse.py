'''
Script to test parsing of single cell data.

'''
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#plt.switch_backend('agg')

root_data = 'C:/Users/medee/spring_rotation_data/data/'


# Functions / classes
def load_data(kind):
    if kind == 'control':
        fn = root_data+'counts_018.tsv'
    elif kind == 'exp':
        fn = root_data+'counts_019.tsv'
    else:
        raise ValueError('Kind of data not understood: {:}'.format(kind))

    data = pd.read_csv(fn, sep='\t', index_col=0)
    return data


def get_names(table):
    return table.index.tolist() #print gene names to new list


def filter_cells_and_genes(data):
    good_cells = ((data >= 3).sum(axis=0) >= 1000) #cells that express at least 1000 genes of at least 3 reads/gene
    good_genes = ((data >= 10).sum(axis=1) >= 10) #genes that are seen in at least 10 cells with at least 10 reads/gene
    data_good = data.loc[good_genes, good_cells] #list of good cells and good gene counts per cell
    coverage = data_good.sum(axis=0)
    data_good = 1.0 * data_good / coverage #normalize expression levels
    return {'data': data_good,
            'coverage': coverage}


def select_top_genes(data, number_of_genes=300):
    top_300 = data.mean(axis=1).nlargest(number_of_genes+5).index[5:] #pull out top 300 expressed genes by averaging expression over all cells, excluding first 4 (spike ins,ribo)
    data_300 = data.loc[top_300] #pull out expression levels of top 300 genes in set of good cells  
    return {'data':data_300,
            'genes': top_300}


def get_gene_names(gene_ids=None, only_good=False):
    
    if gene_ids is not None: #list of geneIDs provided, should only run this once to get file 
        with open(root_data+'gene_names.txt', 'rt') as f:
            gene_names = [line.rstrip('\n') for line in f]
        df = pd.DataFrame([gene_ids, gene_names], index=['id', 'name']).T
        df.to_csv('spring_rotation_data/data/gene_IDs_and_names.tsv', sep='\t', index=False)
        return df
    elif not only_good: 
        return pd.read_csv(root_data+'gene_IDs_and_names.tsv', sep='\t', index_col=0)
    else:
        return pd.read_csv(root_data+'good_geneIDs_and_names.tsv', sep='\t', index_col=0)


def make_tsne(data):
    from sklearn.manifold import TSNE #make 2D tSNE plot 
    model = TSNE(n_components=2)
    tsne = model.fit_transform(data.values.T) 
    tsne = pd.DataFrame(tsne.T, columns=data.columns, index=['tSNE1', 'tSNE2'])
    return tsne 


def plot_tsne(tsne, title='', labels=None, cmap='viridis'):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    if labels is None:
        ax.scatter(tsne.loc['tSNE1'], tsne.loc['tSNE2'], s=15, color='darkred', alpha=0.3, linewidth=0.5, zorder=10)
    else:
        ax.scatter(tsne.loc['tSNE1'], tsne.loc['tSNE2'], s=15, c=labels, alpha=0.3, linewidth=0.5, zorder=10, cmap=cmap)
    ax.set_xlabel('tSNE1')
    ax.set_ylabel('tSNE2')
    ax.grid(True)
    if title:
        ax.set_title(title)
    return {'fig': fig, 'ax': ax}


def kmeans_cluster(data, n_clusters, tsne=None):
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters)
    km_labels = kmeans.fit_predict(data.values.T)
    if tsne is not None:
        plot_tsne(tsne, labels=km_labels)
    return km_labels


def make_pca(data, n_components=2):
    from sklearn.decomposition import PCA 
    model = PCA(n_components=n_components)
    pcs = model.fit_transform(data.values.T)
    pcs = pd.DataFrame(pcs.T, columns=data.columns, index=['PC'+str(i+1) for i in range(n_components)])
    return pcs

def plot_pca(pcs, title = ''):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.scatter(pcs.loc['PC1'], pcs.loc['PC2'], s=15, color='darkred', alpha=0.3, linewidth=0.5, zorder=10)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.grid(True)
    if title:
        ax.set_title(title)
    return {'fig': fig, 'ax': ax}


def hdbscan_cluster(data, tsne=None):
    from sklearn.cluster import DBSCAN
    import hdbscan
    hdb_labels = hdb.labels_
    n_clusters_ = len(set(hdb_labels)) - (1 if -1 in hdb_labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    if tsne is not None:
       plot_tsne(tsne, labels=hdb_labels)
    return hdb_labels

# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parse single cell data')
    parser.add_argument('kind', choices=['control', 'exp'],
                        help='Choose experiment to load')
    parser.add_argument('--filter', action='store_true',
                        help='Filter the data')
    parser.add_argument('--top_genes', type=int, default=0,
                        help='Use only top x genes')
    parser.add_argument('--translatenames', action='store_true',
                        help='Change IDs to names')
    args = parser.parse_args()

    print('Loading data...')
    data = load_data(kind=args.kind)
    print('yay!')

    if args.translatenames:
        good_gene_IDs_and_names = get_gene_names(only_good=True)
        data = data.loc[good_gene_IDs_and_names.index]
        data.index = good_gene_IDs_and_names['name']

    if args.filter:
        # NOTE: this also normalizes
        dic = filter_cells_and_genes(data)
        data = dic['data']

    if args.top_genes>0:
        dic = select_top_genes(data, number_of_genes=args.top_genes)
        data_top = dic['data']

    #table= get_gene_names()
    #dic = Counter(table['name'].values)
    #good_ids = [row['id'] for _, row in table.iterrows() if dic[row['name']] == 1]

    print('Making tsne!')
    tsne = make_tsne(data_top)
    dic_tsne = plot_tsne(tsne)

    print('Making pca!')
    pcs = make_pca(data_top)
    dic_pcs = plot_pca(pcs)

    print('Clustering!')
    labels = kmeans_cluster(data_top, 6, tsne)

    plt.ion()
    plt.show()
    #plt.save()

    data_TRAC = data.loc['TRAC'] #TCR
    df = pd.DataFrame([data_TRAC.values, labels], index=['TRAC', 'label']).T
    df.groupby(['label']).mean()
    plot_tsne(tsne, labels =np.log2(10+ data.loc['TRAC']), cmap = 'Spectral')

    data_CD19 = data.loc['CD19']
    df = pd.DataFrame([data_CD19.values, labels], index=['CD19', 'label']).T
    df.groupby(['label']).mean()
    plot_tsne(tsne, labels =np.log2(10+ data.loc['CD19']), cmap = 'Spectral')

    data_NCAM1 = data.loc['NCAM1'] #CD56
    df = pd.DataFrame([data_NCAM1.values, labels], index=['NCAM1', 'label']).T
    df.groupby(['label']).mean()
    plot_tsne(tsne, labels =np.log2(10+ data.loc['NCAM1']), cmap = 'Spectral')
    
    #data_NCAM1 = data.loc['FCGR3'] #CD16
    #df = pd.DataFrame([data_NCAM1.values, labels], index=['FCGR3', 'label']).T
    #df.groupby(['label']).mean()
    #plot_tsne(tsne, labels =np.log2(10+ data.loc['FCGR3']), cmap = 'Spectral')

    #tsne.T.groupby(['labels']).mean() #cluster center coordinates 
