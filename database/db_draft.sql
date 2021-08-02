use jr429;
DROP TABLE IF EXISTS tornado;

CREATE TABLE tornado (gene_name VARCHAR(30), gene VARCHAR(60), gene_version float, gene_type CHAR(20), HGNC_symbol VARCHAR(25), myAUC float,  avg_diff float,  pwr float,
avg_log2FC_x float, pct_1_auc float,pct_2_auc float,cluster_auc float, avg_log2FC_y float, pct_1_wilcox float, pct_2_wilcox float, pval_adj float, cluster_wilcox float, pval_negbin float, 
avg_log2FC float, pct_1_negbin float, pct_2_negbin float, pval_adj_negbin float,
cluster_negbin float, gene_percent_GC float, d float, PRIMARY KEY (gene_name) );