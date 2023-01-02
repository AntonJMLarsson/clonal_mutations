import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def main():
    parser = argparse.ArgumentParser(description="Plot XCI output")
    parser.add_argument('-i','--input', type=str, help='Input .tsv file')
    parser.add_argument('-s','--sample', type=str, help='Sample name')
    parser.add_argument('-c','--cutoff', type=int, default=50, help='Cutoff for plots')

    args = parser.parse_args()

    vireo_df = pd.read_csv(args.input, sep='\t', index_col=0)
    vireo_df_filtered = vireo_df[vireo_df['n_vars'] > args.cutoff]

    g = sns.barplot(x='donor_id', y='n_vars', data=vireo_df.groupby('donor_id').count().reset_index())
    g.set_xticklabels(['Allele 1', 'Allele 2', 'Unassigned'])
    g.set_title('{} all cells'.format(args.sample))
    g.set_xlabel('Predicted active X allele')
    g.set_ylabel('Counts (cells)')
    plt.savefig('{}_number_of_predicted_active_X_alleles.pdf'.format(args.sample))
    plt.clf()

    plt.hist(vireo_df['n_vars'], bins=50)
    plt.title('{} all cells'.format(args.sample))
    plt.ylabel('Count (cells)')
    plt.xlabel('Observed X-linked variants')
    plt.savefig('{}_X_variants_per_cell.pdf'.format(args.sample))
    plt.clf()

    plt.hist(vireo_df_filtered['n_vars'], bins=50)
    plt.title('{}, cells with > {} variants detected'.format(args.sample, args.cutoff))
    plt.ylabel('Count (cells)')
    plt.xlabel('Observed X-linked variants')
    plt.savefig('{}_X_variants_per_cell_{}_cutoff.pdf'.format(args.sample, args.cutoff))
    plt.clf()

    g = sns.barplot(x='donor_id', y='n_vars', data=vireo_df_filtered.groupby('donor_id').count().reset_index())
    g.set_xticklabels(['Allele 1', 'Allele 2'])
    g.set_title('{}, cells with > {} variants detected'.format(args.sample, args.cutoff))
    g.set_xlabel('Predicted active X allele')
    g.set_ylabel('Counts (cells)')
    plt.savefig('{}_number_of_predicted_active_X_alleles_{}_cutoff.pdf'.format(args.sample, args.cutoff))
    plt.clf()
    

if __name__ == '__main__':
    main()
