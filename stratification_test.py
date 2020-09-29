#!/home/eanderson/Virtual_Envs/General3/bin/python3
from pathlib import Path
import matplotlib.pyplot as plt
import subprocess as sbp
import pandas as pd
import seaborn as sns
import argparse
import sys, os
import errno


def get_arguments():
    parser = argparse.ArgumentParser(description="Stratify SV calls across ga4gh regions")
    parser.add_argument('VCF_Dir', type=str,
                        help="Path to Truvari VCFs")
    parser.add_argument('-t', '--tsv', type=str, help="path to ga4gh tsv file",
                        default="/home-02/eanderson/Sentieon_Binning_and_Stratification/ga4gh_all_coordonly_2.tsv")
    parser.add_argument('-b', '--bed', type=str, help="path to stratification bed dir",
                        default="/home-02/eanderson/Sentieon_Binning_and_Stratification/stratification_regions")
    parser.add_argument('-o', '--outdir', type=str, help="output path",
                        default="Stratification_Results/")
    args = parser.parse_args()
    return args

# get various arguments as Path objects
args = get_arguments()
print(args.VCF_Dir, type(args.VCF_Dir))
vcf_base_path = Path(args.VCF_Dir)
strat_bed_dir = Path(args.bed)
strat_tsv = Path(args.tsv)
outdir = Path(args.outdir)

# Paths of truvari output fp, tp and fn VCFs
VCFs = [vcf_base_path / "fp.vcf", vcf_base_path / 'tp-call.vcf', vcf_base_path / "fn.vcf"]
stratification_map = {}
# map the names and paths of the various stratification bed files
with open(strat_tsv, "r") as fh:
    for line in fh:
        name, strat_bed = line.rstrip().split('\t')
        stratification_map[name] = strat_bed_dir /  strat_bed

# This just counts lines in a file
def count_lines(file_object):
    counter = 0
    for line in file_object:
        if not line.startswith("#"):
            counter += 1

    return counter


# iterate through fp, fn and tp VCFs
results_list = []
for vcf in VCFs:

    results_dict = {}

    # use bcftools to get the vcf results within the bed regions
    for name, strat_path in stratification_map.items():
        output = sbp.run(['bcftools', 'view', '-O', 'v', '-T', strat_path, vcf], stdout=sbp.PIPE)
        print(output.args, file=sys.stderr)
        # count results
        results_dict[name] = count_lines(output.stdout.decode('utf-8').split("\n"))

    # get total count of entries in the vcf
    with open(vcf, "r") as vcf_all:
        results_dict['All'] = count_lines(vcf_all)
    results_list.append(results_dict)

# Create a dataframe of results, transpose it and get totals and fp/fn rates
df = pd.DataFrame(results_list, index=['FP', 'TP', 'FN'])
df2 = df.T
df2['Total_Vars'] = df2['TP'] + df2['FN']
df2['FP_Rate'] = df2['FP'] / df2['Total_Vars']
df2['FN_Rate'] = df2['FN'] /df2['Total_Vars']

# Attempt to make output dir
try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

df2.to_pickle(outdir / "sv_var_results.pkl")

# Plot Totals
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(18,12))

sns.barplot(x=df2.index, y=df2['Total_Vars'], data = df2, ax = ax3)
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=40, ha="right", fontsize='medium', fontweight='normal')

# Plot the FN rate
sns.barplot(x=df2.index, y=df2["FN_Rate"], data=df2, ax = ax2)

# Plot the FP rate
sns.barplot(x=df2.index, y=df2["FP_Rate"], data=df2, ax = ax1)

plt.tight_layout()
plt.savefig(outdir / "results.png")
