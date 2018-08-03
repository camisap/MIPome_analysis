
# coding: utf-8


import numpy as np
import pandas as pd
import sklearn
import seaborn as sns
import matplotlib.pyplot as plt

# Values that I am keeping track of with the code
#     peptide_count_raw ==> number of peptides in raw file. No filtering
#     peptide_count_iso_filtered   ==> number of peptides that that were above the iso_spec_min value
#     peptide_count_iso_intensity_filtered ==> number of peptides that were above the adjusted_intensity_min value
#     peptide_count_iso_intensity_duplicated ==> number of peptides that resulted after duplicated peptides were gone

channel1 = 'Control_A'
channel2 = 'Control_B'
channel3 = 'pH_A'
channel4 = 'pH_B'
channel5 = 'NaCl_A'
channel6 = 'NaCl_B'

iso_spec_min = 0.80
adjusted_intensity_min = 1000
genes_for_norm = ('aup1', 'faf2', 'pnpla2', 'plin2')

# #Note: the adjusted intensities and isolation Specificity columns should be changed from Scientific notation to 
# #numbers with 4 decimal points. then, when saving as an MS-DOS CSV, those columns will be changed to general numbers. 

# User INPUT:
# your .csv file. Mine had the following columns in the following order:
# 'ScanF', 'ScanF Link', 'MS2 ID', 'z', 'PPM', 'XCorr', '&#916;Corr', '&#916;Corr Link', 'Rank/SP', '# Ions', '# Ions Link', 
# 'Reference', 'Reference Link', 'Redun', 'Redun Link', 'Peptide', 'Peptide Link', 'Gene Symbol', 'Annotation', '126 Adjusted Intensity', 
# '127 Adjusted Intensity', '128 Adjusted Intensity', '129 Adjusted Intensity', '130 Adjusted Intensity', '131 Adjusted Intensity',
# 'Isolation Specificity'

filename = 

#this will preface all files/figures that will be saved below
fileToSave = 

# this is a way to keep track of the experiment you're analyzing
experiment_number = 

#clarify what the title of your columns are
col_names = {'Peptide':'peptide', 'Gene Symbol': 'gene_symbol', '126 Adjusted Intensity':channel1, '127 Adjusted Intensity':channel2, 
                       '128 Adjusted Intensity':channel3, '129 Adjusted Intensity':channel4, '130 Adjusted Intensity': channel5, 
                       '131 Adjusted Intensity':channel6, 'Isolation Specificity':'isospec'}

df_raw = pd.read_csv(filename)
FileToSave = fileToSave
Experiment = experiment_number

peptide_count_raw = df_raw.shape

df_analysis = df_raw[[15,17,19,20,21,22,23,24,25]]

df_analysis = df_analysis.rename(columns=col_names)

def parse_Peptide (peptide):
    core_peptide = str(peptide).split('.')
    return core_peptide[1]

def look_for_unmodified_K (peptide):
    for pos in range(len(peptide)):
        if peptide[pos] == 'K':
            if pos==((len(peptide))-1): return False 
            if peptide[(pos+1)] == "#": continue 
            else: return False 
    return True

def termini_modified (peptide):
    if "]" in peptide: return True
    else: return False


df_analysis['All_K_modified']= parse_Peptide(df_analysis["peptide"])
df_analysis['Modified_termini'] = parse_Peptide(df_analysis["peptide"])

for each in range(len(df_analysis)):
    df_analysis.ix[each, 9]= look_for_unmodified_K(parse_Peptide(df_analysis.ix[each, 0]))
    df_analysis.ix[each, 10]= termini_modified(parse_Peptide(df_analysis.ix[each, 0]))

df_analysis_Kmod_peptides = df_analysis[df_analysis["All_K_modified"]==True]
df_analysis_labeled_peptides = df_analysis_Kmod_peptides[df_analysis_Kmod_peptides['Modified_termini']==True]


peptide_count_KModified = df_analysis_Kmod_peptides.shape
peptide_count_labeled = df_analysis_labeled_peptides.shape


col_list= list(df_analysis[[channel1,channel2,channel3,channel4,channel5,channel6]])


df_analysis_labeled_peptides['Raw_intensities_sum'] = df_analysis_labeled_peptides[col_list].sum(axis=1)


df_analysis_labeled_peptides.isnull().sum(axis=0)


df_analysis_labeled_peptides.head()


# I am excluding all IsoSpecificities that are below 0.80. 'iso_spec_min' is set above

df_iso_filtered = df_analysis_labeled_peptides[df_analysis_labeled_peptides["isospec"]>=iso_spec_min]


# Note that I am keeping track of the number of peptides that are being filtered with each step. 

peptide_count_iso_filtered = df_iso_filtered.shape
peptide_count_iso_filtered

hist_isoSpec_filtered = df_iso_filtered["isospec"].hist()
fig1 = hist_isoSpec_filtered.get_figure()
plt.title("isoSpec_filtered")
# fig1.savefig(FileToSave+"_hist1_isoSpec_isoSpecFiltered.eps")

df_iso_intensity_filtered = df_iso_filtered[df_iso_filtered["Raw_intensities_sum"]>adjusted_intensity_min]

peptide_count_iso_intensity_filtered = df_iso_intensity_filtered.shape

hist_iso_intensity_filtered = df_iso_intensity_filtered["isospec"].hist()
fig2 = hist_iso_intensity_filtered.get_figure()
plt.title("isoSpec_intensity_filtered")
#fig2.savefig(FileToSave+'_hist2_isoSpec_intensityFiltered.eps')

peptide_count_iso_intensity_filtered

df_iso_intensity_filtered_sorted=df_iso_intensity_filtered.sort(['peptide','Raw_intensities_sum'], 
                                                                ascending=[True, False])

df_iso_intensity_filtered_sorted.head()

df_filtered = df_iso_intensity_filtered_sorted.drop_duplicates(subset='peptide')

peptide_count_iso_intensity_duplicated = df_filtered.shape

hist_filtered_dupDrpd = df_filtered['isospec'].hist()
fig3 = hist_filtered_dupDrpd.get_figure()
plt.title("DuplicatesDropped_filtered")
# fig3.savefig(FileToSave+'_hist3_isoSpec_duplicatesDropped.eps')

df_filtered.head()

peptide_count_iso_intensity_duplicated


def makeNewColumnLabel(suffix, columnName):
    name = columnName + suffix
    return name

df_relativeIntensities = df_filtered[['peptide','gene_symbol', channel1, channel2, channel3, channel4, 
                        channel5, channel6,'Raw_intensities_sum']]

for x in df_filtered[[channel1, channel2, channel3, channel4, channel5, channel6]]:
    newName1 = makeNewColumnLabel('_Raw_relative_int', str(x))
    df_relativeIntensities[newName1]= df_filtered[x]/df_filtered['Raw_intensities_sum']

df_Raw_proteins = df_relativeIntensities.groupby(by='gene_symbol').mean()

df_Raw_proteins.to_csv(path_or_buf=FileToSave+'_Raw_Proteins.csv', sep=',', na_rep='', 
                         float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', 
                         encoding=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, 
                         tupleize_cols=False, date_format=None, doublequote=True, escapechar=None, decimal='.')

# fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(15,5))

# axes = axes.ravel()

# for idx, col_name in enumerate(df_relativeIntensities.ix[:,9:]):
#     axes[idx].hist(df_relativeIntensities[col_name], bins=20)
#     axes[idx].set_title(col_name)

# fig.savefig(FileToSave+'_preNormalized_relativeIntensities_Histograms.eps')

df_rawRelInt = df_relativeIntensities[[0, 1, 9, 10, 11, 12, 13, 14]]  
df_rawRelInt.to_csv(path_or_buf=FileToSave+'_Raw_filtered_peptides.csv', sep=',', na_rep='', 
                         float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', 
                         encoding=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, 
                         tupleize_cols=False, date_format=None, doublequote=True, escapechar=None, decimal='.')

genes = []
for x in genes_for_norm:
    genes.append(x.upper())

df_normalizing_peptides = df_rawRelInt.loc[df_rawRelInt['gene_symbol'].isin(genes)]

peptide_count_normalizing = df_normalizing_peptides.shape

df_normalizing_peptides.to_csv(FileToSave+'_normalizing_peptides.csv')

df_norm_factors = df_normalizing_peptides.mean(axis=0, skipna=True, numeric_only=True)

df_norm_stdDev = df_normalizing_peptides.std(axis=0,skipna=True, numeric_only=True)

bar_normalizing_factors = df_norm_factors.plot(kind='bar', ax=None, figsize=None, use_index=True, title=None, 
                     grid=None, legend=False, style=None, logx=False, logy=False, 
                     loglog=False, xticks=None, yticks=None, xlim=None, ylim=None, 
                     rot=None, fontsize=None, colormap=None, table=False, yerr=df_norm_stdDev, 
                     xerr=None, label=None, secondary_y=False)

bar_normalizing_factors.set_ylabel('Relative_Intensities')

fig4 = bar_normalizing_factors.get_figure()
plt.title("Normalizing_Peptides")
# fig4.savefig(FileToSave+'_bar_normalizing_peptides.eps')

df_recovery_normalized = df_filtered[['peptide','gene_symbol', channel1, channel2, channel3, channel4, 
                        channel5, channel6]]

counter=0
# this worked! :  df_recovery_normalized[channel1].div(df_norm_factors[0])

for x in df_filtered[[channel1, channel2, channel3, channel4, channel5, channel6]]:
    newName1 = makeNewColumnLabel('_Recov_Norm', str(x))
    df_recovery_normalized[newName1]= df_filtered[x].div(df_norm_factors[counter])
    counter+=1

Norm_col_list= list(df_recovery_normalized[[8,9,10,11,12,13]])


df_recovery_normalized['Norm_intensities_sum'] = df_recovery_normalized[Norm_col_list].sum(axis=1)

df_Norm_relativeIntensities = df_recovery_normalized[[0,1,8,9,10,11,12,13,14]]

for x in df_Norm_relativeIntensities[[2, 3, 4, 5, 6, 7]]:
    newName1 = makeNewColumnLabel('_relative_int', str(x))
    df_Norm_relativeIntensities[newName1]= df_Norm_relativeIntensities[x]/df_Norm_relativeIntensities['Norm_intensities_sum']

df_Norm_RelInt = df_Norm_relativeIntensities[[0, 1, 9, 10, 11, 12, 13, 14]]

df_Norm_RelInt.to_csv(path_or_buf=FileToSave+'_Normalized_peptide_relInt.csv', sep=',', na_rep='', 
                         float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', 
                         encoding=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, 
                         tupleize_cols=False, date_format=None, doublequote=True, escapechar=None, decimal='.')

# fig5, axes5 = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(15,5))

# axes5 = axes5.ravel()

# for idx, col_name in enumerate(df_Norm_relativeIntensities.ix[:,9:15]):
#     axes5[idx].hist(df_Norm_relativeIntensities[col_name], bins=20)
#     axes5[idx].set_title(col_name)

# fig5.savefig(FileToSave+'_postNormalized_relativeIntensities_Histograms.eps')

df_Norm_proteins = df_Norm_RelInt.groupby(by='gene_symbol').mean()

protein_count = df_Norm_proteins.shape

df_Norm_proteins['StDev'] = df_Norm_proteins[[0,1,2,3,4,5]].std(axis=1)

df_Norm_proteins.to_csv(path_or_buf=FileToSave+'_Normalized_Proteins.csv', sep=',', na_rep='', 
                         float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', 
                         encoding=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, 
                         tupleize_cols=False, date_format=None, doublequote=True, escapechar=None, decimal='.')

df_Norm_prot_3values = df_Norm_proteins

df_Norm_prot_3values['Control'] = df_Norm_proteins["Control_A_Recov_Norm_relative_int"]+df_Norm_proteins['Control_B_Recov_Norm_relative_int']
df_Norm_prot_3values['Alkaline_Extraction'] = df_Norm_proteins['pH_A_Recov_Norm_relative_int'] + df_Norm_proteins['pH_B_Recov_Norm_relative_int']
df_Norm_prot_3values['High Salt Concentration'] = df_Norm_proteins['NaCl_A_Recov_Norm_relative_int'] +df_Norm_proteins['NaCl_B_Recov_Norm_relative_int']

df_Norm_prot_Final_3values = df_Norm_prot_3values[[7,8,9]]
df_Norm_proteins.to_csv(path_or_buf=FileToSave+'_Proteins_Final3Values.csv', sep=',', na_rep='', 
                         float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', 
                         encoding=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, 
                         tupleize_cols=False, date_format=None, doublequote=True, escapechar=None, decimal='.')

textToDocument = open(FileToSave+'_Filtering_Specs.csv', 'w')

textToDocument.write("File/Experiment: " + str(Experiment)+'\n')
textToDocument.write("Genes Used to Normalize: " + str(genes_for_norm)+'\n')
textToDocument.write("Isolation Specifity minimum: " + str(iso_spec_min)+'\n')
textToDocument.write("Minimum Summed Intensity: " + str(adjusted_intensity_min)+'\n')
textToDocument.write('\n')
textToDocument.write("Raw peptide count: " + str(peptide_count_raw)+'\n')
textToDocument.write("K modified peptides: " + str(peptide_count_KModified)+'\n')
textToDocument.write("Labeled peptides: " + str(peptide_count_labeled)+'\n')
textToDocument.write("Isolation Specificity filtered count: " + str(peptide_count_iso_filtered)+'\n')
textToDocument.write("Isolation Specificity AND Intensity filtered: " + str(peptide_count_iso_intensity_filtered)+'\n')
textToDocument.write("Duplicates dropped filtered: " + str(peptide_count_iso_intensity_duplicated)+'\n')
textToDocument.write("Number of peptides used to normalize: " + str(peptide_count_normalizing)+'\n')
textToDocument.write("Number of proteins in list after AVERAGING peptides: " + str(protein_count)+'\n')


textToDocument.close()





