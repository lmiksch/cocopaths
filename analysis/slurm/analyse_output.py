import matplotlib.pyplot as plt
import pandas as pd
import os
import ast 
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

def sort_and_print_missing_tsv(file_path, folder_num):
    df = pd.read_csv(file_path, sep='\t')
    df = df[pd.to_numeric(df['Id'], errors='coerce').notnull()]
    df['Id'] = df['Id'].astype(int)

    # Deduplicate by keeping the row with higher 'tries' for each 'acfp'
    df = df.sort_values(by='tries', ascending=False)
    df = df.drop_duplicates(subset='acfp', keep='first')

    print(df.head())

    # Check for missing IDs
    all_ids = set(range(1, max(df['Id']) + 1))
    present_ids = set(df['Id'])
    missing_ids = sorted(list(all_ids - present_ids))

    print(f'{missing_ids = }')

    # Sort the dataframe by 'Id' in ascending order
    sorted_df = df.sort_values(by='Id')

    sorted_df['AFP_list'] = sorted_df['acfp'].apply(lambda x: ast.literal_eval(x))
    sorted_df['refolds'] = sorted_df['AFP_list'].apply(calc_refolds)

    # Save the sorted dataframe to a new TSV file
    sorted_df.to_csv('output_' + folder_num + '/sorted_result.tsv', sep='\t', index=False)



def read_process_tsv(file_paths):
    # Load and concatenate multiple TSV files
    df_list = []
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep='\t')
        df_list.append(df)
    
    # Concatenate all dataframes into a single dataframe
    df = pd.concat(df_list, ignore_index=True)
    
    # Filter out rows where Id starts with a '-'
    df = df[~df['Id'].astype(str).str.startswith('-')]

    # Convert 'nt_success' to numeric, coercing errors to NaN
    df['nt_success'] = pd.to_numeric(df['nt_success'], errors='coerce')

    # Drop rows with NaN values in 'nt_success'
    df = df.dropna(subset=['nt_success'])

    # Determine success at the nucleotide level
    df['nucleotide_success'] = df['nt_success'] > 0.5
    
    print("Total nucleotide successes:", df['nucleotide_success'].sum())
    print("Total nucleotide failures:", (~df['nucleotide_success']).sum())
    print("Total domain successes:", df['domain_success'].sum())
    print("Total domain failures:", (~df['domain_success']).sum())

    print(df.head())
    return df

def plot_out_paper(out_folder,df):

    success_by_domain = pd.crosstab(df['domain_success'], df['nucleotide_success'])

    # Calculate the counts of successes per try
    counts_success_per_try = df[df['nucleotide_success']].groupby('tries').size()

    # Create a figure with two subplots
    fig, axs = plt.subplots(2, 1, figsize=(8.27, 5.85))  # Half of A4 height in inches
    # First plot: Domain Success vs. Nucleotide Success
    success_by_domain.plot(kind='bar', ax=axs[0], color=['lightcoral', 'cornflowerblue'])
    axs[0].set_title('Domain Success vs. Nucleotide Success')
    axs[0].set_xlabel('Domain Success')
    axs[0].set_ylabel('Count')
    axs[0].legend(['Nucleotide Failure', 'Nucleotide Success'])
    axs[0].grid(True)


    axs[0].text(-0.1, 1.1, 'a', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Second plot: Counts of Successes per Try
    counts_success_per_try.plot(kind='bar', ax=axs[1], color='cornflowerblue', edgecolor='black')
    axs[1].set_title('Counts of Successes per Try')
    axs[1].set_xlabel('Number of Tries')
    axs[1].set_ylabel('Count of Successes')
    axs[1].grid(True)

    axs[1].text(-0.1, 1.1, 'b', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the combined figure
    plt.savefig(os.path.join("/home/mescalin/miksch/Documents/figures", 'combined_plot.png'))
    plt.close()


def plot_split_normalized_domain(out_folder, df):

    # Ensure the output directory exists
    os.makedirs(out_folder, exist_ok=True)

    # Compute counts: Crosstab of domain_success vs. nucleotide_success
    success_by_domain = pd.crosstab(df['domain_success'], df['nucleotide_success'])

    # Normalize so the total sum of all values equals 1
    total_sum = success_by_domain.values.sum()
    normalized_counts = success_by_domain / total_sum  # Normalize across all values

    # Rearrange data for grouped bar plotting
    normalized_counts = normalized_counts.reset_index().melt(
        id_vars='domain_success',
        var_name='Nucleotide Outcome',
        value_name='Proportion'
    )

    # Replace binary Nucleotide Outcome with labels
    normalized_counts['Nucleotide Outcome'] = normalized_counts['Nucleotide Outcome'].map({False: 'Failure', True: 'Success'})

    # Calculate medians for Domain Success and Failure
    success_median = df[df['domain_success'] == True]['tries'].median()
    failure_median = df[df['domain_success'] == False]['tries'].median()

    # Print medians to terminal
    print(f"Median number of tries for Domain Success = True: {success_median}")
    print(f"Median number of tries for Domain Success = False: {failure_median}")

    # Plot: Split bars for Domain Success and Nucleotide Outcomes
    fig, ax_main = plt.subplots(figsize=(10, 6), dpi=300)

    # Define bar width and positions
    bar_width = 0.3
    domain_levels = normalized_counts['domain_success'].unique()
    x_positions = np.arange(len(domain_levels))

    # Plot the bars for Nucleotide Failure and Success
    for i, outcome in enumerate(['Failure', 'Success']):
        color = 'lightcoral' if outcome == 'Failure' else 'cornflowerblue'
        subset = normalized_counts[normalized_counts['Nucleotide Outcome'] == outcome]
        plt.bar(
            x_positions + i * bar_width,  # Shift bar positions
            subset['Proportion'],
            width=bar_width,
            color=color,
            edgecolor='black',
            label=f'Nucleotide {outcome}'
        )

    # Label and formatting
    ax_main.set_title('Domain Success vs. Nucleotide Success/Failure (Normalized)', fontsize=14)
    ax_main.set_xlabel('Domain Success', fontsize=12)
    ax_main.set_ylabel('Proportion', fontsize=12)
    ax_main.set_xticks(x_positions + bar_width / 2)
    ax_main.set_xticklabels(['Failure', 'Success'], fontsize=10)
    ax_main.legend(loc='best', fontsize=10)
    ax_main.grid(axis='y', linestyle='--', alpha=0.5)

    # Add an enlarged inset plot
    ax_inset = inset_axes(ax_main, width="80%", height="80%", loc='upper right', borderpad=1.5)  # Larger inset
    inset_data = success_by_domain / success_by_domain.sum().sum()  # Normalize for inset
    inset_data.plot(kind='bar', stacked=True, color=['lightcoral', 'cornflowerblue'], edgecolor='black', ax=ax_inset)
    ax_inset.set_title('Inset: Proportions', fontsize=10)
    ax_inset.set_xlabel('')
    ax_inset.set_ylabel('')
    ax_inset.tick_params(axis='both', labelsize=8)

    # Save the plot
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '04_split_normalized_domain_vs_nucleotide_success_with_inset.png'))
    plt.close()

    
def plot_domain_and_counts(out_folder, df):
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # Ensure the output directory exists
    os.makedirs(out_folder, exist_ok=True)

    # Plot 1: Domain Success vs. Nucleotide Success
    success_by_domain = pd.crosstab(df['domain_success'], df['nucleotide_success'])
    plt.figure(figsize=(10, 6), dpi=300)
    success_by_domain.plot(kind='bar', color=['lightcoral', 'cornflowerblue'], edgecolor='black')
    plt.title('Domain Success vs. nucleotide success', fontsize=16)
    plt.xlabel('Domain success', fontsize=18)
    plt.ylabel('Count', fontsize=18)
    plt.legend(['Nucleotide Failure', 'Nucleotide Success'], loc='best', fontsize=20)
    plt.grid(True)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(0, 1)  # Adjust y-axis to always go from 0 to 1 for comparison
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '04_domain_vs_nucleotide_success.png'))
    plt.close()

    # Compute counts: Crosstab of domain_success vs. nucleotide_success
    total_sum = success_by_domain.values.sum()
    normalized_counts = success_by_domain / total_sum  # Normalize across all values

    # Rearrange data for grouped bar plotting
    normalized_counts = normalized_counts.reset_index().melt(
        id_vars='domain_success',
        var_name='Nucleotide Outcome',
        value_name='Proportion'
    )

    # Replace binary Nucleotide Outcome with labels
    normalized_counts['Nucleotide Outcome'] = normalized_counts['Nucleotide Outcome'].map({False: 'Failure', True: 'Success'})

    # Plot: Split bars for Domain Success and Nucleotide Outcomes
    plt.figure(figsize=(10, 6), dpi=300)
    bar_width = 0.3
    domain_levels = sorted(normalized_counts['domain_success'].unique())
    x_positions = np.arange(1, len(domain_levels) + 1)  # Start at 1

    for i, outcome in enumerate(['Failure', 'Success']):
        color = 'lightcoral' if outcome == 'Failure' else 'cornflowerblue'
        subset = normalized_counts[normalized_counts['Nucleotide Outcome'] == outcome]
        plt.bar(
            x_positions + i * bar_width - bar_width / 2,  # Adjust bar positions
            subset['Proportion'],
            width=bar_width,
            color=color,
            edgecolor='black',
            label=f'Nucleotide {outcome}'
        )

    plt.title('Domain success vs. nucleotide success/failure (normalized)', fontsize=18)
    plt.xlabel('Domain success', fontsize=20)
    plt.ylabel('Proportion', fontsize=20)
    plt.xticks(x_positions, ['Failure', 'Success'], fontsize=20)
    plt.ylim(0, 1)
    plt.yticks(fontsize=16)
    plt.legend(loc='best', fontsize=16)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '04_split_normalized_domain_vs_nucleotide_success.png'))
    plt.close()

    # Compute normalized counts from the full crosstab
    total_sum = success_by_domain.values.sum()
    normalized_counts = success_by_domain / total_sum

    # Rearrange data for grouped bar plotting
    normalized_counts = normalized_counts.reset_index().melt(
        id_vars='domain_success',
        var_name='Nucleotide Outcome',
        value_name='Proportion'
    )
    # Convert binary outcome to labels
    normalized_counts['Nucleotide Outcome'] = normalized_counts['Nucleotide Outcome'].map({False: 'Failure', True: 'Success'})

    # --- New: Plot for Domain Success Only (Split Normalized) ---
    # Filter to only include rows where domain_success is True
    ds_norm = normalized_counts[normalized_counts['domain_success'] == True].copy()

    # Order the outcomes so that 'Failure' comes first, then 'Success'
    order = ['Failure', 'Success']
    ds_norm['Nucleotide Outcome'] = pd.Categorical(ds_norm['Nucleotide Outcome'], categories=order, ordered=True)
    ds_norm = ds_norm.sort_values('Nucleotide Outcome')

    plt.figure(figsize=(10, 6), dpi=300)
    x_positions = np.arange(len(ds_norm))  # This will be 0 and 1 for the two outcomes
    bar_width = 0.5
    # Define colors based on the outcome
    colors = ['lightcoral' if outcome == 'Failure' else 'cornflowerblue' for outcome in ds_norm['Nucleotide Outcome']]

    # Plot the bars for each nucleotide outcome
    plt.bar(x_positions, ds_norm['Proportion'], color=colors, edgecolor='black', width=bar_width)

    # New title and labels for the figure
    plt.title('Normalized Nucleotide Outcomes for Domain Success Only', fontsize=18)
    plt.xlabel('Nucleotide Outcome', fontsize=20)
    plt.ylabel('Proportion', fontsize=20)
    plt.xticks(x_positions, ds_norm['Nucleotide Outcome'], fontsize=20)
    plt.ylim(0, 1)
    plt.yticks(fontsize=16)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '04_split_normalized_domain_success_only.png'))
    plt.close()


    # Filter success and failure data frames for the per try plots
    success_df = df[(df['nucleotide_success']) & (df['domain_success'] == True)]
    failure_df = df[(df['nucleotide_success']) & (df['domain_success'] == False)]

    # Calculate and print medians
    success_median = success_df['tries'].median() if not success_df.empty else None
    failure_median = failure_df['tries'].median() if not failure_df.empty else None
    print(f"Median number of tries for Domain Success = True: {success_median}")
    print(f"Median number of tries for Domain Success = False: {failure_median}")

    # Group data for plotting (per number of tries)
    success_counts = success_df.groupby('tries').size()
    failure_counts = failure_df.groupby('tries').size()

    # Normalize counts for proportion calculation
    success_counts = success_counts / success_counts.sum()
    failure_counts = failure_counts / failure_counts.sum()

    # Pad missing tries to get full range 0–19
    full_try_range = list(range(0, 20))
    success_counts = success_counts.reindex(full_try_range, fill_value=0)
    failure_counts = failure_counts.reindex(full_try_range, fill_value=0)

    # Main Plot: Normalized nucleotide success per try (domain success)
    fig, ax_main = plt.subplots(figsize=(10, 4), dpi=300)
    success_counts.plot(kind='bar', color='#55a868', edgecolor='black', ax=ax_main)
    ax_main.set_facecolor('white')
    ax_main.set_title('Normalized nucleotide success per try (domain success)', fontsize=16)
    ax_main.set_xlabel('Number of tries', fontsize=18)
    ax_main.set_ylabel('Proportion of successes', fontsize=18)
    ax_main.grid(False)
    ax_main.tick_params(axis='both', labelsize=14)
    ax_main.set_ylim(0, 1)
    ax_main.set_xticks(np.arange(len(success_counts)))
    ax_main.set_xticklabels(success_counts.index, rotation=0)

    # Note: Inset section removed. You can create a separate plot for domain failure if needed.
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '04_normalized_counts_success_per_try.png'))
    plt.close()

    # --- New Section: Plotting Crosstab Table in Matplotlib ---
    # Ensure the nt_success_bool column exists (successful if value > 0.5)
    if 'nt_success_bool' not in df.columns:
        df['nt_success_bool'] = df['nucleotide_success'] > 0.5

    # Compute the crosstab with margins to include row/column totals
    crosstab_table = pd.crosstab(df['domain_success'], df['nt_success_bool'], margins=True)

    # Convert table values to a list of lists for matplotlib.table
    data_table = crosstab_table.values.tolist()
    row_labels = list(crosstab_table.index.astype(str))
    col_labels = list(crosstab_table.columns.astype(str))

    # Create a new figure for the table
    plt.figure(figsize=(8, 4))
    ax_table = plt.gca()
    ax_table.axis('tight')
    ax_table.axis('off')
    table = ax_table.table(cellText=data_table, rowLabels=row_labels, colLabels=col_labels, loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.2)
    plt.title("Crosstab: Domain Success vs. Nucleotide Success", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '05_crosstab_domain_vs_nucleotide_table.png'))
    plt.close()



def plot_output(out_folder, df):
    # Plot 1: Distribution of nt_success
    plt.figure(figsize=(10, 6))
    plt.hist(df['nt_success'], bins=20, color='cornflowerblue', edgecolor='black')
    plt.title('Distribution of nt_success')
    plt.xlabel('nt_success')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.savefig(out_folder + 'dist_nt_success.png')
    plt.close()
    
    # Plot 2: Count of nt_success
    plt.figure(figsize=(10, 6))
    df['nucleotide_success'].value_counts().plot(kind='bar', color=['salmon', 'cornflowerblue'])
    plt.title('Count of Nucleotide Success')
    plt.xlabel('Nucleotide Success')
    plt.ylabel('Count')
    plt.grid(True)
    plt.savefig(out_folder + 'count_nt_success.png')
    plt.close()

    # Plot 3 but stacking
    plt.figure(figsize=(10, 6))
    # Plot stacked histogram for domain success and domain failure
    plt.hist(
        [df[df['domain_success'] == True]['nt_success'], df[df['domain_success'] == False]['nt_success']],
        bins=20, 
        color=['cornflowerblue', 'lightcoral'], 
        edgecolor='black', 
        label=['Domain Success', 'Domain Failure'], 
        alpha=0.5,
        stacked=True
    )
    plt.title('Distribution of nt_success by Domain Success')
    plt.xlabel('nt_success')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)
    plt.savefig(out_folder + 'dist_nt_domain_success.png')
    plt.close()
    
    grouped = df.groupby('tries')['nt_success'].mean()

    # Line plot of average nt_success by tries
    plt.figure(figsize=(10, 6))
    grouped.plot(kind='line', marker='o', color='cornflowerblue')
    plt.title('Average nt_success by tries')
    plt.xlabel('tries')
    plt.ylabel('Average nt_success')
    plt.grid(True)
    plt.savefig(out_folder + 'tries_nt_sucess.png')
    plt.close()
    
    # Crosstab of domain success vs nucleotide success
    success_by_domain = pd.crosstab(df['domain_success'], df['nucleotide_success'])

    # Bar plot of domain success vs. nucleotide success
    plt.figure(figsize=(10, 6), dpi=300)  # Set size and DPI explicitly
    success_by_domain.plot(kind='bar', color=['lightcoral', 'cornflowerblue'], edgecolor='black')
    plt.title('Domain Success vs. Nucleotide Success')
    plt.xlabel('Domain Success')
    plt.ylabel('Count')
    plt.legend(['Nucleotide Failure', 'Nucleotide Success'], loc='best')  # Consistent legend placement
    plt.grid(True)
    plt.tight_layout()  # Ensure consistent layout
    plt.savefig(os.path.join(out_folder, 'domain_vs_nucleotide_success.png'), dpi=300)
    plt.close()

    # Counts of nucleotide successes per try, grouped by domain success
    counts_success_per_try_by_domain = df[df['nucleotide_success']].groupby(['tries', 'domain_success']).size().unstack(fill_value=0)

    # Plot the counts of successes per try, distinguishing between domain success and failure
    
    # Counts of nucleotide successes per try, grouped by domain success
    counts_success_per_try_by_domain = df[df['nucleotide_success']].groupby(['tries', 'domain_success']).size().unstack(fill_value=0)

    # Plot the counts of successes per try, distinguishing between domain success and failure
    plt.figure(figsize=(10, 6), dpi=300)  # Explicitly set size and DPI
    counts_success_per_try_by_domain.plot(kind='bar', stacked=True, color=['lightcoral', 'cornflowerblue'], edgecolor='black')
    plt.title('Counts of Successes per Try by Domain Success')
    plt.xlabel('Number of Tries')
    plt.ylabel('Count of Successes')
    plt.legend(['Domain Failure', 'Domain Success'], loc='best')  # Consistent legend placement
    plt.grid(True)
    plt.tight_layout()  # Ensure consistent layout
    plt.savefig(os.path.join(out_folder, 'counts_success_per_try_by_domain.png'), dpi=300)
    plt.close()

    counts_success_per_try = df[df['nucleotide_success']].groupby('tries').size()

    # Plot the counts of successes per try
    plt.figure(figsize=(10, 6), dpi=300)
    counts_success_per_try.plot(kind='bar', color='cornflowerblue', edgecolor='black')
    plt.title('Counts of Successes per Try')
    plt.xlabel('Number of Tries')
    plt.ylabel('Count of Successes')
    plt.grid(True)
    plt.savefig(os.path.join(out_folder, 'counts_success_per_try.png'))
    plt.close()

    plt.figure(figsize=(10, 6), dpi=300)
    success_by_domain_only.plot(kind='bar', color=['lightcoral', 'cornflowerblue'], edgecolor='black')
    plt.title('Only Domain Success: Nucleotide Success vs. Failure', fontsize=16)
    plt.xlabel('Domain Success', fontsize=18)
    plt.ylabel('Count', fontsize=18)
    plt.legend(['Nucleotide Failure', 'Nucleotide Success'], loc='best', fontsize=20)
    plt.grid(True)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '04_domain_success_only.png'))
    plt.close()

    custom_palette = {False: 'lightcoral', True: 'cornflowerblue'}

    plt.figure(figsize=(10, 6))
    sns.histplot(df, x='refolds', hue='nucleotide_success', multiple='stack', palette=custom_palette, edgecolor='black', discrete=True)
    plt.title('Count of Nucleotide Success Compared to Refolds')
    plt.xlabel('Refolds')
    plt.ylabel('Count')
    plt.grid(True)

    # Ensure the output directory exists
    os.makedirs(out_folder, exist_ok=True)

    plt.savefig(os.path.join(out_folder, 'count_nt_success_vs_refolds.png'))
    plt.close()

    # Filter rows where nucleotide_success is True
    df_success = df[df['nucleotide_success'] == True]

    # Plot using seaborn catplot
    plt.figure(figsize=(12, 8))
    sns.catplot(data=df_success, x='tries', kind='count', col='refolds', col_wrap=3, palette='colorblind', edgecolor='black')
    plt.suptitle('Counts of Nucleotide Success per Tries and Refolds', y=1.02)
    plt.xlabel('Number of Tries')
    plt.ylabel('Count of Nucleotide Successes')
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, 'counts_success_per_try_and_refolds_combined.png'))

    # Aggregate counts of nucleotide_success and domain_success per refolds
    nt_success_counts = df.groupby(['refolds', 'nucleotide_success']).size().unstack(fill_value=0)
    domain_success_counts = df.groupby(['refolds', 'domain_success']).size().unstack(fill_value=0)

    # Print the combined counts per refolds for nucleotide_success
    print("Nucleotide Success Counts per Refolds:")
    print(nt_success_counts)

    # Print the combined counts per refolds for domain_success
    print("\nDomain Success Counts per Refolds:")
    print(domain_success_counts)

    # Plotting histogram for domain_success
    plt.figure(figsize=(10, 6))
    sns.histplot(df, x='refolds', hue='domain_success', multiple='stack', palette=custom_palette, edgecolor='black', discrete=True)
    plt.title('Count of Domain Success Compared to Refolds')
    plt.xlabel('Refolds')
    plt.ylabel('Count')
    plt.grid(True)
    plt.savefig(os.path.join(out_folder, 'count_domain_success_vs_refolds.png'))


    filtered_df = df[df['refolds'].isin([0, 1, 2, 3, 4])]
    # Create boxplots for nucleotide_success based on refolds values
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=filtered_df, x='refolds', y='nucleotide_success', palette= "colorblind")
    plt.title('Boxplots of Nucleotide Success Based on Refolds (0-4)')
    plt.xlabel('Refolds')
    plt.ylabel('Nucleotide Success')
    plt.grid(True)

    # Ensure the output directory exists
    os.makedirs(out_folder, exist_ok=True)

    plt.savefig(os.path.join(out_folder, 'boxplot_nt_success_refolds_0_4.png'))
    plt.close()


    df_success = df[df['nucleotide_success'] == True]

    # Create a stacked histogram of successful nucleotide_success attempts by refolds
    plt.figure(figsize=(12, 8))

    # Plot stacked histogram
    plt.hist(
        [df_success[df_success['refolds'] == 1]['tries'],
        df_success[df_success['refolds'] == 2]['tries'],
        df_success[df_success['refolds'] == 3]['tries'],
        df_success[df_success['refolds'] == 4]['tries']],
        bins=range(0, max(df['tries']) + 1),  # Use all tries as bins
        color=['#4c72b0', '#55a868', '#c44e52', '#8172b3'],
        edgecolor='black',
        label=['Refolds 1', 'Refolds 2', 'Refolds 3', 'Refolds 4'],
        alpha=0.7,
        stacked=True
    )

    # Additional refinements
    plt.title('Distribution of Successful Nucleotide Success by Tries and Refolds')
    plt.xlabel('Number of Tries')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)

    plt.savefig(os.path.join(out_folder, 'dist_successful_nt_success_per_try_and_refolds_stacked.png'))
    

        # --- New Section: Plotting Crosstab Table in Matplotlib ---
    # Ensure the nt_success_bool column exists (successful if value > 0.5)
    if 'nt_success_bool' not in df.columns:
        df['nt_success_bool'] = df['nucleotide_success'] > 0.5

    # Compute the crosstab with margins to include row/column totals
    crosstab_table = pd.crosstab(df['domain_success'], df['nt_success_bool'], margins=True)

    # Convert table values to a list of lists for matplotlib.table
    data_table = crosstab_table.values.tolist()
    row_labels = list(crosstab_table.index.astype(str))
    col_labels = list(crosstab_table.columns.astype(str))

    # Create a new figure for the table
    plt.figure(figsize=(8, 4))
    ax_table = plt.gca()
    ax_table.axis('tight')
    ax_table.axis('off')
    table = ax_table.table(cellText=data_table, rowLabels=row_labels, colLabels=col_labels, loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.2)
    plt.title("Crosstab: Domain Success vs. Nucleotide Success", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, '05_crosstab_domain_vs_nucleotide_table.png'))
    plt.close()




def compare_data(folder_1,folder_2):


    df1 = read_process_tsv([folder_1 + "/sorted_result.tsv"])
    df2 = read_process_tsv([folder_2 + "/sorted_result.tsv"])

    
    column_name = "nucleotide_success"

    # Add an index column for merging
    df1['index'] = df1.index
    df2['index'] = df2.index

    # Merge the DataFrames on the index column
    merged_df = df1.merge(df2, on='index', suffixes=('_file1', '_file2'))

    # Filter rows where the specified column differs
    differing_rows = merged_df[merged_df[f'{column_name}_file1'] != merged_df[f'{column_name}_file2']]


    relevant_columns = ['AFP_file1','nt_success_file1', 'nt_success_file2', 'nucleotide_success_file1', 'nucleotide_success_file2',"tries_file1","tries_file2" ]
    # Print the rows with differences
    print(differing_rows[relevant_columns])

    differing_rows = differing_rows[relevant_columns]

    # Write the differing rows to a new TSV file
    differing_rows.to_csv('differences.tsv', sep='\t', index=False)

def calc_refolds(afp): 
    """Based on the afp it calculates the number of refolds and adds them into an extra column in the df. 

    Args:
         afp(list): list of afp
    Returns:
        num_refolds: number of refolds
    """
    num_refolds = 0

    for i,step in enumerate(afp[:-1]):

        bp_dict_1 = get_base_pairings_dict(step)
        bp_dict_2 = get_base_pairings_dict(afp[i+1])

        #print(f'{bp_dict_1 = }')
        #print(f'{bp_dict_2 = }')
        #print(f'{step = }')
        #print(f'{afp[i+1] = }')
        if compare_pairing_dicts(bp_dict_1,bp_dict_2):
            num_refolds += 1
        '''
        if step[-1] == ".":
            print(f'{afp[i+1][-2] = }')
            if afp[i+1][-2] != '(':
               num_refolds += 1 

        else:
            if step != afp[i+1][:len(step)]:
                num_refolds += 1 
        '''

    return num_refolds

def compare_pairing_dicts(bp_dict_1, bp_dict_2):
    # Identify the shared keys in both dictionaries
    shared_keys = set(bp_dict_1.keys()).intersection(bp_dict_2.keys())

    # Collect changes in pairings for shared keys
    changes = {}
    for key in shared_keys:
        if bp_dict_1.get(key) != bp_dict_2.get(key):
            return True

    return False


def get_base_pairings_dict(dot_bracket):
    stack = []
    pairings_dict = {}

    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                opening_index = stack.pop()
                pairings_dict[opening_index] = i
                pairings_dict[i] = opening_index
    
    return pairings_dict

def main():

    folder_num = "paper"
    out_folder = 'output_' + folder_num + '/stat_figures/'

    sort_and_print_missing_tsv('output_' + folder_num + '/6_steps_out.txt',folder_num)

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
        print(f"Folder '{out_folder}' created.")

    print('Making some nice plots')

    df_list = ['output_' + folder_num + '/sorted_result.tsv']
    df = read_process_tsv(df_list)

    df['AFP_list'] = df['acfp'].apply(lambda x: ast.literal_eval(x))
    df['refolds'] = df['AFP_list'].apply(calc_refolds)

    print(df.info())
    print(df.head())
    # Create a new boolean column for nt_success (successful if value > 0.5)
    df['nt_success_bool'] = df['nt_success'] > 0.5

    # Create a cross-tabulation table with row and column sums (margins)
    success_table = pd.crosstab(df['domain_success'], df['nt_success_bool'], margins=True)
    
    print("Cross Tabulation of Domain Success vs Nucleotide Success:")
    print(success_table)

    #plot_out_paper(out_folder,df)
    plot_domain_and_counts(out_folder, df)
    #plot_split_normalized_domain(out_folder, df)

if __name__ == "__main__":



    #calc_refolds(['.', '..', '(.)', '(())', '(().)', '(())()'])
    main()


    #compare_data("output","output_1")

