import matplotlib.pyplot as plt
import pickle
import time
from math import sqrt
import numpy as np
from scipy import linalg

# import all the required functions for sparse coding
from utils import * 

with open("FISTA_v_FLIPS_uniform_eps_seq.pkl", "rb") as file:
    NMSE_seq_FISTA, NMSE_seq_FLIPS, Conv_iter_FISTA, Conv_iter_FLIPS, SNR_seq, reg_param_seq, epsilon_seq, m, N, sigma_noise = pickle.load(file)

def sync_ylim_and_plot(fig, axs, data_dict, x_data, labels, titles, xlabel, ylabel, markers, suptitle, add_crosshairs=False, markersize=5, short_y_labels=False):
    """Helper function to create synchronized plots"""
    # Collect all y-values for synchronization
    all_y_values = []
    for SNR in SNR_seq:
        for key in data_dict:
            all_y_values.extend(data_dict[key][SNR])
    
    y_min, y_max = min(all_y_values), max(all_y_values)
    y_padding = (y_max - y_min) * 0.05
    y_limits = [y_min - y_padding, y_max + y_padding]
    
    # Store minimum NMSE info for printing cumulative values
    min_nmse_info = {}
    
    # Plot data and store line objects for color reference
    line_colors = {}
    for i, (key, data) in enumerate(data_dict.items()):
        for SNR in SNR_seq:
            y_values = np.array(data[SNR])
            x_values = x_data[i] if isinstance(x_data[i], np.ndarray) else x_data[i][SNR]
            line = axs[i].plot(x_values, y_values, marker=markers[i], markersize=markersize, 
                              linestyle='-', label=fr' SNR = ${SNR}$dB')
            
            # Store color for crosshairs
            if add_crosshairs:
                color = line[0].get_color()
                if SNR not in line_colors:
                    line_colors[SNR] = color
                
                # Find minimum NMSE and corresponding x-value
                min_idx = np.argmin(y_values)
                min_nmse = y_values[min_idx]
                min_x = x_values[min_idx] if hasattr(x_values, '__getitem__') else x_values
                
                # Store for later use in cumulative printing
                if SNR not in min_nmse_info:
                    min_nmse_info[SNR] = {}
                min_nmse_info[SNR][key] = {'min_idx': min_idx, 'min_x': min_x}
                
                # Add crosshairs
                axs[i].axhline(y=min_nmse, color=color, linestyle='--', alpha=0.7, linewidth=1)
                axs[i].axvline(x=min_x, color=color, linestyle='--', alpha=0.7, linewidth=1)
        
        # Configure subplot
        axs[i].set_ylim(y_limits)
        axs[i].set_title(titles[i], fontsize=15)
        axs[i].set_xlabel(xlabel[i], fontsize=15)
        axs[i].set_ylabel(ylabel, fontsize=15)
        axs[i].legend()
        
        # Format y-axis labels if requested
        if short_y_labels:
            axs[i].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format_number_short(x)))
    
    plt.tight_layout()
    
    # Return min_nmse_info for use in printing cumulative values
    return min_nmse_info

def format_number_short(num):
    """Format numbers in short form (10k instead of 10000)"""
    if abs(num) >= 1000000:
        return f'{num/1000000:.1f}M'
    elif abs(num) >= 1000:
        return f'{num/1000:.1f}k'
    else:
        return f'{num:.0f}'

# NMSE plots
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

nmse_data = {'FLIPS': NMSE_seq_FLIPS, 'FISTA': NMSE_seq_FISTA}
nmse_x_data = [
    {SNR: epsilon_seq[SNR]/(sigma_noise * np.sqrt(m * N)) for SNR in SNR_seq},
    reg_param_seq
]
nmse_titles = [
    'FLIPS',
    'FISTA'
]
nmse_xlabel = [r'$ \frac{\epsilon}{\sigma_w \sqrt{m N}} $', r'regularization parameter: $\lambda$']

min_nmse_info = sync_ylim_and_plot(fig, axs, nmse_data, nmse_x_data, SNR_seq, nmse_titles, 
                   nmse_xlabel, "NMSE", ['o', 'x'], r'NMSE of reconstruction: $ \frac{\| F_{true} - \widehat{F} (\epsilon) \|_{fro}^2}{\| F_{true} \|_{fro}^2}$', 
                   add_crosshairs=True, markersize = 2)

# Convergence iterations plots
fig1, axs1 = plt.subplots(1, 2, figsize=(14, 4))

conv_data = {'FLIPS': Conv_iter_FLIPS, 'FISTA': Conv_iter_FISTA}
conv_x_data = nmse_x_data  # Same x-axis data as NMSE plots
conv_titles = [
    r'FLIPS',
    r'FISTA'
]
conv_xlabel = [r'$ \frac{\epsilon}{\sigma_w \sqrt{m N}} $', r'regularization parameter: $\lambda$']

sync_ylim_and_plot(fig1, axs1, conv_data, conv_x_data, SNR_seq, conv_titles, 
                   conv_xlabel, "# iterations", ['o', 'x'], 'Iterations required for convergence', markersize = 3)

# Cumulative convergence iterations plots
fig2, axs2 = plt.subplots(1, 2, figsize=(14, 4))

# Calculate cumulative iterations
cumulative_conv_FLIPS = {}
cumulative_conv_FISTA = {}

for SNR in SNR_seq:
    # For FLIPS - simple cumulative sum
    conv_iter_flips = np.array(Conv_iter_FLIPS[SNR])
    cumulative_conv_FLIPS[SNR] = np.cumsum(conv_iter_flips)
    
    # For FISTA - simple cumulative sum
    conv_iter_fista = np.array(Conv_iter_FISTA[SNR])
    cumulative_conv_FISTA[SNR] = np.cumsum(conv_iter_fista)

cumulative_conv_data = {'FLIPS': cumulative_conv_FLIPS, 'FISTA': cumulative_conv_FISTA}
cumulative_conv_titles = [
    r'FLIPS',
    r'FISTA'
]

sync_ylim_and_plot(fig2, axs2, cumulative_conv_data, conv_x_data, SNR_seq, cumulative_conv_titles, 
                   conv_xlabel, "Cumulative # iterations", ['o', 'x'], 'Cumulative iterations required for convergence',
                   short_y_labels=True, markersize = 2)

# Print cumulative iteration values at minimum NMSE points
print("\nCumulative iterations at minimum NMSE points:")
for SNR in SNR_seq:
    print(f"\nSNR = {SNR} dB:")
    for method in ['FLIPS', 'FISTA']:
        if SNR in min_nmse_info and method in min_nmse_info[SNR]:
            min_idx = min_nmse_info[SNR][method]['min_idx']
            if method == 'FLIPS':
                cum_iter = cumulative_conv_FLIPS[SNR][min_idx]
            else:
                cum_iter = cumulative_conv_FISTA[SNR][min_idx]
            print(f"  {method}: {cum_iter} iterations")

plt.show()