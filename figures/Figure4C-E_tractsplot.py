#!/usr/bin/env python

"""
Create a stacked plot with three methods (RFMix, Simplai, Mosaic) in one PDF
Final version with all requested customizations
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from FancyPlot import FancyPlot
import os.path as path
import numpy as np

def create_final_stacked_plot(output_file='stacked_comparison.pdf', 
                             width=3, height=4):
    """
    Create a single PDF with three stacked subplots with all requested customizations
    """
    
    # Configuration for each method
    methods = [
        {
            'name': 'RFMix',
            'input_dir': '../rfmix',
            'color': 'red',
            'display_name': 'RFMix'
        },
        {
            'name': 'Simplai', 
            'input_dir': '../simplai',
            'color': '#e377c2',
            'display_name': 'Simplai'
        },
        {
            'name': 'Mosaic',
            'input_dir': '../mosaic', 
            'color': 'green',
            'display_name': 'Mosaic'
        }
    ]
    
    # Common settings
    common_name = 'one_pulse'
    pop_names = ['HG', 'Farmer']
    markers = ['^', 's']  # triangle for HG, square for Farmer
    
    # Create figure with custom dimensions
    fig, axes = plt.subplots(3, 1, figsize=(width, height))
    
    # Adjust spacing between subplots
    fig.subplots_adjust(hspace=0.25, top=0.95, bottom=0.1, left=0.15, right=0.95)
    
    # Track if any data was successfully loaded
    successful_plots = 0
    
    # Load and plot data for each method
    for i, (method, ax) in enumerate(zip(methods, axes)):
        try:
            print(f"Processing {method['name']}...")
            
            # Set up file paths
            paths = {}
            paths['bins'] = path.join(method['input_dir'], common_name + '_bins')
            paths['dat'] = path.join(method['input_dir'], common_name + '_dat')
            paths['preds'] = [path.join(method['input_dir'], common_name + '_pred')]
            
            # Load data
            fp = FancyPlot.load(paths['bins'], paths['dat'], paths['preds'],
                               pop_names, [common_name])
            
            # Get colors - both populations same color for each method
            colors = [method['color'], method['color']]
            pop_colors, theory_colors = fp.choose_colors(colors=colors)
            
            # CUSTOM DRAW FUNCTION WITH SMALLER MARKERS
            # Instead of using fp.draw(), we'll draw manually to control marker size
            for j, (pop, color, t_colors) in enumerate(zip(fp.populations, pop_colors, theory_colors)):
                # Draw theories (lines and confidence bands)
                for theory, t_color in zip(pop.theories, t_colors):
                    theory.draw(ax, t_color, bin_scale=fp.bin_scale)
                
                # Draw scatter points with SMALLER markers (s=20 instead of default)
                marker = markers[j] if j < len(markers) else 'o'
                ax.scatter(
                    fp.bin_scale * pop.bins[:-1],
                    pop.data[:-1],
                    color=color,
                    marker=marker,
                    s=20,  # SMALLER MARKER SIZE
                    label=(pop.name + " data").replace("_", "-"),
                )
            
            # Set y-scale to log
            ax.set_yscale("log")
            
            # X-AXIS: More ticks and only bottom plot gets x-label
            ax.set_xticks([0, 50, 100, 150, 200, 250])
            if i == len(methods) - 1:  # Bottom plot only
                ax.set_xlabel("tract length (cM)", fontsize=10)
            else:
                ax.tick_params(labelbottom=False)  # Hide x-axis labels for upper plots
            
            # Y-AXIS: Only middle plot (Simplai) gets y-label
            if i == 1:  # Middle plot (Simplai, i=1)
                ax.set_ylabel("number of tracts", fontsize=10)
            
            # METHOD LABELS: No box around the text
            ax.text(0.02, 0.95, method['display_name'], transform=ax.transAxes, 
                    fontsize=9, fontweight='bold', verticalalignment='top')
            
            # Set limits
            ax.set_ylim(bottom=0.92)
            ax.set_xlim(left=0.0, right=250)
            
            # Remove individual legends
            if ax.get_legend():
                ax.get_legend().remove()
                
            successful_plots += 1
            print(f"Successfully plotted {method['name']}")
            
        except Exception as e:
            print(f"Error processing {method['name']}: {e}")
            # Create empty plot with error message
            ax.text(0.5, 0.5, f"Error loading {method['name']}\n{str(e)}", 
                   transform=ax.transAxes, ha='center', va='center')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
    
    if successful_plots > 0:
        # LEGEND: Simple black box, no 3D effects, smaller markers
        legend_elements = [
            plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='black', 
                      markersize=5, label='HG', linestyle='None'),  # Small markers
            plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='black', 
                      markersize=5, label='Farmer', linestyle='None')  # Small markers
        ]
        
        # Add legend to the top subplot - SIMPLE BLACK BOX
        legend = axes[0].legend(handles=legend_elements, loc='upper right', 
                               frameon=True, fancybox=False, shadow=False, 
                               fontsize=8, edgecolor='black', facecolor='white',
                               framealpha=1.0)
    
    # Save the figure
    fig.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
    print(f"Stacked plot saved to: {output_file}")
    print(f"Figure dimensions: {width}\" x {height}\"")
    
    plt.close()
    return successful_plots > 0

def main():
    """Main function"""
    
    success = create_final_stacked_plot(
        output_file='stacked_comparison.pdf',
        width=3,
        height=4
    )
    
    if success:
        print("Plot creation completed successfully!")
        print("✓ Method names without boxes")
        print("✓ Smaller data points") 
        print("✓ Y-axis label only on middle plot (Simplai)")
        print("✓ X-axis ticks: 0, 50, 100, 150, 200, 250")
        print("✓ Simple legend with black box (no 3D effects)")
    else:
        print("No plots were successfully created. Check your file paths.")

if __name__ == "__main__":
    main()
