#!/usr/bin/env python3
"""
Script to extract MIP scores and not found reads fractions from traversal.log files
and create a scatter plot.
"""

import os
import re
import matplotlib.pyplot as plt
import glob
bad_tangle_file = "/gpfs/gsfs11/users/antipovd2/res/TTT_paper/HG002/tangle7_problems.list"

def extract_scores_from_log(log_path):
    """
    Extract normalized MIP score and not found reads fraction from a log file.
    
    Returns:
        tuple: (mip_score, not_found_fraction, directory_name) or None if data not found
    """
    mip_score = None
    not_found_fraction = None
    
    try:
        with open(log_path, 'r') as f:
            content = f.read()
            
        # Extract MIP score (handles both regular and exponential notation)
        mip_pattern = r"Normalized score of the MIP solution:\s*([\d\.]+(?:[eE][+-]?\d+)?)"
        mip_matches = re.findall(mip_pattern, content)
        if mip_matches:
            mip_score = float(mip_matches[-1])
            
        # Extract not found reads fraction (handles both regular and exponential notation)
        not_found_pattern = r"Not found reads fraction\s*:\s*([\d\.]+(?:[eE][+-]?\d+)?)"
        not_found_matches = re.findall(not_found_pattern, content)
        if not_found_matches:
            not_found_fraction = float(not_found_matches[-1])
            
        if mip_score is not None and not_found_fraction is not None:
            # Get directory name for labeling
            dir_name = os.path.basename(os.path.dirname(os.path.dirname(log_path)))
            #tangle_X
            dir_name = dir_name.split('_')[1]
            return (mip_score, not_found_fraction, dir_name)
            
    except Exception as e:
        print(f"Error reading {log_path}: {e}")
        
    return None

def get_problematic_ids(filename):
    res = set()
    with open(filename, 'r') as f:
        for line in f:
            res.add(line.strip())
    return res

def main():
    # Find all traversal.log files
    log_files = glob.glob("**//TTT_8/traversal.log", recursive=True)
    
    print(f"Found {len(log_files)} traversal.log files")
    
    # Extract data from all log files
    data_points = []
    
    for log_file in log_files:
        if  log_file.find("tangle_143")!= -1:
            continue
        if  log_file.find("tangle_rDNA")!= -1:
            continue
        result = extract_scores_from_log(log_file)
        if result:
            mip_score, not_found_fraction, dir_name = result
            data_points.append((mip_score, not_found_fraction, dir_name))
            print(f"{dir_name}: coverage inconsistency={mip_score:.6f}, alignment inconsistency={not_found_fraction:.6f}")
        else:
            print(f"Could not extract data from {log_file}")
    
    if not data_points:
        print("No data points found!")
        return
    
    # Separate X and Y coordinates
    x_values = [point[0] for point in data_points]  # MIP scores
    y_values = [point[1] for point in data_points]  # Not found fractions
    labels = [point[2] for point in data_points]    # Directory names
    red_points = get_problematic_ids(bad_tangle_file)
    # Create the plot
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(x_values, y_values, alpha=0.7, s=20)

    # Highlight problematic points in red
    for red_point in red_points:
        if red_point in labels:
            idx = labels.index(red_point)
            plt.scatter(x_values[idx], y_values[idx], color='red', s=20)

    # Add labels to points
    for i, (x, y, label) in enumerate(zip(x_values, y_values, labels)):
        plt.annotate(label, (x, y), xytext=(5, 5), textcoords='offset points', 
                    fontsize=8, alpha=0.4)
    #scatter = plt.scatter(x_values, y_values, alpha=0.7, s=50)
    plt.xlabel('Multiplicity inconsistency (X)')
    plt.ylabel('Alignments inconsistency (Y)')
    plt.title('Multiplicity vs Alignments inconsistency fractions')
    plt.grid(True, alpha=0.3)
    
    # Add some statistics
    plt.text(0.02, 0.98, f'Total points: {len(data_points)}', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save the plot
    output_file = 'score_analysis_plot.png'
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved as {output_file}")
    
    # Also save data as CSV
    csv_file = 'extracted_scores.csv'
    with open(csv_file, 'w') as f:
        f.write("Directory,inconsistent coverage,inconsistent alignment\n")
        for mip_score, not_found_fraction, dir_name in data_points:
            f.write(f"{dir_name},{mip_score},{not_found_fraction}\n")
    print(f"Data saved as {csv_file}")
    
    plt.show()

if __name__ == "__main__":
    main()
