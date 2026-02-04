import matplotlib.pyplot as plt
import sys

def parse_emr_data(filename):
    # Data structure: dict of lists, keyed by R2
    # {r2_val: {'H': [], 'R': []}}
    data = {}
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if "H(T)" in line or not line.strip():
            continue
        try:
            parts = line.split()
            # H, R, R2, Sigma1, Sigma2
            # Use columns 0, 1, 2
            h = float(parts[0])
            r = float(parts[1])
            r2 = float(parts[2])
            
            if r2 not in data:
                data[r2] = {'H': [], 'R': []}
                
            data[r2]['H'].append(h)
            data[r2]['R'].append(r)
            
        except ValueError:
            continue
            
    return data

def plot_emr(data):
    plt.figure(figsize=(10, 6))
    
    # Sort R2 keys for consistent legend
    r2_keys = sorted(data.keys())
    
    for r2 in r2_keys:
        h_vals = data[r2]['H']
        r_vals = data[r2]['R']
        
        # Sort by H
        # Combine into pairs, sort, unzip
        pairs = sorted(zip(h_vals, r_vals))
        h_sorted = [p[0] for p in pairs]
        r_sorted = [p[1] for p in pairs]
        
        # Calculate MR % (R(H) - R(0)) / R(0) * 100
        # Find R at approx H=0 (min abs H)
        min_abs_h = float('inf')
        r0 = 1.0
        
        for h, r in pairs:
            if abs(h) < min_abs_h:
                min_abs_h = abs(h)
                r0 = r
        
        mr = [(r - r0) / r0 * 100.0 for r in r_sorted]
        
        plt.plot(h_sorted, mr, marker='o', label=f'R2 = {r2:.2e} m')
        
    plt.xlabel('Magnetic Field (T)')
    plt.ylabel('Magnetoresistance (%)')
    plt.title('EMR vs Magnetic Field')
    plt.legend()
    plt.grid(True)
    plt.savefig('../output/EMR_plot.png')
    print("Plot saved to ../output/EMR_plot.png")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        fname = sys.argv[1]
    else:
        fname = "EMRdata.out"
        
    try:
        data = parse_emr_data(fname)
        if not data:
            print("No data found in " + fname)
        else:
            plot_emr(data)
    except Exception as e:
        print(f"Error plotting: {e}")
