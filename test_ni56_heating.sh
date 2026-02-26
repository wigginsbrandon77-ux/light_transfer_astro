#!/bin/bash
#
# test_ni56_heating.sh
# Demonstrate Ni-56 decay heating effects
#

if [ $# -lt 1 ]; then
    echo "Usage: ./test_ni56_heating.sh <vtk_file>"
    echo ""
    echo "Example:"
    echo "  ./test_ni56_heating.sh output_000050.vtk"
    echo ""
    exit 1
fi

VTK_FILE=$1

echo "=========================================="
echo "NI-56 DECAY HEATING DEMONSTRATION"
echo "=========================================="
echo ""
echo "VTK file: $VTK_FILE"
echo ""

# Make sure the code is compiled
if [ ! -f ./radtransfer ]; then
    echo "Building radtransfer..."
    make
    echo ""
fi

# Create output directory
mkdir -p ni56_test

echo "=========================================="
echo "Test 1: Without Ni-56 (adiabatic cooling)"
echo "=========================================="
./radtransfer "$VTK_FILE" ni56_test/spectrum_no_heating.txt \
    -opacity gray -nfreq 100
echo ""

echo "=========================================="
echo "Test 2: With Ni-56 at 10 days"
echo "=========================================="
./radtransfer "$VTK_FILE" ni56_test/spectrum_10days.txt \
    -ni56 0.6 -time 10 -opacity gray -nfreq 100
echo ""

echo "=========================================="
echo "Test 3: With Ni-56 at 20 days (peak)"
echo "=========================================="
./radtransfer "$VTK_FILE" ni56_test/spectrum_20days.txt \
    -ni56 0.6 -time 20 -opacity gray -nfreq 100
echo ""

echo "=========================================="
echo "Test 4: With Ni-56 at 50 days (late)"
echo "=========================================="
./radtransfer "$VTK_FILE" ni56_test/spectrum_50days.txt \
    -ni56 0.6 -time 50 -opacity gray -nfreq 100
echo ""

echo "=========================================="
echo "Test 5: Different Ni-56 masses at 20 days"
echo "=========================================="
for mni in 0.3 0.6 1.0; do
    echo "  Ni-56 mass: ${mni} M_sun"
    ./radtransfer "$VTK_FILE" ni56_test/spectrum_mni${mni}.txt \
        -ni56 $mni -time 20 -opacity gray -nfreq 100
done
echo ""

echo "=========================================="
echo "TESTS COMPLETE!"
echo "=========================================="
echo ""
echo "Output files in ni56_test/:"
ls -lh ni56_test/
echo ""
echo "Compare these spectra to see Ni-56 effects:"
echo ""
echo "1. No heating vs. With heating:"
echo "   - no_heating.txt should be dimmer and redder"
echo "   - With Ni-56: brighter and bluer"
echo ""
echo "2. Time evolution (10d, 20d, 50d):"
echo "   - Luminosity decreases with time"
echo "   - 20d should be brightest (near peak)"
echo "   - 50d shows Co-56 decay dominance"
echo ""
echo "3. Mass dependence (0.3, 0.6, 1.0 M_sun):"
echo "   - More Ni-56 → brighter spectrum"
echo "   - Roughly linear: 2x mass ≈ 2x flux"
echo ""
echo "To plot these:"
echo "  1. Open each file in Excel"
echo "  2. Plot column 1 vs 3 (log-log scale)"
echo "  3. Compare total flux and peak location"
echo ""
echo "Or upload all files to an LLM and ask:"
echo "  'Compare these spectra and show how Ni-56 affects luminosity'"
echo ""

# Create Python plotting script if available
if command -v python3 &> /dev/null; then
    cat > ni56_test/plot_ni56_comparison.py << 'EOF'
#!/usr/bin/env python3
"""
Plot comparison of spectra with different Ni-56 parameters
"""
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

# Find all spectrum files
files = sorted(glob('spectrum_*.txt'))

if len(files) == 0:
    print("No spectrum files found!")
    exit(1)

print(f"Found {len(files)} spectrum files")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Colors for different spectra
colors = plt.cm.viridis(np.linspace(0, 1, len(files)))

for i, file in enumerate(files):
    try:
        data = np.loadtxt(file)
        freq = data[:, 0]
        wavelength = data[:, 1]
        flux_nu = data[:, 2]
        
        # Extract label from filename
        label = file.replace('spectrum_', '').replace('.txt', '')
        
        # Plot
        ax1.loglog(freq, flux_nu, color=colors[i], lw=2, 
                  label=label, alpha=0.8)
        
        # Calculate total luminosity (integrate flux)
        total_flux = np.trapz(flux_nu, freq)
        print(f"{label:20s}: Total flux = {total_flux:.3e}")
        
    except Exception as e:
        print(f"Error reading {file}: {e}")

ax1.set_xlabel('Frequency (Hz)', fontsize=12)
ax1.set_ylabel('Flux Density (arbitrary)', fontsize=12)
ax1.set_title('Ni-56 Decay Heating Effects on Spectrum', fontsize=14)
ax1.legend(fontsize=9, loc='best')
ax1.grid(True, alpha=0.3)

# Plot integrated flux vs time (if time-series files)
time_files = [f for f in files if any(x in f for x in ['10days', '20days', '50days'])]
if len(time_files) >= 2:
    times = []
    fluxes = []
    
    for file in time_files:
        # Extract time from filename
        if '10days' in file:
            t = 10
        elif '20days' in file:
            t = 20
        elif '50days' in file:
            t = 50
        else:
            continue
        
        data = np.loadtxt(file)
        freq = data[:, 0]
        flux_nu = data[:, 2]
        total_flux = np.trapz(flux_nu, freq)
        
        times.append(t)
        fluxes.append(total_flux)
    
    if times:
        ax2.semilogy(times, fluxes, 'o-', lw=2, ms=8, color='red')
        ax2.set_xlabel('Time Since Explosion (days)', fontsize=12)
        ax2.set_ylabel('Integrated Flux (arbitrary)', fontsize=12)
        ax2.set_title('Light Curve: Luminosity vs Time', fontsize=14)
        ax2.grid(True, alpha=0.3)
        
        # Add Ni-56 and Co-56 decay exponentials for reference
        t_plot = np.linspace(5, 100, 100)
        tau_ni = 6.075 / np.log(2)
        tau_co = 77.27 / np.log(2)
        
        # Normalized decay curves
        L_ni = np.exp(-t_plot / tau_ni)
        L_co = 1.5 * (np.exp(-t_plot / tau_ni) - np.exp(-t_plot / tau_co))
        L_total = L_ni + L_co
        
        # Scale to match data
        if fluxes:
            scale = fluxes[0] / L_total[int((times[0]-5)*2)]
            ax2.plot(t_plot, scale * L_total, '--', lw=1.5, 
                    color='blue', alpha=0.7, label='Expected decay')
            ax2.legend(fontsize=10)
else:
    ax2.text(0.5, 0.5, 'Time-series data not found\n(need 10days, 20days, 50days)', 
            ha='center', va='center', transform=ax2.transAxes, fontsize=12)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)

plt.tight_layout()
plt.savefig('ni56_comparison.png', dpi=150)
print("\nPlot saved to ni56_comparison.png")
EOF

    chmod +x ni56_test/plot_ni56_comparison.py
    
    echo "Python detected! Creating comparison plot..."
    cd ni56_test
    python3 plot_ni56_comparison.py
    cd ..
    echo ""
    echo "Comparison plot saved: ni56_test/ni56_comparison.png"
fi

echo ""
echo "=========================================="
