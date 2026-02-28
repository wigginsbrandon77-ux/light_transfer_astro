#!/usr/bin/env python3
"""
make_light_curve.py

Process all VTK files in a directory with radiative transfer
and create a supernova light curve (luminosity vs time)

Usage:
    python make_light_curve.py <vtk_directory> [options]
"""

import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
from glob import glob
import argparse

class LightCurveGenerator:
    def __init__(self, vtk_dir, output_dir='lightcurve_output'):
        self.vtk_dir = vtk_dir
        self.output_dir = output_dir
        Path(output_dir).mkdir(exist_ok=True)
        
        # Storage for light curve data
        self.times = []
        self.luminosities = []
        self.filenames = []
        
    def find_vtk_files(self):
        """Find all VTK files in directory"""
        pattern = os.path.join(self.vtk_dir, '*.vtk')
        files = sorted(glob(pattern))
        print(f"\nFound {len(files)} VTK files in {self.vtk_dir}")
        return files
    
    def extract_timestep(self, filename):
        """
        Extract timestep number from filename
        Assumes format like: output_000050.vtk
        """
        match = re.search(r'_(\d+)\.vtk', filename)
        if match:
            return int(match.group(1))
        else:
            # Try to extract any number from filename
            match = re.search(r'(\d+)', os.path.basename(filename))
            if match:
                return int(match.group(1))
        return 0
    
    def run_radiative_transfer(self, vtk_file, output_spectrum, 
                               ni56_mass=0.0, ni56_radius=0.3, 
                               time_days=10.0, opacity='gray', 
                               opacity_multiplier=1.0,
                               octant_mode=False, nrays_octant=10,
                               nfreq=50, star_mass=8.0, star_radius=7.0,
                               axis='z', doppler=False):
        """
        Run the C radiative transfer code on a VTK file
        """
        cmd = ['./radtransfer', vtk_file, output_spectrum,
               '-opacity', opacity,
               '-nfreq', str(nfreq),
               '-mass', str(star_mass),
               '-radius', str(star_radius),
               '-kappa', str(opacity_multiplier)]
        
        if octant_mode:
            cmd.extend(['-octant', str(nrays_octant)])
        else:
            cmd.extend(['-axis', axis])
        
        if ni56_mass > 0:
            cmd.extend(['-ni56', str(ni56_mass)])
            cmd.extend(['-ni56_radius', str(ni56_radius)])
            cmd.extend(['-time', str(time_days)])
        
        if doppler:
            cmd.append('-doppler')
        
        print(f"  Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, 
                                   timeout=600)  # 10 min timeout
            if result.returncode != 0:
                print(f"  ERROR: {result.stderr}")
                return False
            return True
        except subprocess.TimeoutExpired:
            print(f"  ERROR: Timeout after 10 minutes")
            return False
        except FileNotFoundError:
            print(f"  ERROR: radtransfer executable not found!")
            print(f"  Make sure to run 'make' first")
            return False
    
    def integrate_spectrum(self, spectrum_file):
        """
        Read spectrum file and integrate to get total luminosity
        Returns integrated flux in arbitrary units
        """
        try:
            # Read data (skip comment lines)
            data = np.loadtxt(spectrum_file)
            
            frequency = data[:, 0]  # Hz
            flux_nu = data[:, 2]    # Arbitrary units
            
            # Integrate flux over frequency
            # L = ∫ F_ν dν
            luminosity = np.trapz(flux_nu, frequency)
            
            return luminosity
            
        except Exception as e:
            print(f"  ERROR reading {spectrum_file}: {e}")
            return None
    
    def process_all_vtk_files(self, time_per_step=None, t0=0.0,
                               ni56_mass=0.0, ni56_radius=0.3,
                               opacity='gray', opacity_multiplier=1.0, 
                               octant_mode=False, nrays_octant=10,
                               nfreq=50,
                               star_mass=8.0, star_radius=7.0,
                               axis='z', doppler=False,
                               skip_existing=False):
        """
        Process all VTK files in directory
        
        Parameters:
        -----------
        time_per_step : float
            Time per output step (days). If None, uses step number as time.
        t0 : float
            Time of first snapshot (days)
        ni56_mass : float
            Ni-56 mass in M_sun
        ni56_radius : float
            Ni-56 radius in code units
        opacity : str
            Opacity model ('kramers', 'freefree', 'gray')
        opacity_multiplier : float
            Opacity boost factor (default: 1.0)
        nfreq : int
            Number of frequency bins
        star_mass : float
            Stellar mass in M_sun
        star_radius : float
            Stellar radius in R_sun
        axis : str
            Viewing axis ('x', 'y', 'z')
        doppler : bool
            Include Doppler shifting
        skip_existing : bool
            Skip if spectrum file already exists
        """
        
        vtk_files = self.find_vtk_files()
        
        if len(vtk_files) == 0:
            print(f"No VTK files found in {self.vtk_dir}")
            return
        
        print(f"\n{'='*60}")
        print(f"PROCESSING {len(vtk_files)} VTK FILES")
        print(f"{'='*60}\n")
        
        for i, vtk_file in enumerate(vtk_files):
            basename = os.path.basename(vtk_file)
            print(f"\n[{i+1}/{len(vtk_files)}] Processing {basename}")
            
            # Extract timestep
            step = self.extract_timestep(vtk_file)
            
            # Calculate time
            if time_per_step is not None:
                time_days = t0 + step * time_per_step
            else:
                time_days = t0 + step  # Use step number as time
            
            print(f"  Timestep: {step}, Time: {time_days:.2f} days")
            
            # Output spectrum filename
            spectrum_file = os.path.join(self.output_dir, 
                                        f'spectrum_step{step:06d}.txt')
            
            # Skip if already exists (for restarting)
            if skip_existing and os.path.exists(spectrum_file):
                print(f"  Skipping (already exists): {spectrum_file}")
                luminosity = self.integrate_spectrum(spectrum_file)
                if luminosity is not None:
                    self.times.append(time_days)
                    self.luminosities.append(luminosity)
                    self.filenames.append(basename)
                continue
            
            # Run radiative transfer
            success = self.run_radiative_transfer(
                vtk_file, spectrum_file,
                ni56_mass=ni56_mass,
                ni56_radius=ni56_radius,
                time_days=time_days,
                opacity=opacity,
                opacity_multiplier=opacity_multiplier,
                octant_mode=octant_mode,
                nrays_octant=nrays_octant,
                nfreq=nfreq,
                star_mass=star_mass,
                star_radius=star_radius,
                axis=axis,
                doppler=doppler
            )
            
            if not success:
                print(f"  Failed to process {basename}")
                continue
            
            # Integrate spectrum
            luminosity = self.integrate_spectrum(spectrum_file)
            
            if luminosity is not None:
                self.times.append(time_days)
                self.luminosities.append(luminosity)
                self.filenames.append(basename)
                print(f"  Integrated luminosity: {luminosity:.3e}")
        
        print(f"\n{'='*60}")
        print(f"PROCESSING COMPLETE")
        print(f"Successfully processed {len(self.times)} files")
        print(f"{'='*60}\n")
    
    def save_light_curve(self, filename='lightcurve.txt'):
        """Save light curve data to text file"""
        output_file = os.path.join(self.output_dir, filename)
        
        header = (
            "# Supernova Light Curve\n"
            "# Generated from radiative transfer calculations\n"
            "#\n"
            f"# Number of points: {len(self.times)}\n"
            "#\n"
            "# Columns:\n"
            "#   1. Time (days since explosion)\n"
            "#   2. Integrated luminosity (arbitrary units)\n"
            "#   3. VTK filename\n"
            "#\n"
        )
        
        with open(output_file, 'w') as f:
            f.write(header)
            f.write(f"# {'Time (days)':>15s}  {'Luminosity':>18s}  Filename\n")
            f.write("#\n")
            
            for t, L, fname in zip(self.times, self.luminosities, self.filenames):
                f.write(f"  {t:15.6f}  {L:18.10e}  {fname}\n")
        
        print(f"Light curve data saved to {output_file}")
        return output_file
    
    def plot_light_curve(self, output_file='lightcurve.png',
                        ni56_mass=0.0, show_decay=False):
        """
        Plot the light curve
        
        Parameters:
        -----------
        output_file : str
            Output filename for plot
        ni56_mass : float
            Ni-56 mass (if >0, show expected decay curve)
        show_decay : bool
            Show theoretical Ni-56 decay curve for comparison
        """
        if len(self.times) == 0:
            print("No data to plot!")
            return
        
        # Sort by time
        sorted_indices = np.argsort(self.times)
        times = np.array(self.times)[sorted_indices]
        luminosities = np.array(self.luminosities)[sorted_indices]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot light curve
        ax.plot(times, luminosities, 'o-', color='red', lw=2, ms=8,
               label='Simulation', alpha=0.8)
        
        # Show theoretical decay curve if requested
        if show_decay and ni56_mass > 0:
            t_theory = np.linspace(times.min(), times.max(), 200)
            
            # Ni-56 and Co-56 decay
            tau_ni = 6.075 / np.log(2)  # days
            tau_co = 77.27 / np.log(2)  # days
            
            # Luminosity components (normalized)
            L_ni = np.exp(-t_theory / tau_ni)
            L_co = (tau_co / (tau_co - tau_ni)) * \
                   (np.exp(-t_theory / tau_ni) - np.exp(-t_theory / tau_co))
            L_total = L_ni + L_co
            
            # Scale to match peak of simulation
            scale = luminosities.max() / L_total.max()
            
            ax.plot(t_theory, scale * L_total, '--', color='blue', lw=2,
                   alpha=0.7, label=f'Ni-56 decay (M={ni56_mass:.2f} M☉)')
        
        ax.set_xlabel('Time Since Explosion (days)', fontsize=14)
        ax.set_ylabel('Integrated Luminosity (arbitrary units)', fontsize=14)
        ax.set_title('Supernova Light Curve', fontsize=16, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3, which='both')
        ax.legend(fontsize=12, loc='best')
        
        # Add text with info
        info_text = (
            f"Points: {len(times)}\n"
            f"Time range: {times.min():.1f} - {times.max():.1f} days\n"
            f"Peak: {luminosities.max():.2e} at {times[np.argmax(luminosities)]:.1f} days"
        )
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        plot_path = os.path.join(self.output_dir, output_file)
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Light curve plot saved to {plot_path}")
        plt.close()
    
    def plot_spectra_evolution(self, output_file='spectra_evolution.png',
                              max_spectra=20):
        """
        Plot how spectra evolve over time
        Shows a selection of spectra at different epochs
        """
        spectrum_files = sorted(glob(os.path.join(self.output_dir, 'spectrum_*.txt')))
        
        if len(spectrum_files) == 0:
            print("No spectrum files found for evolution plot")
            return
        
        # Select subset if too many
        if len(spectrum_files) > max_spectra:
            indices = np.linspace(0, len(spectrum_files)-1, max_spectra, dtype=int)
            spectrum_files = [spectrum_files[i] for i in indices]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(spectrum_files)))
        
        for i, spec_file in enumerate(spectrum_files):
            try:
                data = np.loadtxt(spec_file)
                freq = data[:, 0]
                wavelength = data[:, 1]
                flux_nu = data[:, 2]
                
                # Extract timestep from filename
                step = self.extract_timestep(spec_file)
                
                ax1.loglog(freq, flux_nu, color=colors[i], lw=1.5, alpha=0.7)
                ax2.loglog(wavelength, data[:, 3], color=colors[i], lw=1.5, 
                          alpha=0.7, label=f'Step {step}' if i % 3 == 0 else '')
                
            except Exception as e:
                print(f"Error reading {spec_file}: {e}")
        
        ax1.set_xlabel('Frequency (Hz)', fontsize=12)
        ax1.set_ylabel('Flux Density F$_ν$ (arbitrary)', fontsize=12)
        ax1.set_title('Spectral Evolution', fontsize=14)
        ax1.grid(True, alpha=0.3)
        
        ax2.set_xlabel('Wavelength (Å)', fontsize=12)
        ax2.set_ylabel('Flux Density F$_λ$ (arbitrary)', fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.invert_xaxis()
        ax2.legend(fontsize=8, ncol=3, loc='best')
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap='viridis',
                                   norm=plt.Normalize(0, len(spectrum_files)-1))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=[ax1, ax2], label='Time progression')
        
        plt.tight_layout()
        
        plot_path = os.path.join(self.output_dir, output_file)
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Spectral evolution plot saved to {plot_path}")
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate supernova light curve from VTK files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (no Ni-56, gray opacity)
  python make_light_curve.py output/
  
  # With Ni-56 decay, specify time per step
  python make_light_curve.py output/ -ni56 0.6 -dt 0.1
  
  # Use octant geometry (for simulations with core at origin)
  python make_light_curve.py output/ -octant -nrays 10 -kappa 10
  
  # Full octant mode with Ni-56
  python make_light_curve.py output/ -ni56 0.6 -dt 0.1 -t0 5 \\
      -octant -nrays 15 -kappa 50 -opacity gray
  
  # Skip already processed files (for restart)
  python make_light_curve.py output/ -ni56 0.6 -octant -skip
        """
    )
    
    parser.add_argument('vtk_dir', help='Directory containing VTK files')
    parser.add_argument('-output', default='lightcurve_output',
                       help='Output directory (default: lightcurve_output)')
    
    # Time parameters
    parser.add_argument('-dt', '--time_per_step', type=float,
                       help='Time per output step in days (if not set, uses step number)')
    parser.add_argument('-t0', type=float, default=0.0,
                       help='Time of first snapshot in days (default: 0)')
    
    # Ni-56 parameters
    parser.add_argument('-ni56', type=float, default=0.0,
                       help='Ni-56 mass in M_sun (0 = disabled, default: 0)')
    parser.add_argument('-ni56_radius', type=float, default=0.3,
                       help='Ni-56 radius in code units (default: 0.3)')
    
    # Physics parameters
    parser.add_argument('-opacity', default='gray',
                       choices=['kramers', 'freefree', 'gray'],
                       help='Opacity model (default: gray)')
    parser.add_argument('-kappa', '--opacity_multiplier', type=float, default=1.0,
                       help='Opacity multiplier (default: 1.0, try 10-100 for optically thick)')
    parser.add_argument('-octant', '--octant_mode', action='store_true',
                       help='Use octant geometry (radial rays from origin)')
    parser.add_argument('-nrays', '--nrays_octant', type=int, default=10,
                       help='Number of rays per angle in octant mode (default: 10, gives 10×10=100 rays)')
    parser.add_argument('-nfreq', type=int, default=50,
                       help='Number of frequency bins (default: 50)')
    parser.add_argument('-mass', type=float, default=8.0,
                       help='Stellar mass in M_sun (default: 8.0)')
    parser.add_argument('-radius', type=float, default=7.0,
                       help='Stellar radius in R_sun (default: 7.0)')
    parser.add_argument('-axis', default='z', choices=['x', 'y', 'z'],
                       help='Viewing axis (default: z, ignored if -octant is used)')
    parser.add_argument('-doppler', action='store_true',
                       help='Include Doppler shifting')
    
    # Processing options
    parser.add_argument('-skip', '--skip_existing', action='store_true',
                       help='Skip files that already have spectrum output')
    parser.add_argument('-no_plot', action='store_true',
                       help='Do not create plots (only generate data)')
    parser.add_argument('-show_decay', action='store_true',
                       help='Show theoretical Ni-56 decay curve on plot')
    
    args = parser.parse_args()
    
    # Check if radtransfer executable exists
    if not os.path.exists('./radtransfer'):
        print("ERROR: radtransfer executable not found!")
        print("Please run 'make' first to build the code")
        return 1
    
    # Create light curve generator
    lc = LightCurveGenerator(args.vtk_dir, args.output)
    
    # Process all VTK files
    lc.process_all_vtk_files(
        time_per_step=args.time_per_step,
        t0=args.t0,
        ni56_mass=args.ni56,
        ni56_radius=args.ni56_radius,
        opacity=args.opacity,
        opacity_multiplier=args.opacity_multiplier,
        octant_mode=args.octant_mode,
        nrays_octant=args.nrays_octant,
        nfreq=args.nfreq,
        star_mass=args.mass,
        star_radius=args.radius,
        axis=args.axis,
        doppler=args.doppler,
        skip_existing=args.skip_existing
    )
    
    if len(lc.times) == 0:
        print("\nNo data collected. Exiting.")
        return 1
    
    # Save light curve data
    lc.save_light_curve('lightcurve.txt')
    
    # Create plots (unless disabled)
    if not args.no_plot:
        lc.plot_light_curve('lightcurve.png', 
                           ni56_mass=args.ni56,
                           show_decay=args.show_decay)
        lc.plot_spectra_evolution('spectra_evolution.png')
    
    print("\n" + "="*60)
    print("LIGHT CURVE GENERATION COMPLETE!")
    print("="*60)
    print(f"\nResults in {args.output}/:")
    print(f"  - lightcurve.txt        (data file)")
    if not args.no_plot:
        print(f"  - lightcurve.png        (light curve plot)")
        print(f"  - spectra_evolution.png (spectral evolution)")
    print(f"  - spectrum_*.txt        (individual spectra)")
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
