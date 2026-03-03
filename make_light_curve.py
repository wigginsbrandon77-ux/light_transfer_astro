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
from multiprocessing import Pool, cpu_count

# Physical constants (CGS)
H_CGS = 6.62607015e-27
K_CGS = 1.380649e-16
C_CGS = 2.99792458e10

# Band definitions in Hz (UV, Optical, IR)
BANDS = {
    'UV':      (C_CGS / (3000e-8), C_CGS / (100e-8)),
    'Optical': (C_CGS / (7000e-8), C_CGS / (3000e-8)),
    'IR':      (C_CGS / (1e5 * 1e-8), C_CGS / (7000e-8)),
}
BAND_COLORS = {'UV': 'blueviolet', 'Optical': 'gold', 'IR': 'tomato'}


class LightCurveGenerator:
    def __init__(self, vtk_dir, output_dir='lightcurve_output'):
        self.vtk_dir = vtk_dir
        self.output_dir = output_dir
        Path(output_dir).mkdir(exist_ok=True)

        self.times = []
        self.luminosities = []
        self.color_temperatures = []
        self.band_luminosities = {b: [] for b in BANDS}
        self.filenames = []

    def find_vtk_files(self):
        pattern = os.path.join(self.vtk_dir, '*.vtk')
        files = sorted(glob(pattern))
        print(f"\nFound {len(files)} VTK files in {self.vtk_dir}")
        return files

    def extract_timestep(self, filename):
        match = re.search(r'_(\d+)\.vtk', filename)
        if match:
            return int(match.group(1))
        match = re.search(r'(\d+)', os.path.basename(filename))
        if match:
            return int(match.group(1))
        return 0

    def run_radiative_transfer(self, vtk_file, output_spectrum,
                               ni56_mass=0.0, ni56_radius=0.3,
                               time_days=10.0, opacity='kramers',
                               opacity_multiplier=1.0,
                               octant_mode=False, nrays_octant=10,
                               nfreq=100, star_mass=13.0, star_radius=6.41,
                               axis='z', doppler=False):
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

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if result.returncode != 0:
                print(f"  ERROR: {result.stderr[:200]}")
                return False
            return True
        except subprocess.TimeoutExpired:
            print(f"  ERROR: Timeout")
            return False
        except FileNotFoundError:
            print(f"  ERROR: radtransfer not found. Run 'make' first.")
            return False

    def integrate_spectrum(self, spectrum_file):
        """
        Returns (bolometric_lum, T_color, band_dict) or None on failure.

        Color temperature is derived from the Wien peak:
            T_color = h * nu_peak / (2.821 * k)
        """
        try:
            data = np.loadtxt(spectrum_file)
            freq    = data[:, 0]
            flux_nu = data[:, 2]

            if np.all(flux_nu == 0) or len(freq) < 3:
                return None

            luminosity = np.trapz(flux_nu, freq)

            ipeak  = np.argmax(flux_nu)
            nu_peak = freq[ipeak]
            T_color = (H_CGS * nu_peak) / (2.821 * K_CGS) if nu_peak > 0 else np.nan

            bands = {}
            for band_name, (nu_lo, nu_hi) in BANDS.items():
                mask = (freq >= nu_lo) & (freq <= nu_hi)
                bands[band_name] = np.trapz(flux_nu[mask], freq[mask]) if mask.sum() > 1 else 0.0

            return luminosity, T_color, bands

        except Exception as e:
            print(f"  ERROR reading {spectrum_file}: {e}")
            return None

    def _process_single(self, args):
        """Worker: process one VTK file. Must be picklable (no self-mutation)."""
        i, n_total, vtk_file, kwargs = args
        basename  = os.path.basename(vtk_file)
        step      = self.extract_timestep(vtk_file)
        dt        = kwargs['time_per_step']
        t0        = kwargs['t0']
        time_days = t0 + (step * dt if dt is not None else step)

        print(f"  [{i+1}/{n_total}] {basename}  (t = {time_days:.2f} d)")

        spectrum_file = os.path.join(self.output_dir, f'spectrum_step{step:06d}.txt')

        if kwargs['skip_existing'] and os.path.exists(spectrum_file):
            result = self.integrate_spectrum(spectrum_file)
            return (time_days, result, basename) if result else None

        success = self.run_radiative_transfer(
            vtk_file, spectrum_file,
            ni56_mass=kwargs['ni56_mass'], ni56_radius=kwargs['ni56_radius'],
            time_days=time_days, opacity=kwargs['opacity'],
            opacity_multiplier=kwargs['opacity_multiplier'],
            octant_mode=kwargs['octant_mode'], nrays_octant=kwargs['nrays_octant'],
            nfreq=kwargs['nfreq'], star_mass=kwargs['star_mass'],
            star_radius=kwargs['star_radius'], axis=kwargs['axis'],
            doppler=kwargs['doppler'],
        )

        if not success:
            return None

        result = self.integrate_spectrum(spectrum_file)
        return (time_days, result, basename) if result else None

    def process_all_vtk_files(self, time_per_step=None, t0=0.0,
                               ni56_mass=0.0, ni56_radius=0.3,
                               opacity='kramers', opacity_multiplier=1.0,
                               octant_mode=False, nrays_octant=10,
                               nfreq=100, star_mass=13.0, star_radius=6.41,
                               axis='z', doppler=False,
                               skip_existing=False, nworkers=1):

        vtk_files = self.find_vtk_files()
        if not vtk_files:
            return

        print(f"\n{'='*60}")
        print(f"PROCESSING {len(vtk_files)} VTK FILES  (workers: {nworkers})")
        print(f"{'='*60}\n")

        kwargs = dict(
            time_per_step=time_per_step, t0=t0,
            ni56_mass=ni56_mass, ni56_radius=ni56_radius,
            opacity=opacity, opacity_multiplier=opacity_multiplier,
            octant_mode=octant_mode, nrays_octant=nrays_octant,
            nfreq=nfreq, star_mass=star_mass, star_radius=star_radius,
            axis=axis, doppler=doppler, skip_existing=skip_existing,
        )
        task_args = [(i, len(vtk_files), f, kwargs) for i, f in enumerate(vtk_files)]

        if nworkers > 1:
            with Pool(processes=nworkers) as pool:
                results = pool.map(self._process_single, task_args)
        else:
            results = [self._process_single(a) for a in task_args]

        for r in results:
            if r is None:
                continue
            time_days, (lum, T_color, bands), basename = r
            self.times.append(time_days)
            self.luminosities.append(lum)
            self.color_temperatures.append(T_color)
            for b in BANDS:
                self.band_luminosities[b].append(bands.get(b, 0.0))
            self.filenames.append(basename)
            print(f"  t={time_days:.2f} d  L={lum:.3e}  T_color={T_color:.3e} K")

        print(f"\n{'='*60}")
        print(f"PROCESSING COMPLETE — {len(self.times)} files succeeded")
        print(f"{'='*60}\n")

    def save_light_curve(self, filename='lightcurve.txt'):
        output_file = os.path.join(self.output_dir, filename)
        order = np.argsort(self.times)

        with open(output_file, 'w') as f:
            f.write("# Supernova Light Curve\n")
            f.write(f"# Points: {len(self.times)}\n#\n")
            f.write(f"# {'Time(d)':>12s}  {'L_bol':>14s}  {'T_color(K)':>12s}"
                    f"  {'L_UV':>14s}  {'L_Optical':>14s}  {'L_IR':>14s}  Filename\n#\n")
            for idx in order:
                f.write(
                    f"  {self.times[idx]:12.6f}"
                    f"  {self.luminosities[idx]:14.6e}"
                    f"  {self.color_temperatures[idx]:12.3e}"
                    f"  {self.band_luminosities['UV'][idx]:14.6e}"
                    f"  {self.band_luminosities['Optical'][idx]:14.6e}"
                    f"  {self.band_luminosities['IR'][idx]:14.6e}"
                    f"  {self.filenames[idx]}\n"
                )

        print(f"Light curve data saved to {output_file}")
        return output_file

    def plot_light_curve(self, output_file='lightcurve.png',
                         ni56_mass=0.0, show_decay=False):
        if not self.times:
            print("No data to plot!")
            return

        order   = np.argsort(self.times)
        times   = np.array(self.times)[order]
        lum     = np.array(self.luminosities)[order]
        T_color = np.array(self.color_temperatures)[order]

        fig, axes = plt.subplots(3, 1, figsize=(12, 13), sharex=True)
        fig.suptitle('Supernova Light Curve', fontsize=16, fontweight='bold', y=0.98)
        fig.patch.set_facecolor('#1a1a2e')
        for ax in axes:
            ax.set_facecolor('#111111')

        # --- Panel 1: Bolometric ---
        ax = axes[0]
        ax.plot(times, lum, 'o-', color='white', lw=2, ms=6, alpha=0.9, label='Bolometric')

        if show_decay and ni56_mass > 0:
            t_th   = np.linspace(times.min(), times.max(), 300)
            tau_ni = 6.075 / np.log(2)
            tau_co = 77.27 / np.log(2)
            L_th   = np.exp(-t_th / tau_ni) + \
                     (tau_co / (tau_co - tau_ni)) * \
                     (np.exp(-t_th / tau_ni) - np.exp(-t_th / tau_co))
            ax.plot(t_th, lum.max() / L_th.max() * L_th, '--', color='cyan',
                    lw=1.5, alpha=0.7, label=f'Ni-56 decay ({ni56_mass:.2f} M☉)')

        ax.set_ylabel('Bolometric Flux (arb.)', fontsize=12, color='white')
        ax.set_yscale('log')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.2, which='both')
        ax.tick_params(colors='white')

        # --- Panel 2: Color temperature ---
        ax = axes[1]
        valid = np.isfinite(T_color) & (T_color > 0)
        if valid.any():
            ax.plot(times[valid], T_color[valid], 's-', color='orange',
                    lw=2, ms=6, alpha=0.9, label='T$_{color}$ (Wien peak)')
            ax.set_yscale('log')
            ax.legend(fontsize=10)

            # Right axis: Wien wavelength
            ax_r = ax.twinx()
            lam_A = (C_CGS / T_color[valid]) * (2.898e-1 / 1.0) * 1e8  # Wien: lam_peak*T = 0.2898 cm*K
            ax_r.plot(times[valid], lam_A, alpha=0)
            ax_r.set_ylabel('Wien peak λ (Å)', fontsize=10, color='orange')
            ax_r.tick_params(axis='y', labelcolor='orange')
            ax_r.set_yscale('log')
        else:
            ax.text(0.5, 0.5, 'No valid T_color', ha='center', va='center',
                    transform=ax.transAxes, color='gray')

        ax.set_ylabel('Color Temperature (K)', fontsize=12, color='white')
        ax.grid(True, alpha=0.2, which='both')
        ax.tick_params(colors='white')

        # --- Panel 3: Band light curves ---
        ax = axes[2]
        for band_name, color in BAND_COLORS.items():
            band_lum = np.array(self.band_luminosities[band_name])[order]
            if band_lum.max() > 0:
                ax.plot(times, band_lum, 'o-', color=color, lw=2, ms=5,
                        alpha=0.85, label=band_name)

        ax.set_xlabel('Time Since Explosion (days)', fontsize=12, color='white')
        ax.set_ylabel('Band Flux (arb.)', fontsize=12, color='white')
        ax.set_yscale('log')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.2, which='both')
        ax.tick_params(colors='white')

        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plot_path = os.path.join(self.output_dir, output_file)
        plt.savefig(plot_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
        print(f"Light curve plot saved to {plot_path}")
        plt.close()

    def plot_spectra_evolution(self, output_file='spectra_evolution.png',
                               max_spectra=20):
        spectrum_files = sorted(glob(os.path.join(self.output_dir, 'spectrum_*.txt')))
        if not spectrum_files:
            print("No spectrum files found for evolution plot")
            return

        if len(spectrum_files) > max_spectra:
            idx = np.linspace(0, len(spectrum_files)-1, max_spectra, dtype=int)
            spectrum_files = [spectrum_files[i] for i in idx]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        fig.patch.set_facecolor('#1a1a2e')
        for ax in (ax1, ax2):
            ax.set_facecolor('#111111')

        colors = plt.cm.plasma(np.linspace(0, 1, len(spectrum_files)))

        for i, spec_file in enumerate(spectrum_files):
            try:
                data    = np.loadtxt(spec_file)
                freq    = data[:, 0]
                wl      = data[:, 1]
                flux_nu  = data[:, 2]
                flux_lam = data[:, 3]
                step    = self.extract_timestep(spec_file)
                label   = f'Step {step}' if i % max(1, len(spectrum_files)//8) == 0 else ''

                ax1.loglog(freq, flux_nu,  color=colors[i], lw=1.2, alpha=0.8, label=label)
                ax2.loglog(wl,   flux_lam, color=colors[i], lw=1.2, alpha=0.8)

                # Tick the Wien peak on the frequency panel
                ipeak = np.argmax(flux_nu)
                ax1.axvline(freq[ipeak], color=colors[i], lw=0.4, alpha=0.25, ls='--')

            except Exception as e:
                print(f"Error reading {spec_file}: {e}")

        # Band shading on wavelength panel
        for band_name, (nu_hi, nu_lo) in BANDS.items():
            lam_lo = C_CGS / nu_hi * 1e8
            lam_hi = C_CGS / nu_lo * 1e8
            ax2.axvspan(lam_lo, lam_hi, alpha=0.08,
                        color=BAND_COLORS[band_name], label=band_name)

        for ax, xlabel, ylabel in [
            (ax1, 'Frequency (Hz)', 'F$_ν$ (arb.)'),
            (ax2, 'Wavelength (Å)', 'F$_λ$ (arb.)'),
        ]:
            ax.set_xlabel(xlabel, fontsize=12, color='white')
            ax.set_ylabel(ylabel, fontsize=12, color='white')
            ax.grid(True, alpha=0.2)
            ax.tick_params(colors='white')
            ax.legend(fontsize=8, ncol=3)

        ax1.set_title('Spectral Evolution', fontsize=14, color='white')
        ax2.invert_xaxis()

        sm = plt.cm.ScalarMappable(cmap='plasma',
                                   norm=plt.Normalize(0, len(spectrum_files)-1))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=[ax1, ax2], label='Time progression →')
        cbar.ax.yaxis.label.set_color('white')
        cbar.ax.tick_params(colors='white')

        plt.tight_layout()
        plot_path = os.path.join(self.output_dir, output_file)
        plt.savefig(plot_path, dpi=150, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
        print(f"Spectral evolution plot saved to {plot_path}")
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate supernova light curve from VTK files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python make_light_curve.py output/ -mass 13.0 -radius 6.41
  python make_light_curve.py output/ -ni56 0.6 -dt 0.1 -octant -nrays 15
  python make_light_curve.py output/ -octant -nrays 10 -workers 4
  python make_light_curve.py output/ -octant -skip
        """
    )

    parser.add_argument('vtk_dir')
    parser.add_argument('-output', default='lightcurve_output')

    parser.add_argument('-dt', '--time_per_step', type=float)
    parser.add_argument('-t0', type=float, default=0.0)

    parser.add_argument('-ni56', type=float, default=0.0)
    parser.add_argument('-ni56_radius', type=float, default=0.3)

    parser.add_argument('-opacity', default='kramers',
                        choices=['kramers', 'freefree', 'gray'])
    parser.add_argument('-kappa', '--opacity_multiplier', type=float, default=1.0)
    parser.add_argument('-octant', '--octant_mode', action='store_true')
    parser.add_argument('-nrays', '--nrays_octant', type=int, default=10)
    parser.add_argument('-nfreq', type=int, default=100)
    # Corrected defaults for 26 Msun star (1 code unit = 13 Msun, 6.41 Rsun)
    parser.add_argument('-mass',   type=float, default=13.0,
                        help='Code mass unit in Msun (default: 13.0)')
    parser.add_argument('-radius', type=float, default=6.41,
                        help='Code length unit in Rsun (default: 6.41)')
    parser.add_argument('-axis', default='z', choices=['x', 'y', 'z'])
    parser.add_argument('-doppler', action='store_true')

    parser.add_argument('-workers', type=int, default=1,
                        help='Parallel workers (default: 1; -1 = all CPUs)')
    parser.add_argument('-skip', '--skip_existing', action='store_true')
    parser.add_argument('-no_plot', action='store_true')
    parser.add_argument('-show_decay', action='store_true')

    args = parser.parse_args()

    if not os.path.exists('./radtransfer'):
        print("ERROR: radtransfer not found. Run 'make' first.")
        return 1

    nworkers = cpu_count() if args.workers == -1 else args.workers

    lc = LightCurveGenerator(args.vtk_dir, args.output)

    lc.process_all_vtk_files(
        time_per_step=args.time_per_step, t0=args.t0,
        ni56_mass=args.ni56, ni56_radius=args.ni56_radius,
        opacity=args.opacity, opacity_multiplier=args.opacity_multiplier,
        octant_mode=args.octant_mode, nrays_octant=args.nrays_octant,
        nfreq=args.nfreq, star_mass=args.mass, star_radius=args.radius,
        axis=args.axis, doppler=args.doppler,
        skip_existing=args.skip_existing, nworkers=nworkers,
    )

    if not lc.times:
        print("\nNo data collected. Exiting.")
        return 1

    lc.save_light_curve('lightcurve.txt')

    if not args.no_plot:
        lc.plot_light_curve('lightcurve.png',
                            ni56_mass=args.ni56, show_decay=args.show_decay)
        lc.plot_spectra_evolution('spectra_evolution.png')

    print(f"\n{'='*60}\nDONE\n{'='*60}")
    print(f"\nResults in {args.output}/:")
    print(f"  lightcurve.txt        (bol + bands + T_color)")
    if not args.no_plot:
        print(f"  lightcurve.png        (3-panel: bol / T_color / bands)")
        print(f"  spectra_evolution.png")
    print(f"  spectrum_*.txt\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())

