/*
 * radtransfer.c
 * Radiative transfer calculations for stellar explosion simulations
 * 
 * Reads VTK structured points data and computes emergent spectrum
 * via simple ray tracing with Kramers opacity
 */

#include "radtransfer.h"

/*
 * Allocate a 3D array
 */
double*** alloc_3d_array(int nx, int ny, int nz) {
    double ***arr = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i < nx; i++) {
        arr[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j < ny; j++) {
            arr[i][j] = (double*)calloc(nz, sizeof(double));
        }
    }
    return arr;
}

/*
 * Free a 3D array
 */
void free_3d_array(double ***arr, int nx, int ny) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            free(arr[i][j]);
        }
        free(arr[i]);
    }
    free(arr);
}

/*
 * Set up code unit conversions
 */
void setup_units(StellarData *data) {
    /* Length unit based on stellar radius */
    data->L_unit = data->star_radius * R_SUN_CGS;
    
    /* Mass unit based on stellar mass */
    data->M_unit = data->star_mass * M_SUN_CGS;
    
    /* Time unit from G=1 in code units */
    data->T_unit = sqrt(pow(data->L_unit, 3.0) / (G_CGS * data->M_unit));
    
    /* Derived units */
    data->rho_unit = data->M_unit / pow(data->L_unit, 3.0);
    data->P_unit = data->M_unit / (data->L_unit * pow(data->T_unit, 2.0));
    data->v_unit = data->L_unit / data->T_unit;
    
    printf("\n=== Code Unit Conversions ===\n");
    printf("Length unit:   %.3e cm (%.2f R_sun)\n", 
           data->L_unit, data->L_unit / R_SUN_CGS);
    printf("Mass unit:     %.3e g (%.2f M_sun)\n", 
           data->M_unit, data->M_unit / M_SUN_CGS);
    printf("Time unit:     %.3e s (%.2f hours)\n", 
           data->T_unit, data->T_unit / 3600.0);
    printf("Density unit:  %.3e g/cm³\n", data->rho_unit);
    printf("Pressure unit: %.3e erg/cm³\n", data->P_unit);
    printf("Velocity unit: %.3e cm/s (%.1f km/s)\n", 
           data->v_unit, data->v_unit / 1e5);
}

/*
 * Read VTK structured points file
 */
int read_vtk_file(const char *filename, StellarData *data) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file %s\n", filename);
        return -1;
    }
    
    printf("\nReading VTK file: %s\n", filename);
    
    char line[256];
    int reading_density = 0, reading_pressure = 0, reading_velocity = 0;
    int cell_count = 0;
    
    /* Parse header */
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "DIMENSIONS", 10) == 0) {
            sscanf(line, "DIMENSIONS %d %d %d", &data->nx, &data->ny, &data->nz);
            printf("Grid: %d × %d × %d\n", data->nx, data->ny, data->nz);
            
            /* Allocate arrays */
            data->density = alloc_3d_array(data->nx, data->ny, data->nz);
            data->pressure = alloc_3d_array(data->nx, data->ny, data->nz);
            data->vx = alloc_3d_array(data->nx, data->ny, data->nz);
            data->vy = alloc_3d_array(data->nx, data->ny, data->nz);
            data->vz = alloc_3d_array(data->nx, data->ny, data->nz);
            
        } else if (strncmp(line, "ORIGIN", 6) == 0) {
            sscanf(line, "ORIGIN %lf %lf %lf", 
                   &data->xmin, &data->ymin, &data->zmin);
                   
        } else if (strncmp(line, "SPACING", 7) == 0) {
            sscanf(line, "SPACING %lf %lf %lf", 
                   &data->dx, &data->dy, &data->dz);
                   
        } else if (strncmp(line, "SCALARS density", 15) == 0) {
            fgets(line, sizeof(line), fp); /* Skip LOOKUP_TABLE line */
            reading_density = 1;
            reading_pressure = 0;
            reading_velocity = 0;
            cell_count = 0;
            
        } else if (strncmp(line, "SCALARS pressure", 16) == 0) {
            fgets(line, sizeof(line), fp); /* Skip LOOKUP_TABLE line */
            reading_density = 0;
            reading_pressure = 1;
            reading_velocity = 0;
            cell_count = 0;
            
        } else if (strncmp(line, "VECTORS velocity", 16) == 0) {
            reading_density = 0;
            reading_pressure = 0;
            reading_velocity = 1;
            cell_count = 0;
            
        } else {
            /* Read data values */
            if (reading_density) {
                double val;
                if (sscanf(line, "%lf", &val) == 1) {
                    int i = cell_count / (data->ny * data->nz);
                    int j = (cell_count % (data->ny * data->nz)) / data->nz;
                    int k = cell_count % data->nz;
                    data->density[i][j][k] = val;
                    cell_count++;
                }
            } else if (reading_pressure) {
                double val;
                if (sscanf(line, "%lf", &val) == 1) {
                    int i = cell_count / (data->ny * data->nz);
                    int j = (cell_count % (data->ny * data->nz)) / data->nz;
                    int k = cell_count % data->nz;
                    data->pressure[i][j][k] = val;
                    cell_count++;
                }
            } else if (reading_velocity) {
                double vx_val, vy_val, vz_val;
                if (sscanf(line, "%lf %lf %lf", &vx_val, &vy_val, &vz_val) == 3) {
                    int i = cell_count / (data->ny * data->nz);
                    int j = (cell_count % (data->ny * data->nz)) / data->nz;
                    int k = cell_count % data->nz;
                    data->vx[i][j][k] = vx_val;
                    data->vy[i][j][k] = vy_val;
                    data->vz[i][j][k] = vz_val;
                    cell_count++;
                }
            }
        }
    }
    
    fclose(fp);
    printf("VTK file read successfully\n");
    return 0;
}

/*
 * Convert code units to physical (CGS) units
 */
void convert_to_physical_units(StellarData *data) {
    /* Allocate physical arrays */
    data->density_cgs = alloc_3d_array(data->nx, data->ny, data->nz);
    data->pressure_cgs = alloc_3d_array(data->nx, data->ny, data->nz);
    data->temperature = alloc_3d_array(data->nx, data->ny, data->nz);
    data->vx_cgs = alloc_3d_array(data->nx, data->ny, data->nz);
    data->vy_cgs = alloc_3d_array(data->nx, data->ny, data->nz);
    data->vz_cgs = alloc_3d_array(data->nx, data->ny, data->nz);
    
    double T_min = 1e100, T_max = 0;
    double rho_min = 1e100, rho_max = 0;
    
    for (int i = 0; i < data->nx; i++) {
        for (int j = 0; j < data->ny; j++) {
            for (int k = 0; k < data->nz; k++) {
                /* Convert to CGS */
                data->density_cgs[i][j][k] = data->density[i][j][k] * data->rho_unit;
                data->pressure_cgs[i][j][k] = data->pressure[i][j][k] * data->P_unit;
                data->vx_cgs[i][j][k] = data->vx[i][j][k] * data->v_unit;
                data->vy_cgs[i][j][k] = data->vy[i][j][k] * data->v_unit;
                data->vz_cgs[i][j][k] = data->vz[i][j][k] * data->v_unit;
                
                /* Calculate temperature from ideal gas: P = (ρ/μm_p) k T */
                double rho_cgs = data->density_cgs[i][j][k];
                double P_cgs = data->pressure_cgs[i][j][k];
                
                if (rho_cgs > 1e-15) {
                    data->temperature[i][j][k] = (P_cgs * MU * M_P_CGS) / (rho_cgs * K_CGS);
                } else {
                    data->temperature[i][j][k] = 1000.0; /* Floor temperature */
                }
                
                /* Track min/max */
                if (data->temperature[i][j][k] > T_max) T_max = data->temperature[i][j][k];
                if (data->temperature[i][j][k] < T_min) T_min = data->temperature[i][j][k];
                if (rho_cgs > rho_max) rho_max = rho_cgs;
                if (rho_cgs > 1e-15 && rho_cgs < rho_min) rho_min = rho_cgs;
            }
        }
    }
    
    printf("\nTemperature range: %.3e - %.3e K\n", T_min, T_max);
    printf("Density range:     %.3e - %.3e g/cm³\n", rho_min, rho_max);
}

/*
 * Planck function B_ν(T) in erg/s/cm²/Hz/sr
 */
double planck_function(double nu, double T) {
    double x = (H_CGS * nu) / (K_CGS * T);
    
    /* Avoid overflow */
    if (x > 700.0) return 0.0;
    
    double numerator = 2.0 * H_CGS * pow(nu, 3.0) / pow(C_CGS, 2.0);
    double denominator = exp(x) - 1.0;
    
    return numerator / denominator;
}

/*
 * Kramers opacity κ = κ₀ ρ T^(-3.5) in cm²/g
 */
double opacity_kramers(double rho, double T) {
    /* Add floors to avoid issues */
    if (rho < 1e-15) rho = 1e-15;
    if (T < 1000.0) T = 1000.0;
    
    return KAPPA_0 * rho * pow(T, -3.5);
}

/*
 * Get cell value with axis permutation
 */
double get_cell_value_3d(double ***arr, int i, int j, int k, 
                         int nx, int ny, int nz, char axis) {
    if (axis == 'z') {
        if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
            return arr[i][j][k];
    } else if (axis == 'y') {
        if (i >= 0 && i < nx && k >= 0 && k < nz && j >= 0 && j < ny)
            return arr[i][j][k];
    } else if (axis == 'x') {
        if (j >= 0 && j < ny && k >= 0 && k < nz && i >= 0 && i < nx)
            return arr[i][j][k];
    }
    return 0.0;
}

/*
 * Perform ray tracing to compute spectrum
 */
void raytrace_spectrum(StellarData *data, SpectrumData *spec) {
    printf("\n=== Ray Tracing Along %c-Axis ===\n", spec->axis);
    printf("Frequency range: %.2e - %.2e Hz\n", spec->freq_min, spec->freq_max);
    printf("Number of frequencies: %d\n", spec->nfreq);
    
    /* Allocate frequency and spectrum arrays */
    spec->frequencies = (double*)malloc(spec->nfreq * sizeof(double));
    spec->spectrum = (double*)calloc(spec->nfreq, sizeof(double));
    
    /* Create log-spaced frequency grid */
    double log_fmin = log10(spec->freq_min);
    double log_fmax = log10(spec->freq_max);
    for (int ifreq = 0; ifreq < spec->nfreq; ifreq++) {
        double log_f = log_fmin + (log_fmax - log_fmin) * ifreq / (spec->nfreq - 1);
        spec->frequencies[ifreq] = pow(10.0, log_f);
    }
    
    /* Determine ray dimensions and step size */
    int nx_ray, ny_ray, nz_ray;
    double ds;
    
    if (spec->axis == 'z') {
        nx_ray = data->nx;
        ny_ray = data->ny;
        nz_ray = data->nz;
        ds = data->dz * data->L_unit;
    } else if (spec->axis == 'y') {
        nx_ray = data->nx;
        ny_ray = data->nz;
        nz_ray = data->ny;
        ds = data->dy * data->L_unit;
    } else { /* 'x' */
        nx_ray = data->ny;
        ny_ray = data->nz;
        nz_ray = data->nx;
        ds = data->dx * data->L_unit;
    }
    
    /* Loop over frequencies */
    for (int ifreq = 0; ifreq < spec->nfreq; ifreq++) {
        double nu = spec->frequencies[ifreq];
        
        if ((ifreq + 1) % 10 == 0 || ifreq == spec->nfreq - 1) {
            printf("Processing frequency %d/%d (%.2e Hz)\n", 
                   ifreq + 1, spec->nfreq, nu);
        }
        
        /* Loop over spatial pixels */
        for (int i = 0; i < nx_ray; i++) {
            for (int j = 0; j < ny_ray; j++) {
                double I_nu = 0.0; /* Start with zero background */
                
                /* Integrate along ray (from back to front) */
                for (int k = nz_ray - 1; k >= 0; k--) {
                    /* Get cell values */
                    double rho_cell, T_cell, v_cell;
                    
                    if (spec->axis == 'z') {
                        rho_cell = data->density_cgs[i][j][k];
                        T_cell = data->temperature[i][j][k];
                        v_cell = data->vz_cgs[i][j][k];
                    } else if (spec->axis == 'y') {
                        rho_cell = data->density_cgs[i][k][j];
                        T_cell = data->temperature[i][k][j];
                        v_cell = data->vy_cgs[i][k][j];
                    } else { /* 'x' */
                        rho_cell = data->density_cgs[k][i][j];
                        T_cell = data->temperature[k][i][j];
                        v_cell = data->vx_cgs[k][i][j];
                    }
                    
                    /* Skip low-density regions */
                    if (rho_cell < 1e-15) continue;
                    
                    /* Apply Doppler shift if requested */
                    double nu_emit = nu;
                    if (spec->include_doppler) {
                        nu_emit = nu * (1.0 + v_cell / C_CGS);
                    }
                    
                    /* Source function (LTE: S_ν = B_ν) */
                    double B_nu = planck_function(nu_emit, T_cell);
                    
                    /* Opacity */
                    double kappa = opacity_kramers(rho_cell, T_cell);
                    
                    /* Optical depth through cell */
                    double tau = kappa * rho_cell * ds;
                    
                    /* Formal solution: I_out = I_in * exp(-tau) + B_nu * (1 - exp(-tau)) */
                    if (tau < 1e-4) {
                        /* Optically thin limit */
                        I_nu += B_nu * tau;
                    } else {
                        double exp_tau = exp(-tau);
                        I_nu = I_nu * exp_tau + B_nu * (1.0 - exp_tau);
                    }
                }
                
                /* Add to integrated spectrum */
                spec->spectrum[ifreq] += I_nu;
            }
        }
    }
    
    printf("Ray tracing complete!\n");
}

/*
 * Write spectrum to text file with detailed header
 */
void write_spectrum(const char *filename, StellarData *data, SpectrumData *spec) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not open %s for writing\n", filename);
        return;
    }
    
    /* Write comprehensive header */
    fprintf(fp, "# ========================================================================\n");
    fprintf(fp, "# RADIATIVE TRANSFER SPECTRUM OUTPUT\n");
    fprintf(fp, "# ========================================================================\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# This file contains the emergent spectrum from a stellar explosion\n");
    fprintf(fp, "# simulation computed via ray tracing with radiative transfer.\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# PHYSICAL PARAMETERS:\n");
    fprintf(fp, "#   Stellar mass:     %.2f M_sun\n", data->star_mass);
    fprintf(fp, "#   Stellar radius:   %.2f R_sun\n", data->star_radius);
    fprintf(fp, "#   Grid dimensions:  %d × %d × %d\n", data->nx, data->ny, data->nz);
    fprintf(fp, "#   Viewing axis:     %c-axis\n", spec->axis);
    fprintf(fp, "#   Doppler shifting: %s\n", spec->include_doppler ? "Yes" : "No");
    fprintf(fp, "#   Opacity model:    Kramers (bound-free)\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# CODE UNITS:\n");
    fprintf(fp, "#   Length unit:   %.3e cm\n", data->L_unit);
    fprintf(fp, "#   Mass unit:     %.3e g\n", data->M_unit);
    fprintf(fp, "#   Time unit:     %.3e s\n", data->T_unit);
    fprintf(fp, "#   Density unit:  %.3e g/cm³\n", data->rho_unit);
    fprintf(fp, "#   Pressure unit: %.3e erg/cm³\n", data->P_unit);
    fprintf(fp, "#\n");
    fprintf(fp, "# SPECTRUM DATA:\n");
    fprintf(fp, "#   Number of frequency bins: %d\n", spec->nfreq);
    fprintf(fp, "#   Frequency range: %.3e - %.3e Hz\n", spec->freq_min, spec->freq_max);
    fprintf(fp, "#   Wavelength range: %.1f - %.1f Angstroms\n", 
            C_CGS / spec->freq_max * 1e8, C_CGS / spec->freq_min * 1e8);
    fprintf(fp, "#\n");
    fprintf(fp, "# COLUMN DESCRIPTIONS:\n");
    fprintf(fp, "#   1. Frequency (Hz) - Photon frequency in Hertz\n");
    fprintf(fp, "#   2. Wavelength (Angstrom) - Photon wavelength (1 Å = 10^-8 cm)\n");
    fprintf(fp, "#   3. Flux_nu (arbitrary) - Flux density F_ν in arbitrary units\n");
    fprintf(fp, "#                             Integrated over spatial pixels\n");
    fprintf(fp, "#                             Relative normalization only\n");
    fprintf(fp, "#   4. Flux_lambda (arbitrary) - Flux density F_λ in arbitrary units\n");
    fprintf(fp, "#                                 F_λ = F_ν × (c/λ²)\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# PLOTTING SUGGESTIONS:\n");
    fprintf(fp, "#   - For F_ν vs ν:      Plot column 3 vs column 1 (log-log scale)\n");
    fprintf(fp, "#   - For F_λ vs λ:      Plot column 4 vs column 2 (log-log scale)\n");
    fprintf(fp, "#   - Wavelength bands:  UV (< 3000 Å), Optical (3000-7000 Å),\n");
    fprintf(fp, "#                        IR (> 7000 Å)\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# EXAMPLE PLOTS (copy into Excel, Python, or ask an LLM):\n");
    fprintf(fp, "#   Excel: Select columns 1 and 3, insert scatter plot, set log scales\n");
    fprintf(fp, "#   Python: import numpy as np; import matplotlib.pyplot as plt\n");
    fprintf(fp, "#           data = np.loadtxt('spectrum.txt')\n");
    fprintf(fp, "#           plt.loglog(data[:,0], data[:,2])\n");
    fprintf(fp, "#           plt.xlabel('Frequency (Hz)')\n");
    fprintf(fp, "#           plt.ylabel('Flux Density')\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# ========================================================================\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# %18s \t %18s \t %18s \t %18s\n", 
            "Frequency_Hz", "Wavelength_Ang", "Flux_nu", "Flux_lambda");
    fprintf(fp, "#\n");
    
    /* Write data */
    for (int i = 0; i < spec->nfreq; i++) {
        double freq = spec->frequencies[i];
        double wavelength = C_CGS / freq * 1e8; /* Angstroms */
        double flux_nu = spec->spectrum[i];
        double flux_lambda = flux_nu * C_CGS / (wavelength * wavelength * 1e16);
        
        fprintf(fp, "  %18.10e \t %18.10e \t %18.10e \t  %18.10e\n",
                freq, wavelength, flux_nu, flux_lambda);
    }
    
    fclose(fp);
    printf("\nSpectrum saved to %s\n", filename);
}

/*
 * Print summary statistics
 */
void print_summary(StellarData *data) {
    double T_min = 1e100, T_max = 0;
    double rho_min = 1e100, rho_max = 0;
    double v_max = 0;
    
    for (int i = 0; i < data->nx; i++) {
        for (int j = 0; j < data->ny; j++) {
            for (int k = 0; k < data->nz; k++) {
                double T = data->temperature[i][j][k];
                double rho = data->density_cgs[i][j][k];
                double v = sqrt(pow(data->vx_cgs[i][j][k], 2.0) + 
                               pow(data->vy_cgs[i][j][k], 2.0) + 
                               pow(data->vz_cgs[i][j][k], 2.0));
                
                if (T > T_max) T_max = T;
                if (T < T_min) T_min = T;
                if (rho > 1e-15 && rho < rho_min) rho_min = rho;
                if (rho > rho_max) rho_max = rho;
                if (v > v_max) v_max = v;
            }
        }
    }
    
    printf("\n=== Summary Statistics ===\n");
    printf("Temperature: %.3e - %.3e K\n", T_min, T_max);
    printf("Density:     %.3e - %.3e g/cm³\n", rho_min, rho_max);
    printf("Max velocity: %.2f km/s\n", v_max / 1e5);
}
