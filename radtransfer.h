/*
 * radtransfer.h
 * Header file for stellar radiative transfer calculations
 */

#ifndef RADTRANSFER_H
#define RADTRANSFER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Physical constants (CGS units) */
#define C_CGS 2.99792458e10     /* Speed of light (cm/s) */
#define H_CGS 6.62607015e-27    /* Planck constant (erg⋅s) */
#define K_CGS 1.380649e-16      /* Boltzmann constant (erg/K) */
#define M_P_CGS 1.67262192e-24  /* Proton mass (g) */
#define G_CGS 6.67430e-8        /* Gravitational constant (cm³/g/s²) */
#define SIGMA_CGS 5.670374e-5   /* Stefan-Boltzmann (erg/cm²/s/K⁴) */
#define M_SUN_CGS 1.98847e33    /* Solar mass (g) */
#define R_SUN_CGS 6.96e10       /* Solar radius (cm) */

/* Opacity parameters */
#define KAPPA_0 4.0e24          /* Kramers opacity normalization */
#define MU 0.6                  /* Mean molecular weight (ionized gas) */

/* Grid and frequency parameters */
#define MAX_GRID_SIZE 256       /* Maximum grid dimension */
#define MAX_FREQ 200            /* Maximum number of frequency bins */

/* Data structures */
typedef struct {
    int nx, ny, nz;             /* Grid dimensions */
    double xmin, ymin, zmin;    /* Domain origin */
    double dx, dy, dz;          /* Cell spacing */
    
    /* Data arrays (3D) */
    double ***density;          /* Density (code units) */
    double ***pressure;         /* Pressure (code units) */
    double ***vx, ***vy, ***vz; /* Velocity components (code units) */
    
    /* Physical arrays (CGS) */
    double ***density_cgs;      /* Density (g/cm³) */
    double ***pressure_cgs;     /* Pressure (erg/cm³) */
    double ***temperature;      /* Temperature (K) */
    double ***vx_cgs, ***vy_cgs, ***vz_cgs; /* Velocity (cm/s) */
    
    /* Unit conversions */
    double L_unit;              /* Length unit (cm) */
    double M_unit;              /* Mass unit (g) */
    double T_unit;              /* Time unit (s) */
    double rho_unit;            /* Density unit (g/cm³) */
    double P_unit;              /* Pressure unit (erg/cm³) */
    double v_unit;              /* Velocity unit (cm/s) */
    
    /* Stellar parameters */
    double star_mass;           /* Stellar mass (solar masses) */
    double star_radius;         /* Stellar radius (solar radii) */
} StellarData;

typedef struct {
    int nfreq;                  /* Number of frequencies */
    double freq_min;            /* Minimum frequency (Hz) */
    double freq_max;            /* Maximum frequency (Hz) */
    double *frequencies;        /* Frequency array (Hz) */
    double *spectrum;           /* Integrated spectrum (arbitrary units) */
    
    /* Ray tracing parameters */
    char axis;                  /* Viewing axis ('x', 'y', or 'z') */
    int include_doppler;        /* Include Doppler shifting (0 or 1) */
} SpectrumData;

/* Function declarations */

/* Memory management */
double*** alloc_3d_array(int nx, int ny, int nz);
void free_3d_array(double ***arr, int nx, int ny);

/* VTK I/O */
int read_vtk_file(const char *filename, StellarData *data);

/* Unit conversions */
void setup_units(StellarData *data);
void convert_to_physical_units(StellarData *data);

/* Physics */
double planck_function(double nu, double T);
double opacity_kramers(double rho, double T);
double get_cell_value_3d(double ***arr, int i, int j, int k, 
                         int nx, int ny, int nz, char axis);

/* Ray tracing */
void raytrace_spectrum(StellarData *data, SpectrumData *spec);

/* Output */
void write_spectrum(const char *filename, StellarData *data, SpectrumData *spec);
void print_summary(StellarData *data);

#endif /* RADTRANSFER_H */
