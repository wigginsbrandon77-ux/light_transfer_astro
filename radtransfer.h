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

/* Nickel-56 decay parameters */
#define NI56_HALFLIFE 6.075     /* Ni-56 half-life (days) */
#define CO56_HALFLIFE 77.27     /* Co-56 half-life (days) */
#define NI56_ENERGY 1.75e49     /* Energy per solar mass of Ni-56 decay (erg/M_sun) */
#define CO56_ENERGY 3.61e49     /* Energy per solar mass of Co-56 decay (erg/M_sun) */
#define SECONDS_PER_DAY 86400.0 /* Seconds per day */

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
    
    /* Ni-56 decay parameters */
    int enable_ni56;            /* Enable Ni-56 decay heating (0 or 1) */
    double ni56_mass;           /* Initial Ni-56 mass (solar masses) */
    double ni56_radius;         /* Radius containing Ni-56 (code units) */
    double time_since_explosion;/* Time since explosion (seconds) */
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
    char opacity_type[20];      /* Opacity model: "kramers", "freefree", "gray" */
    double opacity_multiplier;  /* Opacity boost factor (default: 1.0) */
    int octant_mode;            /* Use octant geometry (dense core at origin) */
    int nrays_octant;           /* Number of radial rays for octant mode */
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
double opacity_kramers(double rho, double T, double multiplier);
double opacity_freefree(double rho, double T, double nu, double multiplier);
double opacity_gray(double rho, double T, double multiplier);
double get_cell_value_3d(double ***arr, int i, int j, int k, 
                         int nx, int ny, int nz, char axis);

/* Ni-56 decay heating */
double ni56_luminosity(double ni56_mass_msun, double time_days);
double ni56_heating_rate(StellarData *data, int i, int j, int k);
void apply_ni56_heating(StellarData *data);

/* Ray tracing */
void raytrace_spectrum(StellarData *data, SpectrumData *spec);

/* Output */
void write_spectrum(const char *filename, StellarData *data, SpectrumData *spec);
void print_summary(StellarData *data);

#endif /* RADTRANSFER_H */
