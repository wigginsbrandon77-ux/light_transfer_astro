/*
 * main.c
 * Main program for radiative transfer calculation
 * 
 * Usage: ./radtransfer <vtk_file> [output_file] [options]
 * 
 * Options:
 *   -axis <x|y|z>          Viewing axis (default: z)
 *   -nfreq <int>           Number of frequency bins (default: 50)
 *   -mass <float>          Stellar mass in M_sun (default: 8.0)
 *   -radius <float>        Stellar radius in R_sun (default: 7.0)
 *   -doppler               Include Doppler shifting (default: no)
 */

#include "radtransfer.h"

void print_usage(const char *prog_name) {
    printf("\nUsage: %s <vtk_file> [output_file] [options]\n\n", prog_name);
    printf("Options:\n");
    printf("  -axis <x|y|z>       Viewing axis (default: z)\n");
    printf("  -nfreq <int>        Number of frequency bins (default: 50)\n");
    printf("  -mass <float>       Stellar mass in M_sun (default: 8.0)\n");
    printf("  -radius <float>     Stellar radius in R_sun (default: 7.0)\n");
    printf("  -opacity <type>     Opacity model: kramers, freefree, gray (default: kramers)\n");
    printf("  -doppler            Include Doppler shifting (default: off)\n");
    printf("\n");
    printf("Ni-56 decay heating:\n");
    printf("  -ni56 <mass>        Enable Ni-56 heating with given mass in M_sun\n");
    printf("  -ni56_radius <r>    Radius containing Ni-56 in code units (default: 0.3)\n");
    printf("  -time <days>        Time since explosion in days (default: 10)\n");
    printf("\nOpacity models:\n");
    printf("  kramers   - Bound-free absorption, κ ∝ ρ T^-3.5 (frequency independent)\n");
    printf("  freefree  - Free-free absorption, κ ∝ ρ T^-3.5 ν^-3 (shows frequency dependence)\n");
    printf("  gray      - Simple gray opacity (best for seeing blackbody peak)\n");
    printf("\nExample:\n");
    printf("  %s output_000050.vtk spectrum.txt -opacity gray -nfreq 100\n", 
           prog_name);
    printf("  %s output_000050.vtk spectrum.txt -ni56 0.6 -time 20\n\n", 
           prog_name);
}

int main(int argc, char *argv[]) {
    printf("========================================\n");
    printf("STELLAR RADIATIVE TRANSFER CALCULATOR\n");
    printf("========================================\n");
    
    /* Check arguments */
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    /* Default parameters */
    char *vtk_file = argv[1];
    char *output_file = "spectrum.txt";
    char axis = 'z';
    int nfreq = 50;
    double star_mass = 8.0;
    double star_radius = 7.0;
    int include_doppler = 0;
    char opacity_type[20] = "kramers";
    
    /* Ni-56 parameters */
    int enable_ni56 = 0;
    double ni56_mass = 0.0;
    double ni56_radius = 0.3;  /* Code units */
    double time_days = 10.0;   /* Days since explosion */
    
    /* Parse optional arguments */
    if (argc >= 3 && argv[2][0] != '-') {
        output_file = argv[2];
    }
    
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-axis") == 0 && i + 1 < argc) {
            axis = argv[++i][0];
        } else if (strcmp(argv[i], "-nfreq") == 0 && i + 1 < argc) {
            nfreq = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-mass") == 0 && i + 1 < argc) {
            star_mass = atof(argv[++i]);
        } else if (strcmp(argv[i], "-radius") == 0 && i + 1 < argc) {
            star_radius = atof(argv[++i]);
        } else if (strcmp(argv[i], "-opacity") == 0 && i + 1 < argc) {
            strncpy(opacity_type, argv[++i], sizeof(opacity_type) - 1);
        } else if (strcmp(argv[i], "-ni56") == 0 && i + 1 < argc) {
            enable_ni56 = 1;
            ni56_mass = atof(argv[++i]);
        } else if (strcmp(argv[i], "-ni56_radius") == 0 && i + 1 < argc) {
            ni56_radius = atof(argv[++i]);
        } else if (strcmp(argv[i], "-time") == 0 && i + 1 < argc) {
            time_days = atof(argv[++i]);
        } else if (strcmp(argv[i], "-doppler") == 0) {
            include_doppler = 1;
        }
    }
    
    /* Validate parameters */
    if (axis != 'x' && axis != 'y' && axis != 'z') {
        fprintf(stderr, "Error: axis must be 'x', 'y', or 'z'\n");
        return 1;
    }
    
    if (nfreq < 10 || nfreq > MAX_FREQ) {
        fprintf(stderr, "Error: nfreq must be between 10 and %d\n", MAX_FREQ);
        return 1;
    }
    
    if (strcmp(opacity_type, "kramers") != 0 && 
        strcmp(opacity_type, "freefree") != 0 && 
        strcmp(opacity_type, "gray") != 0) {
        fprintf(stderr, "Error: opacity must be 'kramers', 'freefree', or 'gray'\n");
        return 1;
    }
    
    /* Initialize stellar data structure */
    StellarData data;
    data.star_mass = star_mass;
    data.star_radius = star_radius;
    data.enable_ni56 = enable_ni56;
    data.ni56_mass = ni56_mass;
    data.ni56_radius = ni56_radius;
    data.time_since_explosion = time_days * SECONDS_PER_DAY;
    
    /* Set up units */
    setup_units(&data);
    
    /* Read VTK file */
    if (read_vtk_file(vtk_file, &data) != 0) {
        return 1;
    }
    
    /* Convert to physical units */
    convert_to_physical_units(&data);
    
    /* Apply Ni-56 heating if enabled */
    if (enable_ni56) {
        apply_ni56_heating(&data);
    }
    
    /* Print summary */
    print_summary(&data);
    
    /* Initialize spectrum structure */
    SpectrumData spec;
    spec.axis = axis;
    spec.nfreq = nfreq;
    spec.freq_min = 1e14;  /* ~30,000 Angstrom (IR) */
    spec.freq_max = 1e16;  /* ~300 Angstrom (UV) */
    spec.include_doppler = include_doppler;
    strncpy(spec.opacity_type, opacity_type, sizeof(spec.opacity_type) - 1);
    
    /* Perform ray tracing */
    raytrace_spectrum(&data, &spec);
    
    /* Write output */
    write_spectrum(output_file, &data, &spec);
    
    /* Cleanup */
    free_3d_array(data.density, data.nx, data.ny);
    free_3d_array(data.pressure, data.nx, data.ny);
    free_3d_array(data.vx, data.nx, data.ny);
    free_3d_array(data.vy, data.nx, data.ny);
    free_3d_array(data.vz, data.nx, data.ny);
    free_3d_array(data.density_cgs, data.nx, data.ny);
    free_3d_array(data.pressure_cgs, data.nx, data.ny);
    free_3d_array(data.temperature, data.nx, data.ny);
    free_3d_array(data.vx_cgs, data.nx, data.ny);
    free_3d_array(data.vy_cgs, data.nx, data.ny);
    free_3d_array(data.vz_cgs, data.nx, data.ny);
    free(spec.frequencies);
    free(spec.spectrum);
    
    printf("\n========================================\n");
    printf("Calculation complete!\n");
    printf("========================================\n\n");
    
    return 0;
}
