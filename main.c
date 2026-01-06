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
    printf("  -doppler            Include Doppler shifting (default: off)\n");
    printf("\nExample:\n");
    printf("  %s output_000050.vtk spectrum.txt -axis z -nfreq 50 -doppler\n\n", 
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
    
    /* Initialize stellar data structure */
    StellarData data;
    data.star_mass = star_mass;
    data.star_radius = star_radius;
    
    /* Set up units */
    setup_units(&data);
    
    /* Read VTK file */
    if (read_vtk_file(vtk_file, &data) != 0) {
        return 1;
    }
    
    /* Convert to physical units */
    convert_to_physical_units(&data);
    
    /* Print summary */
    print_summary(&data);
    
    /* Initialize spectrum structure */
    SpectrumData spec;
    spec.axis = axis;
    spec.nfreq = nfreq;
    spec.freq_min = 1e14;  /* ~30,000 Angstrom (IR) */
    spec.freq_max = 1e16;  /* ~300 Angstrom (UV) */
    spec.include_doppler = include_doppler;
    
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
