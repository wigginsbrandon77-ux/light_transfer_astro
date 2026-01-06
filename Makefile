# Makefile for stellar radiative transfer code

# Compiler and flags
CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c99
LDFLAGS = -lm

# Target executable
TARGET = radtransfer

# Source files
SOURCES = main.c radtransfer.c
HEADERS = radtransfer.h
OBJECTS = $(SOURCES:.c=.o)

# Default target
all: $(TARGET)
	@echo ""
	@echo "=========================================="
	@echo "Build complete!"
	@echo "=========================================="
	@echo "Run with: ./$(TARGET) <vtk_file>"
	@echo "Example:  ./$(TARGET) output_000050.vtk"
	@echo ""

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile source files to object files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)
	@echo "Cleaned build artifacts"

# Clean and rebuild
rebuild: clean all

# Install (optional - copies to /usr/local/bin)
install: $(TARGET)
	@echo "Installing $(TARGET) to /usr/local/bin..."
	@sudo cp $(TARGET) /usr/local/bin/
	@echo "Installation complete. You can now run '$(TARGET)' from anywhere."

# Uninstall
uninstall:
	@echo "Removing $(TARGET) from /usr/local/bin..."
	@sudo rm -f /usr/local/bin/$(TARGET)
	@echo "Uninstallation complete."

# Help
help:
	@echo "Available targets:"
	@echo "  make         - Build the radiative transfer code"
	@echo "  make clean   - Remove build artifacts"
	@echo "  make rebuild - Clean and rebuild"
	@echo "  make install - Install to /usr/local/bin (requires sudo)"
	@echo "  make help    - Show this help message"
	@echo ""
	@echo "Usage after building:"
	@echo "  ./$(TARGET) <vtk_file> [output_file] [options]"
	@echo ""
	@echo "Options:"
	@echo "  -axis <x|y|z>       Viewing axis (default: z)"
	@echo "  -nfreq <int>        Number of frequency bins (default: 50)"
	@echo "  -mass <float>       Stellar mass in M_sun (default: 8.0)"
	@echo "  -radius <float>     Stellar radius in R_sun (default: 7.0)"
	@echo "  -doppler            Include Doppler shifting"
	@echo ""
	@echo "Example:"
	@echo "  ./$(TARGET) output_000050.vtk spectrum.txt -axis z -nfreq 50 -doppler"

.PHONY: all clean rebuild install uninstall help
