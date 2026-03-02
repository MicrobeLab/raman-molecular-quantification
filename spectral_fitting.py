import numpy as np
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='/path/to/input_spec.txt (two column format)')
parser.add_argument('--ref', type=str, help='/path/to/ref_spec.txt (one compound spectrum per column)')
parser.add_argument('--ref-list', type=str, help='/path/to/ref_list.txt (one compound name per line)')
parser.add_argument('--plot', action='store_true', help='plot original vs fit spectrum')
parser.add_argument('--quant_out', type=str, help='/path/to/quantification_output.txt)', default=None)
parser.add_argument('--residual_out', type=str, help='/path/to/residual_output.txt', default=None)
parser.add_argument('--start_wn', type=int, help='start wn of cell spectrum to cut', default=None)
parser.add_argument('--end_wn', type=int, help='end wn of cell spectrum to cut', default=None)

args = parser.parse_args()

reference_spectra = np.loadtxt(args.ref)

reference_spectra_with_intercept = np.hstack([
    reference_spectra, 
    np.ones((reference_spectra.shape[0], 1))
])

complex_spectrum_full = np.loadtxt(args.input)
if args.start_wn is not None and args.end_wn is not None:
    mask = (complex_spectrum_full[:, 0] >= args.start_wn) & (complex_spectrum_full[:, 0] <= args.end_wn)
    complex_spectrum_cut = complex_spectrum_full[mask]
    wn_values = complex_spectrum_cut[:, 0]
    complex_spectrum = complex_spectrum_cut[:, 1]
else:
    wn_values = complex_spectrum_full[:, 0]
    complex_spectrum = complex_spectrum_full[:, 1]

coefficients, residual_norm = nnls(reference_spectra_with_intercept, complex_spectrum)

fitted_spectrum = np.dot(coefficients[:-1], reference_spectra.T) + coefficients[-1]

residual = complex_spectrum - fitted_spectrum

print("Coefficients (reference spectra):", coefficients[:-1])
print("Intercept (constant term):", coefficients[-1])

if args.quant_out is not None:
    if args.ref_list is not None:
        with open(args.ref_list, 'r') as f:
            compounds = [line.strip() for line in f if line.strip()]
        
        compounds_with_intercept = compounds + ['intercept']
        
        with open(args.quant_out, 'w') as fh:
            fh.write('# Compound\tCoefficient\n')
            for compound, coeff in zip(compounds_with_intercept, coefficients):
                fh.write(f'{compound}\t{coeff}\n')
    else:
        with open(args.quant_out, 'w') as fh:
            for i in coefficients:
                fh.write(str(i) + '\n')

if args.residual_out is not None:
    with open(args.residual_out, 'w') as fh:
        fh.write(f"# Residual components: Original - Fitted spectrum\n")
        fh.write(f"# Original spectrum: {args.input}\n")
        fh.write(f"# Reference spectra: {args.ref}\n")
        fh.write(f"# Wavenumber range: {args.start_wn if args.start_wn else 'min'} to {args.end_wn if args.end_wn else 'max'}\n")
        
        residual_norm_value = np.linalg.norm(residual)
        relative_residual_norm = residual_norm_value / np.linalg.norm(complex_spectrum)
        residual_mean = np.mean(residual)
        residual_std = np.std(residual)
        r_squared = 1 - np.sum(residual**2) / np.sum((complex_spectrum - np.mean(complex_spectrum))**2)
        
        fh.write(f"# Residual norm: {residual_norm_value:.6f}\n")
        fh.write(f"# Relative residual norm: {relative_residual_norm:.6f}\n")
        fh.write(f"# Residual mean: {residual_mean:.6f}\n")
        fh.write(f"# Residual std: {residual_std:.6f}\n")
        fh.write(f"# R-squared: {r_squared:.6f}\n")
        
        fh.write("# Wavenumber\tResidual\tOriginal\tFitted\n")
        
        for i in range(len(wn_values)):
            fh.write(f"{wn_values[i]:.2f}\t{residual[i]:.8f}\t{complex_spectrum[i]:.8f}\t{fitted_spectrum[i]:.8f}\n")
    
    print(f"\nResidual Statistics:")
    print(f"  Residual norm: {residual_norm_value:.6f}")
    print(f"  Relative residual norm: {relative_residual_norm:.6f}")
    print(f"  Residual mean: {residual_mean:.6f}")
    print(f"  Residual std: {residual_std:.6f}")
    print(f"  R-squared: {r_squared:.6f}")
    print(f"  Residual data saved to: {args.residual_out}")

if args.plot:
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 1, 1)
    plt.plot(wn_values, complex_spectrum, label='Complex Spectrum', linewidth=2)
    plt.plot(wn_values, fitted_spectrum, label='Fitted Spectrum', linestyle='--', linewidth=2)
    plt.legend(loc='best')
    plt.title('Complex Spectrum and Fitted Spectrum (with intercept)')
    plt.xlabel('Wavenumber')
    plt.ylabel('Intensity')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 1, 2)
    plt.plot(wn_values, residual, label='Residual (Original - Fitted)', color='red', linewidth=1.5)
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    if args.residual_out is not None:
        plt.text(0.02, 0.95, f'R² = {r_squared:.4f}', transform=plt.gca().transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.text(0.02, 0.85, f'Mean = {residual_mean:.4f}', transform=plt.gca().transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.text(0.02, 0.75, f'Std = {residual_std:.4f}', transform=plt.gca().transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.legend(loc='best')
    plt.title('Residual Components')
    plt.xlabel('Wavenumber')
    plt.ylabel('Residual Intensity')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
