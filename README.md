# IIR-FILTER-DESIGN

# EXP 3 B: DESIGN OF LOW PASS CHEBYSHEV IIR FILTER USING BILINEAR TRANSFORMATION

# AIM: 

To a design of low pass Chebyshev IIR filter using Bilinear Transformation.

# APPARATUS REQUIRED: 
PC installed with SCILAB. 

# PROGRAM: 
```
clc;
close;

// --- 1. GIVEN SPECIFICATIONS (Linear Gain) ---
// We are using the parameters from the previous problem
wp = 0.35 * %pi; // Digital Passband Frequency
ws = 0.7 * %pi;  // Digital Stopband Frequency
Ap_lin = 0.6;    // Passband Gain (Linear)
As_lin = 0.1;    // Stopband Gain (Linear)
T = 1;           // Sampling Time

// Chebyshev formulas require Attenuation in dB, so we convert them:
alphap = -20 * log10(Ap_lin);
alphas = -20 * log10(As_lin);

disp(alphap, 'Passband attenuation (dB) =');
disp(alphas, 'Stopband attenuation (dB) =');

// --- 2. PRE-WARPING (Bilinear Transformation) ---
omegap = (2/T) * tan(wp/2);
omegas = (2/T) * tan(ws/2);

disp(omegap, 'Analog Passband Freq (Omega_p) =');
disp(omegas, 'Analog Stopband Freq (Omega_s) =');

// --- 3. CALCULATE FILTER ORDER (N) ---
num_term = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1)));
den_term = acosh(omegas/omegap);
N_exact = num_term / den_term;

disp(N_exact, 'Calculated Exact Order N =');
N = ceil(N_exact);
disp(N, 'Round off value of N =');

// --- 4. EXPLICIT FORMULAS FROM YOUR NOTES (Epsilon, mu, a, b) ---
Epsilon = sqrt((10^(0.1*alphap)) - 1);
disp(Epsilon, 'Epsilon =');

mu = (1/Epsilon) + sqrt((1/(Epsilon^2)) + 1);
disp(mu, 'mu =');

// Calculate major and minor axes of the ellipse
a = omegap * ((mu^(1/N) - mu^(-1/N)) / 2);
b = omegap * ((mu^(1/N) + mu^(-1/N)) / 2);
disp(a, 'Minor axis (a) =');
disp(b, 'Major axis (b) =');

// --- 5. FIND THE POLES & DENOMINATOR ---
s = poly(0, 's');
den_s = 1; // Initialize denominator

disp('Calculating Poles:');
for k = 1:N
    // Formula for angles from your notes: phi_k = pi/2 + (2k-1)*pi / 2N
    phi_k = (%pi/2) + ((2*k - 1)*%pi) / (2*N);
    
    // Formula for poles: S_k = a*cos(phi) + j*b*sin(phi)
    S_k = a * cos(phi_k) + %i * b * sin(phi_k);
    disp(S_k, 'Pole S_' + string(k) + ' =');
    
    // Multiply out the roots to form the Denominator polynomial
    den_s = den_s * (s - S_k);
end

// Clean up microscopic imaginary residuals from floating-point math
den_s = real(den_s); 

// --- 6. FORM THE NUMERATOR ---
// Using the even/odd N shortcut rule from your handwritten notes!
x = horner(den_s, 0); // Substitute S = 0 into the denominator

if modulo(N, 2) == 0 then
    // If N is even, divide by sqrt(1 + epsilon^2)
    num_s = x / sqrt(1 + Epsilon^2);
else
    // If N is odd, Numerator is just x
    num_s = x;
end

hs = num_s / den_s;
disp('Analog Chebyshev LPF Transfer function H(S) =');
disp(hs);

// --- 7. BILINEAR TRANSFORMATION (S -> Z) ---
z = poly(0, 'z');

// Apply horner() to substitute S = (2/T) * ((z-1)/(z+1))
Hz = horner(hs, (2/T)*((z - 1)/(z + 1))); 
Hz = clean(Hz); // Clean final digital polynomial

disp('Digital LPF Transfer function H(Z) =');
disp(Hz);

// --- 8. PLOT FREQUENCY RESPONSE ---
HW = frmag(Hz, 512); 
w = 0 : %pi/511 : %pi; 
plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency (w / \pi)');
ylabel('Magnitude');
title('Chebyshev IIR LPF (Manual Formula Method)');
```


# OUTPUT:

<img width="1917" height="895" alt="image" src="https://github.com/user-attachments/assets/62380b48-ee99-48d1-9426-28a1f2a2155e" />

# RESULT: 
Thus design of Chebyshev Low pass IIR filter waveforms were plotted and output was verified.
