# computational-photonic-crystals

Develop numerical and visual models of one-dimensional and two-dimensional
photonic crystals, progressing from Bragg scattering and Bloch waves to band
structures, photonic band gaps, defect cavities, and waveguides.

## Implemented Modules

### P01 — One-Dimensional Periodic Dielectric

Construct a one-dimensional periodic dielectric profile consisting of alternating
materials with different refractive indices.

#### Parameters

- Material A refractive index: $n_1 = 1.0$
- Material B refractive index: $n_2 = 3.5$
- Lattice constant: $a = 1.0$
- Number of unit cells: $10$
- Fill fraction: $0.5$

#### Model

The refractive-index profile alternates periodically between two materials within
each unit cell.

#### Output

<p align="center">
  <img src="figures/p01_periodic_dielectric.png" width="650">
</p>

This periodic refractive-index profile serves as the fundamental model for
subsequent photonic-crystal simulations, including Bloch-wave propagation,
band-structure calculations, and transfer-matrix analysis.

### P02 — Bloch Wave Visualization

Visualize the structure of a one-dimensional Bloch wave and compare its periodic
part, plane-wave factor, and complete spatial form.

#### Parameters

- Lattice constant: $a = 1.0$
- Number of unit cells: $8$
- Wave vector: $k = 0.6\pi/a$
- Periodic modulation amplitude: $A = 0.35$

#### Model

A Bloch wave is written as

$$
E_k(x) = u_k(x)e^{ikx},
$$

where the periodic part satisfies

$$
u_k(x+a) = u_k(x).
$$

For visualization, the periodic part is modeled as

$$
u_k(x) = 1 + A\cos(Gx),
$$

with

$$
G = \frac{2\pi}{a}.
$$

The script compares:

1. the periodic part $u_k(x)$,
2. the real part of the plane wave $e^{ikx}$,
3. the real part of the complete Bloch wave $E_k(x)$.

Vertical dashed lines indicate neighboring unit-cell boundaries.

#### Output

<p align="center">
  <img src="figures/p02_bloch_wave_visualization.png" width="650">
</p>

The first figure compares the periodic part, the plane-wave factor, and the
complete Bloch wave.

<p align="center">
  <img src="figures/p02_bloch_translation_verification.png" width="650">
</p>

The second figure verifies the Bloch translation property by comparing
$E_k(x+a)$ with $e^{ika}E_k(x)$. The two curves overlap, confirming that

$$
E_k(x+a) = e^{ika}E_k(x).
$$

### P03 — Bragg Scattering and Standing Waves

Visualize the formation of standing waves produced by Bragg coupling at the
boundary of the first Brillouin zone, and compare their field-intensity
distributions within the periodic dielectric structure.

#### Parameters

- Material A refractive index: $n_1 = 1.0$
- Material B refractive index: $n_2 = 3.5$
- Lattice constant: $a = 1.0$
- Number of unit cells: $6$
- Fill fraction: $0.5$
- Bragg wave vector: $k_{\mathrm{B}} = \pi/a$

#### Model

At the boundary of the first Brillouin zone,

$$
k_{\mathrm{B}} = \frac{\pi}{a},
$$

the reciprocal-lattice vector is

$$
G = \frac{2\pi}{a}.
$$

The periodic dielectric structure couples the forward-propagating wave
$e^{ik_{\mathrm{B}}x}$ to the backward-propagating wave
$e^{-ik_{\mathrm{B}}x}$ because

$$
k_{\mathrm{B}} - G = -k_{\mathrm{B}}.
$$

The forward- and backward-propagating waves can be combined as

$$
e^{ik_{\mathrm{B}}x} + e^{-ik_{\mathrm{B}}x}
= 2\cos(k_{\mathrm{B}}x),
$$

and

$$
e^{ik_{\mathrm{B}}x} - e^{-ik_{\mathrm{B}}x}
= 2i\sin(k_{\mathrm{B}}x).
$$

The symmetric and antisymmetric combinations of these counter-propagating
waves form two standing-wave modes:

$$
E_{\cos}(x) = \cos(k_{\mathrm{B}}x),
$$

$$
E_{\sin}(x) = \sin(k_{\mathrm{B}}x).
$$

Their corresponding field intensities are

$$
|E_{\cos}(x)|^2 = \cos^2(k_{\mathrm{B}}x),
$$

$$
|E_{\sin}(x)|^2 = \sin^2(k_{\mathrm{B}}x).
$$

Although the two modes have the same wave-vector magnitude, their intensity
maxima occur in different parts of the unit cell. As a result, they overlap
differently with the high- and low-index materials and acquire different
eigenfrequencies.

This frequency splitting at the Brillouin-zone boundary illustrates the
physical origin of the photonic band gap.

#### Output

<p align="center">
  <img src="figures/p03_bragg_standing_waves.png" width="650">
</p>

The figure compares the two standing-wave fields and their corresponding
intensity distributions against the periodic dielectric background.

## Project Structure

```text
computational-photonic-crystals/
├── scripts/
│   ├── p01_periodic_dielectric.py
│   ├── p02_bloch_wave_visualization.py
│   └── p03_bragg_standing_waves.py
├── figures/
│   ├── p01_periodic_dielectric.png
│   ├── p02_bloch_wave_visualization.png
│   ├── p02_bloch_translation_verification.png
│   └── p03_bragg_standing_waves.png
├── docs/
├── utils/
├── requirements.txt
├── .gitignore
└── README.md
```

## Status

This project is under development.