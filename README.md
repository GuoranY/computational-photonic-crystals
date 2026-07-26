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
  <img
    src="figures/p01_periodic_dielectric.png"
    alt="One-dimensional periodic dielectric refractive-index profile"
    width="650"
  >
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

```math
E_k(x) = u_k(x)e^{ikx},
```

where the periodic part satisfies

```math
u_k(x+a) = u_k(x).
```

For visualization, the periodic part is modeled as

```math
u_k(x) = 1 + A\cos(Gx),
```

with

```math
G = \frac{2\pi}{a}.
```

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

```math
E_k(x+a) = e^{ika}E_k(x).
```

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

```math
k_{\mathrm{B}} = \frac{\pi}{a},
```

the reciprocal-lattice vector is

```math
G = \frac{2\pi}{a}.
```

The periodic dielectric structure couples the forward-propagating wave
$e^{ik_{\mathrm{B}}x}$ to the backward-propagating wave
$e^{-ik_{\mathrm{B}}x}$ because

```math
k_{\mathrm{B}} - G = -k_{\mathrm{B}}.
```

The forward- and backward-propagating waves can be combined as

```math
e^{ik_{\mathrm{B}}x} + e^{-ik_{\mathrm{B}}x} = 2\cos(k_{\mathrm{B}}x),
```

and

```math
e^{ik_{\mathrm{B}}x} - e^{-ik_{\mathrm{B}}x} = 2i\sin(k_{\mathrm{B}}x).
```

The symmetric and antisymmetric combinations of these counter-propagating
waves form two standing-wave modes:

```math
E_{\cos}(x) = \cos(k_{\mathrm{B}}x),
```

```math
E_{\sin}(x) = \sin(k_{\mathrm{B}}x).
```

Their corresponding field intensities are

```math
|E_{\cos}(x)|^2 = \cos^2(k_{\mathrm{B}}x),
```

```math
|E_{\sin}(x)|^2 = \sin^2(k_{\mathrm{B}}x).
```

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

### P04 — Transfer Matrix and Transmission Spectrum

Calculate the transmission spectrum of a finite one-dimensional photonic crystal
using the transfer-matrix method.

#### Parameters

- Material A refractive index: $n_1 = 1.0$
- Material B refractive index: $n_2 = 3.5$
- Lattice constant: $a = 1.0$
- Fill fraction: $0.5$
- Number of unit cells: $10$

#### Model

Within each homogeneous dielectric layer, the electric field is represented as
a superposition of right- and left-propagating waves:

```math
E(x) = A e^{iqx} + B e^{-iqx},
```

where

```math
q = \frac{n\omega}{c}.
```

Propagation through a layer of thickness $d$ is described by the propagation
matrix

```math
P(n,d,\omega) = \begin{pmatrix} e^{iqd} & 0 \\ 0 & e^{-iqd} \end{pmatrix}.
```

At an interface between two dielectric materials with refractive indices
$n_i$ and $n_j$, the wave amplitudes are related by an interface matrix.

Combining the propagation and interface matrices gives the transfer matrix of
one unit cell. Repeating the unit-cell matrix over multiple periods gives the
total transfer matrix of the finite photonic crystal:

```math
M_{\mathrm{total}} = M_{\mathrm{cell}}^N,
```

where $N$ is the number of unit cells.

For incidence and exit through the same surrounding medium, the power
transmission coefficient is calculated from the transmission amplitude $t$ as

```math
T = |t|^2.
```

#### Output

<p align="center">
  <img src="figures/p04_transmission_spectrum.png" width="650">
</p>

The transmission spectrum contains frequency intervals with high transmission
and intervals in which transmission is strongly suppressed. The low-transmission
regions indicate the formation of photonic stop bands in the finite periodic
structure.

### P05 — Allowed and Forbidden Frequency Bands

Use the transfer matrix of one unit cell to determine the allowed and forbidden
frequency regions of an infinite one-dimensional photonic crystal.

#### Parameters

- Material A refractive index: $n_1 = 1.0$
- Material B refractive index: $n_2 = 3.5$
- Lattice constant: $a = 1.0$
- Fill fraction: $0.5$

#### Model

According to Bloch's theorem, the field amplitudes in neighboring unit cells are
related by

```math
\begin{pmatrix} A(x+a) \\ B(x+a) \end{pmatrix} = e^{ika} \begin{pmatrix} A(x) \\ B(x) \end{pmatrix}.
```

Therefore, the Bloch factor $e^{ika}$ is an eigenvalue of the unit-cell transfer
matrix $M$.

For a lossless unit cell,

```math
\det(M) = 1,
```

and the eigenvalue equation gives the Bloch dispersion relation

```math
\cos(ka) = \frac{1}{2}\mathrm{Tr}(M).
```

Define the Bloch function

```math
F(\omega) = \frac{1}{2}\mathrm{Tr}(M).
```

When

```math
|F(\omega)| \leq 1,
```

a real Bloch wavevector exists, so the corresponding frequency belongs to an
allowed photonic band.

When

```math
|F(\omega)| > 1,
```

the Bloch wavevector becomes complex. The field then decays exponentially
through the crystal, and the corresponding frequency lies inside a forbidden
band or photonic band gap.

#### Output

<p align="center">
  <img src="figures/p05_bloch_function.png" width="650">
</p>

The Bloch function is compared with the boundaries $F=1$ and $F=-1$.
Frequencies for which the curve lies between these boundaries support real
Bloch wavevectors.

<p align="center">
  <img src="figures/p05_allowed_forbidden_bands.png" width="650">
</p>

The second figure classifies the frequency axis into allowed and forbidden
regions using the condition $|F(\omega)| \leq 1$.

### P06 — Complex Bloch Wavevector

Calculate the real and imaginary parts of the Bloch wavevector across both
allowed and forbidden frequency regions.

#### Parameters

- Material A refractive index: $n_1 = 1.0$
- Material B refractive index: $n_2 = 3.5$
- Lattice constant: $a = 1.0$
- Fill fraction: $0.5$

#### Model

The Bloch wavevector is obtained from the dispersion relation

```math
\cos(ka) = F(\omega),
```

where

```math
F(\omega) = \frac{1}{2}\mathrm{Tr}(M).
```

It can therefore be written as

```math
k(\omega) = \frac{1}{a}\arccos\left[F(\omega)\right].
```

Inside an allowed band, $|F(\omega)|\leq 1$, so $k$ is real and the Bloch wave
propagates through the periodic structure.

Inside a forbidden band, $|F(\omega)|>1$, so $k$ becomes complex:

```math
k = k_{\mathrm{r}} + i k_{\mathrm{i}}.
```

The corresponding Bloch factor is

```math
e^{ikx} = e^{ik_{\mathrm{r}}x}e^{-k_{\mathrm{i}}x}.
```

The real part $k_{\mathrm{r}}$ describes the spatial phase variation, while the
imaginary part $k_{\mathrm{i}}$ gives the exponential attenuation rate inside
the photonic band gap.

#### Output

<p align="center">
  <img src="figures/p06_bloch_wavevector.png" width="650">
</p>

The real part of the Bloch wavevector varies across the allowed bands, while a
nonzero imaginary part appears inside the forbidden frequency regions. This
directly shows the transition from propagating Bloch waves to evanescent Bloch
waves.

### P07 — Photonic Band Structure

Construct the photonic band structure of the one-dimensional periodic dielectric
from the real Bloch wavevectors obtained in the allowed frequency regions.

#### Parameters

- Material A refractive index: $n_1 = 1.0$
- Material B refractive index: $n_2 = 3.5$
- Lattice constant: $a = 1.0$
- Fill fraction: $0.5$

#### Model

The photonic band structure represents the relationship between the allowed
frequencies and the corresponding real Bloch wavevectors.

For each normalized frequency, the unit-cell transfer matrix is constructed and
the Bloch relation

```math
\cos(ka) = \frac{1}{2}\mathrm{Tr}(M)
```

is evaluated.

Only frequencies satisfying

```math
\left| \frac{1}{2}\mathrm{Tr}(M) \right| \leq 1
```

are included because these frequencies produce real Bloch wavevectors.

The horizontal axis is the normalized Bloch wavevector

```math
\frac{ka}{\pi},
```

and the vertical axis is the normalized frequency

```math
\frac{\omega a}{2\pi c}.
```

The band structure is displayed along the one-dimensional high-symmetry path

```math
X \rightarrow \Gamma \rightarrow X,
```

where

```math
\Gamma: k=0
```

and

```math
X: k=\pm\frac{\pi}{a}.
```

#### Output

<p align="center">
  <img src="figures/p07_photonic_band_structure.png" width="650">
</p>

Each continuous branch represents an allowed photonic band. The empty frequency
intervals separating neighboring branches are photonic band gaps, in which no
real Bloch wavevector exists.

### P08 — Band-Gap Parameter Study

Investigate how the first photonic band gap of a one-dimensional periodic
dielectric depends on refractive-index contrast, fill fraction, and the number
of unit cells.

#### Base Parameters

- Material A refractive index: $n_1 = 1.0$
- Material B refractive index: $n_2 = 3.5$
- Lattice constant: $a = 1.0$
- Fill fraction: $f = d_1/a = 0.5$
- Number of unit cells: $N = 10$
- Incident-medium refractive index: $n_{\mathrm{incident}} = 1.0$
- Exit-medium refractive index: $n_{\mathrm{exit}} = 1.0$

#### Model

The first photonic band gap is identified from the unit-cell transfer matrix
using the Bloch condition

```math
\cos(ka) = \frac{1}{2}\mathrm{Tr}(M_{\mathrm{cell}}).
```

A normalized frequency belongs to a forbidden band when

```math
\left|
\frac{1}{2}\mathrm{Tr}(M_{\mathrm{cell}})
\right| > 1.
```

For each parameter value, the lower and upper edges of the first band gap are
extracted as

```math
\nu_{\mathrm{lower}},
\qquad
\nu_{\mathrm{upper}},
```

where

```math
\nu = \frac{\omega a}{2\pi c}.
```

The absolute band-gap width is

```math
\Delta\nu =\nu_{\mathrm{upper}} - \nu_{\mathrm{lower}},
```

and the relative band-gap width is

```math
\frac{\Delta\nu}{\nu_{\mathrm{mid}}},
\qquad
\nu_{\mathrm{mid}}
=
\frac{
\nu_{\mathrm{lower}}
+
\nu_{\mathrm{upper}}
}{2}.
```

Three parameter studies are performed.

#### Refractive-Index Contrast

The refractive index of material B is varied over

```math
n_2 = 1.5,\ 2.0,\ 2.5,\ 3.5,\ 4.5,
```

while the physical fill fraction remains fixed.

Increasing the refractive-index contrast shifts the first band gap toward lower
normalized frequencies. The absolute gap width is non-monotonic because both
the refractive-index contrast and the optical thickness of the second layer
change simultaneously. In contrast, the relative gap width increases with
refractive-index contrast.

For example, increasing $n_2$ from $1.5$ to $4.5$ shifts the first band-gap
center from approximately $0.399$ to $0.166$, while the relative gap width
increases from approximately $0.243$ to $0.541$.

<p align="center">
  <img
    src="figures/p08_index_contrast_transmission.png"
    alt="Transmission spectra for different refractive-index contrasts"
    width="650"
  >
</p>

<p align="center">
  <img
    src="figures/p08_index_contrast_gap_edges.png"
    alt="First band-gap edges versus refractive-index contrast"
    width="650"
  >
</p>

<p align="center">
  <img
    src="figures/p08_index_contrast_gap_width.png"
    alt="Absolute first band-gap width versus refractive-index contrast"
    width="650"
  >
</p>

<p align="center">
  <img
    src="figures/p08_index_contrast_relative_gap_width.png"
    alt="Relative first band-gap width versus refractive-index contrast"
    width="650"
  >
</p>

#### Fill Fraction

The fill fraction

```math
f = \frac{d_1}{a}
```

is varied while maintaining

```math
d_1 = fa,
\qquad
d_2 = (1-f)a.
```

The relative band-gap width reaches a maximum near $f \approx 0.8$. This agrees
with the equal-optical-thickness condition

```math
n_1d_1 = n_2d_2,
```

which predicts

```math
f_{\mathrm{opt}}
=
\frac{n_2}{n_1+n_2}
=
\frac{3.5}{4.5}
\approx 0.778.
```

As the fill fraction approaches either extreme, one dielectric layer becomes
very thin and the periodic modulation weakens.

The relative gap width increases from approximately $0.093$ at $f=0.10$ to a
maximum of approximately $0.745$ near $f=0.80$, before decreasing as
$f$ approaches $1$.

<p align="center">
  <img
    src="figures/p08_fill_fraction_transmission.png"
    alt="Transmission spectra for different fill fractions"
    width="650"
  >
</p>

<p align="center">
  <img
    src="figures/p08_fill_fraction_gap_edges.png"
    alt="First band-gap edges versus fill fraction"
    width="650"
  >
</p>

<p align="center">
  <img
    src="figures/p08_fill_fraction_gap_width.png"
    alt="Absolute first band-gap width versus fill fraction"
    width="650"
  >
</p>

<p align="center">
  <img
    src="figures/p08_fill_fraction_relative_gap_width.png"
    alt="Relative first band-gap width versus fill fraction"
    width="650"
  >
</p>

#### Number of Unit Cells

The number of unit cells is varied over

```math
N = 2,\ 4,\ 6,\ 10,\ 20.
```

Changing $N$ does not alter the theoretical band-gap boundaries because the
unit-cell structure remains unchanged. Instead, increasing $N$ strengthens the
suppression of transmission inside the band gap and makes the band edges
sharper.

The transmission at the center of the first band gap decreases approximately
exponentially with crystal length:

```math
T(N) \propto e^{-2\kappa Na},
```

where $\kappa$ is the imaginary part of the Bloch wavevector inside the gap.

At the center of the first band gap, the transmission decreases from
approximately $8.26\times10^{-2}$ for $N=2$ to
$3.43\times10^{-15}$ for $N=20$.

<p align="center">
  <img
    src="figures/p08_cell_number_transmission.png"
    alt="Transmission spectra for different numbers of unit cells"
    width="650"
  >
</p>

<p align="center">
  <img
    src="figures/p08_cell_number_gap_suppression.png"
    alt="Band-gap-center transmission versus number of unit cells"
    width="650"
  >
</p>

The results distinguish changes to the intrinsic band structure from
finite-size effects. Refractive-index contrast and fill fraction modify the
unit cell and therefore change the band gap itself, whereas the number of unit
cells controls how strongly a finite crystal suppresses transmission within
that gap.

### P11 — Two-Dimensional Periodic Dielectric

Construct and visualize a two-dimensional square lattice of dielectric rods
embedded in an air background.

#### Parameters

- Background refractive index: $n_{\mathrm{bg}} = 1.0$
- Rod refractive index: $n_{\mathrm{rod}} = 3.5$
- Lattice constant: $a = 1.0$
- Rod radius: $r = 0.2a$
- Number of unit cells: $6 \times 6$

#### Model

Identical circular dielectric rods are placed at the center of each unit cell
of a two-dimensional square lattice. The refractive-index distribution is
periodic along both the $x$- and $y$-directions:

```math
n(x + a, y) = n(x, y),
```

```math
n(x, y + a) = n(x, y).
```

A point $(x,y)$ lies inside a rod centered at $(x_c,y_c)$ when

```math
(x-x_c)^2 + (y-y_c)^2 \leq r^2.
```

The relative permittivity is obtained from

```math
\varepsilon_r(x,y) = n^2(x,y).
```

#### Output

<p align="center">
  <img
    src="figures/p11_2d_dielectric_structure.png"
    alt="Two-dimensional square lattice of dielectric rods in air"
    width="650"
  >
</p>

The resulting structure serves as the real-space geometric foundation for
two-dimensional reciprocal-lattice, Brillouin-zone, and photonic-band
calculations.

## Project Structure

```text
computational-photonic-crystals/
├── scripts/
│   ├── p01_periodic_dielectric.py
│   ├── p02_bloch_wave_visualization.py
│   ├── p03_bragg_standing_waves.py
│   ├── p04_transfer_matrix.py
│   ├── p05_bloch_band_structure.py
│   ├── p06_bloch_wavevector.py
│   ├── p07_photonic_band_structure.py
│   ├── p08_band_gap_parameter_study.py
│   └── p11_2d_dielectric_structure.py
├── utils/
│   ├── transfer_matrix_utils.py
│   ├── bloch_utils.py
│   └── transmission_utils.py
├── figures/
│   ├── p01_periodic_dielectric.png
│   ├── p02_bloch_wave_visualization.png
│   ├── p02_bloch_translation_verification.png
│   ├── p03_bragg_standing_waves.png
│   ├── p04_transmission_spectrum.png
│   ├── p05_bloch_function.png
│   ├── p05_allowed_forbidden_bands.png
│   ├── p06_bloch_wavevector.png
│   ├── p07_photonic_band_structure.png
│   ├── p08_cell_number_gap_suppression.png
│   ├── p08_cell_number_transmission.png
│   ├── p08_fill_fraction_gap_edges.png
│   ├── p08_fill_fraction_gap_width.png
│   ├── p08_fill_fraction_relative_gap_width.png
│   ├── p08_fill_fraction_transmission.png
│   ├── p08_index_contrast_gap_edges.png
│   ├── p08_index_contrast_gap_width.png
│   ├── p08_index_contrast_relative_gap_width.png
│   ├── p08_index_contrast_transmission.png
│   └── p11_2d_dielectric_structure.png
├── docs/
├── requirements.txt
├── .gitignore
└── README.md
```

## Status

The one-dimensional photonic-crystal foundation is complete, including periodic
dielectric modeling, Bloch waves, Bragg standing waves, transfer-matrix analysis,
photonic band-gap identification, complex Bloch wavevectors, and photonic band
structures.

Development has now progressed to two-dimensional photonic crystals, beginning
with the construction of a square lattice of dielectric rods in air. Future
modules will introduce reciprocal lattices, Brillouin zones, two-dimensional
band structures, defect cavities, and waveguides.
