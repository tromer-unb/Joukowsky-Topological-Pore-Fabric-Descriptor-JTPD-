# JTPD descriptor definition

## Overview

The Joukowsky–Topological Pore-Fabric Descriptor (JTPD) is a scale-aware vector of interpretable quantities computed from a segmented 2-D rock-section field.

For a pore mask `P` observed over physical area `A`, the reference implementation combines abundance, topology, spanning, interface complexity, contour anisotropy, and a nonlinear Joukowsky response.

## 1. Pore-area fraction

`φ_p = N_pore / N_total`

## 2. Digital topology

The implementation uses **8-connectivity** for the pore foreground and **4-connectivity** for the complement when counting enclosed holes.

- `β₀`: connected pore components.
- `β₁`: enclosed holes in the pore phase.
- `χ = β₀ − β₁`: Euler characteristic.

Area-normalized quantities `β₀/A`, `β₁/A`, and `χ/A` are reported in mm⁻².

## 3. Directional spanning

Each pore component is tested for opposite-boundary intersection: left-to-right (X) and top-to-bottom (Y).

## 4. Interface density

`I = L_interface / A`, reported in mm mm⁻². Digital interface length is resolution-sensitive and must be interpreted with effective pixel size.

## 5. Contour anisotropy

For centered contour points with covariance eigenvalues `λ_max` and `λ_min`, the implementation reports approximately `λ_max / λ_min`, with a small numerical regularization.

## 6. Joukowsky nonlinear contour response

Each contour is centered and normalized by its mean radius. The reference implementation evaluates edge-based orientations and averages the transformed geometric measures to reduce dependence on global orientation and contour ordering.

For normalized complex coordinate `u`, mean radius `r`, and Joukowsky strength `α ≥ 0`:

`w = r [u + α²/u]`

Two dimensionless responses are retained:

- `J_A(α) = A_J(α)/A`
- `J_L(α) = L_J(α)/L`

At `α = 0`, both are expected to be approximately 1. The reference grid is `0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35`.

## 7. Physical calibration

The descriptor is not scale-free. Effective pixel size is updated after downsampling:

`effective_pixel_size = source_pixel_size × source_width / target_width`

For the included 4096-pixel fields analyzed at 256 pixels: `0.44 × 4096/256 = 7.04 µm px⁻¹`.

## Recommended reporting

Always report segmentation convention, field dimensions, source/effective pixel size, connectivity convention, minimum retained contour area, Joukowsky alpha grid, and the exact software version or commit.

## Scope

JTPD is a 2-D morphological descriptor. It does not by itself determine 3-D connectivity, permeability, hydraulic conductance, or representative elementary volume.
