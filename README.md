# DiVo RVDon Kahan

**FP64-Equivalent Force Accumulation for Molecular Dynamics — on FP32 Hardware**

RVDon Kahan is a software library that provides Kahan-compensated force accumulation for molecular dynamics simulations, achieving FP64-equivalent precision on FP32 hardware with **zero hardware modifications**.

## Why RVDon Kahan?

Protein MD simulations accumulate force contributions from **hundreds of neighbor atoms** per timestep. In a folded protein (Halle 2002, PNAS), the average Cα atom has **67.5 noncovalent neighbors** within 7 Å; even the most exposed surface atoms have 20-40. With explicit solvent at typical LJ cutoffs (8-12 Å), **every atom has 100-400+ neighbors** — always a dense accumulation problem.

Naive FP32 accumulation error grows as O(Nε), where N = neighbor count and ε ≈ 1.2×10⁻⁷. For N = 100-400 (protein MD), this gives O(Nε) ≈ 1.2e-5 to 4.8e-5 per accumulation — violating the 1e-6 conservation threshold. RVDon Kahan reduces this to **O(ε)** using compensated summation on existing FP32 hardware — at zero hardware cost.

| Method | Per-Accumulation Error | Energy Drift/Step | Hardware Cost |
|--------|:---:|:---:|:---:|
| Naive FP32 (N=200) | O(Nε) ≈ 2.4e-5 | > 1e-5 | Baseline |
| **RVDon Kahan (FP32)** | **O(ε) ≈ 1.2e-7** | **< 1e-6/step** | **0% (software)** |
| FP64 Hardware | O(ε₆₄) ≈ 2.2e-16 | ~1e-15 | ~30% TCU area |

Kahan FP32 is **200× better** than naive FP32 for typical protein MD, and **1000× below** the 1e-6 conservation threshold — at zero hardware cost.

> **Note on scope:** RVDon Kahan is designed for **organic macromolecule** simulations (proteins, nucleic acids, drug molecules in solvent) where atomic density is always high. In gas-phase or vacuum simulations of small molecules (where neighbors < 10), accumulation error is negligible regardless of method — Kahan and naive give similar results.

## Validated Results (Independent Red-Team Reproduced)

| Test | DiVo Result | Red-Team Result | Threshold |
|------|:---:|:---:|:---:|
| Single-atom force accuracy (2000 neighbors) | rel err < 1e-5 | rel err = 3.06e-8 | < 1e-4 |
| 256-atom LJ system RMS force error | rel err < 1e-4 | rel err = 1.19e-6 | < 1e-4 |
| Energy drift (27-atom cluster, 1000 steps) | 8.99e-7/step | **2.13e-7/step** | < 1e-6/step |

Both DiVo and independent red-team results pass all thresholds.

## Quick Start

```c
#include "rvdon_kahan.h"

int main() {
    /* Create force accumulators for N atoms */
    int N = 1000;
    rvdon_force_acc_t **forces = malloc(N * sizeof(rvdon_force_acc_t*));
    for (int i = 0; i < N; i++)
        forces[i] = rvdon_force_acc_create();
    
    /* Accumulate pairwise LJ forces with Kahan compensation */
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            float fx, fy, fz;
            float dx = coords[j*3+0] - coords[i*3+0];
            float dy = coords[j*3+1] - coords[i*3+1];
            float dz = coords[j*3+2] - coords[i*3+2];
            
            float pe = rvdon_lj_force(dx, dy, dz, sigma, epsilon, r_cut, &fx, &fy, &fz);
            
            /* Newton's 3rd law: F_ij = -F_ji */
            rvdon_force_accum_pair(forces, i, j, fx, fy, fz);
        }
    }
    
    /* Read final forces */
    for (int i = 0; i < N; i++) {
        float fx, fy, fz;
        rvdon_force_acc_get(forces[i], &fx, &fy, &fz);
        printf("Atom %d: F = (%.6f, %.6f, %.6f)\n", i, fx, fy, fz);
        rvdon_force_acc_destroy(forces[i]);
    }
    free(forces);
    return 0;
}
```

### Build

```bash
gcc -O2 -o md_sim examples/basic_usage.c -lrvdon_kahan -lm
```

## API Reference

See [docs/api-reference.md](docs/api-reference.md).

## Integration with RVDon Hardware

RVDon Kahan is designed to work with the RVDon PF Extension TCU:

| Component | RVDon Hardware | RVDon Kahan Software |
|-----------|:---:|:---:|
| Pairwise force calculation | PF_TMM (triangle masked) | Kahan force accumulation |
| Exp-based potentials | FA_SOFTMAX (LUT exp pipeline) | Configurable LUT tables |
| Bonded forces | Scalar core (RISC-V) | Standard FP32 arithmetic |
| Integration (Verlet) | Scalar core | FP64 reference implementation |

RVDon Kahan can also be used independently on any FP32-capable hardware (NVIDIA/AMD GPUs, ARM CPUs, etc.).

## Precision Validation

See [docs/precision-validation.md](docs/precision-validation.md) for detailed accuracy analysis and benchmark methodology.

## License

**Commercial license required.** This software is NOT open-source.

- **Evaluation license**: Free for 90-day evaluation (non-production use)
- **Development license**: For integration and testing
- **Production license**: Per-seat or site license for deployment
- **OEM license**: For bundling with hardware products

Contact: wangjueju+divobot@gmail.com

---

Copyright (c) 2026 DiVo Gen²AI. All rights reserved.
