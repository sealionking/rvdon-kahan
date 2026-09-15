/**
 * rvdon_kahan.h — DiVo RVDon Kahan Public API
 *
 * FP64-equivalent force accumulation for molecular dynamics
 * on FP32 hardware. Zero hardware modifications required.
 *
 * Copyright (c) 2026 DiVo Gen²AI. All rights reserved.
 *
 * CONFIDENTIAL AND PROPRIETARY — This software is available only
 * under a commercial license agreement. Unauthorized copying,
 * distribution, reverse engineering, or use is strictly prohibited.
 *
 * Contact: wangjueju+divobot@gmail.com for licensing.
 *
 * Evaluation licenses available for 90-day non-production use.
 *
 * ⚠️ IMPORTANT: Do NOT compile with -ffast-math or equivalent flags.
 * Kahan compensated summation relies on IEEE 754 arithmetic semantics.
 * -ffast-math permits reordering of floating-point operations, which
 * destroys the compensation and degrades accuracy to naive FP32 level.
 */

#ifndef RVDON_KAHAN_H
#define RVDON_KAHAN_H

#ifdef __cplusplus
extern "C" {
#endif

/* ===================================================================
 * Opaque Types — Implementation details hidden in compiled library
 * =================================================================== */

/** Single-dimension Kahan compensated accumulator */
typedef struct rvdon_kahan_s rvdon_kahan_t;

/** Three-dimensional force accumulator with Kahan compensation */
typedef struct rvdon_force_acc_s rvdon_force_acc_t;

/* ===================================================================
 * Scalar Kahan Accumulator
 * =================================================================== */

/** Create a new Kahan accumulator initialized to zero */
rvdon_kahan_t* rvdon_kahan_create(void);

/** Destroy a Kahan accumulator and free memory */
void rvdon_kahan_destroy(rvdon_kahan_t *acc);

/** Kahan compensated accumulation: acc += x */
void rvdon_kahan_add(rvdon_kahan_t *acc, float x);

/** Get the accumulated value */
float rvdon_kahan_get(const rvdon_kahan_t *acc);

/** Reset accumulator to zero */
void rvdon_kahan_reset(rvdon_kahan_t *acc);

/* ===================================================================
 * 3D Force Vector Accumulator
 * =================================================================== */

/** Create a new 3D force accumulator initialized to zero */
rvdon_force_acc_t* rvdon_force_acc_create(void);

/** Destroy a force accumulator and free memory */
void rvdon_force_acc_destroy(rvdon_force_acc_t *acc);

/** Accumulate a 3D force contribution with Kahan compensation */
void rvdon_force_accum(rvdon_force_acc_t *acc, float fx, float fy, float fz);

/** Get accumulated force vector */
void rvdon_force_acc_get(const rvdon_force_acc_t *acc,
                          float *fx, float *fy, float *fz);

/** Reset force accumulator to zero */
void rvdon_force_acc_reset(rvdon_force_acc_t *acc);

/* ===================================================================
 * Batch Operations (N atoms)
 * =================================================================== */

/** Create an array of N force accumulators */
rvdon_force_acc_t** rvdon_force_batch_create(int N);

/** Destroy an array of N force accumulators */
void rvdon_force_batch_destroy(rvdon_force_acc_t **acc, int N);

/** Accumulate pairwise force: atom i gets +f, atom j gets -f (Newton's 3rd) */
void rvdon_force_accum_pair(rvdon_force_acc_t **forces, int i, int j,
                             float fx, float fy, float fz);

/* ===================================================================
 * MD Force Functions
 * =================================================================== */

/**
 * Compute Lennard-Jones force between two atoms.
 *
 * V(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
 * F(r) = 24ε/r [2(σ/r)¹² - (σ/r)⁶]
 *
 * @param dx,dy,dz  Distance vector (r_j - r_i)
 * @param sigma     LJ sigma parameter
 * @param epsilon   LJ epsilon parameter
 * @param r_cut     Cutoff distance (beyond this, force = 0)
 * @param fx,fy,fz  Output: force on atom i.
 *                  符号约定（2026-09-15 三轨交叉验证确认）：F_i = (dV/dr)/r · (dx,dy,dz)，
 *                  即近距排斥时 F_i 指向 −dx（远离 j），远距吸引时指向 +dx。
 *                  实现须满足 F_i = −24ε/r²·(2(σ/r)¹²−(σ/r)⁶)·(dx,dy,dz)。
 * @return          Potential energy at this distance
 */
float rvdon_lj_force(float dx, float dy, float dz,
                      float sigma, float epsilon, float r_cut,
                      float *fx, float *fy, float *fz);

/**
 * Compute all pairwise LJ forces with Kahan-compensated accumulation.
 *
 * Uses triangle symmetry (only j > i) for 50% computation savings.
 * This is the software equivalent of PF_TMM hardware acceleration.
 *
 * @param coords   Atom coordinates [N][3] (x,y,z interleaved)
 * @param N        Number of atoms
 * @param sigma    LJ sigma
 * @param epsilon  LJ epsilon
 * @param r_cut    Cutoff distance
 * @param forces   Output: pre-allocated force accumulator array [N]
 * @return         Total potential energy
 */
float rvdon_lj_forces_kahan(const float *coords, int N,
                              float sigma, float epsilon, float r_cut,
                              rvdon_force_acc_t **forces);

/**
 * Compute Buckingham force between two atoms.
 * V(r) = A·exp(-r/ρ) - C/r⁶
 *
 * The exp(-r/ρ) term maps to the RVDon FA_SOFTMAX LUT exp pipeline
 * when running on RVDon hardware.
 */
float rvdon_buckingham_force(float dx, float dy, float dz,
                              float A, float rho, float C, float r_cut,
                              float *fx, float *fy, float *fz);

/* ===================================================================
 * Version & Build Info
 * =================================================================== */

/** Library version string (e.g., "1.0.0") */
const char* rvdon_kahan_version(void);

/** Build configuration string */
const char* rvdon_kahan_build_info(void);

#ifdef __cplusplus
}
#endif

#endif /* RVDON_KAHAN_H */
