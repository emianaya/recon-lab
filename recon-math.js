/*
 * Core parallel-beam Radon transform math shared by the Recon Lab artifact
 * and its test suite. Kept dependency-free (typed arrays only) so it loads
 * as a plain <script> in the browser and as a CommonJS module in Node.
 */
(function (root, factory) {
  if (typeof module !== "undefined" && module.exports) {
    module.exports = factory();
  } else {
    root.ReconMath = factory();
  }
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  // Builds the sparse system matrix A (as CSR-by-ray: idx/w/rowStart) for an
  // N x N image, nAngles views over nBins each, sampled every `ds` pixels out
  // to +/-rmax along each ray, with bilinear splatting onto the pixel grid.
  function buildSystemMatrix(N, nAngles, nBins, ds, rmax, angles) {
    const cx = (N - 1) / 2, cy = (N - 1) / 2;
    const idxArr = [], wArr = [];
    const nRays = nAngles * nBins;
    const rowStart = new Int32Array(nRays + 1);
    let ptr = 0;
    for (let a = 0; a < nAngles; a++) {
      const th = (angles[a] * Math.PI) / 180;
      const cosT = Math.cos(th), sinT = Math.sin(th);
      for (let b = 0; b < nBins; b++) {
        const r = b - (nBins - 1) / 2;
        const row = a * nBins + b;
        rowStart[row] = ptr;
        const x0 = cx + r * cosT, y0 = cy + r * sinT;
        for (let t = -rmax; t <= rmax; t += ds) {
          const x = x0 - t * sinT, y = y0 + t * cosT;
          const xi = Math.floor(x), yi = Math.floor(y);
          const fx = x - xi, fy = y - yi;
          for (let dx = 0; dx <= 1; dx++) {
            const xx = xi + dx;
            if (xx < 0 || xx >= N) continue;
            const wx = dx ? fx : 1 - fx;
            for (let dy = 0; dy <= 1; dy++) {
              const yy = yi + dy;
              if (yy < 0 || yy >= N) continue;
              const w = wx * (dy ? fy : 1 - fy) * ds;
              if (w < 1e-4) continue;
              idxArr.push(yy * N + xx);
              wArr.push(w);
              ptr++;
            }
          }
        }
      }
    }
    rowStart[nRays] = ptr;
    return { idx: Int32Array.from(idxArr), w: Float32Array.from(wArr), rowStart };
  }

  // Backprojection of an all-ones sinogram: how much total weight lands on
  // each pixel. Used by MLEM to normalize its multiplicative update.
  function computeSensitivity(mat, npix) {
    const s = new Float32Array(npix);
    for (let i = 0; i < mat.idx.length; i++) s[mat.idx[i]] += mat.w[i];
    return s;
  }

  // Forward projection: A * image -> sinogram (one value per ray).
  function forwardProject(mat, image) {
    const nRays = mat.rowStart.length - 1;
    const proj = new Float32Array(nRays);
    for (let row = 0; row < nRays; row++) {
      let s = 0;
      const end = mat.rowStart[row + 1];
      for (let p = mat.rowStart[row]; p < end; p++) s += mat.w[p] * image[mat.idx[p]];
      proj[row] = s;
    }
    return proj;
  }

  // Backprojection: Aᵀ * sinogram -> image. Must stay the exact transpose of
  // forwardProject (same idx/w, same rows) — that's the invariant the
  // adjoint test below is checking.
  function backProject(mat, sino, npix) {
    const nRays = mat.rowStart.length - 1;
    const back = new Float32Array(npix);
    for (let row = 0; row < nRays; row++) {
      const rv = sino[row];
      if (rv === 0) continue;
      const end = mat.rowStart[row + 1];
      for (let p = mat.rowStart[row]; p < end; p++) back[mat.idx[p]] += mat.w[p] * rv;
    }
    return back;
  }

  // ---------------------------------------------------------------------
  // MLEM
  // ---------------------------------------------------------------------

  // One multiplicative MLEM update: estimate_{k+1} = estimate_k * (Aᵀ(measured / A·estimate_k)) / sensitivity
  function mlemStep(mat, measured, estimate, sens, npix, eps) {
    eps = eps === undefined ? 1e-6 : eps;
    const fwd = forwardProject(mat, estimate);
    const nRays = fwd.length;
    const ratio = new Float32Array(nRays);
    for (let i = 0; i < nRays; i++) ratio[i] = measured[i] / (fwd[i] + eps);
    const back = backProject(mat, ratio, npix);
    const next = new Float32Array(npix);
    for (let j = 0; j < npix; j++) next[j] = (estimate[j] * back[j]) / (sens[j] + eps);
    return next;
  }

  // Poisson log-likelihood of `estimate` given `measured` counts, dropping the
  // m*log(m) - log(m!) term (constant w.r.t. the estimate, so it doesn't affect
  // monotonicity). MLEM's EM-algorithm guarantee is that this never decreases
  // from one iteration to the next.
  function poissonLogLikelihood(mat, measured, estimate, eps) {
    eps = eps === undefined ? 1e-4 : eps;
    const fwd = forwardProject(mat, estimate);
    let s = 0;
    for (let i = 0; i < fwd.length; i++) {
      const f = Math.max(fwd[i], eps);
      s += measured[i] * Math.log(f) - f;
    }
    return s;
  }

  // ---------------------------------------------------------------------
  // FBP
  // ---------------------------------------------------------------------

  // Naive O(M^2) DFT/IDFT — fine at M = N_BINS ~ 64, no library needed.
  function dft(p) {
    const M = p.length, re = new Float32Array(M), im = new Float32Array(M);
    for (let k = 0; k < M; k++) {
      let sr = 0, si = 0;
      for (let n = 0; n < M; n++) {
        const ang = (-2 * Math.PI * k * n) / M;
        sr += p[n] * Math.cos(ang);
        si += p[n] * Math.sin(ang);
      }
      re[k] = sr; im[k] = si;
    }
    return { re, im };
  }
  function idftReal(re, im) {
    const M = re.length, out = new Float32Array(M);
    for (let n = 0; n < M; n++) {
      let s = 0;
      for (let k = 0; k < M; k++) {
        const ang = (2 * Math.PI * k * n) / M;
        s += re[k] * Math.cos(ang) - im[k] * Math.sin(ang);
      }
      out[n] = s / M;
    }
    return out;
  }
  // Hann-windowed ramp filter: |freq| up to cutoff*Nyquist, tapered smoothly to
  // 0 at the cutoff (avoids the ringing a hard-edged cutoff would cause).
  function rampFilterHann(M, cutoff) {
    const filt = new Float32Array(M);
    const nyq = M / 2;
    for (let k = 0; k < M; k++) {
      const freq = k <= nyq ? k : k - M;
      const af = Math.abs(freq);
      const limit = Math.max(1, cutoff * nyq);
      if (af > limit) continue;
      const taper = 0.5 * (1 + Math.cos((Math.PI * af) / limit));
      filt[k] = af * taper;
    }
    return filt;
  }
  function computeFBPSpectra(measured, nAngles, nBins) {
    const spectra = [];
    for (let a = 0; a < nAngles; a++) {
      const row = measured.subarray(a * nBins, a * nBins + nBins);
      spectra.push(dft(row));
    }
    return spectra;
  }
  function fbpReconstruct(mat, spectra, cutoff, nAngles, nBins, npix) {
    const filt = rampFilterHann(nBins, cutoff);
    const filtered = new Float32Array(nAngles * nBins);
    for (let a = 0; a < nAngles; a++) {
      const { re, im } = spectra[a];
      const fre = new Float32Array(nBins), fim = new Float32Array(nBins);
      for (let k = 0; k < nBins; k++) { fre[k] = re[k] * filt[k]; fim[k] = im[k] * filt[k]; }
      const row = idftReal(fre, fim);
      filtered.set(row, a * nBins);
    }
    const back = backProject(mat, filtered, npix);
    const norm = Math.PI / nAngles;
    for (let j = 0; j < npix; j++) back[j] *= norm;
    return back;
  }

  // ---------------------------------------------------------------------
  // Noise simulation
  // ---------------------------------------------------------------------

  // mulberry32: tiny, fast, seedable PRNG (returns a function -> float in [0,1)).
  // Deterministic for a given seed, so a sinogram can be reproduced exactly.
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6d2b79f5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Poisson draw via Knuth's method (small lambda) / normal approximation
  // (large lambda). `rng` defaults to Math.random so existing call sites are
  // unaffected; pass a seeded rng (e.g. mulberry32(seed)) for reproducibility.
  function poisson(lambda, rng) {
    rng = rng || Math.random;
    if (lambda <= 0) return 0;
    if (lambda < 35) {
      const L = Math.exp(-lambda);
      let k = 0, p = 1;
      do { k++; p *= rng(); } while (p > L);
      return k - 1;
    }
    const u1 = rng() || 1e-9, u2 = rng();
    const g = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    return Math.max(0, Math.round(g * Math.sqrt(lambda) + lambda));
  }

  return {
    buildSystemMatrix, computeSensitivity, forwardProject, backProject,
    mlemStep, poissonLogLikelihood,
    dft, idftReal, rampFilterHann, computeFBPSpectra, fbpReconstruct,
    mulberry32, poisson,
  };
});
