/*
 * MLEM property tests. Run with: node recon-mlem.test.js
 */
const H = require("./test-helpers.js");
const { ReconMath, mat, sens, NPIX, N_RAYS, phantomShepp } = H;
const check = H.makeCheck();

const phantom = phantomShepp();
const truthProj = ReconMath.forwardProject(mat, phantom);

// ---------------------------------------------------------------------------
// 1. Log-likelihood never decreases. This is MLEM's core mathematical
// guarantee (it's an EM algorithm — every step is guaranteed to not decrease
// the data likelihood) — a dip here means the update rule itself is wrong,
// not just "unlucky noise."
// ---------------------------------------------------------------------------
{
  const scale = 1.2;
  const rng = ReconMath.mulberry32(1);
  const measured = new Float32Array(N_RAYS);
  for (let i = 0; i < N_RAYS; i++) measured[i] = ReconMath.poisson(truthProj[i] * scale, rng);

  let est = new Float32Array(NPIX).fill(0.15);
  let prevLL = ReconMath.poissonLogLikelihood(mat, measured, est);
  let minDelta = Infinity, worstIter = -1;
  const nIter = 60;
  for (let it = 0; it < nIter; it++) {
    est = ReconMath.mlemStep(mat, measured, est, sens, NPIX);
    const ll = ReconMath.poissonLogLikelihood(mat, measured, est);
    const delta = ll - prevLL;
    if (delta < minDelta) { minDelta = delta; worstIter = it; }
    prevLL = ll;
  }
  check(
    "log-likelihood is non-decreasing every iteration",
    minDelta > -1e-4,
    "smallest per-step delta=" + minDelta.toExponential(3) + " (at iter " + worstIter + ")"
  );
}

// ---------------------------------------------------------------------------
// 2. Noise-free data doesn't overfit: with the exact (noiseless) forward
// projection as "measured," RMSE should fall more or less monotonically —
// no U-shape. This is the control that shows the U-curve in the live demo
// comes from Poisson noise, not from a bug in the reconstruction code.
// ---------------------------------------------------------------------------
{
  const measuredClean = truthProj; // no Poisson draw at all
  let est = new Float32Array(NPIX).fill(0.15);
  function rmse(e) {
    let s = 0, n = 0;
    for (let j = 0; j < NPIX; j++) {
      if (sens[j] < 1e-3) continue;
      const d = e[j] - phantom[j];
      s += d * d; n++;
    }
    return Math.sqrt(s / n);
  }
  const nIter = 100;
  const errs = [rmse(est)];
  for (let it = 0; it < nIter; it++) {
    est = ReconMath.mlemStep(mat, measuredClean, est, sens, NPIX);
    errs.push(rmse(est));
  }
  const minErr = Math.min(...errs);
  const finalErr = errs[nIter];
  check(
    "noise-free RMSE settles near its minimum instead of overfitting",
    finalErr <= minErr * 1.15,
    "min=" + minErr.toFixed(5) + " final=" + finalErr.toFixed(5) + " ratio=" + (finalErr / minErr).toFixed(3)
  );
}

// ---------------------------------------------------------------------------
// 3. Non-negativity: MLEM's multiplicative update can't push a pixel below
// zero from a positive start with non-negative data/weights — but a stray
// subtraction or sign error would break that invariant immediately.
// ---------------------------------------------------------------------------
{
  const scale = 4;
  const rng = ReconMath.mulberry32(2);
  const measured = new Float32Array(N_RAYS);
  for (let i = 0; i < N_RAYS; i++) measured[i] = ReconMath.poisson(truthProj[i] * scale, rng);

  let est = new Float32Array(NPIX).fill(0.15);
  let minPixel = Infinity;
  for (let it = 0; it < 40; it++) {
    est = ReconMath.mlemStep(mat, measured, est, sens, NPIX);
    for (let j = 0; j < NPIX; j++) if (est[j] < minPixel) minPixel = est[j];
  }
  check("no pixel ever goes negative", minPixel >= -1e-8, "min pixel value=" + minPixel.toExponential(3));
}

// ---------------------------------------------------------------------------
// 4. Count preservation: sum(A * x_k) should match sum(measured) after every
// iteration (k>=1) — a telescoping identity of the MLEM update, not just an
// asymptotic property, so it should hold tightly from the very first step.
// ---------------------------------------------------------------------------
{
  const scale = 4;
  const rng = ReconMath.mulberry32(3);
  const measured = new Float32Array(N_RAYS);
  for (let i = 0; i < N_RAYS; i++) measured[i] = ReconMath.poisson(truthProj[i] * scale, rng);
  const totalMeasured = H.sum(measured);

  let est = new Float32Array(NPIX).fill(0.15);
  let maxRelErr = 0, worstIter = -1;
  for (let it = 0; it < 40; it++) {
    est = ReconMath.mlemStep(mat, measured, est, sens, NPIX);
    const predictedTotal = H.sum(ReconMath.forwardProject(mat, est));
    const relErr = Math.abs(predictedTotal - totalMeasured) / totalMeasured;
    if (relErr > maxRelErr) { maxRelErr = relErr; worstIter = it; }
  }
  check(
    "sum(A*x_k) tracks sum(measured) after every iteration",
    maxRelErr < 1e-4,
    "worst relative error=" + maxRelErr.toExponential(3) + " (iter " + worstIter + ")"
  );
}

console.log(check.failures() === 0 ? "\nAll tests passed." : "\n" + check.failures() + " test(s) FAILED.");
process.exit(check.failures() === 0 ? 0 : 1);
