/*
 * FBP property tests. Run with: node recon-fbp.test.js
 *
 * Note on scale: FBP's absolute output magnitude depends on this
 * implementation's DFT/ramp-filter normalization convention (it's never
 * calibrated against input units, since the live artifact always
 * renormalizes the image by its own max before display). So the RMSE test
 * below fits a single best scalar (least squares) before comparing to the
 * phantom — that's testing reconstruction *shape* accuracy, which is the
 * property that actually matters, without being tripped up by an arbitrary
 * constant that was never meant to be physically calibrated.
 */
const H = require("./test-helpers.js");
const { ReconMath, mat, NPIX, N_ANGLES, N_BINS, cx, cy, phantomShepp, phantomDisk } = H;
const { mulberry32 } = ReconMath;
const check = H.makeCheck();

// ---------------------------------------------------------------------------
// 1. Noise-free, high cutoff: a clean phantom reconstructs with low
// (scale-fit) RMSE.
// ---------------------------------------------------------------------------
{
  const phantom = phantomShepp();
  const truthProj = ReconMath.forwardProject(mat, phantom);
  const spectra = ReconMath.computeFBPSpectra(truthProj, N_ANGLES, N_BINS);
  const cutoff = 0.9;
  const recon = ReconMath.fbpReconstruct(mat, spectra, cutoff, N_ANGLES, N_BINS, NPIX);

  let num = 0, den = 0;
  for (let j = 0; j < NPIX; j++) { num += recon[j] * phantom[j]; den += recon[j] * recon[j]; }
  const fitScale = num / den;
  let se = 0;
  for (let j = 0; j < NPIX; j++) { const d = fitScale * recon[j] - phantom[j]; se += d * d; }
  const rmse = Math.sqrt(se / NPIX);
  const threshold = 0.15;
  check(
    "noiseless high-cutoff FBP reconstructs with low RMSE",
    rmse < threshold,
    "scale-fit RMSE=" + rmse.toFixed(4) + " threshold=" + threshold
  );
}

// ---------------------------------------------------------------------------
// 2. Uniform disk: the interior should come out flat — no cupping (center
// low) or doming (center high). Checks the ramp filter's normalization is
// applied evenly across spatial frequency, not just "looks roughly right."
// ---------------------------------------------------------------------------
{
  const R = H.N * 0.35;
  const phantom = phantomDisk(0.35, 1);
  const truthProj = ReconMath.forwardProject(mat, phantom);
  const spectra = ReconMath.computeFBPSpectra(truthProj, N_ANGLES, N_BINS);
  const recon = ReconMath.fbpReconstruct(mat, spectra, 0.9, N_ANGLES, N_BINS, NPIX);

  let coreSum = 0, coreN = 0, midSum = 0, midN = 0;
  for (let y = 0; y < H.N; y++) {
    for (let x = 0; x < H.N; x++) {
      const dx = x - cx, dy = y - cy, r = Math.sqrt(dx * dx + dy * dy);
      if (r < R * 0.3) { coreSum += recon[y * H.N + x]; coreN++; }
      else if (r > R * 0.3 && r < R * 0.65) { midSum += recon[y * H.N + x]; midN++; }
    }
  }
  const core = coreSum / coreN, mid = midSum / midN;
  const ratio = core / mid;
  check(
    "uniform disk reconstructs flat (no cupping/doming)",
    Math.abs(ratio - 1) < 0.1,
    "core/mid-band ratio=" + ratio.toFixed(4) + " (1.0 = perfectly flat)"
  );
}

// ---------------------------------------------------------------------------
// 3. Linearity: FBP(a*y1 + b*y2) == a*FBP(y1) + b*FBP(y2). Every stage (DFT,
// filter multiply, IDFT, backprojection) is linear, so this should hold to
// near machine precision — a real failure here would mean something
// nonlinear (a clamp, a normalization by the data itself) snuck in.
// ---------------------------------------------------------------------------
{
  function randomSino(seed) {
    const rng = mulberry32(seed);
    const v = new Float32Array(H.N_RAYS);
    for (let i = 0; i < v.length; i++) v[i] = rng() * 10;
    return v;
  }
  const y1 = randomSino(11), y2 = randomSino(22);
  const a = 1.7, b = -0.6;
  const combined = new Float32Array(y1.length);
  for (let i = 0; i < y1.length; i++) combined[i] = a * y1[i] + b * y2[i];

  const cutoff = 0.6;
  const recon = (y) => ReconMath.fbpReconstruct(mat, ReconMath.computeFBPSpectra(y, N_ANGLES, N_BINS), cutoff, N_ANGLES, N_BINS, NPIX);
  const r1 = recon(y1), r2 = recon(y2), rCombined = recon(combined);

  let maxAbsErr = 0, maxVal = 1e-9;
  for (let j = 0; j < NPIX; j++) {
    maxAbsErr = Math.max(maxAbsErr, Math.abs(a * r1[j] + b * r2[j] - rCombined[j]));
    maxVal = Math.max(maxVal, Math.abs(rCombined[j]));
  }
  const relErr = maxAbsErr / maxVal;
  check(
    "FBP(a*y1 + b*y2) == a*FBP(y1) + b*FBP(y2)",
    relErr < 1e-4,
    "max relative error=" + relErr.toExponential(3)
  );
}

console.log(check.failures() === 0 ? "\nAll tests passed." : "\n" + check.failures() + " test(s) FAILED.");
process.exit(check.failures() === 0 ? 0 : 1);
