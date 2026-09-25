/*
 * Tests for the Radon transform math behind Recon Lab (recon-math.js).
 * Plain Node, no framework — run with: node recon-math.test.js
 */
const { buildSystemMatrix, forwardProject, backProject } = require("./recon-math.js");

// Smaller than the artifact's live 64x64/60x64 config, so the suite runs in
// well under a second, but the geometry (FOV, angle spacing) is the same shape.
const N = 32, N_ANGLES = 30, N_BINS = 32, DS = 1.0, RMAX = 24;
const NPIX = N * N, N_RAYS = N_ANGLES * N_BINS;
const cx = (N - 1) / 2, cy = (N - 1) / 2;
const angles = new Float32Array(N_ANGLES);
for (let a = 0; a < N_ANGLES; a++) angles[a] = (a * 180) / N_ANGLES;

const mat = buildSystemMatrix(N, N_ANGLES, N_BINS, DS, RMAX, angles);

let failures = 0;
function check(name, pass, detail) {
  console.log((pass ? "PASS" : "FAIL") + "  " + name + (detail ? "  (" + detail + ")" : ""));
  if (!pass) failures++;
}

function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}
function randomVec(n, seed) {
  // tiny deterministic PRNG so failures are reproducible
  let s = seed >>> 0;
  const rnd = () => ((s = (s * 1103515245 + 12345) >>> 0) / 4294967296);
  const v = new Float32Array(n);
  for (let i = 0; i < n; i++) v[i] = rnd();
  return v;
}

// ---------------------------------------------------------------------------
// 1. Adjoint test: <Ax, y> == <x, A^T y> for random x, y.
// This is the single most important check — if the forward projector and
// backprojector don't agree on geometry (a flipped angle sign, a transposed
// index, a stray normalization factor in only one of the two), this is what
// catches it. MLEM and FBP will both silently misbehave without erroring if
// this drifts, so it's the one test to keep even if all others are cut.
// ---------------------------------------------------------------------------
{
  const x = randomVec(NPIX, 1);
  const y = randomVec(N_RAYS, 2);
  const Ax = forwardProject(mat, x);
  const Aty = backProject(mat, y, NPIX);
  const lhs = dot(Ax, y);
  const rhs = dot(x, Aty);
  const relErr = Math.abs(lhs - rhs) / Math.max(1e-9, Math.abs(lhs));
  check(
    "adjoint: <Ax,y> ~= <x,A^T y>",
    relErr < 1e-5,
    "<Ax,y>=" + lhs.toFixed(4) + " <x,A^Ty>=" + rhs.toFixed(4) + " relErr=" + relErr.toExponential(2)
  );
}

// ---------------------------------------------------------------------------
// 2. Point source: a single bright pixel should trace r = x*cos(theta) +
// y*sin(theta) across the sinogram as theta sweeps — the textbook sinusoid
// that gives the sinogram its name.
// ---------------------------------------------------------------------------
{
  // an off-center pixel, well inside the FOV, so its offset from center is unambiguous
  const px = cx + 5, py = cy - 3;
  const image = new Float32Array(NPIX);
  image[Math.round(py) * N + Math.round(px)] = 1;
  const proj = forwardProject(mat, image);

  const dx = Math.round(px) - cx, dy = Math.round(py) - cy;
  let maxBinErr = 0, worstAngle = -1;
  for (let a = 0; a < N_ANGLES; a++) {
    const th = (angles[a] * Math.PI) / 180;
    const rPredicted = dx * Math.cos(th) + dy * Math.sin(th);
    const binPredicted = rPredicted + (N_BINS - 1) / 2;

    // peak bin actually produced by the forward projector for this angle
    let peakBin = -1, peakVal = -1;
    for (let b = 0; b < N_BINS; b++) {
      const v = proj[a * N_BINS + b];
      if (v > peakVal) { peakVal = v; peakBin = b; }
    }
    const err = Math.abs(peakBin - binPredicted);
    if (err > maxBinErr) { maxBinErr = err; worstAngle = a; }
  }
  check(
    "point source traces r = x*cos(theta) + y*sin(theta)",
    maxBinErr <= 1.0,
    "worst-case peak-bin error=" + maxBinErr.toFixed(2) + " bins at angle index " + worstAngle
  );
}

// ---------------------------------------------------------------------------
// 3. Mass conservation: for a fixed angle, integrating the projection over
// all bins should recover the image's total mass (every ray at that angle
// tiles the image exactly once). Checked at every angle, not just one.
// ---------------------------------------------------------------------------
{
  // a filled disc well inside the FOV, so no mass is clipped by the edge of the grid
  const image = new Float32Array(NPIX);
  const R = N * 0.3;
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      const dx = x - cx, dy = y - cy;
      if (dx * dx + dy * dy <= R * R) image[y * N + x] = 1;
    }
  }
  let imageTotal = 0;
  for (let i = 0; i < NPIX; i++) imageTotal += image[i];

  const proj = forwardProject(mat, image);
  let maxRelErr = 0, worstAngle = -1;
  for (let a = 0; a < N_ANGLES; a++) {
    let colSum = 0;
    for (let b = 0; b < N_BINS; b++) colSum += proj[a * N_BINS + b];
    const relErr = Math.abs(colSum - imageTotal) / imageTotal;
    if (relErr > maxRelErr) { maxRelErr = relErr; worstAngle = a; }
  }
  check(
    "mass conservation: each angle's column sums to the image total",
    maxRelErr < 0.02,
    "worst-case relative error=" + (maxRelErr * 100).toFixed(2) + "% at angle index " + worstAngle
  );
}

console.log(failures === 0 ? "\nAll tests passed." : "\n" + failures + " test(s) FAILED.");
process.exit(failures === 0 ? 0 : 1);
