/*
 * Shared setup for the recon-math test suite. Uses a smaller grid than the
 * live artifact (32x32 / 30 angles vs. 64x64 / 60) purely for test speed —
 * same geometry shape, sub-second suite runtime.
 */
const ReconMath = require("./recon-math.js");

const N = 32, N_ANGLES = 30, N_BINS = 32, DS = 1.0, RMAX = 24;
const NPIX = N * N, N_RAYS = N_ANGLES * N_BINS;
const cx = (N - 1) / 2, cy = (N - 1) / 2;
const angles = new Float32Array(N_ANGLES);
for (let a = 0; a < N_ANGLES; a++) angles[a] = (a * 180) / N_ANGLES;

const mat = ReconMath.buildSystemMatrix(N, N_ANGLES, N_BINS, DS, RMAX, angles);
const sens = ReconMath.computeSensitivity(mat, NPIX);

// A simplified Shepp-Logan-like phantom (same shape family as the artifact's,
// scaled to this grid) — three overlapping ellipses, non-negative.
function phantomShepp() {
  const g = new Float32Array(NPIX);
  const ellipses = [
    { cx: 0, cy: 0, a: 0.69, b: 0.92, ang: 0, val: 0.55 },
    { cx: 0, cy: -0.0184, a: 0.6624, b: 0.874, ang: 0, val: -0.35 },
    { cx: 0.22, cy: 0, a: 0.11, b: 0.31, ang: -18, val: 0.4 },
  ];
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      const u = (x - cx) / (N * 0.47), v = (y - cy) / (N * 0.47);
      let val = 0;
      for (const e of ellipses) {
        const th = (e.ang * Math.PI) / 180;
        const du = u - e.cx, dv = v - e.cy;
        const ur = du * Math.cos(th) + dv * Math.sin(th);
        const vr = -du * Math.sin(th) + dv * Math.cos(th);
        if ((ur * ur) / (e.a * e.a) + (vr * vr) / (e.b * e.b) <= 1) val += e.val;
      }
      g[y * N + x] = Math.max(0, val);
    }
  }
  return g;
}

// A filled disc of given radius/value, for tests that want a flat, simple object.
function phantomDisk(radiusFrac, val) {
  const g = new Float32Array(NPIX);
  const R = N * (radiusFrac === undefined ? 0.35 : radiusFrac);
  const v = val === undefined ? 1 : val;
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      const dx = x - cx, dy = y - cy;
      if (dx * dx + dy * dy <= R * R) g[y * N + x] = v;
    }
  }
  return g;
}

function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}
function sum(a) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i];
  return s;
}
function randomVec(n, seed) {
  const rng = ReconMath.mulberry32(seed);
  const v = new Float32Array(n);
  for (let i = 0; i < n; i++) v[i] = rng();
  return v;
}

// Shared PASS/FAIL logger. Every test file calls process.exit(failures===0?0:1)
// itself so each can also be run standalone.
function makeCheck() {
  let failures = 0;
  function check(name, pass, detail) {
    console.log((pass ? "PASS" : "FAIL") + "  " + name + (detail ? "  (" + detail + ")" : ""));
    if (!pass) failures++;
  }
  check.failures = () => failures;
  return check;
}

module.exports = {
  ReconMath, N, N_ANGLES, N_BINS, DS, RMAX, NPIX, N_RAYS, cx, cy, angles,
  mat, sens, phantomShepp, phantomDisk, dot, sum, randomVec, makeCheck,
};
