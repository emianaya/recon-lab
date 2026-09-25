/*
 * End-to-end behavior tests: does the full pipeline (forward project -> add
 * Poisson noise -> run MLEM) actually show the dose/overfitting trend the
 * demo claims? Averaged over several seeds so a single unlucky noise draw
 * can't make the suite flaky.
 * Run with: node recon-e2e.test.js
 */
const H = require("./test-helpers.js");
const { ReconMath, mat, sens, NPIX, N_RAYS } = H;
const check = H.makeCheck();

const phantom = H.phantomShepp();
const truthProj = ReconMath.forwardProject(mat, phantom);
const N_ITER = 100;
const SEEDS = [1, 2, 3, 4, 5];

function runMLEM(scale, seed) {
  const rng = ReconMath.mulberry32(seed);
  const measured = new Float32Array(N_RAYS);
  for (let i = 0; i < N_RAYS; i++) measured[i] = ReconMath.poisson(truthProj[i] * scale, rng);

  let est = new Float32Array(NPIX).fill(0.15);
  function rmse(e) {
    let s = 0, n = 0;
    for (let j = 0; j < NPIX; j++) {
      if (sens[j] < 1e-3) continue;
      const d = e[j] / scale - phantom[j];
      s += d * d; n++;
    }
    return Math.sqrt(s / n);
  }
  const errs = [rmse(est)];
  for (let it = 0; it < N_ITER; it++) {
    est = ReconMath.mlemStep(mat, measured, est, sens, NPIX);
    errs.push(rmse(est));
  }
  let best = 0, bestVal = Infinity;
  for (let i = 0; i < errs.length; i++) if (errs[i] < bestVal) { bestVal = errs[i]; best = i; }
  return { best, bestVal, final: errs[N_ITER] };
}

function averageOverSeeds(scale) {
  let bestSum = 0, ratioSum = 0;
  for (const seed of SEEDS) {
    const r = runMLEM(scale, seed);
    bestSum += r.best;
    ratioSum += r.final / r.bestVal;
  }
  return { avgBest: bestSum / SEEDS.length, avgRatio: ratioSum / SEEDS.length };
}

// ---------------------------------------------------------------------------
// 1. Low photons: the best iteration lands well before the last one, and
// final RMSE ends up noticeably above the minimum — i.e. the demo's low-dose
// preset actually overfits, on average, not just in one lucky/unlucky run.
// ---------------------------------------------------------------------------
{
  const LOW_DOSE = 1.2;
  const { avgBest, avgRatio } = averageOverSeeds(LOW_DOSE);
  check(
    "low dose: best iteration is well before the last (avg over " + SEEDS.length + " seeds)",
    avgBest < N_ITER * 0.3,
    "avg best iter=" + avgBest.toFixed(1) + " / " + N_ITER
  );
  check(
    "low dose: final RMSE ends up noticeably above the minimum",
    avgRatio > 1.3,
    "avg final/min ratio=" + avgRatio.toFixed(2)
  );
}

// ---------------------------------------------------------------------------
// 2. More photons pushes the best iteration later — the physical trend the
// demo's dose slider claims to demonstrate.
// ---------------------------------------------------------------------------
{
  const LOW_DOSE = 1.2, HIGH_DOSE = 14;
  const low = averageOverSeeds(LOW_DOSE);
  const high = averageOverSeeds(HIGH_DOSE);
  check(
    "higher dose delays the best iteration",
    high.avgBest > low.avgBest,
    "avg best @ low dose=" + low.avgBest.toFixed(1) + ", @ high dose=" + high.avgBest.toFixed(1)
  );
}

console.log(check.failures() === 0 ? "\nAll tests passed." : "\n" + check.failures() + " test(s) FAILED.");
process.exit(check.failures() === 0 ? 0 : 1);
