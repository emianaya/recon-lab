/*
 * Noise-simulation tests. Run with: node recon-noise.test.js
 */
const H = require("./test-helpers.js");
const { ReconMath } = H;
const check = H.makeCheck();

// ---------------------------------------------------------------------------
// 1. Poisson sanity: over many draws, sample mean and sample variance should
// both land near the expected count (lambda) — the defining property of a
// Poisson distribution. Checked across both code paths in poisson() (Knuth's
// method for small lambda, normal approximation for large).
// ---------------------------------------------------------------------------
{
  const rng = ReconMath.mulberry32(999);
  const nDraws = 20000;
  const lambdas = [3, 15, 200]; // spans both the small-lambda and large-lambda code paths
  let worstMeanErr = 0, worstVarErr = 0, worstLambda = -1;
  for (const lambda of lambdas) {
    let s = 0, s2 = 0;
    for (let i = 0; i < nDraws; i++) {
      const k = ReconMath.poisson(lambda, rng);
      s += k; s2 += k * k;
    }
    const mean = s / nDraws;
    const variance = s2 / nDraws - mean * mean;
    const meanErr = Math.abs(mean - lambda) / lambda;
    const varErr = Math.abs(variance - lambda) / lambda;
    if (meanErr > worstMeanErr) worstMeanErr = meanErr;
    if (varErr > worstVarErr) { worstVarErr = varErr; worstLambda = lambda; }
  }
  check(
    "Poisson sampler: mean ~= variance ~= expected counts",
    worstMeanErr < 0.1 && worstVarErr < 0.1,
    "worst mean error=" + (worstMeanErr * 100).toFixed(2) + "%, worst variance error=" +
      (worstVarErr * 100).toFixed(2) + "% (at lambda=" + worstLambda + ")"
  );
}

// ---------------------------------------------------------------------------
// 2. Seeded RNG: the same seed reproduces the same sinogram exactly, so the
// demo (and any bug report against it) is reproducible.
// ---------------------------------------------------------------------------
{
  const { mat, N_RAYS } = H;
  const phantom = H.phantomShepp();
  const truthProj = ReconMath.forwardProject(mat, phantom);
  const scale = 1.2;

  function sinogramWithSeed(seed) {
    const rng = ReconMath.mulberry32(seed);
    const out = new Float32Array(N_RAYS);
    for (let i = 0; i < N_RAYS; i++) out[i] = ReconMath.poisson(truthProj[i] * scale, rng);
    return out;
  }
  const a = sinogramWithSeed(42);
  const b = sinogramWithSeed(42);
  const c = sinogramWithSeed(43);

  let identical = true;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) { identical = false; break; }
  let anyDifferentFromC = false;
  for (let i = 0; i < a.length; i++) if (a[i] !== c[i]) { anyDifferentFromC = true; break; }

  check("same seed reproduces the exact same sinogram", identical);
  check("a different seed does not (sanity check on the check above)", anyDifferentFromC);
}

console.log(check.failures() === 0 ? "\nAll tests passed." : "\n" + check.failures() + " test(s) FAILED.");
process.exit(check.failures() === 0 ? 0 : 1);
