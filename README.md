# Recon Lab

An interactive PET/CT scanner and image-reconstruction demo: simulate a scan, watch the raw sinogram the scanner actually records, and reconstruct an image from it with either MLEM (iterative) or FBP (filtered backprojection) — including the classic overfitting curve where reconstruction error bottoms out and then rises as an iterative algorithm starts fitting noise instead of signal.

**Live demo:** https://emianaya.github.io/recon-lab/
**Design rationale:** https://claude.ai/artifact/YGwEYkpuiFad93kXWo5hqs

## Files

- `index.html` — the tool itself (self-contained, loads `recon-math.js`)
- `recon-math.js` — the actual Radon transform / MLEM / FBP math, shared between the browser and the test suite (loaded as a plain `<script>` in the browser, `require()`d as a CommonJS module in Node) — this is deliberate: tests exercise the real implementation, not a hand-copied duplicate that could drift out of sync
- `test-helpers.js` — shared test geometry/phantom setup
- `recon-math.test.js`, `recon-mlem.test.js`, `recon-fbp.test.js`, `recon-noise.test.js`, `recon-e2e.test.js` — plain Node test scripts (no framework)
- `recon-ui.test.js` — Playwright end-to-end tests against the actual page

## Running locally

Just open `index.html` in a browser — no build step, no server needed.

## Tests

```
npm test          # math/reconstruction/noise/e2e tests (fast, no browser)
npm install       # only needed once, to pull in Playwright for the UI test
npx playwright install chromium   # one-time browser download
npm run test:ui   # headless Playwright test against index.html
```

What's covered:

- **Adjoint test** — `⟨Ax,y⟩ ≈ ⟨x,Aᵀy⟩` for random x, y. Catches a mismatched forward/backprojector, the most common silent bug in reconstruction code.
- **Point source** — a single bright pixel traces `r = x·cosθ + y·sinθ` across the sinogram.
- **Mass conservation** — each view's projection sums to the image total.
- **MLEM**: log-likelihood is non-decreasing every iteration (its core mathematical guarantee); noise-free data doesn't overfit (control showing the U-curve comes from noise, not a bug); non-negativity; count preservation (`sum(A·x_k) ≈ sum(measured)`, an exact identity of the update rule).
- **FBP**: noiseless high-cutoff RMSE below threshold; a uniform disk reconstructs flat (no cupping/doming); linearity.
- **Noise**: Poisson sampler's mean ≈ variance ≈ expected counts; a seeded RNG reproduces the same sinogram exactly.
- **End-to-end**: low photon counts push the best iteration well before the last (averaged over several seeds); more photons delay it further — the actual physics trend the demo claims, not just unit-level math.
- **UI**: the scan runs automatically on load; dragging the chart scrubs the iteration slider and MLEM image both ways; changing the photon slider triggers a rescan.
