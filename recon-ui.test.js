/*
 * UI tests (Playwright, headless Chromium). Run with: node recon-ui.test.js
 * Loads the artifact directly from disk (file://) — no server needed.
 */
const { chromium } = require("playwright");
const path = require("path");

let failures = 0;
function check(name, pass, detail) {
  console.log((pass ? "PASS" : "FAIL") + "  " + name + (detail ? "  (" + detail + ")" : ""));
  if (!pass) failures++;
}

async function main() {
  const browser = await chromium.launch();
  const page = await browser.newPage();

  const pageErrors = [];
  page.on("pageerror", (err) => pageErrors.push(err.message));
  page.on("console", (msg) => { if (msg.type() === "error") pageErrors.push(msg.text()); });

  const fileUrl = "file://" + path.resolve(__dirname, "recon-lab.html");
  await page.goto(fileUrl);

  // ---------------------------------------------------------------------
  // 1. On load, the scan runs automatically and the chart renders.
  // ---------------------------------------------------------------------
  await page.waitForFunction(() => !document.getElementById("playBtn").disabled, { timeout: 15000 });
  const iterCountText = await page.textContent("#iterCount");
  const errReadoutText = await page.textContent("#errReadout");
  check(
    "scan + reconstruction run automatically on load, no click needed",
    iterCountText.includes("100"),
    "iterCount=" + iterCountText
  );
  check(
    "chart renders real numbers, not the placeholder",
    errReadoutText !== "—" && /\d/.test(errReadoutText),
    "errReadout=" + errReadoutText
  );

  // ---------------------------------------------------------------------
  // 2. Dragging the chart updates the slider + MLEM image, and the reverse.
  // ---------------------------------------------------------------------
  const sliderBefore = await page.$eval("#iterSlider", (el) => el.value);
  const reconBefore = await page.$eval("#reconCanvas", (c) => c.toDataURL());
  const chartLocator = page.locator("#chartCanvas");
  await chartLocator.scrollIntoViewIfNeeded();
  const box = await chartLocator.boundingBox();
  await chartLocator.click({ position: { x: box.width * 0.25, y: box.height * 0.5 } });
  await page.waitForTimeout(150);
  const sliderAfterClick = await page.$eval("#iterSlider", (el) => el.value);
  const reconAfterClick = await page.$eval("#reconCanvas", (c) => c.toDataURL());
  check(
    "clicking the chart moves the iteration slider",
    sliderAfterClick !== sliderBefore,
    "before=" + sliderBefore + " after=" + sliderAfterClick
  );
  check("clicking the chart updates the MLEM image", reconAfterClick !== reconBefore);

  const chartBefore = await page.$eval("#chartCanvas", (c) => c.toDataURL());
  await page.$eval("#iterSlider", (el) => {
    el.value = 50;
    el.dispatchEvent(new Event("input"));
  });
  await page.waitForTimeout(150);
  const chartAfter = await page.$eval("#chartCanvas", (c) => c.toDataURL());
  const iterCountAfterSlider = await page.textContent("#iterCount");
  check("moving the slider (the reverse direction) redraws the chart", chartAfter !== chartBefore);
  check(
    "moving the slider updates the iteration readout to match",
    iterCountAfterSlider.includes("50"),
    "iterCount=" + iterCountAfterSlider
  );

  // ---------------------------------------------------------------------
  // 3. Changing the photon slider triggers a rescan.
  // ---------------------------------------------------------------------
  const meanCountsBefore = await page.textContent("#meanCounts");
  await page.$eval("#doseSlider", (el) => {
    el.value = 3; // "many" end
    el.dispatchEvent(new Event("input"));
  });
  const scannerTagRightAfter = await page.textContent("#scannerTag");
  check(
    "changing the photon slider starts a new scan immediately",
    scannerTagRightAfter !== "idle" && scannerTagRightAfter !== "done",
    "scannerTag=" + scannerTagRightAfter
  );
  await page.waitForFunction(() => !document.getElementById("playBtn").disabled, { timeout: 15000 });
  const meanCountsAfter = await page.textContent("#meanCounts");
  check(
    "rescan completes with a different (higher-dose) mean count",
    meanCountsAfter !== meanCountsBefore,
    "before=" + meanCountsBefore + " after=" + meanCountsAfter
  );

  check("no JS errors were thrown during the whole run", pageErrors.length === 0, pageErrors.join(" | "));

  await browser.close();
  console.log(failures === 0 ? "\nAll tests passed." : "\n" + failures + " test(s) FAILED.");
  process.exit(failures === 0 ? 0 : 1);
}

main();
