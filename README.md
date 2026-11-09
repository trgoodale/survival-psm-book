# Survival Analysis & Propensity Scores in R: a Quarto book

A hands-on tutorial for complete beginners covering **Kaplan–Meier curves**,
**Cox regression** and **propensity score matching**, with interactive apps
and self-check quizzes.

Students use the published website. This README is for whoever maintains
and publishes the book.

## Live R exercises

The exercises in Chapters 1–7 and 9–11 run R in the reader's browser, using the
[quarto-live](https://r-wasm.github.io/quarto-live/) extension (webR).
Students get hints, solutions and automatic checking, and need nothing
installed.

- The extension is in `_extensions/r-wasm/live/` (version 0.1.3-dev),
  committed with the book.
- `_quarto.yml` uses format `live-html`. The `webr:` settings list the
  package repositories, and `live:` switches hints, solutions and grading on.
- Each chapter with exercises has a small front-matter block that installs
  `survival` when the page loads, and an `{{< include ... _knitr.qmd >}}`
  line after its setup chunk. Each exercise is a `{webr}` block with
  matching `setup`, hint, solution and `check` blocks. Heavier packages
  (broom, ggsurvfit, survminer, MatchIt, cobalt) install the first time an
  exercise that needs them is run.
- The exercise setup blocks in Chapters 5–7, 9 and 11 contain a copy of `simulate_cohort()`.
  **If you change `R/simulate_cohort.R`, update those copies too.**
- Readers need internet access, to download R and the packages the first
  time.

**Editing an exercise's checker:** a `check` block receives `.result` (the
last value), `.user_code` (the learner's code as text) and `.envir_result`
(the environment after running it). It returns
`list(correct = TRUE/FALSE, message = "...")`. Test any change by running the
solution and some wrong answers.

## Editing and previewing

1. Install R (≥ 4.2), RStudio and Quarto (bundled with recent RStudio).
2. Install the packages:
   ```r
   install.packages(c("survival", "ggsurvfit", "survminer", "MatchIt", "cobalt",
                      "broom", "dplyr", "ggplot2", "tidyr", "patchwork", "shiny"))
   ```
3. Open `survival-psm-book.Rproj` and run `quarto preview` in the terminal.
   This opens the book in your browser and refreshes it as you edit.

The interactive apps and quizzes (Observable JS) only work when the book
is viewed through a web server, such as `quarto preview` or GitHub Pages.
They won't run if you double-click an HTML file.

## Publishing on GitHub Pages

The site is rebuilt automatically by GitHub Actions whenever you push to
the `main` branch.

### How it works

The book uses `freeze: auto`. **R code runs only on your own computer.**
When you render locally, Quarto stores the results of every R chunk in
`_freeze/`. GitHub then builds the pages from those stored results, so the
workflow needs Quarto but not R. That keeps builds fast and stops package
updates from changing numbers that the text quotes.

The one rule: **after changing any R code, run `quarto render` locally and
commit the updated `_freeze/` folder with the chapter.** If you forget, the
build fails (there's no R on GitHub to run the code) and the live site
stays as it was. Text-only edits don't need a local render.

The workflow uses **Quarto 1.9.38** (set in
`.github/workflows/publish.yml`). Use the same version on your computer.

### One-time setup

**1. Install the tools** (once per computer)

- [Git](https://git-scm.com/downloads). On Windows, the default options are fine.
- [Quarto 1.9.38](https://github.com/quarto-dev/quarto-cli/releases/tag/v1.9.38):
  download the installer for your system from the release page's
  *Assets* list (on Windows, the `.msi` file).
  `quarto --version` in a terminal should then print `1.9.38`.
- R and the packages listed in the book's "Before you start" chapter, so
  you can render chapters with R code.
- A [GitHub](https://github.com) account.

**2. Start from a clean copy of the book**

Unzip `survival-psm-book.zip` into a new, empty folder. Check that the
hidden `.github/workflows/publish.yml` file and the `_freeze/` folder are
there.

**3. Check it builds on your computer**

Open a terminal in the folder (in RStudio: open `survival-psm-book.Rproj`,
then the *Terminal* tab) and run:

```
quarto preview
```

The book opens in your browser. Close the preview with Ctrl+C.

**4. Create an empty repository on GitHub**

On GitHub, click **New repository**. Give it a name (for example
`survival-psm-book`), make it **Public** (GitHub Pages needs a public
repository on a free account), and **don't** add a README, licence or
`.gitignore`: the book already has them.

**5. Push the book to GitHub**

In the terminal, in the book's folder, run these one at a time, replacing
`YOUR-USERNAME` and the repository name with yours:

```
git init -b main
git add .
git commit -m "First version of the book"
git remote add origin https://github.com/YOUR-USERNAME/survival-psm-book.git
git push -u origin main
```

The first push opens a browser window to sign in to GitHub.

**6. Create the `gh-pages` branch**

Still in the terminal, run:

```
quarto publish gh-pages
```

Answer **Y** when it asks to publish. This renders the book, creates a
`gh-pages` branch holding the website, and prints the site's address.

**7. Turn on GitHub Pages**

On GitHub, open the repository's **Settings ▸ Pages**. Under *Build and
deployment*, set *Source* to **Deploy from a branch**, choose the
**gh-pages** branch and the **/ (root)** folder, and click **Save**.

Then check **Settings ▸ Actions ▸ General ▸ Workflow permissions** is set
to **Read and write permissions**, so the workflow can update the
`gh-pages` branch.

**8. Check the published site**

After a minute or two, the book is at
`https://YOUR-USERNAME.github.io/survival-psm-book/`. Try a live exercise
in a few chapters: R takes a little longer to start on GitHub Pages than
on your computer.

### Making changes

1. Edit the chapters.
2. If you changed any **R code**, run `quarto render` (or
   `quarto render 03-cox-regression.qmd` for one chapter).
3. Commit and push, including `_freeze/`:

   ```
   git add -A
   git commit -m "Describe your change"
   git push
   ```

4. Watch the build in the repository's **Actions** tab. It takes a couple
   of minutes, and the site updates when it finishes.

### If something goes wrong

- **The Actions build fails while rendering a chapter.** You probably
  changed R code without rendering locally. Run `quarto render`, commit
  `_freeze/` and push again.
- **The build fails with a permissions error.** Check the Workflow
  permissions setting in step 7.
- **The site shows a 404 page.** Wait a few minutes after the first
  publish, and check the Pages settings in step 7.
- **Upgrading Quarto.** Install the new version locally, run `quarto
  render` and check the book, including the live exercises. If all is well,
  change `version:` in `.github/workflows/publish.yml`, commit and push.

## Files students download

These are published with the site and linked from the "Before you start"
page and the chapters. Students download them, so the book never refers to
a local folder:

| Source file | Downloaded as |
|---|---|
| `R/simulate_cohort.R` | `simulate_cohort.R` |
| `data/cohort.csv` | `cohort.csv` |
| `R/pairwise_logrank.R` | `pairwise_logrank.R` |
| `apps/survival-explorer/app.R` | `survival-explorer.R` |
| `apps/psm-explorer/app.R` | `psm-explorer.R` |

If you change `R/simulate_cohort.R`, regenerate the CSV:
`source("R/simulate_cohort.R"); write.csv(simulate_cohort(), "data/cohort.csv", row.names = FALSE)`.

## Contents

| File | What |
|---|---|
| `index.qmd` | Welcome, data, suggested route |
| `00-before-you-start.qmd` | Assumed knowledge with refreshers, R basics, setup and downloads, guide to the boxes |
| `01-survival-data.qmd` | Censoring, `Surv()`. App: censoring explorer |
| `02-kaplan-meier.qmd` | KM, log-rank, pairwise tests. Apps: KM by hand, two-arm simulator |
| `03-cox-regression.qmd` | Cox models, HRs, comparing groups, trend tests, PH checks and what to do if PH fails. App: PH explorer |
| `04-traps.qmd` | Informative censoring, competing risks, immortal time bias. Apps: informative censoring, competing risks, immortal time |
| `05-propensity-scores.qmd` | Confounding, the simulated cohort, propensity scores, matching by hand, calipers. App: matching simulator |
| `06-matching-in-r.qmd` | Matching with MatchIt, checking balance with cobalt |
| `07-putting-it-together.qmd` | Complete matched survival analysis |
| `part-going-further.qmd` | Landing page for the Going further part |
| `08-shiny-apps.qmd` | Guide to the two Shiny apps |
| `09-weighting.qmd` | IPTW, WeightIt, matching with replacement, exact and full matching. App: matching vs weighting |
| `10-non-ph.qmd` | Non-proportional hazards: period-specific HRs, `tt()`, RMST, AFT models. App: RMST time horizon |
| `11-missing-data.qmd` | Missing data: complete-case bias, multiple imputation with mice, MatchThem |
| `a1-cheatsheet.qmd`, `a2-glossary.qmd`, `references.qmd` | Appendices: cheat sheet, glossary, references |
| `js/survtools.js` and `js/extras.js` | KM, log-rank, Cox, matching and quizzes in JavaScript, for the in-book apps |
| `R/`, `data/`, `apps/` | Files students download |

## Licence

© 2026 Tom Goodale. The text, figures and teaching material are licensed
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), and the
code (including `R/`, `apps/`, `js/` and the simulated data in `data/`)
under the MIT licence. See [LICENSE.md](LICENSE.md) for details, including
third-party components.

