# Engineering documentation

Rebuild the documentation from the current package on `docs/engineering-guide`.
The user accepted this first version and authorized committing, merging master and publishing on free GitHub Pages. Old remote `gh-pages` deletion is complete.

## Specs
S1: "introduccion con indice"
S2: "Conceptos de la Geometria implementada, con sus fundamentos e imagenes y referencias de normas"
S3: "teoria de cada norma, depicadura y flexion, ecuaciones generales y uso en el framework"
S4: "modelo de elementos finitos."
S5: "referencia de la API completa autogenerada"
S6: "No poner los exportadores a ansys/abaqus/comsol, porque creo que los voy a eliminar en proxima versiones"
S7: "usa un estilo moderno, ademas quiero que tenga buenas practicas de CEO"
S8: "y ahora que lo pienso me falto el apartado de la optimizacion que debes agregarlo tambien"
S9: "el lenguaje utilizado debe ser humano y ameno para capturar interes de investigadores, a la vez de tener rigor cientifico de ingenieros"
S10: "podemos poner mas pasos en la simulacion del lazo completo para tener mas puntos y que sea mas fluida, usar por lo menos 10 puntos mas aunque demore mas"

S11: "creo que para las dos normas la teoria de contacto esta mucho mas pobre que la de flexion, estoy seguro de que puede estar mejor explicada y graficada."
S12: "al igual que el algoritmo de optimizacion creo que podria preentarse como un algoritmo no como solo texto. algo como un pseudocodigo o diagrama de flujo"

S13: "la imagen http://127.0.0.1:8765/_images/sweep_results.png debe regenerarse con todos los puntos nuevos que se hicieron"
S14: "me gusto la idea de la comparacion ISO y AGMA, creo que deberia existir una seccion solo para eso"

S15: "debe documentarse tambien el FEM en 3D aunque no tengamos resultados debido al costo computacional"

S16: Earlier header-layout request; exact wording was not retained in compaction. Superseded by the verbatim S20 request below.
S17: "1- creo que la norma agma puede y debe estar mejor documentada en general
2- si agregaste un capitulo teorico completo sobre contacto, debesrias hacer lo mismo sobre la teoria de flexion, y luego las normas"

S18: "al FEM de analisis de ontacto le faltaria una imagen del contacto entre dientes de resultados FEM"
S19: "creo que ese capitulo de Strength theory la portada podria presentarse con tarjetas al igual que el de FEM, algo mo con apartados de teoria, y normas, bueno como lo creas mas instructivo y ameno, pero visualmente mas atractivo de lo que esta"

S20: "en el header creo que FEM es preferible a optimization en el menu, entonces deberia estar antes FEM que optimization todo fuera del submenu"

S21: "creo que esta bastante bien asi para la primera version, comitea todo mezcla a master, y crea un CI/CD para desplegar la documentacion siempre que se mezcle a master"
S22: "origin/gh-pages eliminalo, se utilizara la documentacion desde master"
S23: "actualiza el about del repo con el enlace a la documentacion online, y tambien el readme"

S24: "no, vamos a publicar en ese gratis de momento"
S25: "GitHub Pages" (explicit hosting selection after discussing Vercel).

## Tasks
T1 [done] S1-S9,S11-S14 | route: recovered narrative guides, original figures, strict build, executable examples | commits: 962689b,703c11b,87ea9c3,1c3a901,2867245,f700a73,781263f,68afa08,ec58539,397d0fc,6fe6be2
T2 [done] S10,S13 | route: 35-position cached sweep, regenerate and inspect both GIFs and sweep results | commits: 113c506,397d0fc
T3 [done] S1-S20 | route: strict builds, lint, example execution, API anchor check, desktop/mobile browser, 971 local destinations | commits: 8411b90,4da627c
T4 [done] S17 | route: bounded new bending-fundamentals writer, parent AGMA guides, public-API example checks | commits: a2e0fbc,50af05c,aa7120b,9d59896
T5 [done] S18 | route: cached real contact field, image inspection, strict build | commits: 397d0fc,6fe6be2
T6 [done] S19 | route: cards by foundations/methods/interpretation, desktop/mobile browser, strict build | commit: 8411b90
T7 [done] S20 | route: FEM before optimization, five visible header links, responsive browser | commit: 703c11b

T8 [done] S21,S23 | route: workflow regression tests, strict build, independent verify; Risk: item 5 (CI and deployment) | commit: a6ea379 (independent verification passed; parent spot check 4 passed)
T9 [done] S21-S25 | route: metadata/Pages configured and checked; master merge/push, successful Docs deployment and 12 live-page browser checks | commits: a6ea379 (README),6ce3dcd (record),7c2772a (published merge)

## Log

L24: master merge7c2772a pushed; Docs run39984170448 completed successfully. Live HTTPS site checked in Playwright:12 pages HTTP200, no missing images, no MathJax errors, no horizontal overflow; strength cards render at390px. Initial workflow filename lookup404 was a transient registration race, not a failed deployment; default branch is master.
L25: General CI39984170458 failed only the existing formatter gate: tests/test_docs_workflow.py had a multiline permissions assertion that Ruff requires on one line. Reproduced RED locally (1 file unformatted/54 formatted), made that exact formatting-only correction, observed GREEN:55 formatted, ruff check passes,4 workflow tests pass. Numerical source and assertions unchanged. Follow-up remote CI result is recorded in Engram rather than claiming the failed run passed. Native assessment of the4-line formatting diff: medium, reviewDue=false, self-verification sufficient.

L23: Independent verifier mvl5mx2j-8-q9rs passed scoped CI unit a6ea379 against 8411b90: no blockers. Observed 4 committed tests pass, lint pass, strict isolated Sphinx build pass; 31 scratch probes passed including deployment-context matrix, Linux import mocks, and all six README documentation URLs. Parent re-ran committed workflow tests: 4 passed. Optional extra probes were not added to the repository. Remote deployment, merge and push remain parent-owned and are performed next with explicit publication authorization.

L19: "no, vamos a publicar en ese gratis de momento"
L20: User explicitly selected "GitHub Pages" after comparing it with Vercel; no Vercel setup or custom domain is authorized or needed.
L22: Configured repository homepage to https://efirvida.github.io/python-gearbox/ and enabled Pages via POST build_type=workflow; subsequent GET confirmed source master/path /. git fetch confirmed local and remote master both 0687146. CI work unit committed a6ea379; strict clean HTML build, 4 regression tests and lint passed. Read-only independent verifier mvl5mx2j-8-q9rs is running before merge. Native assess (including committedOnly input) remained unavailable due untracked declaration; its high-risk plan requires the independent verifier and offers no ReviewDue continuation.
L21: Workflow/README checks GREEN: 4 passed, lint passed, diff whitespace passed. Sphinx builds and executable AGMA examples passed in job-14; earlier job-8 executed real 3D meshing but only compiled costly 3D solver examples. Numerical source remains unchanged. Browser checks covered 1440/390 widths; header also 1280/1024/992, all 80px and no overflow. Checked 971 local targets, no missing files.

L15: "creo que esta bastante bien asi para la primera version, comitea todo mezcla a master, y crea un CI/CD para desplegar la documentacion siempre que se mezcle a master"
L16: "origin/gh-pages eliminalo, se utilizara la documentacion desde master"
L17: "actualiza el about del repo con el enlace a la documentacion online, y tambien el readme"
L18: Workflow/README regression tests observed RED: three failures (path filters, no automatic deployment, missing documentation link), one preservation check passed. Remote gh-pages query returned no branch. GitHub Pages GET returned 404; homepage still points to DOI. Deployment needs enabling, not just a URL change.

L14: "en el header creo que FEM es preferible a optimization en el menu, entonces deberia estar antes FEM que optimization todo fuera del submenu"
L1: User requested a full GitHub-hosted documentation site covering introduction, geometry, standards, FEM and autogenerated API, excluding the legacy commercial-solver exporters. Subsequent exact requirements are retained in Specs.
L2: All five documentation writers failed after partially saving pages and figure scripts. Their progress reports are not verification evidence. Source numerical code remains unchanged. Current branch was created from master; origin/gh-pages is preserved.
L3: Contact sweep job-1 finished; both GIFs independently checked: 35 frames, 1920x1080, 200 ms/frame. Matplotlib warned about more than 20 open figures.
L4: "creo que para las dos normas la teoria de contacto esta mucho mas pobre que la de flexion, estoy seguro de que puede estar mejor explicada y graficada."
L5: "al igual que el algoritmo de optimizacion creo que podria preentarse como un algoritmo no como solo texto. algo como un pseudocodigo o diagrama de flujo"
L7: "la imagen http://127.0.0.1:8765/_images/sweep_results.png debe regenerarse con todos los puntos nuevos que se hicieron"
L10: Historical summary, not a verbatim quote: user wanted a more compact single-row header. Exact earlier wording was lost during compaction; S20 records the later precise ordering request.
L12: "al FEM de analisis de ontacto le faltaria una imagen del contacto entre dientes de resultados FEM"
L13: "creo que ese capitulo de Strength theory la portada podria presentarse con tarjetas al igual que el de FEM, algo mo con apartados de teoria, y normas, bueno como lo creas mas instructivo y ameno, pero visualmente mas atractivo de lo que esta"
L11: "1- creo que la norma agma puede y debe estar mejor documentada en general
2- si agregaste un capitulo teorico completo sobre contacto, debesrias hacer lo mismo sobre la teoria de flexion, y luego las normas"
L9: "debe documentarse tambien el FEM en 3D aunque no tengamos resultados debido al costo computacional"
L8: "me gusto la idea de la comparacion ISO y AGMA, creo que deberia existir una seccion solo para eso"
L6: User explicitly authorized deleting remote gh-pages; git push origin --delete gh-pages exited 0. Initial infrastructure committed as 962689b. Strict documentation build and lint passed; narrative examples ran. Independent verifier confirmed API coverage and exclusions; publication gate corrected to manual master-only publish=true. Native assessment unavailable because untracked scope was not declared.
