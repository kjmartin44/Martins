@echo off
REM Every analysis this package has, including docking, in one run.
REM
REM Order matters: docking runs FIRST, because the main analysis folds its
REM result in through --docking-results. A docking run that refuses - no
REM structure, or a pocket that cannot recover its own known ligands -
REM writes nothing, and the analysis then proceeds without it rather than
REM failing. That is the intended behaviour, not a problem to work around.
REM
REM Usage:  run_everything.bat            (prompts for everything)
REM         run_everything.bat GENE DIR   (gene and up/down given up front)
setlocal enabledelayedexpansion
cd /d "%~dp0.."

REM --- the interpreter -------------------------------------------------
REM A virtual environment in the checkout wins. Without this the script
REM runs against whatever "python" is first on PATH, which on a machine
REM with several installs is rarely the one the package was installed to.
set "PY=python"
if exist ".venv\Scripts\python.exe" set "PY=.venv\Scripts\python.exe"
if exist "venv\Scripts\python.exe" set "PY=venv\Scripts\python.exe"
echo Interpreter: %PY%

echo.
echo === Checking dependencies ===
"%PY%" -c "import numpy" 2>nul || ("%PY%" -m pip install numpy)
"%PY%" -c "import rdkit, sklearn" 2>nul || ("%PY%" -m pip install rdkit scikit-learn)

REM --- motif file (required) -------------------------------------------
if not exist data\jaspar2024.txt (
    echo === Downloading JASPAR motifs ^(~284 KB^) ===
    curl.exe -fL -o data\jaspar2024.txt "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
)
REM Checked by CONTENT, not by size. A captive portal or an error page can
REM be any length, and a 2 MB HTML apology passes a size test.
findstr /m /c:">MA" data\jaspar2024.txt >nul 2>&1
if errorlevel 1 (
    echo JASPAR file is not a PFM file - delete data\jaspar2024.txt and retry.
    pause
    exit /b 1
)

REM --- regulon: ChIP-derived if built, else curated ---------------------
if not exist data\trrust.tsv (
    echo === Downloading TRRUST regulons ^(optional^) ===
    curl.exe -fL -o data\trrust.tsv "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
)
REM --- dosage tolerance: is there room to move this gene, and is
REM --- moving it safe? Neither file is required and neither is large
REM --- enough to think about, except gnomAD, which is asked for.
if not exist data\ClinGen_gene_curation_list_GRCh38.tsv (
    echo === Downloading ClinGen dosage curation ^(~250 KB^) ===
    curl.exe -fL -o data\ClinGen_gene_curation_list_GRCh38.tsv "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
)
REM Content, not size: the header is a comment line, so a gene column
REM in it is the thing that proves this is the curation list.
findstr /m /c:"Haploinsufficiency Score" data\ClinGen_gene_curation_list_GRCh38.tsv >nul 2>&1
if errorlevel 1 (
    echo   ClinGen file is not the curation list - deleting it.
    del data\ClinGen_gene_curation_list_GRCh38.tsv >nul 2>&1
)

if not exist data\gnomad.v4.1.constraint_metrics.tsv (
    echo.
    echo gnomAD constraint says whether this gene's dose is under
    echo selection. It is 91 MB and takes about a minute. Without it the
    echo report says dosage was not assessed, which is not the same as
    echo saying the gene tolerates any dose.
    set "GETGNOMAD="
    set /p GETGNOMAD="Download it? Y/n: "
    if /i not "!GETGNOMAD!"=="n" (
        echo === Downloading gnomAD v4.1 constraint ^(91 MB^) ===
        curl.exe -fL -o data\gnomad.v4.1.constraint_metrics.tsv "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
        findstr /m /c:"lof.oe_ci.upper" data\gnomad.v4.1.constraint_metrics.tsv >nul 2>&1
        if errorlevel 1 (
            echo   That is not the constraint table - deleting it.
            del data\gnomad.v4.1.constraint_metrics.tsv >nul 2>&1
        )
    )
)

REM --- route of administration: the FDA Orange Book. About 1 MB, read
REM --- from the zip directly. Without it, route is known for roughly
REM --- seventy hand-entered drugs and unknown for the rest, which the
REM --- report says rather than treating unknown as systemic.
set "HAVEOB="
if exist data\EOBZIP.zip set "HAVEOB=1"
if exist data\products.txt set "HAVEOB=1"
if not defined HAVEOB (
    echo.
    echo The FDA Orange Book gives each approved drug's route of
    echo administration ^(oral, injection, topical ...^). It is about 1 MB.
    echo Without it the report cannot tell a topical-only drug from a
    echo systemic one, and shows the route as unknown.
    set "GETOB="
    set /p GETOB="Download it? Y/n: "
    if /i not "!GETOB!"=="n" (
        echo === Downloading the FDA Orange Book ^(~1 MB^) ===
        curl.exe -fL -A "Mozilla/5.0 promoter-tf-drugs" -o data\EOBZIP.zip "https://www.fda.gov/media/76860/download?attachment"
        "%PY%" -c "import sys,zipfile; sys.exit(0 if any(n.lower().endswith('products.txt') for n in zipfile.ZipFile(sys.argv[1]).namelist()) else 1)" data\EOBZIP.zip >nul 2>&1
        if errorlevel 1 (
            echo   That is not the Orange Book - deleting it.
            del data\EOBZIP.zip >nul 2>&1
        )
    )
)

set "REGULON="
if exist data\trrust.tsv set "REGULON=--regulon-file data\trrust.tsv"
REM A ChIP-derived table covers far more factors. It is occupancy rather
REM than curation and the report says exactly that, which is why it is
REM still preferred: more factors testable, with the weaker claim stated.
if exist data\chip_regulon.tsv set "REGULON=--regulon-file data\chip_regulon.tsv"
if "!REGULON!"=="" echo   ^(no regulon table - using bundled tiers^)

REM --- inputs ----------------------------------------------------------
set "GENE=%~1"
if "!GENE!"=="" set /p GENE="Gene symbol: "
if "!GENE!"=="" (
    echo No gene entered.
    pause
    exit /b 1
)
REM Letters, digits, hyphen, underscore. Anything else is a typo, or an
REM attempt to put a second command on the line.
REM
REM The value is echoed QUOTED, and the pattern matches the quotes. That
REM is not decoration: unquoted, `echo !GENE!| findstr ...` with a gene
REM of `CCR5 & del *` runs as `echo CCR5` and then `del *| findstr ...`,
REM so findstr reads del's output, errorlevel comes back 0, and the
REM value is accepted because it contained a separator. Delayed
REM expansion does not protect a pipeline. Quoting does.
echo "!GENE!"| findstr /r /c:"^\"[A-Za-z0-9_-][A-Za-z0-9_-]*\"$" >nul
if errorlevel 1 (
    echo "!GENE!" is not a gene symbol.
    pause
    exit /b 1
)

set "DIR=%~2"
if "!DIR!"=="" set /p DIR="Direction - up or down [up]: "
if "!DIR!"=="" set "DIR=up"
if /i not "!DIR!"=="up" if /i not "!DIR!"=="down" (
    echo Direction must be up or down, not "!DIR!".
    pause
    exit /b 1
)

set "PROLIF=unknown"
set /p PROLIF="Tissue state - cycling / post-mitotic / mixed / unknown [unknown]: "
if "!PROLIF!"=="" set "PROLIF=unknown"

set "TISS="
set /p CNS="Is the target tissue in the CNS? y/N: "
if /i "!CNS!"=="y" set "TISS=--tissue cns"

REM --- which tissue to infer factor activity in ------------------------
REM Spelled as the expression table spells it: "liver", "lymph node",
REM "skeletal muscle". Left blank, --tf-activity is not passed at all,
REM because the CLI refuses that combination up front rather than failing
REM a quarter of an hour into the run.
set "ACTIVITY="
set /p ATISSUE="Tissue to infer factor activity in, blank to skip: "
if not "!ATISSUE!"=="" set "ACTIVITY=--tf-activity --activity-tissue "!ATISSUE!""

REM --- the shape of the dosage problem ---------------------------------
REM This fixes what counts as enough, and the ranking below cannot tell
REM you it: a list of two hundred routes looks the same whether the plan
REM needs a doubling or a halving. Declared rather than inferred, because
REM nothing in this package reads a variant - guessing here would be
REM inventing the premise of the answer.
REM
REM Skipping is a real option and the default. The report then ranks
REM routes without stating a requirement, which is what it did before
REM this prompt existed.
echo.
echo What is the shape of the problem for !GENE!?
echo    1  haploinsufficiency  - one copy broken, raise the other
echo    2  hypomorph           - both copies weak
echo    3  toxic gain          - lower a harmful product
echo    4  null                - no working product at all
echo    5  silenced allele     - an intact copy switched off
echo    6  modifier            - no mutation; moving a normal gene helps
echo    0  skip                - rank routes, state no requirement
set "SNUM="
set /p SNUM="Choose [0]: "

set "SCEN="
if "!SNUM!"=="1" set "SCEN=haploinsufficiency"
if "!SNUM!"=="2" set "SCEN=hypomorph"
if "!SNUM!"=="3" set "SCEN=toxic_gain"
if "!SNUM!"=="4" set "SCEN=null"
if "!SNUM!"=="5" set "SCEN=silenced_allele"
if "!SNUM!"=="6" set "SCEN=modifier"
if not "!SNUM!"=="" if not "!SNUM!"=="0" if "!SCEN!"=="" (
    echo "!SNUM!" is not one of the choices.
    pause
    exit /b 1
)

REM A null allele is the one case where the answer is known before the
REM run: there is no working copy to raise, so no amount of
REM transcription produces a working product. The report says so and
REM the routes become context rather than candidates. Worth saying here
REM too, rather than after half an hour of docking and scanning.
if "!SCEN!"=="null" (
    echo.
    echo   With no working copy, raising transcription cannot produce a
    echo   working product. The report will say so and treat every route
    echo   below it as context rather than as a candidate.
    set "GOON="
    set /p GOON="Run anyway? y/N: "
    if /i not "!GOON!"=="y" exit /b 0
)

REM A scenario fixes which way transcription has to move, so a direction
REM that disagrees is refused here rather than after the run: the CLI
REM refuses the same pair, but only once it is started.
set "NEEDDIR="
if "!SCEN!"=="haploinsufficiency" set "NEEDDIR=up"
if "!SCEN!"=="hypomorph" set "NEEDDIR=up"
if "!SCEN!"=="silenced_allele" set "NEEDDIR=up"
if "!SCEN!"=="toxic_gain" set "NEEDDIR=down"
if defined NEEDDIR if /i not "!DIR!"=="!NEEDDIR!" (
    echo.
    echo The scenario "!SCEN!" needs direction "!NEEDDIR!", but "!DIR!" was chosen.
    echo The requirement and the ranking would point opposite ways. Start again
    echo with the direction that matches, or choose another scenario.
    pause
    exit /b 1
)

REM What the broken copy makes. This matters because transcription is
REM not allele-specific: raising the gene raises BOTH copies, so a
REM truncated product gets amplified alongside the working one, and the
REM requirement above stops being the whole specification.
set "BROKEN="
set "ASKBROKEN="
if "!SCEN!"=="haploinsufficiency" set "ASKBROKEN=1"
if "!SCEN!"=="silenced_allele" set "ASKBROKEN=1"
if defined ASKBROKEN (
    echo.
    echo What does the broken copy make?
    echo    1  nothing - no transcript or no protein
    echo    2  a truncated or altered product
    echo    3  unknown
    set "BNUM="
    set /p BNUM="Choose [3]: "
    set "BMAKES=unknown"
    if "!BNUM!"=="1" set "BMAKES=null"
    if "!BNUM!"=="2" set "BMAKES=truncated"
    set "BROKEN=--broken-allele-makes !BMAKES!"
)

REM For a lowering plan, how much you want left. Half is the usual ask
REM and is the default the package uses.
set "FRACTION="
if "!SCEN!"=="toxic_gain" (
    set "TF="
    set /p TF="Share of the current level you want left, e.g. 0.5 [0.5]: "
    if not "!TF!"=="" (
        echo "!TF!"| findstr /r /c:"^\"0\.[0-9]*[1-9][0-9]*\"$" >nul
        if errorlevel 1 (
            echo "!TF!" is not a fraction above 0 and below 1, such as 0.5.
            pause
            exit /b 1
        )
        set "FRACTION=--target-fraction !TF!"
    )
)

REM For a gene with no mutation in it, what is being treated. Letters,
REM digits, spaces and hyphens only - this goes onto a command line.
set "CONDITION="
if "!SCEN!"=="modifier" (
    set "COND="
    set /p COND="Condition being treated, blank to skip: "
    if not "!COND!"=="" (
        echo "!COND!"| findstr /r /c:"^\"[A-Za-z0-9][A-Za-z0-9 _-]*\"$" >nul
        if errorlevel 1 (
            echo "!COND!" has characters that cannot go on a command line.
            pause
            exit /b 1
        )
        set "CONDITION=--condition "!COND!""
    )
)

set "DOSAGE="
if not "!SCEN!"=="" set "DOSAGE=--scenario !SCEN! !BROKEN! !FRACTION! !CONDITION!"

set "OUT=results\!GENE!_!DIR!"

echo.
echo === Every method, for !GENE! ^(!DIR!regulation^) ===
echo Output: !OUT!
echo.

REM --- local data, each gated on actually being there ------------------
set "CHIP=--chip-api"
if exist data\remap2022.bed.gz set "CHIP=--chip-bed data\remap2022.bed.gz"
if "!CHIP!"=="--chip-api" echo   ^(no local ReMap BED - falling back to the API^)

set "EXPR="
if exist data\rna_tissue_consensus.tsv set "EXPR=--expression-file data\rna_tissue_consensus.tsv"
if "!EXPR!"=="" echo   ^(no local HPA table - using the Open Targets API^)

set "GTF="
if exist data\gencode.basic.gtf.gz set "GTF=--gtf data\gencode.basic.gtf.gz"

set "LINCS=--lincs"
if exist data\lincs_chem_pert_up.gmt if exist data\lincs_chem_pert_down.gmt set "LINCS=--lincs --lincs-gmt data\lincs_chem_pert_up.gmt data\lincs_chem_pert_down.gmt"
if exist data\lincs_geneinfo.txt set "LINCS=!LINCS! --lincs-geneinfo data\lincs_geneinfo.txt"

set "LIBRARY="
if exist data\approved_library.tsv set "LIBRARY=--library data\approved_library.tsv"
if exist data\repurposing_drugs.txt set "LIBRARY=--library data\repurposing_drugs.txt"

set "SENS="
if exist data\prism_lfc.csv set "SENS=--sensitivity-scan"
if exist data\Repurposing_Public_24Q2_LFC_COLLAPSED.csv set "SENS=--sensitivity-scan"

REM ---------------------------------------------------------------------
REM 1. Docking, against the gene's own protein.
REM ---------------------------------------------------------------------
REM First, because the analysis below folds the result in. It refuses more
REM often than it succeeds and that is correct: no structure carrying a
REM co-crystallised ligand, or a pocket whose docking cannot recover its
REM own known ligands, means the scores would carry no information.
REM
REM --state any because the conformation classifier reads coactivator and
REM corepressor peptides, which only nuclear receptors carry. For a GPCR,
REM a kinase or an enzyme it cannot be read at all, and the result then
REM records the conformation as unverified so no direction is claimed.
REM Several compounds at once beats one compound on many cores, by about
REM 2.5x: a docking program parallelises a single job with heavy
REM diminishing returns. Measured on a 40-core machine, eight compounds
REM on four cores each was 2.0s per compound against 11.4s for one on
REM sixteen, with identical scores.
REM
REM Sized to leave the machine with headroom rather than to fill it:
REM jobs x cores lands near 80%% of the cores, because oversubscribing
REM makes every job slower without finishing any sooner. Memory is not
REM the constraint - each docking process is around a gigabyte - and on
REM a GPU neither is the card, at roughly 900 MB per process.
set "JOBS=2"
set "CORES=2"
if defined NUMBER_OF_PROCESSORS (
    set /a JOBS=%NUMBER_OF_PROCESSORS% / 5
    if %NUMBER_OF_PROCESSORS% GEQ 16 (set "CORES=4") else (
        if %NUMBER_OF_PROCESSORS% GEQ 8 (set "CORES=2") else set "CORES=1")
)
if !JOBS! LSS 1 set "JOBS=1"
if !CORES! LSS 1 set "CORES=1"

echo === 1 of 4: docking the !GENE! protein ===
echo   ^(needs gnina or smina. Skipped cleanly when neither is installed.^)
echo   !JOBS! compound^(s^) at a time on !CORES! core^(s^) each; 600 compounds
echo   by default, which is roughly twenty minutes once a pocket validates.
if not exist "!OUT!" mkdir "!OUT!"
set "DOCKDIR=!OUT!\docking"
if not exist "!DOCKDIR!" mkdir "!DOCKDIR!"
set "DOCKJSON=!DOCKDIR!\product.json"
if exist "!DOCKJSON!" del "!DOCKJSON!"

"%PY%" -m promoter_tf_drugs.docking --product --gene !GENE! ^
    --direction !DIR! --state any !LIBRARY! ^
    --jobs !JOBS! --cpu !CORES! --top 40 ^
    --work-dir "!OUT!\docking_work" --out "!DOCKJSON!" -v
if errorlevel 1 (
    echo.
    echo   Nothing from the product dock. The run continues without it;
    echo   the message above says which of the refusals it was.
    echo.
)

REM ---------------------------------------------------------------------
REM 2. Docking the transcription factors.
REM ---------------------------------------------------------------------
REM This needs to know WHICH factors, and that comes from the analysis -
REM which is the step that wants the docking. So a fast pass first, with
REM every expensive stage off: the motif scan, the ChIP intersection and
REM the annotation, which is a couple of minutes rather than thirty. Its
REM only job is to name the factors.
REM
REM Unlike the product dock, the conformation here is NOT "any". A
REM factor's role and your direction decide which pocket you want
REM stabilised - to lower a gene through an activator you want that
REM activator inactive - and the docking derives it. A factor whose role
REM is unresolved is left out by top_factors.py rather than docked
REM against a guessed conformation.
set "NTF=3"
set /p NTF="How many transcription factors to dock [3, 0 to skip]: "
if "!NTF!"=="" set "NTF=3"

set "TFLIST=!DOCKDIR!\factors.txt"
if exist "!TFLIST!" del "!TFLIST!"

REM Last run's per-factor results, which are NOT overwritten the way
REM product.json is. They are named <FACTOR>.json, the set of factors
REM changes with the scan, and merge_docking.py globs the directory -
REM so a stale STAT3.json joins today's SP1 and NFKB1 and is handed to
REM the analysis as part of one result.
REM
REM Which would be harmless if a docking run meant the same thing every
REM time, and it does not. The conformation comes from the direction
REM and the factor's role: lowering a gene through an activator docks
REM that activator's INACTIVE pocket. Yesterday's "up" run leaves a
REM screen of the opposite pocket under exactly the filename today's
REM run would write.
REM
REM Nothing is lost. That run's merged docking_all.json is still there.
for %%F in ("!DOCKDIR!\*.json") do (
    if /i not "%%~nxF"=="product.json" del "%%F" >nul 2>&1
)

if not "!NTF!"=="0" (
    echo.
    echo === 2 of 4: finding the factors worth docking ===
    echo   A fast scan, expensive stages off. Two to three minutes.
    "%PY%" -m promoter_tf_drugs --gene !GENE! --direction !DIR! ^
        --pfm-file data\jaspar2024.txt !REGULON! ^
        !CHIP! --upstream 1000 --downstream 200 --pvalue 1e-5 !GTF! ^
        --no-grade --no-outcome --no-improvements --no-mechanisms ^
        --no-direct --no-stabilizers ^
        --out "!DOCKDIR!\scan" >nul 2>&1
    if exist "!DOCKDIR!\scan\results.json" (
        REM --out rather than a redirect: cmd and PowerShell disagree
        REM about byte-order marks, and a BOM on the first line makes
        REM the first symbol one that no gene lookup resolves.
        "%PY%" scripts\top_factors.py "!DOCKDIR!\scan\results.json" ^
            --top !NTF! --verbose --out "!TFLIST!"
    ) else (
        echo   The scan produced nothing, so no factor can be chosen.
    )
)

REM Each line is "SYMBOL role". The role is passed explicitly because
REM the docking CLI resolves it from local annotation, which comes back
REM unknown for many factors the scan knows the role of - and without a
REM role the conformation is not filtered at all, which is the whole
REM reason to dock a factor rather than a product. With it, lowering a
REM gene through an activator docks the INACTIVE conformation.
if exist "!TFLIST!" (
    for /f "usebackq tokens=1,2" %%T in ("!TFLIST!") do (
        echo.
        echo   --- docking %%T ^(%%U^) ---
        set "TFROLE="
        if not "%%U"=="unknown" set "TFROLE=--role %%U"
        "%PY%" -m promoter_tf_drugs.docking --gene %%T !TFROLE! ^
            --direction !DIR! !LIBRARY! ^
            --jobs !JOBS! --cpu !CORES! --top 40 ^
            --work-dir "!OUT!\docking_work" ^
            --out "!DOCKDIR!\%%T.json" -v
        if errorlevel 1 echo       %%T produced nothing; continuing.
    )
)

REM One file for the analysis to read. Refusals are kept rather than
REM dropped: a pocket that could not be validated is a result, and the
REM report has somewhere to say so.
REM Written outside !DOCKDIR! on purpose: a merged file inside the
REM pattern would be folded back into itself on a second run.
set "DOCKING="
set "DOCKALL=!OUT!\docking_all.json"
if exist "!DOCKALL!" del "!DOCKALL!"
"%PY%" scripts\merge_docking.py "!DOCKALL!" "!DOCKDIR!\*.json"
if exist "!DOCKALL!" set "DOCKING=--docking-results "!DOCKALL!""

REM ---------------------------------------------------------------------
REM 3. The analysis itself.
REM ---------------------------------------------------------------------
REM A 6 kb window at p=1e-4 across ~880 matrices expects ~740 chance hits.
REM Real promoters do not contain 900 genuine sites, so that combination
REM returns a list that is ~83%% noise whichever gene you run. A core
REM promoter window at a stricter threshold is what makes the scan
REM interpretable. Widen deliberately, and tighten --pvalue if you do.
set "WINDOW=--upstream 1000 --downstream 200 --pvalue 1e-5"

REM --tf-weight 0.7 puts most of each grade on the promoter call, which is
REM only justified once ChIP evidence is present. Without it the weight
REM sits at 0.5 so a motif guess cannot carry a grade.
set "TFW=--tf-weight 0.5"
if not "!CHIP!"=="--chip-api" set "TFW=--tf-weight 0.7"

echo === 3 of 4: analysing !GENE! ===
echo   Ten to thirty minutes on a first run. Checkpoints print as it goes.
echo.

"%PY%" -m promoter_tf_drugs --gene !GENE! --direction !DIR! ^
    --pfm-file data\jaspar2024.txt !REGULON! ^
    !CHIP! !WINDOW! --all-tss ^
    !TFW! --proliferation !PROLIF! --cell-type-check !TISS! !EXPR! !GTF! ^
    --predict-direction --upstream-network --upstream-depth 3 ^
    !LINCS! --genetics --screen --design-screen --design-modality ^
    --promoter-dna --screen-product --find-analogs ^
    !ACTIVITY! !SENS! !DOCKING! !DOSAGE! ^
    --out "!OUT!" -v

if errorlevel 1 (
    echo.
    echo Run failed - see the error above.
    echo.
    echo Common causes:
    echo   * Ensembl temporarily down ^(HTTP 500^) - wait and rerun.
    echo     Everything already fetched is cached, so a rerun resumes quickly.
    echo   * Gene symbol not recognised - check it at https://www.genenames.org
    echo.
    pause
    exit /b 1
)

REM ---------------------------------------------------------------------
REM 4. The narrowing funnel, with a control.
REM ---------------------------------------------------------------------
REM ACTB is a housekeeping gene with no specific pharmacology. A funnel
REM that keeps as many routes for it as for your gene is not
REM discriminating, and without a control an empty shortlist and a
REM misconfigured filter look identical.
echo.
echo === 4 of 4: narrowing, against a control ===
echo   Puts !GENE! and ACTB through the same funnel.
echo.

REM The scan settings are the ones stage 3 used, and they have to be.
REM The funnel narrows a route list; built from a different window it
REM narrows a different population, and the number it reports can no
REM longer be cross-checked against `shortlist --preset specific` on
REM the same gene - which is the check stage 4 exists to support, since
REM stage 4's regulon stage and that preset apply the same criterion.
REM
REM Left out of the defaults this would scan 6 kb at p=1e-4, which the
REM comment above stage 3 explains is ~83%% noise whichever gene you
REM run. Observed: 241 sites against 187 expected by chance.
REM
REM --proliferation and --tissue describe the target tissue rather than
REM either gene, so the control gets them too; that is what makes it a
REM control rather than a second analysis.
"%PY%" -m promoter_tf_drugs.panel --genes !GENE! --control-gene ACTB ^
    --direction !DIR! ^
    --pfm-file data\jaspar2024.txt !REGULON! ^
    !CHIP! !WINDOW! --all-tss !TFW! ^
    --proliferation !PROLIF! !TISS! !EXPR! !GTF! --cell-type-check ^
    !LINCS! --genetics --upstream-network ^
    --promoter-dna --screen-product !ACTIVITY! !SENS! !DOCKING! ^
    --max-regulon 600 --max-burden moderate ^
    --out "!OUT!_panel" -v
if errorlevel 1 (
    echo   The panel failed. The analysis above is unaffected.
)

echo.
echo === Done ===
echo Report:   !OUT!\report.md
echo Summary:  !OUT!\summary.md
echo Table:    !OUT!\tf_summary.tsv
echo Browser:  !OUT!\sites.bed
if exist "!DOCKALL!" echo Docking:  !DOCKALL!
if exist "!TFLIST!" echo Factors docked: see !TFLIST!
if exist "!OUT!_panel\panel.md" echo Funnel:   !OUT!_panel\panel.md
echo.
echo Read the evidence table at the top of report.md before the ranking
echo underneath it. A source that failed is recorded as failed and never
echo as zero, and the verdict says whether a pillar of the run is missing.
pause
