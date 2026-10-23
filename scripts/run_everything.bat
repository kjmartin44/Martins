@echo off
REM Every analysis this package has, including docking, in one run.
REM
REM Order matters: docking runs FIRST, because the main analysis folds its
REM result in through --docking-results. A docking run that refuses - no
REM structure, or a pocket that cannot recover its own known ligands -
REM writes nothing, and the analysis then proceeds without it rather than
REM failing. That is the intended behaviour, not a problem to work around.
REM
REM After the five stages, two optional steps read the saved results: a dose
REM check (what blood level a person reaches against the concentration the
REM evidence was seen at), which is on unless declined, and an exploratory dock
REM into a pocket copied onto the AlphaFold model from a related structure that
REM the script finds by itself. The exploratory dock is OFF unless you answer y
REM to its question, because it is unvalidated and takes hours. When the gene's
REM own pocket is refused and you said no, the command to run it later is
REM printed instead.
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

REM --- reference files -------------------------------------------------
REM The files below are git-ignored, so a fresh checkout has none of them.
REM The offers further down fetch a few; scripts\setup_data.bat fetches all
REM of them, and leaves any that are already here alone.
set "NEEDSETUP="
for %%P in (data\jaspar2024.txt data\remap2022.bed.gz data\chip_regulon.tsv data\collectri.tsv data\gencode.basic.gtf.gz data\rna_tissue_consensus.tsv data\lincs_chem_pert_up.gmt data\lincs_chem_pert_down.gmt data\lincs_geneinfo.txt data\approved_library.tsv data\transcript_sites.tsv) do if not exist %%P set "NEEDSETUP=1"
if defined NEEDSETUP (
    echo.
    echo Some reference files are not in data\. Without them routes are skipped
    echo or fall back to web services, and the report says which. The setup
    echo script fetches every missing file and is safe to rerun.
    set "RUNSETUP="
    set /p RUNSETUP="Run scripts\setup_data.bat now? y/N: "
    if /i "!RUNSETUP!"=="y" (
        set "SETUP_NOPAUSE=1"
        call scripts\setup_data.bat --yes
    )
)

REM --- motif file (required) -------------------------------------------
if not exist data\jaspar2024.txt (
    echo === Downloading JASPAR motifs ^(~284 KB^) ===
    REM To a .part, and curl's exit code is read. A dropped connection
    REM leaves a partial file at whatever name curl was given, and the
    REM content check below cannot see that: findstr matches the FIRST
    REM matrix header and a truncated PFM still has one. So a short
    REM motif set was kept for ever under the real name, and every
    REM later run scanned with it and reported fewer sites.
    curl.exe -fL --retry 3 -o data\jaspar2024.txt.part "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    set "JOK="
    if not errorlevel 1 (
        findstr /m /c:">MA" data\jaspar2024.txt.part >nul 2>&1
        if not errorlevel 1 set "JOK=1"
    )
    if defined JOK (
        move /y data\jaspar2024.txt.part data\jaspar2024.txt >nul
    ) else (
        echo   The JASPAR download was incomplete or is not a PFM file.
        del data\jaspar2024.txt.part >nul 2>&1
    )
)
REM Checked by CONTENT, not by size. A captive portal or an error page can
REM be any length, and a 2 MB HTML apology passes a size test.
findstr /m /c:">MA" data\jaspar2024.txt >nul 2>&1
if errorlevel 1 (
    echo The motif file is missing or is not a PFM file. It is required.
    echo Delete data\jaspar2024.txt if it is there, then run this again.
    pause
    exit /b 1
)

REM --- regulon: ChIP-derived if built, else curated ---------------------
if not exist data\trrust.tsv (
    echo === Downloading TRRUST regulons ^(optional^) ===
    REM The same .part treatment. This one had no check at all, and it
    REM is wired in below as --regulon-file: an error page there reads
    REM as a regulon table with nothing in it.
    curl.exe -fL --retry 3 -o data\trrust.tsv.part "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
    set "TROK="
    if not errorlevel 1 (
        findstr /m /c:"Activation" /c:"Repression" data\trrust.tsv.part >nul 2>&1
        if not errorlevel 1 set "TROK=1"
    )
    if defined TROK (
        move /y data\trrust.tsv.part data\trrust.tsv >nul
    ) else (
        echo   TRRUST could not be fetched. Skipped - it is optional.
        del data\trrust.tsv.part >nul 2>&1
    )
)
REM --- dosage tolerance: is there room to move this gene, and is
REM --- moving it safe? Neither file is required and neither is large
REM --- enough to think about, except gnomAD, which is asked for.
if not exist data\ClinGen_gene_curation_list_GRCh38.tsv (
    echo === Downloading ClinGen dosage curation ^(~250 KB^) ===
    REM Content, not size: the header is a comment line, so a gene
    REM column in it is what proves this is the curation list. Checked
    REM on the .part, so a bad download never reaches the real name.
    curl.exe -fL --retry 3 -o data\ClinGen_gene_curation_list_GRCh38.tsv.part "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
    set "CGOK="
    if not errorlevel 1 (
        findstr /m /c:"Haploinsufficiency Score" data\ClinGen_gene_curation_list_GRCh38.tsv.part >nul 2>&1
        if not errorlevel 1 set "CGOK=1"
    )
    if defined CGOK (
        move /y data\ClinGen_gene_curation_list_GRCh38.tsv.part data\ClinGen_gene_curation_list_GRCh38.tsv >nul
    ) else (
        echo   The ClinGen download is not the curation list. Skipped.
        del data\ClinGen_gene_curation_list_GRCh38.tsv.part >nul 2>&1
    )
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
        REM .part, and curl's exit code is read. The content check is
        REM on the HEADER, which a truncated 91 MB download still has,
        REM so on its own it cannot tell a short file from a whole one.
        curl.exe -fL --retry 3 -o data\gnomad.v4.1.constraint_metrics.tsv.part "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
        set "GNOK="
        if not errorlevel 1 (
            findstr /m /c:"lof.oe_ci.upper" data\gnomad.v4.1.constraint_metrics.tsv.part >nul 2>&1
            if not errorlevel 1 set "GNOK=1"
        )
        if defined GNOK (
            move /y data\gnomad.v4.1.constraint_metrics.tsv.part data\gnomad.v4.1.constraint_metrics.tsv >nul
        ) else (
            echo   The constraint table could not be fetched. Skipped.
            del data\gnomad.v4.1.constraint_metrics.tsv.part >nul 2>&1
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
        curl.exe -fL -A "Mozilla/5.0 generoutes" -o data\EOBZIP.zip "https://www.fda.gov/media/76860/download?attachment"
        "%PY%" -c "import sys,zipfile; sys.exit(0 if any(n.lower().endswith('products.txt') for n in zipfile.ZipFile(sys.argv[1]).namelist()) else 1)" data\EOBZIP.zip >nul 2>&1
        if errorlevel 1 (
            echo   That is not the Orange Book - deleting it.
            del data\EOBZIP.zip >nul 2>&1
        )
    )
)

REM --- GWAS Catalog: which traits variants near the gene are linked to. About
REM --- 70 MB from EBI as a zip, unpacked to the table the package reads.
REM --- Without it genetics rests on Open Targets alone.
set "HAVEGWAS="
for %%P in (data\gwas_catalog_associations.tsv data\gwas_catalog.tsv) do if exist %%P set "HAVEGWAS=1"
if not defined HAVEGWAS (
    echo.
    echo The GWAS Catalog lists the traits that variants near a gene are linked
    echo to. It is about 70 MB. Without it the genetics section uses Open
    echo Targets alone and says so.
    set "GETGWAS="
    set /p GETGWAS="Download it? Y/n: "
    if /i not "!GETGWAS!"=="n" (
        echo === Downloading the GWAS Catalog associations ^(~70 MB^) ===
        curl.exe -fL -o data\gwas_catalog_associations.zip "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
        if not errorlevel 1 (
            "%PY%" scripts\unpack_gwas.py data\gwas_catalog_associations.zip data\gwas_catalog_associations.tsv
            if not errorlevel 1 set "HAVEGWAS=1"
        )
        if not defined HAVEGWAS echo   The GWAS Catalog could not be fetched or read. Genetics uses Open Targets alone.
        if exist data\gwas_catalog_associations.zip del data\gwas_catalog_associations.zip >nul 2>&1
    )
)

REM --- DepMap and PRISM, from figshare, where the Broad posts each release as
REM --- open data (CC BY 4.0). The DepMap portal itself asks for a browser check,
REM --- figshare does not. Each file is downloaded under a .part name and kept
REM --- only if its published MD5 matches, because a dropped connection on a
REM --- 400 MB file otherwise leaves something that looks like data.
REM ---
REM --- Both are needed for the measured phenocopy route. DepMap alone also
REM --- scores how much each candidate factor is needed by a cell.
set "HAVEDM0="
for %%P in (data\CRISPRGeneEffect.csv data\crispr_gene_effect.csv data\depmap_gene_effect.csv) do if exist %%P set "HAVEDM0=1"
if not defined HAVEDM0 (
    echo.
    echo DepMap's CRISPR knockout scores say what a cell loses without each
    echo candidate factor. The 24Q2 release is about 420 MB ^(plus a 1 MB table of
    echo cell lines^) from figshare.
    set "GETDM="
    set /p GETDM="Download it? Y/n: "
    if /i not "!GETDM!"=="n" (
        echo === Downloading DepMap 24Q2 CRISPRGeneEffect.csv ^(~420 MB^) ===
        curl.exe -fL -o data\CRISPRGeneEffect.csv.part "https://ndownloader.figshare.com/files/46489063"
        set "DMOK="
        if not errorlevel 1 (
            "%PY%" scripts\check_download.py data\CRISPRGeneEffect.csv.part d10d73692672380b3792231549e97d5a
            if not errorlevel 1 set "DMOK=1"
        )
        if defined DMOK (
            move /y data\CRISPRGeneEffect.csv.part data\CRISPRGeneEffect.csv >nul
            if not exist data\Model.csv curl.exe -fL -o data\Model.csv "https://ndownloader.figshare.com/files/46489732"
        ) else (
            echo   The DepMap table could not be fetched or did not match its checksum. Skipped.
            del data\CRISPRGeneEffect.csv.part >nul 2>&1
        )
    )
)

set "HAVEPR0="
for %%P in (data\Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv data\Repurposing_Public_24Q2_LFC_COLLAPSED.csv data\Repurposing_Public_LFC.csv data\prism_lfc.csv data\PRISM_Repurposing_Public.csv) do if exist %%P set "HAVEPR0=1"
if not defined HAVEPR0 (
    echo.
    echo PRISM measured how about 6,000 drugs change the growth of the same cell
    echo lines. It is what finds drugs that behave as though they remove a factor.
    echo The 24Q2 matrix is about 70 MB from figshare.
    set "GETPR="
    set /p GETPR="Download it? Y/n: "
    if /i not "!GETPR!"=="n" (
        echo === Downloading PRISM 24Q2 Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv ^(~70 MB^) ===
        curl.exe -fL -o data\Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv.part "https://ndownloader.figshare.com/files/46630984"
        set "PROK="
        if not errorlevel 1 (
            "%PY%" scripts\check_download.py data\Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv.part 042b15d3ebcc138a80a571b9286451ae
            if not errorlevel 1 set "PROK=1"
        )
        if defined PROK (
            move /y data\Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv.part data\Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv >nul
            REM The compound list turns the BRD ids in that table into drug names.
            if not exist data\Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv curl.exe -fL -o data\Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv "https://ndownloader.figshare.com/files/46630981"
        ) else (
            echo   The PRISM table could not be fetched or did not match its checksum. Skipped.
            del data\Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv.part >nul 2>&1
        )
    )
)

REM --- Factor knockout signatures, from Enrichr. They show whether your gene
REM --- changes when a transcription factor is knocked out, which replaces a
REM --- fixed placeholder for "responsiveness". About 19 MB for both. Each file
REM --- is downloaded under a .part name and kept only if it reads as a gene
REM --- set library, so an error page is not mistaken for data.
set "HAVEKO="
for %%P in (data\knockout_lincs_crispr_ko.gmt data\knockout_tf_perturbations.gmt) do if exist %%P set "HAVEKO=1"
if not defined HAVEKO (
    echo.
    echo Knockout signatures say whether your gene moves when a factor is lost.
    echo Without them the report cannot say, for any factor. About 19 MB.
    set "GETKO="
    set /p GETKO="Download them? Y/n: "
    if /i not "!GETKO!"=="n" (
        echo === Downloading the Enrichr knockout libraries ^(~19 MB^) ===
        curl.exe -fL -o data\knockout_lincs_crispr_ko.gmt.part "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LINCS_L1000_CRISPR_KO_Consensus_Sigs"
        if not errorlevel 1 "%PY%" scripts\check_gmt.py data\knockout_lincs_crispr_ko.gmt.part 1000
        if not errorlevel 1 (
            move /y data\knockout_lincs_crispr_ko.gmt.part data\knockout_lincs_crispr_ko.gmt >nul
        ) else (
            echo   The CRISPR knockout library could not be fetched or read. Skipped.
            del data\knockout_lincs_crispr_ko.gmt.part >nul 2>&1
        )
        curl.exe -fL -o data\knockout_tf_perturbations.gmt.part "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=TF_Perturbations_Followed_by_Expression"
        if not errorlevel 1 "%PY%" scripts\check_gmt.py data\knockout_tf_perturbations.gmt.part 500
        if not errorlevel 1 (
            move /y data\knockout_tf_perturbations.gmt.part data\knockout_tf_perturbations.gmt >nul
        ) else (
            echo   The TF perturbation library could not be fetched or read. Skipped.
            del data\knockout_tf_perturbations.gmt.part >nul 2>&1
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
set /p CNS="Is the target tissue in the CNS (c), in skin (s), or neither? c/s/N: "
if /i "!CNS!"=="c" set "TISS=--tissue cns"
if /i "!CNS!"=="y" set "TISS=--tissue cns"
if /i "!CNS!"=="s" set "TISS=--tissue skin"

REM The dose check is cheap and reads only the finished results, so it is
REM on unless declined.
set "DOSECHK=1"
set /p DOSECHK_ANS="Run the dose check on the candidates at the end? Y/n: "
if /i "!DOSECHK_ANS!"=="n" set "DOSECHK="

REM The upstream search looks for drugs at the proteins that act on the gene's
REM factors, and by default it covers the 8 best-ranked factors. A larger number
REM adds drugs that reach the gene only through the network, all of them tier 4,
REM and a few more minutes. Left blank, nothing is passed and the default stands.
set "UPFACT="
set "UPF_ANS="
set /p UPF_ANS="Factors to search for upstream drugs, a number, blank for the default 8: "
if not "!UPF_ANS!"=="" (
    echo "!UPF_ANS!"| findstr /r /c:"^\"[1-9][0-9]*\"$" >nul
    if errorlevel 1 (
        echo   "!UPF_ANS!" is not a whole number above 0 - the default 8 is used.
    ) else (
        set "UPFACT=--upstream-factors !UPF_ANS!"
    )
)

REM The exploratory dock is the opposite: off unless asked for. It docks the
REM library into a pocket copied onto an AlphaFold model, the result is not
REM validated, and it can take hours, so it is not something to start by
REM default. It only applies when the gene's own pocket is refused.
set "XDOCK_ANS="
set /p XDOCK_ANS="If the gene's own pocket is refused, also run the exploratory AlphaFold dock? Unvalidated, hours. y/N: "

REM --- which tissue to infer factor activity in ------------------------
REM Spelled as the expression table spells it: "liver", "lymph node",
REM "skeletal muscle". Left blank, --tf-activity is not passed at all,
REM because the CLI refuses that combination up front rather than failing
REM a quarter of an hour into the run.
set "ACTIVITY="
set "ATISSUE="
set "APICK="
set "TISSUETABLE="
if exist data\rna_tissue_consensus.tsv set "TISSUETABLE=data\rna_tissue_consensus.tsv"
if defined TISSUETABLE (
    echo.
    echo Tissues in the expression table:
    "%PY%" scripts\tissues.py !TISSUETABLE!
)
set /p ATISSUE="Tissue to infer factor activity in, a number or a name, blank to skip: "
if not "!ATISSUE!"=="" (
    echo "!ATISSUE!"| findstr /r /c:"^\"[A-Za-z0-9][A-Za-z0-9 -]*\"$" >nul
    if errorlevel 1 (
        echo   "!ATISSUE!" is not a tissue name - letters, digits, spaces and hyphens only.
        set "ATISSUE="
    )
)
if not "!ATISSUE!"=="" (
    if defined TISSUETABLE (
        REM Resolved against the table, so the name that reaches the command is
        REM spelled as the table spells it. A name the table lacks used to log a
        REM warning and leave the route silently empty.
        for /f "usebackq delims=" %%T in (`%PY% scripts\tissues.py !TISSUETABLE! --pick "!ATISSUE!"`) do set "APICK=%%T"
        if not defined APICK echo   Not one tissue in the table - factor activity is skipped.
    ) else (
        set "APICK=!ATISSUE!"
    )
)
REM 2000 permutations build an empirical null in place of the analytic one,
REM which assumes targets are independent and co-expression makes false.
REM The CLI's own help says to use 2000 or more before relying on a number.
if defined APICK set "ACTIVITY=--tf-activity --activity-tissue "!APICK!" --activity-permutations 2000"

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
REM Asked for the recessive scenarios too. A truncated product is what
REM turns the designed screen from a promoter reporter into a protein
REM screen, and that is as true of two impaired copies as of one.
set "BROKEN="
set "ASKBROKEN="
if "!SCEN!"=="haploinsufficiency" set "ASKBROKEN=1"
if "!SCEN!"=="silenced_allele" set "ASKBROKEN=1"
if "!SCEN!"=="hypomorph" set "ASKBROKEN=1"
if "!SCEN!"=="biallelic_hypomorph" set "ASKBROKEN=1"
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

REM A premature stop is worth its own question rather than being inferred
REM from "truncated". It decides two things: the designed screen reads
REM protein instead of a reporter, and readthrough agents become a
REM modality worth listing. Answer no if you do not know.
set "PTC="
if defined ASKBROKEN (
    echo.
    echo Does the transcript carry a premature stop codon ^(nonsense or
    echo frameshift^)? This makes the designed screen measure protein
    echo rather than a promoter reporter.
    set "HASPTC="
    set /p HASPTC="y/N: "
    if /i "!HASPTC!"=="y" set "PTC=--has-ptc"
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

REM Whether a DNA-damaging route is acceptable. Asked rather than assumed:
REM many of the chemicals the literature records as RAISING a gene are
REM carcinogens, tumour promoters or irritants, and for a patient already
REM at raised cancer risk those are not candidates at any score.
echo.
echo Should DNA-damaging compounds be demoted? Say yes for a patient with
echo a cancer predisposition, an existing tumour, or who is a child.
set "GENOTOX="
set "AVOIDG="
set /p AVOIDG="Y/n: "
if /i not "!AVOIDG!"=="n" set "GENOTOX=--avoid-genotoxic"

set "DOSAGE="
if not "!SCEN!"=="" set "DOSAGE=--scenario !SCEN! !BROKEN! !PTC! !FRACTION! !CONDITION!"
if "!SCEN!"=="" if defined PTC set "DOSAGE=!PTC!"

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

REM The measured phenocopy route needs BOTH the PRISM matrix and the DepMap
REM knockout matrix, and the factor dependency cost needs DepMap. Neither is
REM downloaded here (DepMap is about 500 MB and its release names change), so
REM say what is skipped instead of letting the report simply lack the section.
set "HAVEPRISM="
for %%P in (data\Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv data\Repurposing_Public_24Q2_LFC_COLLAPSED.csv data\Repurposing_Public_LFC.csv data\prism_lfc.csv data\PRISM_Repurposing_Public.csv) do if exist %%P set "HAVEPRISM=1"
set "HAVEDEPMAP="
for %%P in (data\CRISPRGeneEffect.csv data\crispr_gene_effect.csv data\depmap_gene_effect.csv) do if exist %%P set "HAVEDEPMAP=1"
set "SENS="
if defined HAVEPRISM if defined HAVEDEPMAP set "SENS=--sensitivity-scan"
if not defined HAVEPRISM echo   ^(no PRISM matrix in data\ - the measured phenocopy route is skipped^)
if not defined HAVEDEPMAP echo   ^(no DepMap CRISPRGeneEffect.csv in data\ - the phenocopy route is skipped and factor dependency is not scored^)
set "NEEDDM="
if not defined HAVEPRISM set "NEEDDM=1"
if not defined HAVEDEPMAP set "NEEDDM=1"
if defined NEEDDM (
    echo   The offers above fetch both from figshare. By hand, DepMap 24Q2 is at
    echo   https://figshare.com/articles/dataset/DepMap_Public_24Q2/25880521 and
    echo   PRISM at https://figshare.com/articles/dataset/Repurposing_Public_24Q2/25917643
    echo   ^(CRISPRGeneEffect.csv and Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv in data\^).
)
if not defined HAVEGWAS echo   ^(no GWAS Catalog table in data\ - genetics uses Open Targets only^)
if not defined ACTIVITY echo   ^(no tissue given for factor activity - that route is skipped^)

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

REM The whole approved library rather than the first 600, and exhaustiveness
REM 16 rather than 8, which the docking help describes as slower and steadier.
REM Both only cost time once a pocket validates; a pocket that cannot be
REM validated is refused before the library is touched.
set "LIBLIMIT=100000"
set "EXH=16"

echo === 1 of 5: docking the !GENE! protein ===
echo   ^(needs gnina or smina. Skipped cleanly when neither is installed.^)
echo   !JOBS! compound^(s^) at a time on !CORES! core^(s^) each; the whole
echo   approved library ^(about 3,300 compounds^) at exhaustiveness !EXH!, which
echo   is a couple of hours once a pocket validates and seconds if it refuses.
if not exist "!OUT!" mkdir "!OUT!"
set "DOCKDIR=!OUT!\docking"
if not exist "!DOCKDIR!" mkdir "!DOCKDIR!"
set "DOCKJSON=!DOCKDIR!\product.json"
if exist "!DOCKJSON!" del "!DOCKJSON!"

"%PY%" -m generoutes.docking --product --gene !GENE! ^
    --direction !DIR! --state any !LIBRARY! ^
    --jobs !JOBS! --cpu !CORES! --top 40 ^
    --library-limit !LIBLIMIT! --exhaustiveness !EXH! ^
    --work-dir "!OUT!\docking_work" --out "!DOCKJSON!" -v
if errorlevel 1 (
    echo.
    echo   Nothing from the product dock. The run continues without it;
    echo   the message above says which of the refusals it was.
    echo.
)

REM Whether the pocket was refused decides whether the exploratory dock is
REM available, and the answer to the question near the top decides whether it
REM runs. When it runs it is at the end, with no further input: it searches the
REM PDB with the AlphaFold model's sequence, takes the closest structure whose
REM ligand transplants onto the model, and labels the result unvalidated.
set "EXPLORE="
set "XUNI=-"
set "PSTATE=refused"
for /f "tokens=1,2" %%A in ('%PY% scripts\product_status.py !DOCKJSON!') do (
    set "PSTATE=%%A"
    set "XUNI=%%B"
)
if /i "!PSTATE!"=="refused" if not "!XUNI!"=="-" (
    echo.
    echo The pocket of the !GENE! protein was refused.
    if /i "!XDOCK_ANS!"=="y" (
        echo An exploratory dock will run at the end instead. It finds a related
        echo structure that holds a ligand by itself and labels the result unvalidated.
        set "EXPLORE=1"
    ) else (
        echo An exploratory dock on an AlphaFold model is available and is not run
        echo by default. To run it later, from the same folder:
        echo   python -m generoutes.exploratory --gene !GENE! --uniprot !XUNI! --reference --out "!OUT!_exploratory"
    )
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
set "NTF=5"
set /p NTF="How many transcription factors to dock [5, 0 to skip]: "
if "!NTF!"=="" set "NTF=5"

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
    echo === 2 of 5: finding the factors worth docking ===
    echo   A fast scan, expensive stages off. Two to three minutes.
    "%PY%" -m generoutes --gene !GENE! --direction !DIR! ^
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
        "%PY%" -m generoutes.docking --gene %%T !TFROLE! ^
            --direction !DIR! !LIBRARY! ^
            --jobs !JOBS! --cpu !CORES! --top 40 ^
            --library-limit !LIBLIMIT! --exhaustiveness !EXH! ^
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

echo === 3 of 5: analysing !GENE! ===
echo   Thirty minutes to a couple of hours on a first run. Checkpoints print as it goes.
echo.

"%PY%" -m generoutes --gene !GENE! --direction !DIR! ^
    --pfm-file data\jaspar2024.txt !REGULON! ^
    !CHIP! !WINDOW! --all-tss ^
    !TFW! --proliferation !PROLIF! --cell-type-check !TISS! !EXPR! !GTF! ^
    --predict-direction --upstream-network --upstream-depth 3 !UPFACT! ^
    --regulator-route --regulator-steps 2 ^
    !LINCS! --genetics --screen --design-screen --design-modality ^
    --promoter-dna --screen-product --find-analogs ^
    !ACTIVITY! !SENS! !DOCKING! !DOSAGE! !GENOTOX! ^
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
echo === 4 of 5: narrowing, against a control ===
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
"%PY%" -m generoutes.panel --genes !GENE! --control-gene ACTB ^
    --direction !DIR! ^
    --pfm-file data\jaspar2024.txt !REGULON! ^
    !CHIP! !WINDOW! --all-tss !TFW! ^
    --proliferation !PROLIF! !TISS! !EXPR! !GTF! --cell-type-check ^
    !LINCS! --genetics --upstream-network !UPFACT! --regulator-route --regulator-steps 2 ^
    --promoter-dna --screen-product !ACTIVITY! !SENS! !DOCKING! ^
    --max-regulon 600 --max-burden moderate ^
    --out "!OUT!_panel" -v
if errorlevel 1 (
    echo   The panel failed. The analysis above is unaffected.
)

REM ---------------------------------------------------------------------
REM 5. One ranked list from docking and the rest of the analysis.
REM ---------------------------------------------------------------------
REM Docking hits against the gene's own protein enter the main ranking as
REM tier 4 and sit below every curated factor route however well the pocket
REM validated. This fuses the rankings instead of comparing their tiers: the
REM panel's surviving routes, and every docking result whose validation gate
REM passed (the product dock and each factor dock). A refused dock is listed
REM with its reason and not used.
echo.
echo === 5 of 5: one ranked list from docking and the rest ===
echo.
set "DOCKARGS="
for %%F in ("!DOCKDIR!\*.json") do set "DOCKARGS=!DOCKARGS! --docking "%%F""
set "COMBPANEL="
if exist "!OUT!_panel\panel.json" set "COMBPANEL=--panel "!OUT!_panel\panel.json""
"%PY%" -m generoutes.combined --results "!OUT!\results.json" ^
    !COMBPANEL! !DOCKARGS! --gene !GENE! --top 40 ^
    --out "!OUT!_combined"
if errorlevel 1 (
    echo   The combined list failed. The analysis above is unaffected.
)

REM ---------------------------------------------------------------------
REM Optional: the dose check, and the exploratory dock.
REM ---------------------------------------------------------------------
REM Both read what is already written. A failure in either leaves every file
REM above as it was.
if defined DOSECHK (
    echo.
    echo === Dose check ===
    set "CMAXARG="
    if exist data\cmax_table.tsv set "CMAXARG=--cmax-table data\cmax_table.tsv"
    if not defined CMAXARG echo   ^(no data\cmax_table.tsv - using blood levels from FDA labels and ChEMBL^)
    REM Reads the FDA label of each approved candidate for a human peak level and
    REM adds what it can stand behind to data\fda_label_cmax.tsv, which the dose
    REM check reads. Most labels give no usable number, and those are listed in
    REM data\fda_label_cmax_review.tsv. A failure here only means fewer values.
    "%PY%" scripts\build_fda_cmax.py --from-results "!OUT!\results.json"
    if errorlevel 1 echo   The FDA label table could not be built. Continuing without it.
    "%PY%" -m generoutes.exposure "!OUT!\results.json" !CMAXARG! !LIBRARY! ^
        --out "!OUT!_exposure"
    if errorlevel 1 echo   The dose check failed. The analysis above is unaffected.
)
if defined EXPLORE (
    echo.
    echo === Exploratory dock for !GENE! ===
    echo   Not validated. The template structure is chosen by searching the PDB.
    echo   Up to three distinct ligand sites are docked, each costing as much time as
    echo   the first. Add --max-sites 1 to the command below to dock only the best.
    echo   For a kinase the reference dock doubles the time and is what the class
    echo   check is calibrated on.
    "%PY%" -m generoutes.exploratory --gene !GENE! --uniprot !XUNI! ^
        !LIBRARY! --library-limit !LIBLIMIT! --exhaustiveness !EXH! ^
        --jobs !JOBS! --cpu !CORES! --reference ^
        --out "!OUT!_exploratory"
    if errorlevel 1 echo   The exploratory dock failed. The analysis above is unaffected.
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
if exist "!OUT!_combined\combined.md" echo Combined: !OUT!_combined\combined.md
if defined DOSECHK if exist "!OUT!_exposure" echo Dose:     !OUT!_exposure
if defined EXPLORE if exist "!OUT!_exploratory" echo Explore:  !OUT!_exploratory ^(not validated^)
echo.
echo Read the evidence table at the top of report.md before the ranking
echo underneath it. A source that failed is recorded as failed and never
echo as zero, and the verdict says whether a pillar of the run is missing.
pause
