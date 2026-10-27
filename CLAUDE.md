# CircuSNV

## Doel
Een betere, efficiëntere en gebruiksvriendelijkere versie van CircuSNV maken, met exact dezelfde
functionaliteit, in bash.

## Structuur
- `CircuSNV.sh` is het startscript (bash, met `#SBATCH`-kop): instellingen laden, argumenten
  valideren, controles vooraf, daarna PoN en samples. De code staat in `lib/`, één module per stap:
  `config.sh` (modules, resourcepaden, vaste waarden), `args.sh` (opties, help, validatie),
  `common.sh` (logging, modules, venv), `checks.sh`, `callers.sh` (VarDict, Mutect2, LoFreq,
  SiNVICT, retry-wrapper), `filters.sh` (UNMET), `annotate.sh` (OpenCRAVAT, matched control),
  `mnv.sh` (MNVista), `pon.sh`, `samples.sh` (per sample, volgorde, resume).
- Volgorde per sample: 4 callers parallel → samenvoegen met `bcftools isec` → PoN- en
  UNMET-filtering → optioneel MNV-calling → annotatie met OpenCRAVAT → post-filtering → tweede
  annotatie.
- Python-hulpscripts in `scripts/`, aangeroepen via de venv `env/` van de installatie:
  - `sinvict_to_vcf2.py`: SiNVICT-uitvoer naar VCF
  - `isec_to_vcf.py`: `bcftools isec`-sites naar VCF
  - `finish_pon.py`: PoN aanvullen met varianten uit één tool
  - `merge_variants.py`: metrics per tool samenvoegen; modi `no_pon`, `pon`, `phased`
  - `post_filtering.py`: filteren op functioneel effect en populatie-AF
- `scripts/matched_filtering.py` staat in de repo maar wordt niet aangeroepen.
- De twee annotatiestappen zijn bewust: de eerste haalt de gegevens op waarop gefilterd wordt, maar
  het filteren verliest de structuur en layout van de annotatie, dus wordt daarna opnieuw geannoteerd.
  Dit kan later misschien beter; verwijder of combineer deze stappen niet zonder overleg met Nick.

## Draaien
- Alleen op het HPC-cluster: vereist SLURM en `module load` (namen en versies in `lib/config.sh`,
  aan te passen met `--config`). Kan niet lokaal draaien. Er is geen lokale testsuite;
  vergelijkende tests tegen een referentierun draait Nick op Nibbler volgens `docs/nibbler-test.md`
  (vergelijkingsscript in `tests/nibbler/`).
- Het script zoekt zijn installatiemap zelf (ook onder `sbatch`); `env/`, `MNVista`,
  `MNV_PoN_No_overlap_new.txt`, de UNMET-bed en de PoN-werkmap `PoN/` staan daar. Relatieve
  invoerpaden en `-O` gelden vanaf de map waar de pipeline gestart wordt.
- Niet in de repo: `env/`, `MNVista`, `MNV_PoN_No_overlap_new.txt`, de vooraf gebouwde tools
  (LoFreq, SiNVICT, bam-readcount).
- Genoombuild hg38 is hardcoded.
- Bij verschillen tussen README en code geldt voorlopig het gedrag van de code als baseline.
  Openstaande punten staan in `docs/known-issues.md`.

## Werkafspraken
- Geen gedragswijziging zonder vergelijking met de baseline-output. Refactors en gedragswijzigingen
  in aparte commits.
- Bij grotere wijzigingen eerst een plan voorleggen.
- Nooit op `main` werken: feature-branch, kleine commits, niet pushen zonder te vragen.
- Pipeline niet lokaal draaien. Lokaal alleen syntaxis en lint (`bash -n`, `shellcheck` indien
  beschikbaar). Echte runs doet Nick op Nibbler via `sbatch`.
- Geen patiëntdata of echte sequencingdata lezen, maken of committen; alleen kleine testdata.
- Bash-stijl: `set -euo pipefail`, variabelen quoten, geen hardcoded paden, LF-regeleinden.
- Blijf bij bash, tenzij een andere taal duidelijk beter is; leg dat eerst voor.
- Uitleg aan Nick in het Nederlands; code, commentaar en commitberichten in het Engels (openbare repo).

## Conventies
- Scripts met LF-regeleinden (`.gitattributes`).
- Nooit data in git (BAM/VCF/FASTQ, `output/`, `results/`); de repo is openbaar, dus nooit patiëntdata.
