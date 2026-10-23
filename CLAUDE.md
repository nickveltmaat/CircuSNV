# CircuSNV

## Doel
Een betere, efficiëntere en gebruiksvriendelijkere versie van CircuSNV maken, met exact dezelfde
functionaliteit, in bash.

## Structuur
- `CircuSNV.sh` is het hoofdscript (bash, SLURM). Volgorde: argumenten parsen → optioneel Panel of
  Normals genereren → per sample 4 callers parallel (VarDict, Mutect2, LoFreq, SiNVICT) → samenvoegen
  met `bcftools isec` → PoN- en UNMET-filtering → optioneel MNV-calling (MNVista) → annotatie met
  OpenCRAVAT → post-filtering → tweede annotatie.
- Python-hulpscripts in de root, aangeroepen via de venv `./env`:
  - `sinvict_to_vcf2.py`: SiNVICT-uitvoer naar VCF
  - `isec_to_vcf.py`: `bcftools isec`-sites naar VCF
  - `finish_pon.py`: PoN aanvullen met varianten uit één tool
  - `merge_variants.py`: metrics per tool samenvoegen; modi `no_pon`, `pon`, `phased`
  - `post_filtering.py`: filteren op functioneel effect en populatie-AF
- `matched_filtering.py` staat in de repo maar wordt niet aangeroepen.

## Draaien
- Alleen op het HPC-cluster: vereist SLURM en `module load` (versies staan in `CircuSNV.sh`).
  Kan niet lokaal draaien. Er is geen testsuite.
- Het script doet `cd` naar een hardcoded clusterpad; relatieve paden (`./env`, `./MNVista`, `./PoN`,
  `./*.py`, `-O`) gelden ten opzichte van die map.
- Niet in de repo: `./env`, `./MNVista`, `./MNV_PoN_No_overlap_new.txt`, de vooraf gebouwde tools
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
