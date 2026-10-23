# Bekende afwijkingen code vs README

Elk punt vraagt een beslissing over welk gedrag de bedoeling is. Tot die beslissing volgt de baseline
het gedrag van de code.

## `-S` (strand bias): standaardwaarde
- Code: 0.9
- Usage-tekst en README: 0.1
- Status: te beslissen door Nick

## UNMET-drempel
- Code: varianten met een UNMET-score < 0.90 blijven over
- README: filtert scores > 0.93
- Status: te beslissen door Nick

## `-O` (outputmap)
- Code: optioneel, standaard `./output` (relatief aan de hardcoded werkmap)
- README: vermeld als verplicht
- Status: te beslissen door Nick

## Python-dependencies
- Code: `post_filtering.py` gebruikt `openpyxl`
- README: installeert `xlrd` (niet gebruikt) en `glob` (zit in de standaardbibliotheek; `pip install` faalt)
- Status: te beslissen door Nick

## `-T` (matched control)
- Code: voegt alleen read-tellingen uit de control toe aan de varianten; `matched_filtering.py` wordt niet aangeroepen
- README: beschrijft filtering met `matched_filtering.py` (Fisher-test, VAF-verschil, VAF-ratio)
- Status: te beslissen door Nick

## PoN-percentage voor varianten uit één tool
- Code: 20% van de normalen
- Commentaar in `CircuSNV.sh`: 50%
- Status: te beslissen door Nick

## Mappaden met of zonder afsluitende slash
- Code: sommige functies bouwen paden als `$2${base_name}` (verwacht een afsluitende slash), andere als `$2/${base_name}`
- README: niet beschreven
- Status: te beslissen door Nick
