# 🗺️ Stegplan: «Alle-labber»-prosjektet

> **Kodemal:** `analysis_alle_labber.Rmd` (samme mappe). Hvert steg her har en egen seksjon der.
> **Teori:** `oppgavemal.md` · **Skriving:** `skriveplan.md`
> Den gamle `analysis_mal.Rmd` er ikke endret, så dere kan fortsatt bruke den.

---

## ❓ Forskningsspørsmål

**Hovedspørsmål:**
> *Are the gene co-expression modules associated with Braak stage in the fusiform gyrus the same as those associated with Alzheimer's disease diagnosis, and what functions do they have?*

| | Delspørsmål | Labber |
|---|---|---|
| **RQ1** | Which **genes and modules** are associated with Braak stage and with diagnosis? | 2, 3, 4, 6 |
| **RQ2** | Which **biological processes, pathways and regulators** characterise these modules? | 6, 7, 8 |
| **RQ3** | Do the modules contain information about AD **beyond age and sex**? | 5 |

---

## 🔗 Oversikt

```
 Steg 0–1  Lab 1–2  data → filtrering → VST
                │
 Steg 2    Lab 4    PCA + clustering + heatmap (kvalitet)
                │
        ┌───────┴────────┐
 Steg 3 Lab 3 DESeq2   Steg 4–5 Lab 6 WGCNA → module–trait
        │ log2FC             │
        │                    ├─► Steg 6  Lab 6  hubs + AD-risikogener
        │                    ├─► Steg 7  Lab 7  Enrichr GO/KEGG
        └────────────────────┼─► Steg 8  Lab 7  pathview (farget med log2FC)
                             ├─► Steg 9  Lab 8  TF-regulatorer
                             └─► Steg 10 Lab 5  random forest (bonus)
```

---

## 🚦 Prioritet: kutt fra bunnen hvis tiden blir knapp

| | Steg | Hvis dere må kutte |
|---|---|---|
| 🔴 **Må** | 0, 1, 2, 4, 5, 7 | kan ikke kuttes, dette **er** hovedspørsmålet |
| 🟠 **Bør** | 6 (risikogener) + 9 (TF) · **lab 8: omics integration** 🌟 | én tabell med topp-TF og risikogen-anrikning |
| 🟡 **Fint** | 3 (DESeq2), 8 (pathview) | én figur eller én setning hver |
| 🟢 **Bonus** | 10 (ML) | kutt først |

💡 **Hvorfor er lab 8 høyt?** Lærerens slides kaller lab 8 *«ambitious projects!»*. GWAS-gener (**genetikk**) og ChIP-seq-mål (**regulering**) er to ekstra datatyper oppå genuttrykket. Det er **integrasjon**, og det gir mer enn ML.

---

## 📋 Stegene

> **Hvem:** **K** = Kachina · **C** = Charlene (bare et forslag)
> **Ord** = omtrent hvor mye det får i rapporten

### Steg 0 — Last inn data · Lab 1 · **K** · ons 18.11
- **Mål:** `counts` + `traits` med Braak som tall (1–6), diagnose som 0/1, kjønn som 0/1
- **Pass på:** 2 prøver mangler Braak, 1 mangler alder og 1 mangler kjønn
- ✅ **Ferdig når:** `table(traits$braak, traits$braak_num)` viser riktig koding

### Steg 1 — Filtrering + VST · Lab 2 · **K** · ons 18.11
- **Mål:** fjerne lavt uttrykte gener og gjøre variansen stabil
- **Regel:** ≥ 10 counts i ≥ 20 % av prøvene
- ✅ **Ferdig når:** `vsd_mat` har ca. 21 500 gener × 289 prøver
- 📝 **Methods:** filterregel + antall gener igjen (~40 ord)

### Steg 2 — PCA, clustering, heatmap · Lab 4 · **K** · tor 19.11
- **Mål:** se hovedmønstrene og finne outliers
- ✅ **Ferdig når:** PCA farget etter Braak, diagnose og kjønn + prøvetre uten tydelige outliers
- 📝 **Results:** 1 figur (PCA) + 1–2 setninger

### Steg 3 — DESeq2 · Lab 3 · **C** · fre 20.11
- **Mål:** finne **enkeltgener** som endrer seg med Braak
- **Modell:** `~ age + sex + braak_num` → koeffisienten for Braak = log2FC **per Braak-stadium**
- ⏱️ `DESeq()` på 289 prøver tar **5–15 min**. Bruk `cache = TRUE`
- ✅ **Ferdig når:** antall DEGs (padj < 0.05) + volcano plot
- 📝 **Results:** 1 setning + volcano (eller tabell over topp 10)

### Steg 4 — WGCNA: nettverk og moduler · Lab 6 · **K** · lør 21. – søn 22.11
- **Mål:** moduler av gener som varierer sammen
- **Valg å skrive ned:** antall gener (5000), signed, power (R² ≥ 0.8), `minModuleSize`, `mergeCutHeight`
- ✅ **Ferdig når:** dendrogram med modulfarger + tabell med antall gener per modul
- 🚩 **Sjekkpunkt:** har dere ikke moduler innen **søndag 22.11**, be om hjelp med en gang

### Steg 5 — Module ↔ trait · Lab 6 · **K** · man 23.11
- **Mål:** svare på **RQ1** og **hovedspørsmålet**
- **Delsteg:**
  1. Heatmap: moduler × (Braak, diagnose, alder, kjønn), Spearman + FDR
  2. 🌟 **Hovedfigur:** scatter, én prikk per modul, x = cor(Braak), y = cor(diagnose)
  3. `lm(ME ~ braak_num + ad + age + male)`: sier diagnose noe **ekstra**?
- ✅ **Ferdig når:** dere kan si hvilke moduler som følger Braak, diagnose eller begge
- 📝 **Results:** 2 figurer + 1 avsnitt (det viktigste i rapporten)

### Steg 6 — Hub-gener, DEG-sjekk og AD-risikogener · Lab 3 + 6 · **K** · man 23.11
- **Mål:** hvem er de viktigste genene, og stemmer WGCNA med DESeq2?
- **Delsteg:**
  1. Topp 10 kME i hver interessant modul → gensymboler
  2. Hvor mange av hub-genene er DEGs (fra steg 3)?
  3. Er kjente **AD-risikogener** (GWAS Catalog) overrepresentert i modulene? `phyper`, som i lab 6
  4. Ligger *FNDC3A, MED23, NRIP1, PKN2* (Ribeiro-dos-Santos 2023) i de samme modulene?
- 📝 **Results:** 1 tabell (hubs) + 1–2 setninger

### Steg 7 — GO/KEGG med Enrichr · Lab 7 · **C** · tir 24.11
- **Mål:** svare på **RQ2**: hva gjør modulene?
- **Slik:** gener med kME > 0.7 → tekstfil → https://maayanlab.cloud/Enrichr/
- **Biblioteker:** *GO Biological Process 2023* + *KEGG 2021 Human*. **Skriv ned datoen!**
- ✅ **Ferdig når:** topp 5 termer per modul i en tabell eller et søylediagram

### Steg 8 — Pathview · Lab 7 · **C** · ons 25.11
- **Mål:** vise **hvor** i en pathway endringene skjer
- **Farge:** log2FC per Braak-stadium fra DESeq2 (steg 3)
- **Pathway:** velg en fra steg 7 (f.eks. synaptisk vesikkelsyklus `04721` eller Alzheimer `05010`)
- 📝 **Results:** 1 figur

### Steg 9 — Transkripsjonsfaktorer · Lab 8 · **C** · tor 26.11 · 🟠 integrasjon
- **Mål:** har modulene **felles regulatorer**?
- **Slik:** samme hypergeometriske test som lab 8, men med menneske-data fra ChIP-seq (Enrichr *ChEA 2022*)
- **Koble til steg 6:** én felles tabell per modul med **AD-risikogener (GWAS)** + **topp-TF (ChIP-seq)**
- ✅ **Ferdig når:** dere kan si om Braak-modulene er styrt av kjente regulatorer og/eller har flere risikogener enn forventet
- 📝 **Results:** 1 tabell + 2 setninger · **Discussion:** 1–2 setninger om hva integrasjonen viser

### Steg 10 — Random forest · Lab 5 · **K** · fre 27.11 · 🟢 bonus
- **Mål:** svare på **RQ3**
- **Slik:** modell A = alder + kjønn, modell B = alder + kjønn + alle MEs. Samme 5-fold CV for begge
- **Mål for hvor godt:** ROC-AUC + **balanced accuracy** (fordi gruppene er skjeve: 219 AD og 70 kontroll)
- ✅ **Ferdig når:** to ROC-kurver i ett plott + én tabell
- ⚠️ **Gjør bare dette** hvis steg 0–8 er ferdige

### Skriving · **K + C** · lør 28. – man 30.11
- **K:** Results · **C:** Discussion + Introduction
- **Sammen:** Abstract, knit, ordtelling, korrektur
- 📤 **Lever tirsdag 1.12** · 🎤 Øv presentasjonen ons 2. – tor 3.12

---

## 📏 Hvor dypt skal hver analyse være?

| Steg | Dybde | Hva holder |
|---|---|---|
| 1–2 | 🟢 grunt | standard, én figur |
| 3 | 🟡 middels | én modell, ikke prøv mange design |
| 4–5 | 🔴 **dypt** | dette er hovedanalysen: begrunn valgene og sjekk robusthet |
| 6–7 | 🟡 middels | bare de 2–4 mest interessante modulene |
| 8 | 🟢 grunt | én pathway |
| 6.5 + 9 (lab 8) | 🟡 middels | én integrert tabell: GWAS-gener + TF per modul |
| 10 | 🟢 grunt | to modeller, ingen tuning-jakt |

---

## 📊 Figurbudsjett (maks ca. 6 i rapporten)

| # | Figur | Steg |
|---|---|---|
| 1 | PCA farget etter Braak | 2 |
| 2 | Module–trait heatmap | 5 |
| 3 | 🌟 cor(Braak) mot cor(diagnose) | 5 |
| 4 | GO-søylediagram | 7 |
| 5 | 🌟 Integrert tabell: hubs + GWAS-gener + topp-TF per modul | 6 + 9 |
| 6 | Pathview **eller** volcano **eller** ROC | 8 / 3 / 10 |

Resten (dendrogram, soft threshold, prøvetre) kan ligge i et vedlegg eller i presentasjonen.
