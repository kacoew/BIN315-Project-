# ✍️ Skriveplan: Introduction + Methods

> Dette er en **plan**, ikke ferdig tekst. Dere skriver alle setningene selv, på engelsk.
> Fakta er hentet fra de to PDF-ene i PROSJEKT-mappen. Sidetall = side i PDF-en.

**Ordbudsjett:** Introduction ~300 · Methods ~300 (av ~1500 totalt)

---

# Introduction (~300 ord, 4 avsnitt)

```
 ① AD generelt  ─►  ② Braak  ─►  ③ Fusiform gyrus  ─►  ④ Gap + spørsmål
   (bredt)                                                 (smalt)
```

🔻 **Trakt:** start **bredt** (sykdommen) og slutt **smalt** (deres spørsmål).

### ① Alzheimer's disease (~70 ord)
**Svarer på:** hva er AD, og hvorfor er den viktig?

| Fakta å bruke | Kilde |
|---|---|
| Irreversibel nevrodegenerativ sykdom, fører til demens, mest hos eldre | Ribeiro-dos-Santos 2023, s. 1 (siterer Alzheimer's Association 2022) |
| Kjennetegn: **β-amyloid-plakk** og **tau-floker** | Ribeiro-dos-Santos 2023, s. 1 (siterer Jack et al. 2018) |
| Molekylære endringer starter **lenge før** symptomene | Ribeiro-dos-Santos 2023, s. 1 (siterer Aisen et al. 2017) |

### ② Braak stage (~70 ord)
**Svarer på:** hvordan måler man hvor langt sykdommen har kommet?

| Fakta å bruke | Kilde |
|---|---|
| Tau-flokene sprer seg i et **fast mønster**, I → VI | Braak & Braak 1991 (originalkilden, siteres i Srinivasan 2020) |
| Temporal cortex rammes **tidligere** enn frontal cortex | Srinivasan 2020, s. 5 |
| Frontal cortex rammes ved **Braak V–VI**, omtrent når demensen starter | Srinivasan 2020, s. 3 |

💡 **Poeng å få fram:** Braak er en **gradert** skala, så den gir mer informasjon enn bare AD/kontroll.

### ③ Fusiform gyrus (~70 ord)
**Svarer på:** hvorfor akkurat denne hjerneregionen?

| Fakta å bruke | Kilde |
|---|---|
| Viktig for **ansikts- og objektgjenkjenning** (og lesing) | Srinivasan 2020, s. 5 · Ribeiro-dos-Santos 2023, s. 2 |
| Endret **funksjonell konnektivitet** ved mild kognitiv svikt (MCI) | Ribeiro-dos-Santos 2023, s. 2 (siterer Ma et al. 2020) |
| **Atrofi** henger sammen med amyloid i hippocampus | Ribeiro-dos-Santos 2023, s. 2 (siterer Chang et al. 2016) |
| Mekanismene i fusiform gyrus er **lite utforsket** | Ribeiro-dos-Santos 2023, s. 2 |

### ④ Gap + forskningsspørsmål (~90 ord)
**Svarer på:** hva mangler, og hva gjør dere?

- **Gap:** tidligere studier på dette datasettet sammenlignet mest **AD mot kontroll**. Få har sett på **Braak som en gradert skala** i et co-expression-nettverk.
- **Metode i én setning:** WGCNA → moduler → korrelasjon med Braak.
- **Spørsmålet** (én setning, deres egne ord).
- *(Valgfritt)* én setning om **ekstra-analysene**.

> ⚠️ **Sjekk gapet selv:** les Results-delene i begge artiklene og se om de faktisk brukte Braak. Srinivasan 2020 brukte Braak for gensett-scorer (s. 5), men ikke WGCNA-moduler.

---

# Methods (~300 ord)

> 💡 Skriv **Data nå**. Resten skriver dere **mens dere gjør analysen**: fyll inn parametrene dere faktisk brukte.

| Underoverskrift | Hva skal med | Når | ~Ord |
|---|---|---|---|
| **Data** | GEO-ID, Genentech, post-mortem fusiform gyrus, Illumina HiSeq 2500, 289 prøver (219 AD / 70 kontroll), traits, NCBI raw counts (GRCh38.p13), 89 personer overlapper med GSE95587 | **nå** | 70 |
| **Preprocessing** | filterregel (≥ 10 counts i ≥ 20 %), antall gener igjen, VST (DESeq2, `blind = TRUE`) | steg 1 | 40 |
| **Overview** | topp 2000 variable gener, PCA, hierarkisk clustering (avstand + linkage), outlier-sjekk | steg 2 | 40 |
| **Network analysis** | antall gener (5000?), signed/unsigned, valgt power + hvorfor (R² ≥ 0.8), `minModuleSize`, `mergeCutHeight` | steg 3 | 60 |
| **Module–trait** | Spearman, Braak 1–6, FDR (BH), `lm` med alder + kjønn, kME, hub-gener | steg 4 | 40 |
| **Enrichment + pathway** | kME > 0.7, Enrichr (hvilke biblioteker + dato), pathview (hvilken pathway, hvilken verdi) | steg 5–6 | 30 |
| **Software** | R-versjon + pakker med versjoner (fra `sessionInfo()`) | til slutt | 20 |

**Regler for Methods:**
- **Fortid:** *"Genes were filtered …"*, ikke *"We will filter …"*
- **Nok detalj til at noen kan gjenta det**, men ikke forklar teorien (den hører ikke hjemme her)
- Henvis til pakkene: DESeq2, WGCNA, pathview, org.Hs.eg.db

---

# References (samle underveis)

| Kilde | Til |
|---|---|
| Srinivasan K et al. (2020) *Cell Reports*. PMID 32610143 | datasettet + Braak + fusiform |
| Ribeiro-dos-Santos A et al. (2023) *Front. Aging Neurosci.* 15:1138336 | fusiform + tidligere nettverksanalyse |
| Braak H & Braak E (1991) *Acta Neuropathologica* | Braak stage |
| Langfelder P & Horvath S (2008) *BMC Bioinformatics* | WGCNA |
| Love MI et al. (2014) *Genome Biology* | DESeq2 |
| Luo W & Brouwer C (2013) *Bioinformatics* | pathview |
| GEO GSE125583 | data |

💡 **Sekundærkilder:** står det *«(Ma et al. 2020)»* i Ribeiro-dos-Santos, så har **de** lest Ma, ikke dere. Enten finn og les originalen, eller skriv *«as reviewed in Ribeiro-dos-Santos et al. 2023»*.

---

# ✅ Sjekk når et avsnitt er ferdig
```
☐ Én hovedidé per avsnitt
☐ Første setning sier hva avsnittet handler om
☐ Hver påstand har en kilde
☐ Egne ord (ikke omskrevet setning for setning fra artikkelen)
☐ Innenfor ordbudsjettet
```

Skriv `sjekk teksten min` + lim inn avsnittet, så får dere tilbakemelding.
