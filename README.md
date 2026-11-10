# Akashi

[![DUB](https://img.shields.io/badge/DUB-package-red)](https://code.dlang.org/packages/akashi) [![License](https://img.shields.io/badge/License-Apache--2.0-blue)](LICENSE.txt) [![CI](https://github.com/cetio/akashi/actions/workflows/ci.yml/badge.svg)](https://github.com/cetio/akashi/actions/workflows/ci.yml)

Akashi is a D library for chemistry-aware retrieval, source resolution, and text parsing across public scientific and knowledge sources.
It was designed to work with Intuit, a D AI library for easy, deterministic connectors to AI models, but it is fully functional on its own.

## Features

- PubChem Compounds
- PubChem Assays
- PubChem Proteins
- PubChem 3D Conformers
- Wikipedia
- PsychonautWiki
- Erowid (via PsychonautWiki experience reports)
- PubMed
- PubMed Central (PMC)

## Installation

```sh
dub add akashi
```

## Quick Start

### Compounds

```d
import akashi.pubchem : getProperties, similaritySearch;

auto c = getProperties("caffeine");
writeln(c.name, c.properties.formula, c.properties.weight);
writeln(c.synonyms());
writeln(c.description());

auto similar = similaritySearch(c.cid, 90, 10);
```

### Assays and proteins

```d
import akashi.pubchem : getProperties, getAssay, getProtein;

auto c = getProperties("aspirin");
auto results = c.assayResults;  // AssayResult[] with outcome, activity value, target
auto assays  = c.assays;        // unique Assay[] derived from results

auto assay = assays[0];
writeln(assay.name, assay.description(), assay.targets());
auto proteins = assay.proteins;

auto protein = getProtein("P23219");
writeln(protein.name, protein.taxonomy, protein.geneSymbol);
auto proteinAssays = protein.assays;
```

### 3D conformers

```d
import akashi.pubchem : getProperties;

auto conformer = getProperties("caffeine").conformer3D;
foreach (atom; conformer.atoms)
    writeln(atom.element, atom.x, atom.y, atom.z);
foreach (bond; conformer.bonds)
    writeln(bond.atomA, bond.atomB, bond.order);
```

### Wikipedia

```d
import akashi.wikipedia : getPages, resolvePage;
import akashi.pubchem : getProperties;

auto pages = getPages!"wikipedia"("serotonin", 5);
auto page  = resolvePage(getProperties("serotonin"));
writeln(page.preamble);
writeln(page.fulltext);
```

### PsychonautWiki

```d
import akashi.psychonaut : getDosage, resolvePage, getPages;
import akashi.pubchem : getProperties;

auto compound = getProperties("psilocybin");
auto page     = resolvePage(compound);

auto dosage = getDosage(compound);
foreach (d; dosage.dosages)
    writeln(d.route, d.common, d.strong);
```

### Erowid experience reports

```d
import akashi.psychonaut : resolvePage, getReports;
import akashi.pubchem : getProperties;

auto page    = resolvePage(getProperties("mdma"));
auto reports = getReports(page, 20);
foreach (r; reports)
    writeln(r.title, r.url);
```

### PubMed and PMC

```d
import akashi.entrez.pubmed : getPages;
import akashi.entrez.pmc   : getPages, getDOI;

auto abstracts = getPages!"pubmed"("serotonin transporter", 10);
writeln(abstracts[0].fulltext);

auto articles = getPages!"pmc"("CRISPR off-target", 5);
writeln(getDOI(articles[0]));
```

### Text parsing

```d
import akashi.page : Page;

Page page = Page.fromRaw("title", "wikipedia", someWikitextString);
auto sections = page.document.sections;
writeln(page.preamble);
writeln(page.fulltext);
```

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) and [TESTING.md](TESTING.md).

## License

Akashi is licensed under [Apache 2.0](LICENSE.txt).
