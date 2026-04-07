# Akashi

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

## Install

### DUB

Add Akashi to your `dub.json`:

```json
{
  "dependencies": {
    "akashi": "*"
  }
}
```

For local development as a path dependency:

```json
{
  "dependencies": {
    "akashi": {
      "path": "../akashi"
    }
  }
}
```

## Usage

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

Page page     = Page.fromRaw("title", "wikipedia", someWikitextString);
auto sections = page.sections;
writeln(page.preamble);
writeln(page.fulltext);
```

## Modules

- `akashi.pubchem`
  - Look up compounds, assays, proteins, synonyms, descriptions, 3D conformers, and similarity matches from PubChem.
- `akashi.page`
  - Work with a shared `Page` type that lazily fetches content and exposes `raw()`, `document()`, `sections()`, `preamble()`, and `fulltext()`.
- `akashi.wikipedia`
  - Search Wikipedia, fetch page content, and resolve likely compound pages.
- `akashi.psychonaut`
  - Search PsychonautWiki, resolve compound pages, extract dosage information, and find experience reports.
- `akashi.entrez`
  - Query NCBI Entrez, retrieve PubMed abstracts, retrieve PMC full text, and work with article metadata such as DOI values.
- `akashi.text`
  - Parse wikitext and XML into a compact AST for further analysis, rendering, or plain-text extraction.
- `akashi.composer`
  - Reuse Akashi's HTTP and rate-limiting orchestration primitives when building additional source adapters.

## License

Akashi is licensed under the [AGPL-3.0 license](LICENSE.txt).
