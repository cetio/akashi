# Testing

Akashi has one offline unittest plus version-gated live integration tests.

## Prerequisites

Install DUB and a D compiler. The commands here use LDC (`ldc2`). Akashi has no external dependencies.

## Library Build

```sh
dub build --compiler=ldc2 --config=library
```

## Offline Tests

```sh
dub test --compiler=ldc2 --config=library
```

This runs the inline unittests (URL building in `source/akashi/request.d`) and makes no network requests. Use it for routine verification.

## Live Integration Tests

Files under `source/tests/integration/` are guarded by D version flags. Each flag enables real HTTP requests to the matching public service. No credentials are needed, but the tests require network access and depend on each service's availability and rate limits.

| Flag | Service |
| --- | --- |
| `pubchem` | PubChem |
| `wikipedia` | Wikipedia |
| `psychonaut` | PsychonautWiki |
| `pubmed` | PubMed |
| `pmc` | PubMed Central |

```sh
dub test --compiler=ldc2 --d-version=pubchem
```

Combine flags to run several suites at once, for example `--d-version=pubmed --d-version=pmc`. A failure may reflect service state rather than a library defect.
