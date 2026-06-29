<!--
SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC0-1.0
-->

# FPGAlign [![build status][1]][2]
<!-- [![codecov][3]][4] -->
<!--
    Above uses reference-style links with numbers.
    See also https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#links.

    For example, `[![build status][1]][2]` evaluates to the following:
        `[link_text][2]`
        `[2]` is a reference to a link, i.e. `[link_text](https://...)`

        `[link_text]` = `[![build status][1]]`
        `[1]` is once again a reference to a link - this time an image, i.e. `[![build status](https://...)]
        `![build status]` is the text that should be displayed if the linked resource (`[1]`) is not available

    `[![build status][1]][2]` hence means:
    Show the picture linked under `[1]`. In case it cannot be displayed, show the text "build status" instead.
    The picture, or alternative text, should link to `[2]`.
-->

<!--
    This is the CI badge image:
        `https://img.shields.io/github/workflow/status/` - we do not use GitHub's badges as they are not customisable.
        `/seqan/FPGAlign/` - owner/repository
        `CI%20on%20Linux` - name of the workflow as encoded URL (e.g., whitespace = %20)
        `main` - branch to show
        `?style=flat&logo=github` - use a GitHub-style badge
        `&label=FPGAlign%20CI` - text on the badge
        `"Open GitHub actions page"` - this text will be shown on hover
-->
[1]: https://img.shields.io/github/actions/workflow/status/seqan/FPGAlign/ci_linux.yml?branch=main&style=flat&logo=github&label=FPGAlign%20CI "Open GitHub actions page"
<!--
    This is the CI badge link:
        `https://github.com/seqan/FPGAlign/actions` - actions page of owner/repository
        `?query=branch%3Amain` - only show actions that ran on the mater branch
-->
[2]: https://github.com/seqan/FPGAlign/actions?query=branch%3Amain
<!--
    This is the Codecov badge image:
        Codecov offers badges: https://app.codecov.io/gh/seqan/FPGAlign/settings/badge
        While being logged in into Codecov, navigate to Settings->Badge and copy the markdown badge.
        Copy the image part of the markdown badge here.
    `"Open Codecov page"` - this text will be shown on hover
-->
[3]: https://codecov.io/gh/seqan/FPGAlign/branch/main/graph/badge.svg "Open Codecov page"
<!--
    This is the Codecov badge link:
        Codecov offers badges: https://app.codecov.io/gh/seqan/FPGAlign/settings/badge
        While being logged in into Codecov, navigate to Settings->Badge and copy the markdown badge.
        Copy the URL part of the markdown badge here.
-->
[4]: https://codecov.io/gh/seqan/FPGAlign

This repository implements a CPU-based read-mapper that combines an Interleaved Bloom Filter (IBF) and per-bin FM-indexes
to quickly map sequencing reads against a collection of reference genomes organized in user-defined bins.

The implementation is split into two main subcommands: `build` and `search`.

## Background

- The pipeline uses an IBF as a fast probabilistic prefilter to assign reads to candidate bins.
- Candidate bins are searched with per-bin FM-indexes to determine exact match positions.
- Final alignments are computed with SeqAn3 pairwise alignment and exported as a SAM file.
- To avoid I/O and synchronization bottlenecks the pipeline uses a specialized shopping-cart queue (SCQ) to batch and
    asynchronously pass results between pipeline stages.

## Features

- `build`: construct an IBF and per-bin FM-indexes from a file of per-bin reference paths.
- `search`: query the constructed index with FASTA/FASTQ queries and produce a SAM output.
- Parallelized IBF membership, FM-index searching and alignment stages with configurable thread counts.

## Quickstart

### Build from source (example)

```bash
mkdir build
cd build
cmake ..
make -j
```

### Example: construct indexes from a list of bins

```bash
./bin/FPGAlign build \
    --input /path/to/bin_list.txt \
    --output /path/to/output_prefix \
    --kmer 20 --window 20 --hash 2 --fpr 0.05 --threads 4
```

- `bin_list.txt` is a whitespace-separated list where each non-empty line defines one user bin and contains one or more file paths (reference sequences) belonging to that bin.

### Example: search with query reads

```bash
./bin/FPGAlign search \
    --input /path/to/output_prefix \
    --query queries.fasta \
    --output results.sam \
    --threads 4 --errors 2 --queue-capacity 8
```

- `--input` points to the prefix used when building the index (the tool will look for `*.ibf`, `*.meta`, `*.fmindex`, `*.ref` files using that prefix).

## Input format

- Build input: a text file where each non-empty line contains one or more reference file paths (FASTA/FASTQ) that constitute a single bin.
- Search input: query sequences in FASTA/FASTQ format (the code uses SeqAn3 sequence file input for parsing).

## Output files produced by `build`

- `<output_prefix>.ibf` — serialized Interleaved Bloom Filter.
- `<output_prefix>.meta` — metadata describing k-mer/window sizes, number of bins and reference IDs.
- `<output_prefix>.<id>.fmindex` — per-bin serialized FM-index (one file per bin id).
- `<output_prefix>.<id>.ref` — per-bin serialized reference sequences (stored for alignment retrieval).

## Runtime behavior

- The `search` pipeline consists of three asynchronous stages connected by SCQs:
    1. IBF membership agent (prefilter) — produces candidate bin hits.
    2. FM-index lookup per bin — locates exact reference positions for candidate reads.
    3. Pairwise alignment — computes final CIGARs and writes SAM records.

## Key options

- `--kmer`, `--window`: k-mer and window sizes used for minimiser hashing (default `k=20`).
- `--hash`, `--fpr`: number of hash functions and target false-positive rate for the IBF.
- `--errors`: maximum allowed errors for FM-index search.
- `--threads`: number of threads for parallel stages.
- `--queue-capacity`: batching capacity of the shopping-cart queues (SCQ).

## Development & testing

- Build and run tests (from repository root):

```bash
mkdir -p build && cd build
cmake ..
make -j
ctest -j --output-on-failure
```

<br>
<p align="center">
  <a target="_blank" rel="noopener noreferrer" href="#">
  <img src=".github/BMFTR_logo_foerderung.svg" alt="Funded by the BMFTR" width="280">
</p>
