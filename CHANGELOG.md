# Changelog

This changelog starts at 2.0.

It records what changes **for you as a user** — the commands and functions you call, the
structures that come out, and the things you have to do differently. Internal restructuring is
not listed, however large, unless you can see it from outside.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — 2.0.0

DODO 2.0 is a rewrite. **It breaks the 1.x API deliberately**, so read the migration notes below
before upgrading a script.

### Fixed

- Rebuilt-region/folded-domain peptide seams are constrained during CA generation and closed with
  exact bond geometry. The previous post-processing fallback left visible 2.6–4.3 Å C–N gaps;
  paired testing over dnmt3a, arf19 and p300 at three seeds reduced 60 strained seams to zero.
