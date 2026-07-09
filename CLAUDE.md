# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Pono is a SMT-based model checker for verifying properties of transition systems. It supports hardware verification via BTOR2, SMV, and VMT input formats, and implements algorithms including IC3, BMC, k-induction, and interpolant-based model checking. Written in C++17, it uses smt-switch for solver-agnostic SMT operations.

## Build Commands

```bash
# First-time setup (order matters)
./contrib/setup-smt-switch.sh     # Build smt-switch with Bitwuzla
./contrib/setup-btor2tools.sh     # Build btor2tools
./configure.sh                    # Generate build dir (add --debug for debug build)

# Build
ninja -C build

# Run all tests
ninja -C build check

# Run specific test
ctest --test-dir build -R test_ts               # Match by name pattern
build/test_ts --gtest_filter="*specific_test*"  # GTest filter

# Python bindings
./contrib/setup-smt-switch.sh --python
./configure.sh --python
ninja -C build
pip install 'build/python[test]'
pytest tests
```

Key configure options: `--debug`, `--python`, `--with-msat`, `--with-yices2`, `--with-z3`, `--with-btor`, `--docs`, `--static`.

## Architecture

**Data flow:** Parse input → TransitionSystem → Apply modifiers → Choose prover engine → Unroll + send to SMT solver → Produce witness or proof.

### Core (`/core/`)
- **TransitionSystem** (`ts.h`): Base class with `init_` (initial states) and `trans_` (transition relation). Two subclasses:
  - **RelationalTransitionSystem** (`rts.h`): Explicit transition relations
  - **FunctionalTransitionSystem** (`fts.h`): Functional state update assignments via `assign_next()`
- **Unroller** (`unroller.h`): Creates time-step formulas for bounded analysis
- **Prop** (`prop.h`): Property representation

### Engines (`/engines/`)
Prover hierarchy rooted at `Prover` → `SafetyProver` / `LivenessProver`:
- **Bmc**: Bounded model checking (linear and binary search)
- **KInduction**: K-induction proofs
- **IC3Base**: Abstract parameterized IC3; concrete variants: IC3, IC3Bits, IC3SA, IC3IA, MbIC3
- **InterpolantMC / InterpSeqMC**: Interpolant-based model checking
- **SyGusPDR**: Synthesized predicates for IC3

### Frontends (`/frontends/`)
Input parsers: BTOR2Encoder, SmvEncoder (Flex/Bison-generated), VmtEncoder, CoreIREncoder (optional).

### Modifiers (`/modifiers/`)
TransitionSystem transformations: ArrayAbstractor (arrays→UFs), OpsAbstractor, StaticConeOfInfluence, ControlSignals (clock/reset), HistoryModifier, LivenessToSafetyTranslator.

### Refiners (`/refiners/`)
CEGAR refinement: ArrayAxiomEnumerator for abstraction-refinement loops.

### Main entry point
`pono.cpp` — CLI that parses options, reads input, applies modifiers, runs prover, outputs witness (BTOR2 or VCD format).

## Code Style

- Clang-format: Google-based style, 2-space indent, pointer alignment middle, no bin-packing args
- Format check runs in CI; use `clang-format` before committing
- File headers use Doxygen-style `/*! \file ... \brief ... */` blocks
- Tests use GoogleTest framework (`tests/test_*.cpp`)

## Gotchas

- **Solver-term binding:** Nodes are tied to a solver. Use `solver.transfer_term(t)` when moving terms between solvers.
- **Boolector rewriting:** Boolector treats bools and BV1 equivalently and does on-the-fly rewriting. Asserting `x = 0` can rewrite `x = 1` elsewhere to `false`, breaking untimed/timed formula expectations.
- **Boolector bool/BV1:** Transferring terms from Boolector to other solvers can fail due to bool↔BV1 conflation; handle with lazy casting.
