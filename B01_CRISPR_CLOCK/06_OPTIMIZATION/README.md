# Optimization layer

Implementation: `b01_clock/solver/optimize.py`.

Multi-start L-BFGS-B + differential evolution over internal ages, with rates profiled globally in locked bounds. Returns timing certificates with solver diagnostics.

Use exact grid solver for ≤3 free internals; use this layer for larger trees. Topology-ensemble unions live in `08_TOPOLOGY/`.
