# This is a collection of tricky cases in pono

## Solver Issues

* Nodes are tied to a solver. Thus, if you create a new solver, you must transfer all terms that you would like to use to the new solver via `solver.transfer_term(t)`.
* Boolector treats bools and bitvectors of size one equivalently, therefore you can easily run into problems tranferring terms __from__ boolector to a different solver. This can be handled by casting lazily.
* Boolector does on the fly rewriting as you assert things. Therefore, you can run into strange cases where nodes are not what you expect. For example, say that the initial state is `x = 0`. Somewhere else in your formula you have `x = 1`. If you assert the 'untimed' initial state, you might expect that it should not affect your final result, because the rest of the assertions are timed. However, once you assert `x = 0`, then `x=1` is rewritten to `false`. Therefore, when you unroll `x = 1`, expecting to get `x@0 = 1`, you instead get `false`
