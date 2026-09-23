# Runahead decoding

Background notes, written for an agent picking this up cold. This describes an
**already-evaluated design** — it is context to reason *within*, not a proposal
to critique. Source: design discussion with the author (2026-09-23).

Current open work is designing decoders that fit this paradigm; see
`src/decoder/bunchaluts/DESIGN.md` for one candidate.

## The core idea

A real-time decoder is normally asked to be both fast (to meet the reaction-time
requirement) and accurate (to hit a target logical error rate). These two
demands fight each other in a single decoder.

Runahead decoding **decouples them across two decoders**:

| | latency | accuracy | what it determines |
|---|---|---|---|
| **L1** | must be fast | may be poor | the system's **reaction time** |
| **L2** | may be slow | must be excellent | the system's **logical error rate** |

L1 runs at the head and drives program execution — the quantum computer acts on
L1's Pauli frame updates, so reaction time equals L1's latency. L2 lags behind
and **verifies** whether L1 decoded the dependent detection events correctly.
Because the committed answer is always L2's, system LER tracks L2's LER, not
L1's.

The name is by analogy to runahead execution (Mutlu et al., HPCA 2003): a
speculative engine keeps the pipeline moving while the architecturally-correct
state resolves behind it.

## Misprediction recovery: the retirement FIFO

Instructions that have been executed by the quantum computer under L1's Pauli
frames are held in a **retirement FIFO** (the reorder buffer of this analogy).
An instruction stays there until L2 has verified the L1 decisions it depends on.

- **L2 confirms L1** → the instruction retires.
- **L2 finds an L1 decoding error** → the FIFO is **completely drained**. Every
  instruction in it is **uncomputed** (apply its inverse), then **recomputed**
  in its corrected form.

### What actually has to be recomputed

Not everything in the FIFO is a liability:

- A **Clifford at the head of the FIFO retires immediately** — no harm, no foul.
  If the FIFO holds only Cliffords, there is nothing at risk.
- The instructions that matter are Cliffords **sandwiched between
  non-Cliffords**. Those are the ones that must be recomputed.

So the live contents of the FIFO are effectively "instructions since the oldest
unverified non-Clifford," and the squash unit is the span between non-Clifford
operations.

### Why frame errors don't just get absorbed

A tempting-but-wrong intuition (raised and rejected in discussion): a
mispredicted frame bit on a T-teleportation applies `S†` where `S` was correct,
and since desired = `S²`·actual = `Z`·actual, the discrepancy is a **Pauli** and
should fold harmlessly into the frame.

This fails in practice. **A heavily optimized program leaves you with
operations of alternating bases** — a T gate is followed by a non-Clifford in
the X basis, and so on. The discrepancy then meets a gate that conjugates it out
of the Pauli group, so it does not stay absorbable. CNOTs are the partial
exception; those can be delayed somewhat, and that is how they are handled in
the implementation.

**Consequence: effectively every non-Clifford is a commit point.**

### Measurements are a drain point, not a hazard

Another rejected concern: logical measurements have no inverse, so a measurement
inside the speculative shadow would be un-undoable.

This does not arise. At a logical measurement the **reaction-time requirement
disappears** — decoder latency no longer matters and the measurement can simply
be decoded offline at full accuracy. Measurements therefore never sit inside a
live shadow that would need inverting.

## Cost model

Runahead costs **throughput, not correctness**. Slowdown is driven by two terms,
roughly multiplicatively:

```
slowdown  ~  P(L1 decoding error)  ×  2 × (L2 latency × issue rate)
                    |                            |
        how often we squash            how much we undo + redo per squash
```

- **High L1 error rate** → frequent mispredict-and-recompute → little forward
  progress on the program.
- **High L2 latency** → more instructions in flight per squash → more to
  uncompute and recompute.

There is a genuine balance between the two. Iso-slowdown curves are hyperbolas,
so a 10× more accurate L1 buys tolerance for a 10× slower L2.

**Storage** rides on the same axis: FIFO depth = L2 latency × issue rate, i.e.
the same term as the recompute cost. L2 latency therefore sets the sizing of the
runahead microarchitectural structures too. This is a secondary consideration,
but it is another reason to keep L2 latency as small as possible.

## Implications for L1 decoder design

Falls out of the above; useful when evaluating candidate L1 decoders.

1. **L1's figure of merit is its logical error rate, in the coset sense.** L1
   and L2 producing *different corrections* costs nothing if they are
   homologically equivalent — only a different logical coset triggers a squash.
   (When L2 is near-truth, P(L1 disagrees with L2) ≈ L1's LER, so these two
   framings coincide.) This is a weaker bar than "be a good decoder" and is what
   makes a cheap L1 plausible at all.

2. **A detected failure is cheaper than a silent one.** An L1 that can signal "I
   don't know" lets the commit point *stall* rather than speculate into it:

   - stall costs ≈ L2 latency
   - squash costs ≈ 2 × FIFO occupancy, plus re-distilled magic states consumed
     on the wrong path

   So expected cost ≈ `P(silently wrong)·squash + P(knows it doesn't know)·stall`,
   and an L1 with a *worse* raw error rate but high self-detection can win. A
   LUT-ladder decoder gets this for free: ladder exhaustion with lit syndrome
   bits remaining is an explicit unknown signal.

3. **Speculative candidate (not yet validated by the author): a
   problem-reducing L1.** Rather than balancing the two cost terms, an L1 whose
   work also *reduces L2's latency* would improve both terms at once — e.g. an
   L1 that consumes the syndrome bits it confidently explains and hands L2 a
   sparser residual graph with hot spots pre-localized. Flagged as an idea from
   discussion, not as established design — and note this is *not* what
   "transformational" refers to; see the section below.

## Why "transformational"

The paradigm is *transformational* because it **transforms how real-time
decoding has to be tackled at all**. It is a statement about the design
objective, not about any one decoder.

Conventionally, a real-time decoder lives or dies by its latency. Under
runahead, that is no longer true for L2:

> **Slow decoders are fine, as long as they have the throughput to back it up.**

What changes is which constraint binds:

- **Throughput was always fundamental.** This is the backlog problem: consume
  syndrome data more slowly than the machine generates it and the backlog grows
  without bound. Nothing about runahead changes this, and it was never specific
  to runahead.
- **Latency stops being the qualifying constraint for L2.** It shows up only as
  slowdown and storage (see the cost model above): more instructions in flight
  per squash, deeper retirement FIFO. Expensive, but bounded and tunable.

So the challenge is **high-latency, high-throughput** decoders. In principle
throughput is always obtainable with sufficient parallelism — but that is
exactly where the difficulty hides.

### The real difficulty: parallelism scales with latency

By Little's law, sustaining throughput `λ` with per-problem latency `L` requires

```
N  =  λ × L        decoding problems in flight
```

`λ` is fixed by the machine's syndrome generation rate, so **raising latency
raises the required degree of parallelism proportionally** — and for a genuinely
high-latency decoder, `N` can get very large.

That makes the operative design goal:

> **Design decoders such that parallelism is cheap.**

Not "design low-latency decoders," and not merely "design parallelizable
decoders" — the cost that matters is the marginal cost of the `N`-th parallel
instance. Things that make parallelism cheap:

- small per-instance state, so `N` copies stay affordable
- no contended global structure serializing the instances
- read-only state **shared** across instances rather than replicated (a lookup
  table is the clean case: one table, many cheap compute units against it)
- independent sub-problems, so no inter-instance communication

A decoder that is fast but whose parallel instances are individually expensive
can lose to a slow decoder that replicates almost for free.
