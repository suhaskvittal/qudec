<!-- Author: OpenAI Codex. -->

# MPI groups for concatenated MWPM

Build with `-DENABLE_MPI=ON`. Add `--concat-parallel` to `qudec concat_mwpm`
to distribute its matching problems across a group of MPI ranks.
`--concat-cores N` chooses the group size (default 3 when parallelism is enabled).
The normal serial decoder uses one rank per independently sampled shot stream.

For example:

```sh
mpirun -n 64 qudec concat_mwpm 7 --code cc666 -p 0.0005 \
    --concat-ensemble 20 --concat-flag-renormalization \
    --concat-parallel --concat-cores 16
```

This creates four independent sampling groups: ranks 0–15, 16–31, 32–47 and
48–63. Their sampling roots are 0, 16, 32 and 48. Set the launcher's CPU
binding policy as appropriate for the machine; the decoder groups MPI ranks,
and does not itself set operating-system CPU affinity.

World size must be divisible by the requested group size. Incomplete groups
are rejected. Parallel decoding with more than one rank requires an MPI build
and initialized MPI. The chain executable retains its single-rank contract.

## Matching work

Each process owns its matching graphs. Within a decoder group:

1. The root supplies the detector syndrome, broadcast to the group.
2. First-stage matchings for colors, ensemble members and optional boundary
   sectors are distributed round-robin across ranks.
3. Their selected projection IDs are collected in the serial candidate order.
   Identical virtual syndromes are coalesced using the existing ensemble rule.
4. Only the remaining unique second-stage problems are distributed to ranks.
5. Candidates are collected and the root returns the minimum-weight correction,
   preserving the serial candidate order and tie tolerance.

The baseline, seeded ensembling, boundary-gap extension and flag-conditioned
graph updates use the same models as serial decoding. There are no worker
threads or nested MPI launches per shot. More ranks than available matching
tasks leave some ranks idle during that phase.

## Sampler and API contract

`ConcatMWPM::config_type::decoding_cores` selects the group size in C++.
`ConcatMWPM::decoding_core_count()` reports it. The shared helper
`dec::decoding_core_count(decoder)` returns that value when a decoder exposes
the method, and otherwise returns 1 for existing serial decoders.

Both Monte Carlo and rare-event sampling divide work by **decoder groups**.
Only group roots generate shots, compare predictions with truth, invoke error
callbacks, print shot diagnostics and accumulate decoder shot/latency statistics.
Helpers participate in collective decoding of their group's shots. Global
sample limits, error stopping and rare-event conditional-level counts count a
group's shot once. Uneven batches and batches smaller than the number of groups
are supported.

The random seed is the configured seed plus the group index. Thus a run with
four serial ranks and a run with four multi-rank decoder groups use the same
four random streams. Changing the number of groups changes the sharding.

When invoking the parallel C++ decoder directly, every group member must call
`decode` in the same order. The input syndrome and returned observable are
meaningful on the group's root; helpers receive a placeholder result. Only the
root needs the sampled syndrome and observable truth. Group members must use
the same circuit, decoder configuration and ensemble seed.

This implementation replicates graph preparation and flag updates on each
rank, and communicates between matching stages. Whether grouping improves
latency or total sampling throughput depends on problem size, ensemble size
and MPI costs; parallelism is opt-in.
