# Performance Portability Blocks

Here we can discuss EKAT, Kokkos, and all of the highly-technical non-scientific
stuff that makes our heads hurt.

## Kokkos (View, device, parallel stuff)
Briefly go over what kokkos is, and what we use from it
Also, describe how some aspects are furtherly hidden by ekat (e.g., || reduce)

## Hierarchical parallelism
This could be in the previous section, but it is the key block for scream perf,
so maybe it deserves a small subsection? We can then explain how kokkos implements it,
and how EKAT/SCREAM use kokkos for it (e.g., ekat::TeamUtils).

## Vectorization: Packs
Describe vectorization issue, point to https://gmd.copernicus.org/articles/12/1423/2019/
for a study of packs and pack length on CPU performance.
Explain how/why this does not compromise GPU perf.
