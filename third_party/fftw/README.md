# FFTW 3.3.11

Unmodified upstream archive: https://fftw.org/fftw-3.3.11.tar.gz

SHA256: 5630c24cdeb33b131612f7eb4b1a9934234754f9f388ff8617458d0be6f239a1

The default build uses static double-precision FFTW 1D transforms. WCNS supplies its own MPI slab transposes using MPI_Alltoallv; it does not require or call libfftw3_mpi. Set WCNS_ENABLE_FFTW=OFF for the built-in radix-2 backend with the same distributed algorithm.

Upstream COPYRIGHT and COPYING are preserved. FFTW is GPL version 2 or later; this notice does not relicense WCNS. Distribution of a linked executable requires considering FFTW licensing; the source package includes the complete unmodified library archive.
