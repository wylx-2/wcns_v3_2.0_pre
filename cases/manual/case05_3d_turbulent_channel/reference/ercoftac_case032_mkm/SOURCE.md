# ERCOFTAC Case 032 channel DNS reference

- ERCOFTAC case page supplied for this comparison:
  <http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case032>
- Primary `Re_tau=180` data file linked by that page:
  <http://cfd.mace.manchester.ac.uk/ercoftac/lib/exe/fetch.php?media=cdata:case032:simul1.dat>
- Retrieved: 2026-09-17

The primary comparison in this directory uses the ERCOFTAC `simul1.dat` file
verbatim. Its embedded header gives `Re_delta=3250`, `Re_tau=180` and
`Re_theta=282`, and states that the data were compiled from:

> J. Kim, P. Moin and R. Moser, "Turbulence statistics in fully developed
> channel flow at low Reynolds number," Journal of Fluid Mechanics 177,
> 133--166 (1987).

All quantities are normalized by `u_tau` and the kinematic viscosity unless stated
otherwise; `delta` is the channel half-width. The first 65-row block contains
`y/delta`, `y+`, `U+`, `uu+`, `vv+`, `ww+` and `uv+` and is used here.

For provenance and cross-checking, the retrieved HTML page and the later MKM
`chan180.means` / `chan180.reystress` files from the Oden Institute mirror are also
retained, but the plots and error metrics use only ERCOFTAC `simul1.dat`.
