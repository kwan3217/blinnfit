# SPK - trajectory kernels
  1. `jup???.bsp` - Trajectories of the moons
     of Jupiter. Should use the version that
     corresponds to the Voyager trajectory
     to be used.
  2. `sat???.bsp` - Trajectories of the moons
     of Saturn. Same caveat for versions.
  3. `vgr1_jup230.bsp` - trajectory consistent
     with jup230.bsp
  4. `vgr2_ura083.bsp` - trajectory consistent with ura083, 2007. There is a newer
     ura111 documented in [Jacobson 2014](http://dx.doi.org/10.1088/0004-6256/148/5/76)
     but no revised Voyager trajectory with it. I haven't found the paper that says that
     it describes ura083 but [AAS 07-319](https://nextcloud.kwansystems.org/index.php/f/106018) is about the right time and
     correct author.
  5. `nep081.bsp` - 801 Triton, 802 Neried, 
     808 Proteus (N1), documented in [Jacobson 2009](http://dx.doi.org/10.1088/0004-6256/137/5/4322)
  6. `vgr2_nep081.bsp` - Voyager trajectory consistent with nep081.
     Documented as being an update to nep076 solution, which
     itself is documented in [Jacobson 2008](https://doi.org/10.2514/6.2008-7372)
  7. `nep095.bsp` - Triton, Neried, 803 Naiad (N6),
     804 Thalassa (N5), 805 Despina (N3), 806 Galatea (N4),
     807 Larissa (N2), Proteus. 
  8. `nep097.bsp` - Triton. Documented in [Brozovic 2022](https://doi.org/10.1016/j.icarus.2019.113462)
  9. `vgr2_nep097.bsp` - Documented in Brozovic 2022 
 10. `de440.bsp` - JPL DE/LE 440, 2020. Covers AD1550-2650. It has
     a subset `de440s.bsp` with identical positions but shorter coverage
     1850-2150 for a smaller file size. The corresponding `de441.bsp`
     long-term ephemeris covers 13,200 BC to 17,191 AD but with larger file
     size and lower accuracy over historical times. Use `de440s.bsp` for
     the Voyager mission.
11. `voyager_1.ST+1991_a54418u.merged.bsp` and `voyager_2.ST+1992_m05208u.merged.bsp`
    These are supertrajectories -- they are not guaranteed to be the highest quality
    in between encounters but are "good enough." They even show some TCMs. Note -- for
    the super-trajectories, there are newer versions on the NAIF web site, but the 
    comments in those kernels specifically say that they are regressions to the 
    versions used here.