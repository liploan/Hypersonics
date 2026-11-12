# Verified Technical Bibliography (`REFERENCES.md`)

Every reference in this bibliography has been verified against primary publisher records, DOI registries, the NASA Technical Reports Server (NTRS), and ISBN catalogs. An interactive, searchable HTML version with 1-click BibTeX copy buttons is hosted at [`docs/references.html`](../docs/references.html) (**https://liploan.github.io/Hypersonics/references.html**).

---

## A. Supersonic & Hypersonic Gas Dynamics & ODE Integration

1. **Ames Research Staff (1953)** — *"Equations, Tables, and Charts for Compressible Flow,"* NACA Report 1135, Ames Aeronautical Laboratory, Moffett Field, CA, 1953. NASA NTRS ID: [`19930091059`](https://ntrs.nasa.gov/citations/19930091059).  
   *Implemented in*: `+waverider/+flow/obliqueShock.m`, `+waverider/+flow/taylorMaccoll.m`, `docs/physics.js`, `tests/test_physics_validation.m` (Test A4 conical shock chart benchmark).

2. **Taylor, G. I., and Maccoll, J. W. (1933)** — *"The Air Pressure on a Cone Moving at High Speeds—I & II,"* *Proceedings of the Royal Society of London. Series A*, Vol. 139, No. 838, Feb. 1933, pp. 278–311. DOIs: [`10.1098/rspa.1933.0017`](https://doi.org/10.1098/rspa.1933.0017) and [`10.1098/rspa.1933.0018`](https://doi.org/10.1098/rspa.1933.0018).  
   *Implemented in*: `+waverider/+flow/taylorMaccoll.m`, `+waverider/+flow/traceStreamline.m`, `docs/physics.js`.

3. **Anderson, J. D., Jr. (2006)** — *Hypersonic and High-Temperature Gas Dynamics*, 2nd ed., AIAA Education Series, American Institute of Aeronautics and Astronautics, Reston, VA, 2006. ISBN `978-1-56347-780-5`. DOI: [`10.2514/4.861956`](https://doi.org/10.2514/4.861956).  
   *Implemented in*: `+waverider/WaveriderBody.m`, `+waverider/+flow/obliqueShock.m`, `+waverider/+geometry/generateCSTBody.m`, `docs/physics.js`.

4. **Dormand, J. R., and Prince, P. J. (1980)** — *"A family of embedded Runge-Kutta formulae,"* *Journal of Computational and Applied Mathematics*, Vol. 6, No. 1, March 1980, pp. 19–26. DOI: [`10.1016/0771-050X(80)90013-3`](https://doi.org/10.1016/0771-050X(80)90013-3).  
   *Implemented in*: `+waverider/+flow/taylorMaccoll.m` (`ode45` RK5(4)7M integration).

5. **Shampine, L. F., and Reichelt, M. W. (1997)** — *"The MATLAB ODE Suite,"* *SIAM Journal on Scientific Computing*, Vol. 18, No. 1, Jan. 1997, pp. 1–22. DOI: [`10.1137/S1064827594276424`](https://doi.org/10.1137/S1064827594276424).  
   *Implemented in*: `+waverider/+flow/taylorMaccoll.m` (zero-crossing event detection `u_theta = 0`).

---

## B. Waverider Lofting Methodologies & Hypersonic Aerodynamics

6. **Nonweiler, T. R. F. (1959 & 1963)** — *"Aerodynamic Problems of Manned Space Vehicles,"* *Journal of the Royal Aeronautical Society*, Vol. 63, No. 585, Sept. 1959, pp. 521–528; and **Nonweiler, T. (1963)** — *"Delta Wings of Shapes Amenable to Exact Shock-Wave Theory,"* *Journal of the Royal Aeronautical Society*, Vol. 67, No. 625, Jan. 1963, pp. 39–40. DOI (1963): [`10.1017/S0368393100090027`](https://doi.org/10.1017/S0368393100090027).  
   *Implemented in*: `+waverider/+geometry/generateCaretWedge.m`, `docs/physics.js`, `tests/test_physics_validation.m` (Test E1 exact $L/D = \cot\delta$).

7. **Jones, J. G., Moore, K. C., Pike, J., and Roe, P. L. (1968)** — *"A method for designing lifting configurations for high supersonic speeds, using axisymmetric flow fields,"* *Ingenieur-Archiv* (*Archive of Applied Mechanics*), Vol. 37, No. 1, 1968, pp. 56–72. DOI: [`10.1007/BF00532683`](https://doi.org/10.1007/BF00532683).  
   *Implemented in*: `+waverider/+geometry/generateOsculatingCone.m` (`'conical'`), `docs/physics.js`.

8. **Rasmussen, M. L. (1980)** — *"Waverider Configurations Derived from Inclined Circular and Elliptic Cones,"* *Journal of Spacecraft and Rockets*, Vol. 17, No. 6, Nov.–Dec. 1980, pp. 537–545. DOI: [`10.2514/3.57771`](https://doi.org/10.2514/3.57771).  
   *Implemented in*: `+waverider/+geometry/generateOsculatingCone.m` (`'conical'`), `docs/physics.js`.

9. **Bowcutt, K. G., Anderson, J. D., Jr., and Capriotti, D. P. (1987)** — *"Viscous Optimized Hypersonic Waveriders,"* AIAA Paper 87-0272, *25th AIAA Aerospace Sciences Meeting*, Reno, NV, Jan. 1987. DOI: [`10.2514/6.1987-272`](https://doi.org/10.2514/6.1987-272).  
   *Implemented in*: `+waverider/WaveriderConfig.m`, `+waverider/WaveriderBody.m`, `+waverider/+viz/plotTradeSpace.m`, `examples/demo_parametric_sweep.m`.

10. **Sobieczky, H., Dougherty, F. C., and Jones, K. D. (1990)** — *"Hypersonic Waverider Design from Given Shock Waves,"* *Proceedings of the First International Hypersonic Waverider Symposium*, University of Maryland, College Park, MD, Oct. 17–19, 1990.  
    *Implemented in*: `+waverider/+geometry/generateOsculatingCone.m` (`'osculating_cone'`), `docs/physics.js`.

11. **Kulfan, B. M. (2008)** — *"Universal Parametric Geometry Representation Method,"* *Journal of Aircraft*, Vol. 45, No. 1, Jan.–Feb. 2008, pp. 142–158. DOI: [`10.2514/1.29958`](https://doi.org/10.2514/1.29958).  
    *Implemented in*: `+waverider/+geometry/generateCSTBody.m`, `docs/physics.js`.

12. **Ding, F., Liu, J., Shen, C.-b., Liu, Z., Chen, S., and Huang, W. (2018)** — *"An overview of waverider design concept in airframe/inlet integration methodology for air-breathing hypersonic vehicles,"* *Acta Astronautica*, Vol. 152, Nov. 2018, pp. 639–656. DOI: [`10.1016/j.actaastro.2018.09.002`](https://doi.org/10.1016/j.actaastro.2018.09.002).  
    *Implemented in*: `+waverider/+geometry/generateOsculatingCone.m`, `+waverider/+geometry/appendConformalFrustum.m`.

13. **Lees, L. (1956)** — *"Laminar Heat Transfer Over Blunt-Nosed Bodies at Hypersonic Flight Speeds,"* *Jet Propulsion*, Vol. 26, No. 4, April 1956, pp. 259–269. DOI: [`10.2514/8.6977`](https://doi.org/10.2514/8.6977).  
    *Implemented in*: `+waverider/WaveriderBody.m`, `+waverider/+geometry/generateCSTBody.m`, `docs/physics.js` (modified Newtonian $C_{p,\max}\sin^2\theta$).

14. **Eckert, E. R. G. (1955)** — *"Engineering Relations for Friction and Heat Transfer to Surfaces in High Velocity Flow,"* *Journal of the Aeronautical Sciences*, Vol. 22, No. 8, Aug. 1955, pp. 585–587.  
    *Implemented in*: `+waverider/WaveriderBody.m`, `docs/physics.js` (reference-temperature $T^*$ compressible skin friction).

---

## C. High-Frequency Electromagnetics & Concave Cavity RCS

15. **Gordon, W. B. (1975)** — *"Far-Field Approximations to the Kirchhoff-Helmholtz Representations of Scattered Fields,"* *IEEE Transactions on Antennas and Propagation*, Vol. 23, No. 4, July 1975, pp. 590–592. DOI: [`10.1109/TAP.1975.1141105`](https://doi.org/10.1109/TAP.1975.1141105).  
    *Implemented in*: `+waverider/+rcs/facetPhaseIntegral.m`, `+waverider/+rcs/doubleBounceGOPO.m`, `+waverider/+rcs/computeRCS.m`, `docs/physics.js`.

16. **Ling, H., Chou, R.-C., and Lee, S.-W. (1989)** — *"Shooting and Bouncing Rays: Calculating the RCS of an Arbitrarily Shaped Cavity,"* *IEEE Transactions on Antennas and Propagation*, Vol. 37, No. 2, Feb. 1989, pp. 194–205. DOI: [`10.1109/8.18706`](https://doi.org/10.1109/8.18706).  
    *Implemented in*: `+waverider/+rcs/doubleBounceGOPO.m`, `+waverider/+rcs/buildConcavePairList.m`, `docs/physics.js`.

17. **Knott, E. F., Shaeffer, J. F., and Tuley, M. T. (2004)** — *Radar Cross Section*, 2nd ed., SciTech Publishing, Raleigh, NC, 2004. ISBN `978-1-891121-25-8`.  
    *Implemented in*: `+waverider/+rcs/computeRCS.m`, `+waverider/+rcs/doubleBounceGOPO.m`, `tests/test_physics_validation.m` (Tests C1–C5 canonical PEC plate, sphere & $90^\circ$ dihedral benchmarks).

18. **Ufimtsev, P. Ya. (2014)** — *Fundamentals of the Physical Theory of Diffraction*, 2nd ed., Wiley-IEEE Press, Hoboken, NJ, 2014. ISBN `978-1-118-75366-8`. DOI: [`10.1002/9781118753767`](https://doi.org/10.1002/9781118753767).  
    *Implemented in*: `+waverider/+rcs/computeRCS.m` (leading/trailing-edge wedge diffraction).

---

## D. Wideband Radar Signal Processing, HRRP & RTI Sinograms

19. **Harris, F. J. (1978)** — *"On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform,"* *Proceedings of the IEEE*, Vol. 66, No. 1, Jan. 1978, pp. 51–83. DOI: [`10.1109/PROC.1978.10837`](https://doi.org/10.1109/PROC.1978.10837).  
    *Implemented in*: `+waverider/+rcs/computeRangeProfile.m`, `docs/physics.js`, `tests/test_physics_validation.m` (Tests D1–D3 rectangular, Hann & Hamming PSFs).

20. **Wehner, D. R. (1995)** — *High-Resolution Radar*, 2nd ed., Artech House, Norwood, MA, 1995. ISBN `978-0-89006-727-7`.  
    *Implemented in*: `+waverider/+rcs/computeRangeProfile.m`, `+waverider/+viz/plotRangeAngleIntensity.m`, `docs/physics.js`.

---

## E. Computational Geometry, Polygon Clipping & 2-Manifold Topology

21. **Sutherland, I. E., and Hodgman, G. W. (1974)** — *"Reentrant Polygon Clipping,"* *Communications of the ACM*, Vol. 17, No. 1, Jan. 1974, pp. 32–42. DOI: [`10.1145/360767.360802`](https://doi.org/10.1145/360767.360802).  
    *Implemented in*: `+waverider/+rcs/doubleBounceGOPO.m`, `docs/physics.js`.

22. **van Oosterom, A., and Strackee, J. (1983)** — *"The Solid Angle of a Plane Triangle,"* *IEEE Transactions on Biomedical Engineering*, Vol. BME-30, No. 2, Feb. 1983, pp. 125–126. DOI: [`10.1109/TBME.1983.325207`](https://doi.org/10.1109/TBME.1983.325207).  
    *Implemented in*: `tests/test_physics_validation.m`, `docs/physics.js`, `tests/web/validate_physics.mjs`.

23. **Jacobson, A., Kavan, L., and Sorkine-Hornung, O. (2013)** — *"Robust Inside-Outside Segmentation using Generalized Winding Numbers,"* *ACM Transactions on Graphics*, Vol. 32, No. 4 (SIGGRAPH 2013), July 2013, Article 33. DOI: [`10.1145/2461912.2461916`](https://doi.org/10.1145/2461912.2461916).  
    *Implemented in*: `tests/test_physics_validation.m`, `docs/physics.js`, `tests/web/validate_physics.mjs` (Test B2 self-intersection audit).

24. **Zhang, C., and Chen, T. (2001)** — *"Efficient Feature Extraction for 2D/3D Objects in Mesh Representation,"* *Proceedings 2001 International Conference on Image Processing (ICIP)*, Vol. 3, Oct. 2001, pp. 935–938. DOI: [`10.1109/ICIP.2001.958278`](https://doi.org/10.1109/ICIP.2001.958278).  
    *Implemented in*: `+waverider/+mesh/verifyWatertight.m`, `docs/physics.js` (signed tetrahedral enclosed volume & volumetric centroid).
