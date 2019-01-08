# GriD-SSE
Fortran module for Grid-based Determination of Slow Slip Events.

## Dependency
- Lapack
- DC3D by Okada (1992) (http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html)
    - Change REAL\*4 to REAL\*8  
    - Specify the common blocks as threadprivate (e.g., !$omp threadprivate(/c0/) ) for OpenMP parallelization  

## Reference
Takagi, R., Naoki, U., and Obara, K. (2019). Along-strike variation and migration of long-term slow slip events in the western Nankai subduction zone, Japan. Journal of Geophysical Research.
