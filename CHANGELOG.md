# Changelog

## unreleased 
### Added 
- add an option to control the minimum depth used in the vertical strain-rates fit in utils.computeStrainRates.[362ab82](https://github.com/ldeo-glaciology/xapres/pull/99/commits/362ab827688f613bea457c4cc99af7fe02314bdf)
- add computation of the variance in the estimates of strain rate and surface intercept when performing linear fits to the vertical velocity profiles. 

### Changed
- the name of the parameter to control the maximum depth of the data used in strain-rates fit is now max_depth_for_ezz_fit. [lower_limit_on_fit still works, but we print a warning that people should use max_depth_for_ezz_fit instead. [362ab82](https://github.com/ldeo-glaciology/xapres/pull/99/commits/362ab827688f613bea457c4cc99af7fe02314bdf)
- the least-squares fits to the vertical velocity profiles are now weighted by the inverse of the variance of the estimates of the vertical velocity. 
- Errors in how we computed the variance in the phase, displacement and vertical velocity were corrected (https://github.com/ldeo-glaciology/xapres/issues/98). [c438bfa](https://github.com/ldeo-glaciology/xapres/pull/99/commits/c438bfa2627e39a92687d36af2705bf7c685c3d3)

## 0.5.3 - 2025-01-27
### Fixed 

- merge changes in upstream fork of bas-apres (commit in that sub-module: [59dbd62](https://github.com/jkingslake/bas-apres/pull/2/commits/59dbd6211274b2260cba6a57020358f57e6972b6))

## 0.4.0 - 2025-01-16 

### Changed

- Update load.py to using bas-apres instead of our own loading code

### Added

- Add many new tests

### Removed

- Remove ApRES_plot.py

### Fixed

- move chirp_time dtype correction [8c9a5fe](https://github.com/ldeo-glaciology/xapres/pull/67/commits/8c9a5fe49852fac55f0b94622101ce66fd20941b)
- deal with case when nattenuators is different in different bursts in attended lod [bf3028b](https://github.com/ldeo-glaciology/xapres/pull/67/commits/bf3028b8dd80477c26bef68321dc840245e424c1)
- fix bug in attended load using from_dats.load() [8c9a5fe](https://github.com/ldeo-glaciology/xapres/pull/67/commits/8c9a5fe49852fac55f0b94622101ce66fd20941b)
- fix bug aligning profiles in computeProfile [44487cc](https://github.com/ldeo-glaciology/xapres/pull/67/commits/44487cc7ec1513ebf77f0e8fa160056ff235e136)

## 0.3.0 - 2025-01-15

### Added

- Add the displacement_timeseries method in utils
- Add the compute_displacement method in utils
- Add a number of methods that displacement_timeseries and compute_displacement call in utils
- Add the computeProfile method in utils
- Add documentation hosted on github pages
- Add material to usingXApRES.ipynb (now superseded by the documentation)
