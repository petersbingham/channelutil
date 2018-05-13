# channelutil
Python package to convert between asymptotic wavenumbers and energies for multichannel, inelastic collisions.

## Installation

Clone the repository and install with the following commands:

    git clone https://github.com/petersbingham/channelutil.git
    cd channelutil
    python setup.py install
    
## Dependencies
Author Libraries (these will have their own dependencies):
 - pynumwrap https://github.com/petersbingham/pynumwrap

## Usage
As in the introduction. The various calculations supported are mostly simple and can be ascertained by looking in the channelutil/\_\_init\_\_.py file.

### Constructor
The channels are described via the parameters supplied in the constructor to ```AsymCalc``` class in the package scope. It looks like:
```python
__init__(self, units, ls=None, thresholds=None, sign_sel=signs_pos, signs=None)
```
`units` specifies the units to be used in the calculations; either `rydbergs` or `hartrees`, both of which are in package scope.

If the class is constructed using the default parameters then it will provide calculations for a single channel reaction. 

If only `ls` is supplied then `AsymCalc` will provide calculations for a multichannel elastic reaction, otherwise thresholds can be specified using the `thresholds` parameter.

If only `thresholds` is supplied then `AsymCalc` will provide calculations for zero angular momenta in all channels.

`ls` contains the angular momenta for each for the channels. It should be of equal length to the `thresholds` parameter (if supplied) with the indices of the two lists describing the same channels. If `ls` is not supplied then the calculation will be for zero angular momenta in all channels.

`sign_sel` and `signs` parameters. See next section.

### Sign Selection
The package was developed as a means to provide selection of signs when converting from electron energies to wavenumbers that arise from the plus and minus signs in front of the square root energy wavenumber relationship. See "S.A. Rakityansky P.O.G. Ogunbade. S-matrix parametrization as a way of locating quantum resonances and bound states:multichannel case, 2010" for details. Some selections are provided, see below for further details. We plan to provide more elaborate selections in the future. All of the selections below exist in the package scope.
 - `signs_pos` : Uses positive signs always.
 - `signs_specified` : Uses signs as specified with the constructor `signs` parameter.
 - `signs_bndandres` : Chooses signs to place on the unphysical sheet above threshold and physical sheet below threshold.
 - `signs_ana_over_axis` : Chooses signs to attempt to retain analyticity when crossing the real axis. 
 - `signs_ana_over_thres` : Chooses signs to attempt to retain analyticity when crossing the thresholds.

### Example
```python
>>> import channelutil as cu
>>> calc1 = cu.AsymCalc(thresholds=[0.0,1.0], sign_sel=cu.signs_ana_over_axis)
>>> calc2 = cu.AsymCalc(thresholds=[0.0,1.0], sign_sel=cu.signs_ana_over_thres)
>>> calc1.k(0,1.0+0.1j)
(1.0012461141278126+0.04993777183700243j)
>>> calc1.k(0,1.0-0.1j)
(1.0012461141278126-0.04993777183700243j)  # No discontinuity
>>> calc2.k(0,1.0+0.1j)
(-1.0012461141278126-0.04993777183700243j)
>>> calc2.k(0,1.0-0.1j)
(1.0012461141278126-0.04993777183700243j)  # Discontinuity
```
