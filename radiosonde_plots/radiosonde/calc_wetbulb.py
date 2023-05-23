import numpy as np
import scipy.integrate as si
import scipy.optimize as so
import constants.thermodynamics as therm_consts


def vapor_pressure(pressure, mixing_ratio):
    """Calculate water vapor (partial) pressure.
    Given total ``pressure`` and water vapor ``mixing_ratio``, calculates the
    partial pressure of water vapor.

    Parameters
    ----------
    pressure : Total atmospheric pressure
    mixing_ratio : Dimensionless mass mixing ratio

    Returns
    -------
    Ambient water vapor (partial) pressure in the same units as ``pressure``

    Notes
    -----
    This function is a straightforward implementation of the equation given
    in many places, such as [Hobbs1977]_ pg.71:
    .. math:: e = p \frac{r}{r + \epsilon}
    """
    return pressure * mixing_ratio / (therm_consts.CONST_EPSILON + mixing_ratio)


def saturation_vapor_pressure(temperature):
    """Calculate the saturation water vapor (partial) pressure.

    Parameters
    ----------
    temperature : Air temperature

    Returns
    -------
     Saturation water vapor (partial) pressure

    Notes
    -----
    Instead of temperature, dewpoint may be used in order to calculate
    the actual (ambient) water vapor (partial) pressure.
    The formula used is that from [Bolton1980]_ for T in degrees Celsius:
    .. math:: 6.112 e^\frac{17.67T}{T + 243.5}
    """
    # Converted from original in terms of C to use kelvin. Using raw absolute values
    # of C in a formula plays havoc with units support.
    return (
        therm_consts.CONST_ES0
        * 100
        * np.exp(
            17.67 * (temperature - therm_consts.CONST_KELVIN) / (temperature - 29.65)
        )
    )


def mixing_ratio(
    partial_press, total_press, molecular_weight_ratio=therm_consts.CONST_EPSILON
):
    """Calculate the mixing ratio of a gas.
    This calculates mixing ratio given its partial pressure and the total pressure of
    the air. There are no required units for the input arrays, other than that
    they have the same units.

    Parameters
    ----------
    partial_press :  Partial pressure of the constituent gas
    total_press : Total air pressure
    molecular_weight_ratio :
        The ratio of the molecular weight of the constituent gas to that assumed
        for air. Defaults to the ratio for water vapor to dry air
        (:math:`\epsilon\approx0.622`).

    Returns
    -------
    The (mass) mixing ratio, dimensionless (e.g. Kg/Kg or g/g)

    Notes
    -----
    This function is a straightforward implementation of the equation given
    in many places, such as [Hobbs1977]_ pg.73:
    .. math:: r = \epsilon \frac{e}{p - e}
    .. versionchanged:: 1.0
       Renamed ``part_press``, ``tot_press`` parameters to ``partial_press``,
       ``total_press``
    """
    return molecular_weight_ratio * partial_press / (total_press - partial_press)


def saturation_mixing_ratio(total_press, temperature):
    """Calculate the saturation mixing ratio of water vapor.
    This calculation is given total atmospheric pressure and air temperature.

    Parameters
    ----------
    total_press: Total atmospheric pressure
    temperature: Air temperature

    Returns
    -------
    Saturation mixing ratio, dimensionless

    Notes
    -----
    This function is a straightforward implementation of the equation given
    in many places, such as [Hobbs1977]_ pg.73:
    .. math:: r_s = \epsilon \frac{e_s}{p - e_s}
    """
    return mixing_ratio(saturation_vapor_pressure(temperature), total_press)


def dewpoint(vapor_pressure):
    """Calculate the ambient dewpoint given the vapor pressure.

    Parameters
    ----------
    e : Water vapor partial pressure

    Returns
    -------
    Dewpoint temperature

    See Also
    --------
    dewpoint_from_relative_humidity, saturation_vapor_pressure, vapor_pressure

    Notes
    -----
    This function inverts the [Bolton1980]_ formula for saturation vapor
    pressure to instead calculate the temperature. This yield the following
    formula for dewpoint in degrees Celsius:

    .. math:: T = \frac{243.5 log(e / 6.112)}{17.67 - log(e / 6.112)}
    """
    val = np.log(vapor_pressure / (therm_consts.CONST_ES0 * 100))
    return 243.5 * val / (17.67 - val)


def lcl(pressure, temperature, dewpoint, max_iters=50, eps=1e-5):
    """Calculate the lifted condensation level (LCL) from the starting point.
    The starting state for the parcel is defined by `temperature`, `dewpoint`,
    and `pressure`. If these are arrays, this function will return a LCL
    for every index. This function does work with surface grids as a result.

    Parameters
    ----------
    pressure : Starting atmospheric pressure
    temperature : Starting temperature
    dewpoint : Starting dewpoint

    Returns
    -------
    LCL pressure, LCL temperature

    Other Parameters
    ----------------
    max_iters : int, optional
        The maximum number of iterations to use in calculation, defaults to 50.
    eps : float, optional
        The desired relative error in the calculated value, defaults to 1e-5.

    Notes
    -----
    This function is implemented using an iterative approach to solve for the
    LCL. The basic algorithm is:
    1. Find the dewpoint from the LCL pressure and starting mixing ratio
    2. Find the LCL pressure from the starting temperature and dewpoint
    3. Iterate until convergence
    The function is guaranteed to finish by virtue of the `max_iters` counter.
    Only functions on 1D profiles (not higher-dimension vertical cross sections
    or grids). Since this function returns scalar values when given a profile,
    this will return Pint
    """

    def _lcl_iter(p, p0, w, t):
        nonlocal nan_mask
        td = globals()["dewpoint"](vapor_pressure(p, w)) + therm_consts.CONST_KELVIN
        p_new = p0 * (td / t) ** (
            therm_consts.CONST_CP_AIR / therm_consts.CONST_GAS_CONST_AIR
        )  # .m
        nan_mask = nan_mask | np.isnan(p_new)
        return np.where(np.isnan(p_new), p, p_new)

    # Handle nans by creating a mask that gets set by our _lcl_iter function if it
    # ever encounters a nan, at which point pressure is set to p, stopping iteration.
    nan_mask = False
    w = mixing_ratio(saturation_vapor_pressure(dewpoint), pressure)
    lcl_p = so.fixed_point(
        _lcl_iter,
        np.abs(pressure),
        args=(np.abs(pressure), w, temperature),
        xtol=eps,
        maxiter=max_iters,
    )
    lcl_p = np.where(nan_mask, np.nan, lcl_p)

    # np.isclose needed if sfc is LCL due to precision error with np.log in dewpoint.
    # Causes issues with parcel_profile_with_lcl if removed. Issue #1187
    lcl_p = np.where(np.isclose(lcl_p, np.abs(pressure)), np.abs(pressure), lcl_p)

    return (
        lcl_p,
        globals()["dewpoint"](vapor_pressure(lcl_p, w)) + therm_consts.CONST_KELVIN,
    )


def moist_lapse(pressure, temperature, reference_pressure=None):
    """Calculate the temperature at a level assuming liquid saturation processes.
    This function lifts a parcel starting at `temperature`. The starting pressure can
    be given by `reference_pressure`. Essentially, this function is calculating moist
    pseudo-adiabats.
    Parameters
    ----------
    pressure : Atmospheric pressure level(s) of interest
    temperature : Starting temperature
    reference_pressure : Reference pressure; if not given, it defaults to the first
                         element of the pressure array.
    Returns
    -------
    Resulting parcel temperature at levels given by `pressure`
    """

    def _greater_or_close(a, value, **kwargs):
        """Compare values for greater or close to boolean masks.
        Returns a boolean mask for values greater than or equal to a target
        within a specified absolute or relative tolerance (as in :func:`numpy.isclose`).
        Parameters
        ----------
        a : array-like
            Array of values to be compared
        value : float
            Comparison value
        Returns
        -------
        array-like
            Boolean array where values are greater than or nearly equal to value.
        """
        return (a > value) | np.isclose(a, value, **kwargs)

    def dt(t, p):
        # t = units.Quantity(t, temperature.units)
        # p = units.Quantity(p, pressure.units)
        rs = saturation_mixing_ratio(p, t)
        frac = (
            therm_consts.CONST_GAS_CONST_AIR * t
            + therm_consts.CONST_LATENT_HEAT_VAP_WATER * rs
        ) / (
            therm_consts.CONST_CP_AIR
            + (
                therm_consts.CONST_LATENT_HEAT_VAP_WATER**2
                * rs
                * therm_consts.CONST_EPSILON
                / (therm_consts.CONST_GAS_CONST_AIR * t * t)
            )
        )

        return np.abs(frac / p)

    pressure = np.atleast_1d(pressure)
    if reference_pressure is None:
        reference_pressure = pressure[0]

    if np.isnan(reference_pressure):
        return np.full(pressure.shape, np.nan)

    pressure = pressure  # hPa
    reference_pressure = reference_pressure  # hPa
    # org_units = temperature.units
    temperature = np.atleast_1d(temperature)  # .to('kelvin')

    side = "left"

    pres_decreasing = pressure[0] > pressure[-1]
    if pres_decreasing:
        # Everything is easier if pressures are in increasing order
        pressure = pressure[::-1]
        side = "right"

    ref_pres_idx = np.searchsorted(
        np.abs(pressure), np.abs(reference_pressure), side=side
    )

    ret_temperatures = np.empty((0, temperature.shape[0]))

    if _greater_or_close(reference_pressure, pressure.min()):
        # Integrate downward in pressure
        pres_down = np.append(
            np.abs(reference_pressure), np.abs(pressure[(ref_pres_idx - 1) :: -1])
        )
        trace_down = si.odeint(dt, np.abs(temperature).squeeze(), pres_down.squeeze())
        ret_temperatures = np.concatenate((ret_temperatures, trace_down[:0:-1]))

    if reference_pressure < pressure.max():
        # Integrate upward in pressure
        pres_up = np.append(np.abs(reference_pressure), np.abs(pressure[ref_pres_idx:]))
        trace_up = si.odeint(dt, np.abs(temperature).squeeze(), pres_up.squeeze())
        ret_temperatures = np.concatenate((ret_temperatures, trace_up[1:]))

    if pres_decreasing:
        ret_temperatures = ret_temperatures[::-1]

    return ret_temperatures.T.squeeze()


def wet_bulb_temperature(pressure, temperature, dewpoint):
    """Calculate the wet-bulb temperature using Normand's rule.
    This function calculates the wet-bulb temperature using the Normand method.
    The LCL is computed, and that parcel brought down to the starting pressure
    along a moist adiabat.
    The Normand method (and others) are described and compared by [Knox2017]_.
    Parameters
    ----------
    pressure : Initial atmospheric pressure
    temperature : Initial atmospheric temperature
    dewpoint : Initial atmospheric dewpoint
    Returns
    -------
        Wet-bulb temperature
    Notes
    -----
    Since this function iteratively applies a parcel calculation, it should be used with
    caution on large arrays.
    """
    # if not hasattr(pressure, 'shape'):
    #     pressure = np.atleast_1d(pressure)
    #     temperature = np.atleast_1d(temperature)
    #     dewpoint = np.atleast_1d(dewpoint)

    lcl_press, lcl_temp = lcl(pressure, temperature, dewpoint)

    it = np.nditer(
        [np.abs(pressure), np.abs(lcl_press), np.abs(lcl_temp), None],
        op_dtypes=["float", "float", "float", "float"],
        flags=["buffered"],
    )

    for press, lpress, ltemp, ret in it:
        moist_adiabat_temperatures = moist_lapse(press, ltemp, lpress)
        ret[...] = moist_adiabat_temperatures  # .m_as(temperature.units)

    # If we started with a scalar, return a scalar
    ret = it.operands[3]
    if ret.size == 1:
        ret = ret[0]
    return ret
