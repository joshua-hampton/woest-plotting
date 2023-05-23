import numpy as np
from matplotlib.transforms import Transform

from constants.thermodynamics import (
    CONST_KELVIN,
    CONST_TO,
    CONST_CP_AIR,
    CONST_P0,
    CONST_GAS_CONST_AIR,
    CONST_EPSILON,
    CONST_GAS_CONST_VAP,
    CONST_LATENT_HEAT_VAP_WATER,
    CONST_ES0,
)


def convert_pressure_temperature_to_pressure_theta(pressure, temperature):
    """
    Transform pressure and temperature into pressure and potential temperature.
    :param pressure: pressure in hPa
    :type pressure: list-like
    :param temperature: temperature in degC
    :type temperature: list-like
    :return: Pressure in hPa, potential temperature in degC
    :rtype: tuple
    """
    pressure, temperature = np.asarray(pressure), np.asarray(temperature)

    temperature_kelvin = temperature + CONST_KELVIN

    theta_kelvin = temperature_kelvin * (CONST_P0 / pressure) ** (
        CONST_GAS_CONST_AIR / CONST_CP_AIR
    )
    theta = theta_kelvin - CONST_KELVIN
    return pressure, theta


def convert_temperature_theta_to_temperature_pressure(temperature, theta):
    """
    Transform temperature and potential temperature into temperature and pressure
    :param temperature: temperature in degC
    :type temperature: list-like
    :param theta: potential temperature in degC
    :type theta: list-like
    :return: Pressure in hPa, potential temperature in degC
    :rtype: tuple
    """
    temperature, theta = np.asarray(temperature), np.asarray(theta)

    temperature_kelvin = temperature + CONST_KELVIN
    theta_kelvin = theta + CONST_KELVIN

    pressure = CONST_P0 * (temperature_kelvin / theta_kelvin) ** (
        CONST_CP_AIR / CONST_GAS_CONST_AIR
    )
    return temperature, pressure


def convert_pressure_theta_to_temperature_theta(pressure, theta):
    """
    Transform pressure and potential temperature into pressure and temperature.
    :param pressure: pressure in hPa
    :type pressure: list-like
    :param theta: potential temperature in degC
    :type theta: list-like
    :return: Pressure in hPa, temperature in degC
    :rtype: tuple
    """
    pressure, theta = np.asarray(pressure), np.asarray(theta)

    theta_kelvin = theta + CONST_KELVIN

    temperature_kelvin = theta_kelvin * (CONST_P0 / pressure) ** (
        -CONST_GAS_CONST_AIR / CONST_CP_AIR
    )
    temperature = temperature_kelvin - CONST_KELVIN
    return temperature, theta


def convert_pressure_mixing_ratio_to_temperature(pressure, mixing_ratio):
    mixing_ratio_unitless = mixing_ratio / 1000
    es = (mixing_ratio_unitless * pressure) / (mixing_ratio_unitless + CONST_EPSILON)
    temperature_kelvin = 1.0 / (
        (1.0 / CONST_KELVIN)
        - CONST_GAS_CONST_VAP / CONST_LATENT_HEAT_VAP_WATER * np.log(es / CONST_ES0)
    )

    temperature = temperature_kelvin - CONST_KELVIN

    return temperature


def convert_temperature_theta_to_xy(temperature, theta):
    """
    Converts temperature, potential temperature (theta) coordinates into x,y grid
    coordinates as per Appendix B of Thermal Physics of the Atmosphere (2nd ed)
    :param temperature: dry-bulb temperature in degC
    :type temperature: list-like
    :param theta: potential temperature in degC
    :type theta: list-like
    :return: grid coordinates x,y
    :rtype: tuple
    """

    temperature, theta = np.asarray(temperature), np.asarray(theta)
    # theta = np.clip(theta, 1, 1e10)

    temperature_kelvin = temperature + CONST_KELVIN
    theta_kelvin = theta + CONST_KELVIN

    x_coords = np.log(theta_kelvin / CONST_TO) + temperature_kelvin / CONST_TO
    y_coords = np.log(theta_kelvin / CONST_TO) - temperature_kelvin / CONST_TO

    return x_coords, y_coords


def convert_xy_to_temperature_theta(x_coords, y_coords):
    """
    Inverse conversion of convert_temperature_theta_to_xy()
    :param x_coords: x coordinates on grid
    :type x_coords: list-like
    :param y_coords: y coordinates on grid
    :type y_coords: list-like
    :return: temperature, potential temperature (theta)
    :rtype: tuple
    """
    x_coords = np.asarray(x_coords)
    y_coords = np.asarray(y_coords)

    temperature_kelvin = CONST_TO * (x_coords - y_coords) / 2
    theta_kelvin = CONST_TO * np.exp((x_coords + y_coords) / 2)

    temperature = temperature_kelvin - CONST_KELVIN
    theta = theta_kelvin - CONST_KELVIN

    return temperature, theta


class TephigramTransform(Transform):
    """
    Inherits matplotlib.transforms.Transform as a custom transformation class to
    convert temperature and potential temperature to x,y plotting grid coordinates
    """

    # Override attributes
    input_dims = 2
    output_dims = 2
    is_separable = False
    has_inverse = True

    def transform_non_affine(self, values):
        return np.concatenate(
            convert_temperature_theta_to_xy(values[:, 0:1], values[:, 1:2]), axis=1
        )

    def inverted(self):
        return TephigramTransformInverse()


class TephigramTransformInverse(Transform):
    """Inverse transformation of TephigramTransform class"""

    # Override attributes
    input_dims = 2
    output_dims = 2
    is_separable = False
    has_inverse = True

    def transform_non_affine(self, values):
        return np.concatenate(
            convert_xy_to_temperature_theta(values[:, 0:1], values[:, 1:2]), axis=1
        )

    def inverted(self):
        return TephigramTransform()


if __name__ == "__main__":
    print(convert_temperature_theta_to_xy(0, 0))
    print(convert_temperature_theta_to_xy(10, 10))
    print(convert_temperature_theta_to_xy(20, 20))
    print(convert_temperature_theta_to_xy(30, 30))
    print(
        np.array(
            convert_xy_to_temperature_theta(
                convert_temperature_theta_to_xy([-10, 40], [10, 30])[0],
                convert_temperature_theta_to_xy([-10, 40], [10, 30])[1],
            )
        )
    )
