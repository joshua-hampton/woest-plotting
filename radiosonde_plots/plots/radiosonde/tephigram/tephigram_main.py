import numpy as np
from functools import partial
import os.path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from plots.radiosonde.tephigram.tephigram_transforms import (
    TephigramTransform,
    convert_pressure_temperature_to_pressure_theta,
    convert_temperature_theta_to_temperature_pressure,
)
import plots.radiosonde.tephigram.isopleths as isopleths
import plots.radiosonde.tephigram.labels as labels


class _PlotGroup:
    """
    Container for a related group of tephigram isopleths.
    Manages the creation and plotting of all isopleths within the group.
    """

    def __init__(
        self,
        axes,
        plot_func,
        levels,
        text_kwargs=None,
    ):
        self.axes = axes
        self.text_kwargs = text_kwargs
        self.levels = levels

        for level in self.levels:
            plot_func(level)


class _PlotLabel:
    def __init__(
        self,
        axes,
        plot_func,
        levels,
        text_kwargs=None,
    ):
        self.axes = axes
        self.text_kwargs = text_kwargs
        self.levels = levels

        for level in self.levels:
            plot_func(level)


class Tephigram:
    """Generates a tephigram of one or more pressure and tempereature datasets."""

    def __init__(self):
        # Create figure
        # fig = plt.figure()
        # ax = fig.add_subplot()
        #
        # ax.spines['bottom'].set_color('red')
        # ax.spines['top'].set_color('red')
        #
        # plt.show()

        self.figure = plt.figure(figsize=(7.8565, 11.1055))

        # Tephigram transformation
        self.tephi_transform = TephigramTransform()

        # Intialise subplot
        self.axes = self.figure.add_subplot()
        self.axes.axis("off")
        self.axes.set_axis_off()

        # No borders
        plt.subplots_adjust(
            left=0, bottom=0, right=1.0, top=1.0, wspace=None, hspace=None
        )

        # Configure edge axes properties
        # self.axes.axis["top"].toggle(all=False)
        # self.axes.axis["left"].toggle(all=False)
        # self.axes.axis["bottom"].toggle(all=False)
        # self.axes.axis["right"].toggle(all=False)
        # self.axes.gridlines.set_linestyle("solid")

        # self.axes.patch.set_visible(False)
        # self.figure.patch.set_visible(False)

        # Drawing
        self.transform = self.tephi_transform + self.axes.transData
        self.axes.tephi_transform = self.tephi_transform
        self.axes.tephi_inverse = self.tephi_transform.inverted()

        # Draw isotherms
        isotherms_func = partial(
            isopleths.isotherm,
            50,
            1050,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.08},
        )
        _PlotGroup(self.axes, isotherms_func, np.arange(-90, 70, 1))
        isotherms_func = partial(
            isopleths.isotherm,
            50,
            1050,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.30},
        )
        _PlotGroup(self.axes, isotherms_func, np.arange(-90, 70, 10))

        # Draw isentropes
        isentropes_func = partial(
            isopleths.isentropes,
            -90,
            70,
            50,
            1050,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.08},
        )
        _PlotGroup(self.axes, isentropes_func, np.arange(-90, 250, 10))

        # Draw isobars
        isobars_func = partial(
            isopleths.isobar,
            -90,
            70,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.08},
        )
        _PlotGroup(self.axes, isobars_func, np.arange(50, 1051, 10))
        isobars_func = partial(
            isopleths.isobar,
            -90,
            70,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.30},
        )
        _PlotGroup(self.axes, isobars_func, np.arange(100, 1051, 100))

        # Draw moist adiabats
        moist_adiabats_func = partial(
            isopleths.moist_adiabat,
            -50,
            1050,
            1000,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.30},
        )
        _PlotGroup(self.axes, moist_adiabats_func, np.arange(-40, 70, 10))
        moist_adiabats_func = partial(
            isopleths.moist_adiabat,
            -50,
            1050,
            1000,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.08},
        )
        _PlotGroup(self.axes, moist_adiabats_func, np.arange(-40, 70, 2))

        # Draw mixing ratios
        mixing_ratios_func = partial(
            isopleths.mixing_ratio,
            -50,
            50,
            1050,
            self.axes,
            self.transform,
            {"color": "#23CE1F", "linewidth": 0.25, "linestyle": "--"},
        )
        _PlotGroup(
            self.axes,
            mixing_ratios_func,
            np.array(
                [
                    0.10,
                    0.15,
                    0.20,
                    0.30,
                    0.40,
                    0.50,
                    0.60,
                    0.80,
                    1,
                    1.5,
                    2,
                    2.5,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    12,
                    14,
                    16,
                    18,
                    20,
                    24,
                    28,
                    32,
                    36,
                    40,
                    44,
                    48,
                    52,
                    56,
                    60,
                    68,
                    80,
                ]
            ),
        )

        # Isotherm Labels
        isotherm_label_list = np.arange(-40, 70, 10)
        isotherm_label_func = partial(
            labels.isotherm_label, 1000, self.axes, self.transform
        )
        _PlotLabel(self.axes, isotherm_label_func, isotherm_label_list)

        isotherm_label_list = np.arange(-80, 0, 10)
        isotherm_label_func = partial(
            labels.isotherm_label, 190, self.axes, self.transform
        )
        _PlotLabel(self.axes, isotherm_label_func, isotherm_label_list)

        # Isentrope Labels
        isentrope_label_list = np.arange(-40, 220, 20)
        isentrope_label_func = partial(
            labels.isentrope_label, -45, self.axes, self.transform
        )
        _PlotLabel(self.axes, isentrope_label_func, isentrope_label_list)

        # Isobar Labels
        isobar_label_list = np.array([50, 60, 70, 80, 90, 100, 150, 200])
        isobar_label_func = partial(
            labels.isobar_label, "isotherm", -90, self.axes, self.transform
        )
        _PlotLabel(self.axes, isobar_label_func, isobar_label_list)

        isobar_label_list = np.array([300])
        isobar_label_func = partial(
            labels.isobar_label, "isentrope", 0, self.axes, self.transform
        )
        _PlotLabel(self.axes, isobar_label_func, isobar_label_list)

        isobar_label_list = np.arange(400, 1050, 100)
        isobar_label_func = partial(
            labels.isobar_label, "isentrope", -10, self.axes, self.transform
        )
        _PlotLabel(self.axes, isobar_label_func, isobar_label_list)

        # Mixing Ratio Labels
        mixing_ratio_label_list = np.array(
            [
                0.10,
                0.15,
                0.20,
                0.30,
                0.40,
                0.50,
                0.60,
                0.80,
                1,
                1.5,
                2,
                2.5,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                12,
                14,
                16,
                18,
                20,
                24,
                28,
                32,
                36,
                40,
                44,
                48,
                52,
                56,
                60,
                68,
                80,
            ]
        )
        mixing_ratio_label_func = partial(
            labels.mixing_ratio_label, 1054, "right", self.axes, self.transform
        )
        _PlotLabel(self.axes, mixing_ratio_label_func, mixing_ratio_label_list)
        mixing_ratio_label_func = partial(
            labels.mixing_ratio_label, 496, "left", self.axes, self.transform
        )
        _PlotLabel(self.axes, mixing_ratio_label_func, mixing_ratio_label_list)

        # Retain aspect ratio
        self.axes.set_aspect(1.0)

        # Limits
        self.axes.set_xlim(0.48, 1.26)
        self.axes.set_ylim(-1.06, -0.18)

        # Initialise blank profile lists
        self._profiles = []
        self.axes.tephigram_profiles = self._profiles

    def plot_profile(self, pressures, temperatures, **kwargs):
        profile = Profile(pressures, temperatures, self.axes)
        profile.plot(**kwargs)
        self._profiles.append(profile)
        return profile

    def plot_barbs(self, pressures, wind_speeds, wind_directions, **kwargs):
        barbs = Barbs(pressures, wind_speeds, wind_directions, self.axes)
        barbs.plot(**kwargs)

    def plot_main_title(self, metadata, **kwargs):
        title = Title(metadata, self.axes)
        title.plot_main_title()

    def plot_dorset_title(self, metadata):
        title = DorsetTitle(metadata, self.axes)
        title.plot_main_title()

    def read_metadata(
        self, metadata, sonde_lookup_filepath="../../../lookup/SondeStations.txt"
    ):
        self.station_id = (
            f"{metadata.loc['WMO_BLCK_NMBR', 'info']:02.0f}"
            f"{metadata.loc['WMO_STTN_NMBR', 'info']:03.0f}"
        )
        self.year = f"{metadata.loc['YEAR', 'info']:.0f}"
        self.month = f"{metadata.loc['MONTH', 'info']:02.0f}"
        self.day = f"{metadata.loc['DAY', 'info']:02.0f}"
        self.hour = f"{metadata.loc['HOUR', 'info']:02.0f}"
        self.minute = f"{metadata.loc['MINT', 'info']:02.0f}"
        if os.path.isfile(sonde_lookup_filepath):
            self.sonde_lookup_table = pd.read_csv(
                sonde_lookup_filepath, sep=" ", header=None, index_col=0, dtype="string"
            )
            self.sonde_lookup_table.columns = ["name"]
            self.sonde_lookup_table["name"] = self.sonde_lookup_table[
                "name"
            ].str.replace("_", " ")
            self.station_name = self.sonde_lookup_table.loc[self.station_id]["name"]

    def read_metadata_dorset(self, metadata):
        self.station_id = None
        self.year = f"{metadata.loc['YEAR', 'info']:.0f}"
        self.month = f"{metadata.loc['MONTH', 'info']:02.0f}"
        self.day = f"{metadata.loc['DAY', 'info']:02.0f}"
        self.hour = f"{metadata.loc['HOUR', 'info']:02.0f}"
        self.minute = f"{metadata.loc['MINT', 'info']:02.0f}"

    def save_tephi(self, output_dir):
        plt.savefig(
            f"{output_dir}/"
            f"{self.station_name}_{self.station_id}_"
            f"{self.year}{self.month}{self.day}{self.hour}{self.minute}Z.pdf",
            bbox_inches="tight",
            pad_inches=0,
            backend="pgf",
        )
        plt.savefig(
            f"{output_dir}/"
            f"{self.station_name}_{self.station_id}_"
            f"{self.year}{self.month}{self.day}{self.hour}{self.minute}Z.png",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0,
        )

    def save_tephi_manual(self, output_path, **kwargs):
        plt.savefig(output_path, bbox_inches="tight", pad_inches=0, **kwargs)


class Title:
    """Generate a title for the tephigram"""

    def __init__(
        self, metadata, axes, sonde_lookup_filepath="../../../lookup/SondeStations.txt"
    ):
        os.chdir(os.path.dirname(__file__))
        self.axes = axes
        self.station_id = (
            f"{metadata.loc['WMO_BLCK_NMBR', 'info']:02.0f}"
            f"{metadata.loc['WMO_STTN_NMBR', 'info']:03.0f}"
        )
        self.year = f"{metadata.loc['YEAR', 'info']:.0f}"
        self.month = f"{metadata.loc['MONTH', 'info']:02.0f}"
        self.day = f"{metadata.loc['DAY', 'info']:02.0f}"
        self.hour = f"{metadata.loc['HOUR', 'info']:02.0f}"
        self.minute = f"{metadata.loc['MINT', 'info']:02.0f}"
        if os.path.isfile(sonde_lookup_filepath):
            self.sonde_lookup_table = pd.read_csv(
                sonde_lookup_filepath, sep=" ", header=None, index_col=0, dtype="string"
            )
            self.sonde_lookup_table.columns = ["name"]
            self.sonde_lookup_table["name"] = self.sonde_lookup_table[
                "name"
            ].str.replace("_", " ")
            self.station_name = self.sonde_lookup_table.loc[self.station_id]["name"]
            self.plot_title = (
                f"{self.station_name} {self.station_id}\n"
                f"{self.year}-{self.month}-{self.day} {self.hour}{self.minute}Z"
            )
        else:
            self.plot_title = (
                f"{self.station_id} "
                f"{self.year}-{self.month}-{self.day} {self.hour}{self.minute}Z"
            )

    def plot_main_title(self, **kwargs):
        self.axes.annotate(
            self.plot_title,
            xy=(0.02, 0.92),
            xytext=(0, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            fontsize=20,
        )


class DorsetTitle(object):
    def __init__(self, metadata, axes):
        self.axes = axes
        self.year = f"{metadata.loc['YEAR', 'info']:.0f}"
        self.month = f"{metadata.loc['MONTH', 'info']:02.0f}"
        self.day = f"{metadata.loc['DAY', 'info']:02.0f}"
        self.hour = f"{metadata.loc['HOUR', 'info']:02.0f}"
        self.minute = f"{metadata.loc['MINT', 'info']:02.0f}"
        self.plot_title = (
            f"{self.year}-{self.month}-{self.day} {self.hour}{self.minute}Z"
        )

    def plot_main_title(self, **kwargs):
        self.axes.annotate(
            self.plot_title,
            xy=(0.02, 0.92),
            xytext=(0, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            fontsize=20,
        )


class Profile:
    """Generate a temperature profile on a tephigram"""

    def __init__(self, pressures, temperatures, axes):
        pressures, temperatures = np.asarray(pressures), np.asarray(temperatures)
        assert pressures.shape == temperatures.shape
        self.axes = axes
        self._transform = self.axes.tephi_transform + self.axes.transData
        self.pressures = pressures
        self.temperatures = temperatures
        _, self.thetas = convert_pressure_temperature_to_pressure_theta(
            self.pressures, self.temperatures
        )
        # self.line = None

    def plot(self, **kwargs):
        if "zorder" not in kwargs:
            kwargs["zorder"] = 10

        (self.line,) = self.axes.plot(
            self.temperatures, self.thetas, transform=self._transform, **kwargs
        )
        return self.line


class Barbs:
    """Generate wind barbs on tephigram"""

    def __init__(self, pressures, wind_speeds, wind_directions, axes):
        self.axes = axes
        self.pressures = pressures
        self.wind_speeds = wind_speeds
        self.wind_directions = wind_directions
        self.barbs = None
        self._gutter = None
        self._transform = self.axes.tephi_transform + self.axes.transData
        self._kwargs = None
        self._custom_kwargs = None
        self._custom = dict(
            color=["barbcolor", "color", "edgecolor", "facecolor"],
            linewidth=["lw", "linewidth"],
            linestyle=["ls", "linestyle"],
        )

    def _uv(self, magnitude, angle):
        u = magnitude * np.sin(np.deg2rad(angle))
        v = magnitude * np.cos(np.deg2rad(angle))
        return u, v

    def _make_barb(self, temperature, theta, speed, angle):
        u, v = self._uv(speed, angle)
        barb = plt.barbs(
            temperature, theta, u, v, transform=self._transform, **self._kwargs
        )
        return barb

    def _calculate_barb_positions(self, pressures):
        axesfrac_y_points = np.linspace(0, 1, 1001)
        axesfrac_x_points = np.ones_like(axesfrac_y_points) * self._gutter
        axes_frac_xy = np.column_stack((axesfrac_x_points, axesfrac_y_points))
        transAxes = self.axes.transLimits.inverted() + self.axes.tephi_inverse
        temperature_theta_points = transAxes.transform(axes_frac_xy)
        temperature, theta = (
            temperature_theta_points[:, 0],
            temperature_theta_points[:, 1],
        )
        _, pressure_points = convert_temperature_theta_to_temperature_pressure(
            temperature, theta
        )
        interp_func = interp1d(pressure_points, temperature, fill_value="extrapolate")
        temperatures = interp_func(pressures)
        _, thetas = convert_pressure_temperature_to_pressure_theta(
            pressures, temperatures
        )

        return temperatures, thetas

    def plot(self, **kwargs):
        self._gutter = kwargs.pop("gutter", 0.9)
        self._kwargs = dict(length=5, linewidth=0.2, zorder=10)
        self._kwargs.update(kwargs)
        self._custom_kwargs = dict(
            color=None, linewidth=1.0, zorder=self._kwargs["zorder"]
        )
        if hasattr(self.pressures, "__next__"):
            self.pressures = list(self.pressures)
        if hasattr(self.wind_speeds, "__next__"):
            self.wind_speeds = list(self.wind_speeds)
        if hasattr(self.wind_directions, "__next__"):
            self.wind_directions = list(self.wind_directions)
        self.pressures, self.wind_speeds, self.wind_directions = (
            np.asarray(self.pressures),
            np.asarray(self.wind_speeds),
            np.asarray(self.wind_directions),
        )

        temperatures, thetas = self._calculate_barb_positions(self.pressures)
        self._make_barb(temperatures, thetas, self.wind_speeds, self.wind_directions)


if __name__ == "__main__":
    tpg = Tephigram()
    plt.savefig("test.png", dpi=300, pad_inches=0)
