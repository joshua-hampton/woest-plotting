from plots.radiosonde.tephigram.tephigram_transforms import (
    convert_pressure_temperature_to_pressure_theta,
    convert_pressure_theta_to_temperature_theta,
    convert_pressure_mixing_ratio_to_temperature,
)


def isotherm_label(pressure, axes, transform, temperature):
    _, theta = convert_pressure_temperature_to_pressure_theta(pressure, temperature)
    annotation = axes.annotate(
        xy=(temperature, theta),
        xycoords=transform,
        xytext=(-1.0, -2.0),
        textcoords="offset points",
        text=f"{temperature}",
        ha="right",
        va="top",
        fontsize=7,
        color="#23CE1F",
    )

    return annotation


def isentrope_label(temperature, axes, transform, theta):
    annotation = axes.annotate(
        xy=(temperature, theta),
        xycoords=transform,
        xytext=(0, 0),
        textcoords="offset points",
        text=f"{theta}",
        ha="center",
        va="center",
        fontsize=7,
        color="#23CE1F",
    )

    return annotation


def isobar_label(along, isopleth_val, axes, transform, pressure):
    if along == "isotherm":
        _, theta = convert_pressure_temperature_to_pressure_theta(
            pressure, isopleth_val
        )
        annotation = axes.annotate(
            xy=(isopleth_val, theta),
            xycoords=transform,
            xytext=(0, 0),
            textcoords="offset points",
            text=f"{pressure}hPa",
            ha="right",
            va="bottom",
            fontsize=8,
            color="#23CE1F",
        )
    elif along == "isentrope":
        temperature, _ = convert_pressure_theta_to_temperature_theta(
            pressure, isopleth_val
        )
        annotation = axes.annotate(
            xy=(temperature, isopleth_val),
            xycoords=transform,
            xytext=(0, 0),
            textcoords="offset points",
            text=f"{pressure}hPa",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#23CE1F",
        )
    else:
        raise ValueError("Your along option does not exist!")

    return annotation


def mixing_ratio_label(pressure, horizontal_align, axes, transform, mixing_ratio):
    temperature = convert_pressure_mixing_ratio_to_temperature(pressure, mixing_ratio)
    _, theta = convert_pressure_temperature_to_pressure_theta(pressure, temperature)
    annotation = axes.annotate(
        xy=(temperature, theta),
        xycoords=transform,
        xytext=(0, 0),
        textcoords="offset points",
        text=f"{mixing_ratio:g}",
        ha=horizontal_align,
        va="center",
        fontsize=5,
        rotation=55,
        rotation_mode="anchor",
        annotation_clip=True,
        color="#23CE1F",
    )

    return annotation
