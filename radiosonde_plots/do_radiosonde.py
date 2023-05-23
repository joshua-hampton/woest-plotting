import polars as pl
import datetime as dt
import os
from plot_skewt import plot_skewt
import plot_tephigram_woest


def do_radiosondes(file_name, outdir):
    radiosonde_metadata = {}
    radiosonde_metadata["date"] = file_name.split("/")[-1].split("_")[1]
    radiosonde_metadata["time"] = file_name.split("/")[-1].split("_")[2].split(".")[0]

    with open(file_name, encoding="charmap") as f:
        line_skip = 0
        units_line_next = False
        while True:
            data = f.readline()
            if units_line_next:
                # data_units = [i.strip() for i in data.split("\t")]
                break
            elif data.startswith("Elapsed time"):
                units_line_next = True
            elif "System trademark and model" in data:
                radiosonde_metadata["model"] = data.split("System trademark and model")[
                    1
                ].strip()
                line_skip += 1
            elif "Station name" in data:
                radiosonde_metadata["location"] = "_".join(
                    data.split("Station name")[1].strip().split(" ")
                )
                line_skip += 1
            elif "Balloon release date and time" in data:
                start_time = data.split("Balloon release date and time")[1].strip()
                radiosonde_metadata["start_time_dt"] = dt.datetime.strptime(
                    start_time, "%Y-%m-%dT%H:%M:%S"
                )
                line_skip += 1
            else:
                line_skip += 1

    df = pl.read_csv(
        file_name,
        encoding="charmap",
        skip_rows=line_skip,
        separator="\t",
        skip_rows_after_header=1,
        ignore_errors=True,
    )

    for column_name in df.columns:
        df = df.rename({column_name: column_name.strip()})

    plot_skewt(
        df,
        radiosonde_metadata["date"],
        radiosonde_metadata["time"],
        radiosonde_metadata["location"],
        outdir,
    )

    df_small = df.select(
        [
            "Elapsed time",
            "HeightMSL",
            "RH",
            "Lat",
            "Lon",
            "P",
            "Temp",
            "Dir",
            "Speed",
        ]
    )

    df_small_renamed = df_small.rename(
        {
            "Elapsed time": "DataSrvTime",
            "HeightMSL": "Height",
            "RH": "Humidity",
            "Lat": "Latitude",
            "Lon": "Longitude",
            "P": "Pressure",
            "Temp": "Temperature",
            "Dir": "WindDir",
            "Speed": "WindSpeed",
        }
    )

    temp_k = pl.Series([i + 273.15 for i in list(df_small_renamed["Temperature"])])
    df_small_renamed.replace("Temperature", temp_k)

    csv_filename = (
        f"{radiosonde_metadata['model']}-{radiosonde_metadata['location']}-{radiosonde_metadata['date']}-"
        f"{radiosonde_metadata['time']}-radiosonde.csv"
    )
    df_small_renamed.write_csv(file=f"{os.path.split(file_name)[0]}/{csv_filename}")

    plot_tephigram_woest.plot_woest_tephigram(
        f"{os.path.split(file_name)[0]}/{csv_filename}",
        f"{outdir}/{radiosonde_metadata['date']}/tephigram_{radiosonde_metadata['date']}T{radiosonde_metadata['time']}",
        theta_w=True,
        sfc_parcel=True,
    )


if __name__ == "__main__":
    import sys

    file_name = sys.argv[1]
    outdir = sys.argv[2]
    do_radiosondes(file_name, outdir)
