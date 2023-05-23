import os
import pandas as pd

from radiosonde.calc_wetbulb import dewpoint, saturation_vapor_pressure


class DorsetRadiosonde(object):
    def __init__(self, infile):
        self.filename = os.path.split(infile)[1]
        self.df = pd.read_csv(infile)

        self.sonde_model = None
        self.location = None
        self.date_str = None
        self.time_str = None

    def calc_dewpoint(self):
        self.df["Dewpoint"] = dewpoint(
            self.df["Humidity"]
            / 100
            * saturation_vapor_pressure(self.df["Temperature"])
        )

    def get_metadata(self):
        long_filename = self.filename
        if self.filename[-4:] == ".csv":
            long_filename = self.filename[:-4]
        filename_segments = long_filename.split("-", 5)
        self.sonde_model = filename_segments[0]
        self.location = filename_segments[1]
        self.date_str = filename_segments[2]
        self.time_str = filename_segments[3]

        meta_df = pd.DataFrame(
            {
                "field": ["LOCATION", "YEAR", "MONTH", "DAY", "HOUR", "MINT"],
                "info": [
                    self.location,
                    int(self.date_str[0:4]),
                    int(self.date_str[4:6]),
                    int(self.date_str[6:8]),
                    int(self.time_str[0:2]),
                    int(self.time_str[2:4]),
                ],
            }
        )
        meta_df = meta_df.set_index("field")
        return meta_df

    def prune_data(self, prune_pressure_list):
        pruned_profile = pd.DataFrame()
        for prune_pressure in prune_pressure_list:
            idx = self.df["Pressure"].sub(prune_pressure).abs().argmin()
            pruned_profile = pruned_profile.append(self.df.iloc[idx], ignore_index=True)
        pruned_profile.drop_duplicates(subset=["Pressure"])
        return pruned_profile


if __name__ == "__main__":
    sonde = DorsetRadiosonde(
        "/Users/brianlo/Downloads/MW41-DURLSON-20191026-083151-SynchronizedSoundingData.csv"
    )
    sonde.calc_dewpoint()
    pass
