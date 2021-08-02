"""Read in data.

python src/read_in_data.py
"""
import os
from typing import Tuple
import numpy as np
import xarray as xr
import pint
from src.constants import SRC_PATH


class ReadData:
    """An object ot reat the data files, 
        and convert them to a more usual format."""

    def __init__(self, folder_name="output_init") -> None:
        """Initialse the data structure and read in the paramters.

        Args:
            folder_name (str, optional): Folder name. Defaults to "output_init".
        """
        self.output_dir = os.path.join(SRC_PATH, folder_name)
        self.ureg = pint.UnitRegistry()
        self.param = self.get_param()
        self.zgraph = self.get_tab_array(name="zgraph" + ".out", delim=" ")[0]
        self.rgraph = self.get_tab_array(name="rgraph" + ".out", delim=" ")[0]

    def get_param(self, name: str = "hurr.in") -> dict:
        """Get parameters.

        Args:
            name (str, optional): Parameter file. Defaults to "hurr.in".

        Returns:
            dict: dictionary of processed.
        """
        param_d = {}
        with open(os.path.join(self.output_dir, name), "r") as file:
            for line in file:
                line = line.replace("\t", "   ")
                el = [
                    x.strip("\n")
                    for x in [x.strip(" ") for x in line.split("    ")]
                    if x != ""
                ]
                if len(el) >= 4:
                    print(el)
                    param_d[el[2]] = el[0]

        param_d["PDS"] = float(param_d["PDS"]) * self.ureg.millibar
        param_d["CD"] = float(param_d["CD"]) * 1e-3
        param_d["CD1"] = float(param_d["CD1"]) * 1e-5 * self.ureg.second / self.ureg.meter
        param_d["CE"] = float(param_d["CE"]) * 1e-3
        param_d["FC"] = float(param_d["FC"]) * 1e-5 / self.ureg.seconds
        param_d["PDEP"] = str(param_d["PDEP"])
        param_d["DISS"] = float(param_d["DISS"])
        param_d["ZB"] = float(param_d["ZB"]) * self.ureg.kilometers
        param_d["RB"] = float(param_d["RB"]) * self.ureg.kilometers
        param_d["RO"] = float(param_d["RO"]) * self.ureg.kilometers
        param_d["VMAX"] = float(param_d["VMAX"]) * self.ureg.meter / self.ureg.seconds
        param_d["ROG"] = float(param_d["ROG"]) * self.ureg.kilometers
        param_d["RSST"] = float(param_d["RSST"]) * self.ureg.kilometers
        param_d["RADMAX"] = float(param_d["RADMAX"]) * self.ureg.delta_celsius / self.ureg.days
        # 'delta_degC/second'
        param_d["RMAX"] = float(param_d["RMAX"]) * self.ureg.kilometers
        param_d["ZOG"] = float(param_d["ZOG"]) * self.ureg.kilometers
        param_d["NS"] = float(param_d["NS"])
        param_d["XHL"] = float(param_d["XHL"]) * self.ureg.meter
        param_d["XVL"] = float(param_d["XVL"]) * self.ureg.meter
        param_d["DT"] = float(param_d["DT"]) * self.ureg.seconds
        param_d["PLTIME"] = float(param_d["PLTIME"]) * self.ureg.days
        param_d["TIMEPL"] = float(param_d["TIMEPL"]) * self.ureg.days
        param_d["ETIME"] = float(param_d["ETIME"]) * self.ureg.days
        time = 5 * self.ureg.days
        router = 4 * self.ureg.kilometers
        param_d["ROUTER"] = min([param_d["RB"], router])
        param_d["DINT"] = param_d["PLTIME"]
        param_d["ULAST"] = param_d["ETIME"] // param_d["DINT"] - 1
        # param_d["DAYS"] = [x*param_d["DINT"] for x in range(0,  int(param_d["ULAST"]) + 1)]
        param_d["TIME"] = time 
        param_d["TIMMAX"] = float(param_d["TIMMAX"]) * self.ureg.hours
        param_d["TAVE"] = float(param_d["TAVE"]) * self.ureg.days
        param_d["TMID"] = float(param_d["TMID"]) * self.ureg.kelvin
        param_d["TMIN"] = float(param_d["TMIN"]) * self.ureg.kelvin
        param_d["TAUR"] = float(param_d["TAUR"]) * self.ureg.hours
        param_d["VTERM"] = float(param_d["VTERM"]) * self.ureg.meter / self.ureg.second
        param_d["VTERMSNOW"] = float(param_d["VTERMSNOW"]) * self.ureg.meter / self.ureg.second
        param_d["VCAP"] = float(param_d["VCAP"]) * self.ureg.meter / self.ureg.second
        return param_d

    def get_days(self) -> np.array:
        """Get days list.

        Returns:
            np.array: numpy array.
        """
        return np.array(
            [x * self.param["DINT"] for x in range(0, int(self.param["ULAST"]) + 1)]
        )

    def get_s(self, name: str = "s.in") -> xr.Dataset:
        """Get the s profiles as a dataset.

        Args:
            name (str, optional): [description]. Defaults to "s.in".

        Returns:
            xr.Dataset: [description]
        """
        param_d_loc = {}
        s_l = []
        with open(os.path.join(self.output_dir, name), "r") as file:
            for line in file:
                el = [
                    x
                    for x in [x.strip(" ").strip("\n") for x in line.split(" ")]
                    if x != ""
                ]
                if len(el) >= 3:
                    param_d_loc[el[0]] = [el[1], el[2]]
                    s_l.append(el)
        units = []
        for i in [s_l[1][1], s_l[1][3], s_l[1][5]]:
            units.append(i.strip("(").strip(")"))
        names = []
        for i in [s_l[1][0], s_l[1][2], s_l[1][4]]:
            names.append(i)
        s_l.pop(0)
        s_l.pop(0)
        s_l.pop(0)
        s_npa = np.array(s_l)
        ds = xr.Dataset(
            data_vars={
                names[1]: (names[0], s_npa[:, 1].astype("float32")),
                names[2]: (names[0], s_npa[:, 2].astype("float32")),
            },
            coords={names[0]: s_npa[:, 0].astype("float32")},
            attrs={"description": name},
        )
        for i in range(len(names)):
            ds[names[i]].attrs["units"] = units[i]
        return ds

    def get_tab_array(self, name="time", delim="\t") -> np.ndarray:
        """Get the arrays at a particular level.

        Args:
            name (str, optional): Name. Defaults to "time".
            delim (str, optional): Delimiter. Defaults to "\t".

        Returns:
            np.ndarray: numpy.
        """
        out_lol = []
        with open(os.path.join(self.output_dir, name), "r") as file:
            for line in file:
                # print(line)
                el = [
                    float(x)
                    for x in [x.strip(" ").strip("\n") for x in line.split(delim)]
                    if x != ""
                ]
                out_lol.append(el)
        return np.array(out_lol)

    def get_outputs(self, tim: str = "05") -> Tuple[dict, dict]:
        """Get the outputs.

        Args:
            tim (str, optional): Get all the outputs. Defaults to "05".

        Returns:
            Tuple[dict, dict]: Dictionary.
        """
        di_a = {}
        for key in [
            "radius",
            "rgraph",
            "zgraph",
            "time",
            "vmaxz",
            "vthov",
            "uthov",
            "vbhov",
            "ubhov",
            "thehov",
            "xkcon",
        ]:
            di_a[key] = self.get_tab_array(name=key + ".out", delim=" ")
        ti_a = {}
        for key in [
            "ucon",
            "vcon",
            "wcon",
            "pcon",
            "tcon",
            "tfcon",
            "tescon",
            "tecon",
            "qcon",
            "liqcon",
        ]:
            ti_a[key] = self.get_tab_array(name=key + tim + ".out", delim=" ")

        return di_a, ti_a

    def get_time_datasets(self) -> Tuple[xr.Dataset, xr.Dataset]:
        """Tdsts.

        Returns:
            Tuple[xr.Dataset, xr.Dataset]: xarray.
        """

        def make_datasets(d2) -> Tuple[xr.Dataset, xr.Dataset]:
            da_out = []
            da_out2 = []
            for key in d2:
                ds = xr.Dataset(
                    data_vars={
                        key: (["y", "x"], d2[key].astype("float32")),
                    },
                )
                ds[key].attrs["long_name"] = key

                def get_range(npa):
                    return min(npa), max(npa)

                if np.all((0, 265) == get_range(ds[key]["x"].values)):
                    da = ds[key]
                    da = da.assign_coords({"x": ("x", self.rgraph)})
                    da = da.assign_coords({"y": ("y", self.zgraph)})
                    da.y.attrs["units"] = "km"
                    da.x.attrs["long_name"] = "r"
                    da.x.attrs["units"] = "km"
                    da.y.attrs["long_name"] = "z"
                    da_out.append(da)

                else:
                    da_out2.append(ds[key].rename({"x": "y", "y": "x"}))
            return xr.merge(da_out), xr.merge(da_out2)

        tim_l = ["0" + str(i) for i in range(0, 10)]
        tim_l[0] = ""
        ds1_l = []
        ds2_l = []
        for tim in tim_l:
            ti_a = {}
            for key in [
                "ucon",
                "vcon",
                "wcon",
                #"pcon",
                #"tfcon",
                # "qcon",
                #"tescon",
                #"tecon",
                #"liqcon",
            ]:
                # qcon not at 0
                ti_a[key] = self.get_tab_array(name=key + tim + ".out", delim=" ")
            ds1, ds2 = make_datasets(ti_a)
            ds1_l.append(ds1)
            ds2_l.append(ds2)

        ds1 = xr.concat(ds1_l, "T")
        ds2 = xr.concat(ds2_l, "T")
        times = [float(self.param["PLTIME"].magnitude) * (x) for x in range(ds1.sizes["T"])]
        ds1 = ds1.assign_coords({"T": ("T", times)})
        ds1["T"].attrs["units"] = "days"
        ds2 = ds2.assign_coords({"T": ("T", times)})
        ds2["T"].attrs["units"] = "days"
        return ds1, ds2
        # return ds1_l, ds2_l


"""
print(get_tab_array(name="radius" + ".out", delim=" "))
print(get_tab_array(name="rgraph" + ".out", delim=" "))
print(get_tab_array(name="zgraph" + ".out", delim=" "))
print(get_tab_array(name="ucon" + tim + ".out", delim=" "))
print(get_tab_array(name="vcon" + tim + ".out", delim=" "))
print(get_tab_array(name="wcon" + tim + ".out", delim=" "))
print(get_tab_array(name="pcon" + tim + ".out", delim=" "))
print(get_tab_array(name="liqcon" + tim + ".out", delim=" "))
"""

if __name__ == "__main__":
    # python src/read_in_data.py
    print("ok")
    read = ReadData()

"""
    di_a, ti_a = get_outputs()
    print(di_a, ti_a)
    param_d = get_param()
    sin_ds = get_s(name="s.in")
    sout_ds = get_s(name="s.out")
    print(sin_npa)
    print(sout_npa)
    print(param_d)
    time = 5
    router = 4
    router = min([float(param_d["RB"]), router])
    dint = float(param_d["TIMEPL"])
    # print(dint)
    ulast =  (float(param_d["ETIME"]) // dint)  -1
    # print(dint, ulast)
    days = [x*dint for x in range(0,  int(ulast) + 1)]
    # print(days)
"""
