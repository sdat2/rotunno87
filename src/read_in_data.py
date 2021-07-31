"""Read in data.

python src/read_in_data.py
"""
import os
from typing import Tuple
import numpy as np
import pandas as pd
import xarray as xr
from src.constants import SRC_PATH


class ReadData():
    def __init__(self, folder_name="output_init") -> None:
        """Initialse the data structure and read in the paramters.

        Args:
            folder_name (str, optional): [description]. Defaults to "output_init".
        """
        self.output_dir = os.path.join(SRC_PATH, folder_name)
        self.param = self.get_param()

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
                el = [x.strip("\n") for x in [x.strip(" ") for x in line.split("\t")] if x != ""]
                if len(el) >= 4:
                    print(el)
                    param_d[el[2]] = el[0]

        param_d["PDS"] = float(param_d["PDS"])
        param_d["CD"] = float(param_d["CD"]) * 1e-3
        param_d["CD1"] = float(param_d["CD1"]) * 1e-5
        param_d["CE"] = float(param_d["CE"]) * 1e-3
        param_d["FC"] = float(param_d["FC"]) * 1e-5
        param_d["PDEP"] = str(param_d["PDEP"])
        param_d["DISS"] = float(param_d["DISS"])
        param_d["RB"] = float(param_d["RB"])
        param_d["TIMEPL"] = float(param_d["TIMEPL"])
        param_d["ETIME"] = float(param_d["ETIME"])
        time = 5
        router = 4
        param_d["ROUTER"] = min([param_d["RB"], router])
        param_d["DINT"] = param_d["TIMEPL"]
        param_d["ULAST"] =  param_d["ETIME"] // param_d["DINT"] - 1
        # param_d["DAYS"] = [x*param_d["DINT"] for x in range(0,  int(param_d["ULAST"]) + 1)]
        param_d["TIME"] = time
        return param_d

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
                el = [x for x in [x.strip(" ").strip("\n") for x in line.split(" ")] if x != ""]
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
                el = [float(x) for x in [x.strip(" ").strip("\n") for x in line.split(delim)] if x != ""]
                out_lol.append(el)
        return np.array(out_lol)


    def get_outputs(self, tim: str="05") -> Tuple[dict, dict]:
        di_a = {}
        for key in ["radius", "rgraph", "zgraph", "time", "vmaxz", "vthov", "uthov", "vbhov", "ubhov", "thehov"]:
            di_a[key] = self.get_tab_array(name=key + ".out", delim=" ")
        ti_a = {}
        for key in ["ucon", "vcon", "wcon", "pcon", "tfcon", "qcon", "tescon", "tecon", "liqcon"]:
            ti_a[key] = self.get_tab_array(name=key + tim + ".out", delim=" ")

        return di_a, ti_a


    def tdsts(self) -> Tuple[xr.Dataset, xr.Dataset]:
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
                ds[key].attrs["long_name"] =  key
                def get_range(npa):
                    return min(npa), max(npa)

                if np.all((0, 265) == get_range(ds[key]["x"].values)):
                    da_out.append(ds[key])
                else:
                    da_out2.append(ds[key].rename({"x":"y", "y": "x"}))
            return xr.merge(da_out), xr.merge(da_out2)

        tim_l = ["0"+ str(i) for i in range(1, 10)]
        # tim_l[0] = ""
        ds1_l = []
        ds2_l = []
        for tim in tim_l:
            ti_a = {}
            for key in ["ucon", "vcon", "wcon", "pcon", "tfcon", "qcon", "tescon", "tecon", "liqcon"]:
                ti_a[key] = self.get_tab_array(name=key + tim + ".out", delim=" ")
            ds1, ds2 = make_datasets(ti_a)
            ds1_l.append(ds1)
            ds2_l.append(ds2)

        return xr.concat(ds1_l, "T"), xr.concat(ds2_l, "T")
        # return ds1_l, ds2_l


""""
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
