"""Read in data.

python src/read_in_data.py
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
from src.constants import SRC_PATH


print(pd.__version__)
output_dir = os.path.join(SRC_PATH, "output_init")

def get_param(name="hurr.in"):
    param_d = {}
    with open(os.path.join(output_dir, name), "r") as file:
        for line in file:
            el = [x for x in [x.strip(" ") for x in line.split("\t")] if x != ""]
            if len(el) >= 4:
                print(el)
                param_d[el[2]] = el[0]

    param_d["PDS"] = float(param_d["PDS"])
    param_d["CD"] = float(param_d["CD"]) * 1e-3
    param_d["CE"] = float(param_d["CE"]) * 1e-3
    param_d["FC"] = float(param_d["FC"]) * 1e-5
    param_d["PDEP"] = str(param_d["PDEP"])
    param_d["DISS"] = float(param_d["DISS"])
    param_d["RB"] = float(param_d["RB"])
    param_d["TIMEPL"] = float(param_d["TIMEPL"])
    param_d["ETIME"] = float(param_d["ETIME"])
    return param_d

def get_s(name="s.in"):
    param_d_loc = {}
    s_l = []
    with open(os.path.join(output_dir, name), "r") as file:
        for line in file:
            el = [x for x in [x.strip(" ").strip("\n") for x in line.split(" ")] if x != ""]
            if len(el) >= 3:
                print(el)
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


def get_tab_array(name="time", delim="\t"):
    out_lol = []
    with open(os.path.join(output_dir, name), "r") as file:
        for line in file:
            # print(line)
            el = [float(x) for x in [x.strip(" ").strip("\n") for x in line.split(delim)] if x != ""]
            out_lol.append(el)
    return np.array(out_lol)

def get_outputs(tim="05"):
    di_a = {}
    for key in ["radius", "rgraph", "zgraph", "time", "vmaxz", "vthov", "uthov", "vbhov", "ubhov", "thehov"]:
        di_a[key] = get_tab_array(name=key + ".out", delim=" ")
    ti_a = {}
    for key in ["ucon", "vcon", "wcon", "pcon", "tfcon", "qcon", "tescon", "tecon", "liqcon"]:
        ti_a[key] = get_tab_array(name=key + tim + ".out", delim=" ")

    return di_a, ti_a


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
    time = 5
    router = 4
    di_a, ti_a = get_outputs()
    print(di_a, ti_a)
    param_d = get_param()
    sin_npa, names, units = get_s(name="s.in")
    sout_npa, names, units = get_s(name="s.out")
    print(sin_npa)
    print(sout_npa)
    print(param_d)
    router = min([float(param_d["RB"]), router])
    dint = float(param_d["TIMEPL"])
    # print(dint)
    ulast =  (float(param_d["ETIME"]) // dint)  -1
    # print(dint, ulast)
    days = [x*dint for x in range(0,  int(ulast) + 1)]
    # print(days)
