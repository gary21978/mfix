from datetime import datetime
from pathlib import Path
from typing import List, Mapping, Optional, Tuple

import rtl2.params
from rtl2.test_util import git_commit


def read_historic_pressure(fname: Path) -> List[Mapping[str, float]]:
    """Read *.dat file with data from mfix-classic
    Returns:  list of dicts"""

    historic: List[Mapping[str, float]] = []
    with open(fname) as f:
        for line in f:
            t, pg, *_extra = line.split()
            historic.append({"t": float(t), "pg": float(pg)})
    return historic


def read_historic_velocity(fname: Path, scale: float) -> List[Mapping[str, float]]:
    """Read *.dat file with data from mfix-classic
    Returns:  list of dicts"""

    historic: List[Mapping[str, float]] = []
    with open(fname) as f:
        for line in f:
            t, pg, *_extra = line.split()
            historic.append({"t": float(t), "velocity": float(pg) * scale})
    return historic


def read_historic_np(fname: Path) -> List[Mapping[str, float]]:
    """Read *.dat file with data from mfix-classic
    Returns:  list of dicts"""

    historic: List[Mapping[str, float]] = []
    with open(fname) as f:
        for line in f:
            t, np, *_extra = line.split()
            historic.append({"t": float(t), "normalized_num_of_particles": float(np)})
    return historic


def read_two_pressures(pg1_fname: Path, pg2_fname: Path) -> List[Mapping[str, float]]:
    """Read pressure output of current run
    Returns: list of dict (values)"""

    data: List[Mapping[str, float]] = []
    with open(pg1_fname) as pg1_f, open(pg2_fname) as pg2_f:
        next(pg1_f)  # skip header
        next(pg2_f)  # skip header
        for line, line2 in zip(pg1_f, pg2_f):
            t1, pg1, vol1 = line.split()
            t2, pg2, vol2 = line2.split()
            assert t1 == t2
            data.append(
                {
                    "t": float(t1),
                    "pg1": float(pg1),
                    "vol1": float(vol1),
                    "pg2": float(pg2),
                    "vol2": float(vol2),
                }
            )
    return data


def read_avg_values(refdata_fname: Path) -> Tuple[List[int], List[float]]:
    """Returns: list of data points in refdata/runningavg.dat"""

    with open(refdata_fname) as run_avg:
        ys = []
        for line in run_avg:
            val, *_ = line.split()
            if val != "N/A":
                ys.append(float(val))
        xs = list(range(-len(ys), 0))
        return (xs, ys)


def avg_values_within_tolerance(refdata_fname: Path, tolerance: Optional[float]) -> bool:
    """Returns: whether the last data point is within tolerance of last 10 data points"""

    _, ys = read_avg_values(refdata_fname)
    if not ys:
        return True
    latest = ys[-1]
    recent = ys[-10:-1]
    if not len(recent):
        return True
    avg = sum(recent) / len(recent)
    if not bool(avg + latest):
        return True
    return abs((avg - latest) / (avg + latest)) < (0.1 if tolerance is None else tolerance)


def append_avg_dp_value(
    pg_data: List[Mapping[str, float]],
    vel_data: List[Mapping[str, float]],
    runningave_fname: Path,
    weight_over_area_per_particle: float,
) -> None:
    """Append new weighted average value for drop in pressure to end of refdata/runningavg.dat
    Returns: list of data points in (updated) refdata/runningavg.dat"""

    if rtl2.suite.get_suite().post_only:
        return

    WoA = woa(vel_data, weight_over_area_per_particle)
    dp_rate = dp_per_t(pg_data)
    avg_dp = "{:16.8}".format(dp_rate / WoA) if WoA is not None and dp_rate is not None else "N/A"
    branch, sha = git_commit(rtl2.suite.get_suite().sourceDir)
    time = datetime.now().astimezone().isoformat()
    assert runningave_fname.is_file()
    with open(runningave_fname, "a") as run_avg:
        run_avg.write(f"{avg_dp}\t{time}\t{sha}\t{branch}\n")


def read_velocity(vel_p_fname: Path) -> List[Mapping[str, float]]:
    """Read velocity file from current test run
    Returns:  Weighted Average value for current test run"""

    data: List[Mapping[str, float]] = []
    with open(vel_p_fname) as f:
        next(f)  # skip header
        for line in f:
            t, np, up, vp, wp, ke = line.split()
            data.append(
                {
                    "t": float(t),
                    "np": float(np),
                    "up": float(up),
                    "vp": float(vp),
                    "wp": float(wp),
                    "ke": float(ke),
                }
            )
    return data


def woa(data: List[Mapping[str, float]], weight_over_area_per_particle: float) -> Optional[float]:
    """Uses data from the velocity file for number of particles
    Returns:  Weighted Average value for current test run"""

    sumT = 0.0
    sumNp = 0.0
    WoA = None
    for ii in range(1, len(data) - 1):
        if data[ii]["t"] >= 0.5:
            dt = data[ii + 1]["t"] - data[ii]["t"]
            assert dt > 0
            sumNp += dt * 0.5 * (data[ii + 1]["np"] + data[ii]["np"])
            sumT += dt
            aveNp = sumNp / sumT
            WoA = weight_over_area_per_particle * aveNp
    return WoA


def dp_per_t(data: List[Mapping[str, float]]) -> Optional[float]:
    """Compute pressure from list of values"""
    T = 0.0
    total_dp = 0.0
    dp_rate = None
    for ii in range(0, len(data) - 1):
        if data[ii]["t"] >= 0.2:
            dt = data[ii + 1]["t"] - data[ii]["t"]
            assert dt > 0
            dpiL = 1.5 * data[ii]["pg1"] - 0.5 * data[ii]["pg2"]
            dpiR = 1.5 * data[ii + 1]["pg1"] - 0.5 * data[ii + 1]["pg2"]
            T += dt
            total_dp += dt * 0.5 * (dpiL + dpiR)
            dp_rate = total_dp / T
    return dp_rate


def get_mean_vp(vel_p_fname: Path) -> float:
    data = []
    with open(vel_p_fname) as f:
        next(f)  # skip header
        for line in f:
            t, np, up, vp, wp, ke = line.split()
            data.append(
                {
                    "t": float(t),
                    "np": float(np),
                    "up": float(up),
                    "vp": float(vp),
                    "wp": float(wp),
                    "ke": float(ke),
                }
            )

    mean_vp = 0.0
    sumT = 0.0
    sumV = 0.0
    for ii in range(1, len(data) - 1):
        if data[ii]["t"] >= 0.5:
            dt = data[ii + 1]["t"] - data[ii]["t"]
            assert dt > 0
            sumT += dt
            sumV += dt * 0.5 * (data[ii + 1]["up"] + data[ii]["up"])
            mean_vp = sumV / sumT
    return mean_vp
