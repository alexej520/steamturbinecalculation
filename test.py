#!/usr/bin/python3
# -*- coding: utf-8 -*-

# import pandas as pd
from os import system
import numpy as np
from sympy import Point, Line, Segment
from iapws import IAPWS97 as IAPWS97_ORIG
import matplotlib.pyplot as plt
import matplotlib
import plots as diagram
from math import log10
from sympy import S, solve, symbols


def intersect_segment_polyline(segment, x, y):
    c_s = 0  # current segment
    m_s = min(len(x), len(y)) - 1  # last segment
    point = []
    while len(point) < 1 and c_s < m_s:
        point = segment.intersection(
            Segment(Point(x[c_s], y[c_s]), Point(x[c_s + 1], y[c_s + 1])))
        c_s += 1
    return point


def plot_segment(*args, **kwargs):
    point1 = args[0]
    point2 = args[1]
    args = ((point1.x, point2.x), (point1.y, point2.y)) + args[2:len(args)]
    plt.plot(*args, **kwargs)
    return


def plot_point(*args, **kwargs):
    point = args[0]
    args = (point.x, point.y) + args[1:len(args)]
    plt.plot(*args, **kwargs)
    return


class IAPWS97(IAPWS97_ORIG):

    def __init__(self, **kwargs):
        if "t" in kwargs:
            kwargs["T"] = kwargs["t"] + 273.15
        if "p" in kwargs:
            kwargs["P"] = kwargs["p"] * 0.1
        IAPWS97_ORIG.__init__(self, **kwargs)

    @property
    def t(self):
        return self.T - 273.15

    @property
    def p(self):
        return self.P * 10


class MyData:

    def __init__(self, attr_list):
        for attr in attr_list:
            setattr(self, attr[0].decode("UTF-8"), attr[1])


vc = 3  # Value column
sc = 1  # Symbol column
# df = pd.read_csv(
#     "source.csv",
#     names=["Symbol", "Value"],
#     usecols=[sc, vc],
#     index_col=0,
#     skiprows=[0, 1, 2, 3])
# _ = df.Value
# ndt = pd.read_csv("source.csv", nrows=2, index_col=sc)

_ = MyData(np.genfromtxt(
    "source.csv",
    delimiter=",",
    usecols=[sc, vc],
    dtype=None,
    skip_header=4))
Temps = np.genfromtxt(
    "source.csv",
    delimiter=",",
    skip_header=1,
    max_rows=1)[vc:].tolist()
Time = np.genfromtxt(
    "source.csv",
    delimiter=",",
    skip_header=2,
    max_rows=1)[vc:].tolist()

# ########## Построение графиков нагрузки ##########
# Style
font = {
    "family": "DejaVu Sans",
    "weight": "normal"
}
matplotlib.rc("font", **font)
matplotlib.rcParams["mathtext.default"] = "regular"  # "it"

t_os = 70
t_ps_ = 150
p_spm = 0.9 * _.p_Tm
t_s = IAPWS97(p=p_spm, x=0).t
delta_tsp = 5  # 5..7
t_psm = t_s - delta_tsp
delta_ts = t_psm - t_os
alpha_TPS_calc = (t_psm - t_os) / (t_ps_ - t_os)
Q_TTPS = _.Q_Tm / alpha_TPS_calc
Q_y = 0.15 * Q_TTPS
# Temps = ndt.T.Temp[vc - 1:].tolist()
# Time = ndt.T.Time[vc - 1:].tolist()
Temps.append(8)
Time.append(365 - sum(Time))
Hours = [0]
Hours[0] = Time[0] * 24

for i in range(1, len(Time)):
    Hours.insert(i, Hours[i - 1] + Time[i] * 24)

Q_T = [Q_y + (18 - t) / (18 - _.t_calc) * (Q_TTPS - Q_y) for t in Temps]
Q_T[-1] = Q_y
Q_T.append(Q_y)

point_Te = Point(18, 18)  # Thermodynamic equilibrium
line_Q_Tm = Line(Point(20, _.Q_Tm), Point(Temps[0], _.Q_Tm))
line_2 = Line(Point(18, Q_y), Point(Temps[1], Q_TTPS))
line_3 = Line(point_Te, Point(Temps[1], 150))
line_4 = Line(point_Te, Point(Temps[1], t_os))
point_M = line_Q_Tm.intersection(line_2)[0]
point_M1 = line_3.intersection(
    Line(Point(point_M.x, 0), Point(point_M.x, 150)))[0]
point_M11 = line_4.intersection(
    Line(Point(point_M.x, 0), Point(point_M.x, 150)))[0]
point_2 = line_2.intersection(Line(Point(8, 0), Point(8, Q_y)))[0]
delta_t = point_M1.distance(point_M11)
t_ps = t_os + delta_t

if t_psm > t_ps:
    t_psm = t_ps
    p_T = IAPWS97(t=t_psm + delta_tsp, x=0).p / 0.9
    Q = _.Q_Tm
    point_N1 = Point(Temps[0], t_ps)
else:
    p_T = _.p_Tm
    Q = _.Q_Tm * delta_ts / delta_t
    point_N1 = Line(Point(20, t_psm), Point(Temps[0], t_psm)).intersection(
        Line(point_M1, Point(Temps[0], t_ps)))[0]

point_N = Point(point_N1.x, _.Q_Tm)
point_A = line_3.intersection(Line(Point(20, t_os), Point(Temps[0], t_os)))[0]
point_B = line_4.intersection(Line(point_A, Point(point_A.x, 0)))[0]

# plot graphs
solid = {"linestyle": "-", "marker": "o", "color": "black", "linewidth": 2.0}
solid_ = {"linestyle": "-", "color": "black", "linewidth": 2.0}
dashed = {"linestyle": "--", "color": "black"}

fig = plt.figure("Графики тепловых нагрузок")

plt.subplot(2, 2, 4)
plt.plot(Hours, Temps, **solid_)  # 1
ca = plt.gca()
ca.set_xlim(0, 8760)
ca.set_ylim(Temps[0], 20)
ca.invert_yaxis()
ca.set_xlabel(r"$\tau$, ч/год")
ca.set_ylabel(r"$t_{вз}$, $^\circ$C")
ca.grid(color="gray", linestyle="-")

plt.subplot(2, 2, 2)
Hours[-1] = Hours[-2]
Hours.append(8760)
point_M2 = intersect_segment_polyline(
    Segment(Point(0, _.Q_Tm), Point(8760, _.Q_Tm)), Hours, Q_T)[0]
point_N2 = Point(intersect_segment_polyline(
    Segment(Point(0, point_N1.x), Point(8760, point_N1.x)),
    Hours, Temps)[0].x, _.Q_Tm)
plt.plot(Hours, Q_T, **solid_)  # 5
plot_segment(point_M2, point_N2, **solid)  # M2 N2
plot_segment(point_N2, Point(0, Q), **solid)  # -> N2
plot_segment(Point(Hours[-2], Q_T[-3]), Point(Hours[-1], Q_T[-3]), **dashed)
ca = plt.gca()
ca.set_xlim(0, 8760)
ca.set_ylim(0, Q_TTPS)
ca.set_xlabel(r"$\tau$, ч/год")
ca.set_ylabel(r"$Q_T$, МВт")
ca.grid(color="gray", linestyle="-")

plt.subplot(2, 2, 3)
plot_segment(point_Te, point_B, **dashed)  # -> B
plot_segment(point_Te, point_A, **dashed)  # -> A
plot_segment(point_B, Point(Temps[0], t_os), **solid)  # B t_os
plot_segment(point_A, point_M1, **solid)  # A M1
plot_segment(point_M1, Point(Temps[0], 150), **dashed)  # M1 ->
plot_segment(point_M1, Point(Temps[0], t_ps), **solid)  # M1 N1
plot_segment(point_N1, Point(Temps[0], t_psm), **solid)  # N1 t_psm
plot_segment(Point(20, t_os), point_A, **solid)  # -> A
plot_segment(Point(20, point_B.y), point_B, **solid)  # -> B
ca = plt.gca()
ca.set_xlim(Temps[0], 20)
ca.set_ylim(0, 160)
ca.invert_xaxis()
ca.set_xlabel(r"$t_{вз}$, $^\circ$C")
ca.set_ylabel(r"$t_{ПС}$, $t_{ОС}$, $^\circ$C")
ca.grid(color="gray", linestyle="-")

plt.subplot(2, 2, 1)
plot_segment(Point(20, _.Q_Tm), point_M, **dashed)
plot_segment(Point(18, Q_y), point_2, **dashed)  # -> 2
plot_segment(Point(8, Q_y), point_2, **solid)  # -> 2
plot_segment(Point(8, Q_y), Point(Temps[0], Q_y), **dashed)  #
plot_segment(point_M, point_N, **solid)  # M N
plot_segment(point_2, Point(Temps[0], Q_TTPS), **solid)  # 2 M
plot_segment(Point(20, Q_y), Point(8, Q_y), **solid)  #
plot_segment(point_N, Point(Temps[0], Q), **solid)  # N Q
ca = plt.gca()
ca.set_xlim(Temps[0], 20)
ca.set_ylim(0, Q_TTPS)
ca.invert_xaxis()
ca.set_xlabel(r"$t_{вз}$, $^\circ$C")
ca.set_ylabel(r"$Q_T$, МВт")
ca.grid(color="gray", linestyle="-")

# ########## Построение диаграмм ##########

_0_ = IAPWS97(p=_.p_0, t=_.t_0)  # . 0
_0 = IAPWS97(p=_0_.p * 0.95, h=_0_.h)  # . 0*
withSH = not np.isnan(_.t_SH)

_SH = IAPWS97(p=_0.p / 6 if np.isnan(_.p_SH) else _.p_SH, t=_.t_SH)  # . SH Superheater
p_HP = 1.1 * _SH.p  # pressure after HP

_1_s = IAPWS97(p=p_HP, s=_0.s)  # . 1_s
eff_oi0 = 0.85  # for HP
alpha = 2.6e-4  # for HP
eff_oiHP = eff_oi0 - alpha / _0.v
h_1 = _0.h - (_0.h - _1_s.h) * eff_oiHP
_1 = IAPWS97(p=p_HP, h=h_1)  # . 1

_2_s = IAPWS97(p=p_T, s=_SH.s)  # . 2_s
eff_oi0 = 0.82  # for IP
alpha = 1.8e-3  # for IP
eff_oiIP = eff_oi0 - alpha / _SH.v
h_2 = _SH.h - (_SH.h - _2_s.h) * eff_oiIP
_2 = IAPWS97(p=p_T, h=h_2)  # . 2

_C_s = IAPWS97(p=_.p_C, s=_2.s)  # . C_s Condenser s
eff_oi0 = 0.8  # for LP
alpha = 1.8e-3  # for LP
eff_oiLP = eff_oi0 - alpha / _2.v
h_C = _2.h - (_2.h - _C_s.h) * eff_oiLP
_C = IAPWS97(p=_.p_C, h=h_C)  # . C

_C_ = IAPWS97(p=_.p_C, x=0)  # . C'
_D = IAPWS97(p=_.p_D_, x=0)  # . D Deaerator
p_D = 1.2 * (_.p_D_ + 2)  # отбор на деаэратор
_FP_ = IAPWS97(p=1.35 * _0.p, h=_D.h)
# . FP' Feedwater Pump pressure 1.3..1.4 * _.p_0
eff_FP = 0.8  # Feedwater Pump efficiency
delta_h_FP = _FP_.P * 1e-3 * 1e3 / eff_FP  # p_fp * v / eff_fp
_FP = IAPWS97(p=_FP_.p, h=_D.h + delta_h_FP)  # . FP
# steam generator pressure p_0 + 10..15 bar  p_SG t_SG have standard
_SG = IAPWS97(P=13.8 if _0.P + 1 <= 13.8 else 25, t=_0.t + 5)
_CP = IAPWS97(p=11, h=_C_.h)  # . CP Condensate pump 11 bar by default


def to_sp_point(iapws97_point):  # p in MPa
    return Point(iapws97_point.s, iapws97_point.P)


def to_sh_point(iapws97_point):
    return Point(iapws97_point.s, iapws97_point.h)


solid = {"linestyle": "-", "marker": "o", "color": "black", "linewidth": 2.0}
solid_ = {"linestyle": "-", "color": "black", "linewidth": 1.5}
dashed = {"linestyle": "--", "color": "black", "linewidth": 1.5}

fig1 = plt.figure("Цикл энергоблока с промежуточным перегревом пара")
diagram.plot(x="s", y="P", show=False)
plot_segment(to_sp_point(_0_), to_sp_point(_0), **solid)  # 0' 0
plot_segment(to_sp_point(_SH), to_sp_point(_2), **solid)  # SH 2
plot_segment(to_sp_point(_SH), to_sp_point(_2_s), **dashed)  # SH 2_s
plot_segment(to_sp_point(_2), to_sp_point(_C), **solid)  # 2 C
plot_segment(to_sp_point(_2), to_sp_point(_C_s), **dashed)  # 2 C_s
plot_segment(to_sp_point(_C), to_sp_point(_C_), **solid)  # C C'
plot_segment(to_sp_point(_D), to_sp_point(_FP_), **dashed)  # D FP'
plot_segment(to_sp_point(_D), to_sp_point(_FP), **solid)  # D FP
plot_segment(to_sp_point(_FP_), to_sp_point(_SG), **solid)  # FP' SG
plot_segment(to_sp_point(_SG), to_sp_point(_0_), **solid)  # FP' SG
plot_segment(to_sp_point(_C_), to_sp_point(_CP), **solid)  # C' CP
plot_segment(to_sp_point(_CP), to_sp_point(_D), **solid)  # CP D
plot_segment(to_sp_point(_0), to_sp_point(_1), **solid)  # 0 1
plot_segment(to_sp_point(_0), to_sp_point(_1_s), **dashed)  # 0 1_s
plot_segment(to_sp_point(_1), to_sp_point(_SH), **solid)  # 1 SH

# ## Распределение подогрева воды по ступеням ##

# ## Сетевые подгреватели ##

t_ps = t_psm
t_os = t_os
t_sv = (t_ps + t_os) * 0.5

# равномерное распределение по ступеням
n_LP = int(_.n_LP)
n_HP = int(_.n_HP)
delta_t_FP = _FP.p * _FP.v / (eff_FP * _FP.cp)
delta_t_HP = (_.t_FW - _D.t - delta_t_FP) / n_HP
delta_t_D = 20  # >=20
delta_t_OE_PU = 5  # 4..5 подогрев в сальниковом подогревателе
t_LP0 = _C_.t + delta_t_OE_PU
delta_t_LP = (t_sv - t_LP0) / 2

dt_LP = 4  # недогрев
dt_HP = 2
t_h = []  # heater temps
t_d = []  # drainage temps

t_h.append(t_LP0 + delta_t_LP)  # 1
t_d.append(t_h[-1] + dt_LP)
t_h.append(t_sv)  # 2
t_d.append(t_h[-1] + delta_tsp)
t_h.append(t_ps)  # 3
t_d.append(t_h[-1] + delta_tsp)
t_LP0 = t_ps
delta_t_LP = (_D.t - delta_t_D - t_ps) / (n_LP - 3)
for i in np.arange(1, n_LP + 1 - 3):
    t_h.append(t_LP0 + i * delta_t_LP)
    t_d.append(t_h[-1] + dt_LP)

t_h.append(_D.t)
t_d.append(_D.t)

for i in np.arange(n_HP):
    t_h.append(_D.t + i * delta_t_HP + delta_t_HP)
    t_d.append(t_h[-1] + dt_HP)

_h = [IAPWS97(x=0, t=t) for t in t_h]  # water
_d = [IAPWS97(x=0, t=t) for t in t_d]  # drainage
_s = []  # steam

P_s = [1.1 * d.P for d in _d]  # MPa
P_s[n_LP] = p_D * 0.1  # bar to MPa


def intersection_log10(line, segments, **kwargs):
    res_intrsn = []
    log_axis = kwargs["log"] if "log" in kwargs else "y"
    scratch = kwargs["scratch"] if "scratch" in kwargs else "x"
    for segment in segments:
        point0 = segment.points[0]
        point1 = segment.points[1]
        maxy = max(point0.y, point1.y)
        miny = min(point0.y, point1.y)
        maxx = max(point0.x, point1.x)
        minx = min(point0.x, point1.x)
        intrsn = line.intersection(segment)
        if log_axis is "y" and scratch is "x":
            for intrsn_p in intrsn:
                pr = log10(intrsn_p.y / miny) / log10(maxy / miny)
                res_intrsn.append(
                    Point(maxx - pr * (maxx - minx), intrsn_p.y))
        if log_axis is "y" and scratch is "y":
            for intrsn_p in intrsn:
                pr = (intrsn_p.x - minx) / (maxx - minx)
                res_intrsn.append(
                    Point(intrsn_p.x, maxy / 10 ** (pr * log10(maxy / miny))))
    return res_intrsn


_2_C = Segment(to_sp_point(_2), to_sp_point(_C))
_SH_2 = Segment(to_sp_point(_SH), to_sp_point(_2))
_0_1 = Segment(to_sp_point(_0), to_sp_point(_1))
for i in np.arange(len(_d)):
    if i == len(_d) - 2:
        plot_segment(to_sp_point(_d[i]), to_sp_point(_1), **solid_)
        _s.append(_1)
        continue
    p_line = Line(Point(0, P_s[i]), Point(1, P_s[i]))
    its = intersection_log10(p_line, [_2_C, _SH_2, _0_1], log="y", scratch="x")
    if len(its) > 0:
        plot_segment(to_sp_point(_d[i]), its[0], **solid_)
        plot_point(its[0], **solid)
        _s.append(IAPWS97(s=its[0].x, P=its[0].y))
    else:
        print("cannot draw heater %d" % (len(_d) - i))

_CP_D = Segment(to_sp_point(_CP), to_sp_point(_D))
_FP_SG = Segment(to_sp_point(_FP), to_sp_point(_SG))
for i in np.arange(len(_h)):
    p_line = Line(Point(_h[i].s, 0), Point(_h[i].s, 1))
    its = intersection_log10(p_line, [_CP_D, _FP_SG], log="y", scratch="y")
    if len(its) > 0:
        plot_segment(to_sp_point(_h[i]), its[0], **solid_)
        plot_point(its[0], **solid)
    else:
        print("cannot draw heater %d" % (len(_h) - i))

fig2 = plt.figure("Процесс расширения пара в турбине")
diagram.plot(x="s", y="h", show=False)
plot_segment(to_sh_point(_0_), to_sh_point(_0), **solid)  # 0' 0
plot_segment(to_sh_point(_SH), to_sh_point(_2), **solid)  # SH 2
plot_segment(to_sh_point(_SH), to_sh_point(_2_s), **dashed)  # SH 2_s
plot_segment(to_sh_point(_2), to_sh_point(_C), **solid)  # 2 C
plot_segment(to_sh_point(_2), to_sh_point(_C_s), **dashed)  # 2 C_s
plot_segment(to_sh_point(_C), to_sh_point(_C_), **solid)  # C C'
plot_segment(to_sh_point(_D), to_sh_point(_FP_), **dashed)  # D FP'
plot_segment(to_sh_point(_D), to_sh_point(_FP), **solid)  # D FP
plot_segment(to_sh_point(_FP_), to_sh_point(_SG), **solid)  # FP' SG
plot_segment(to_sh_point(_SG), to_sh_point(_0_), **solid)  # FP' SG
plot_segment(to_sh_point(_C_), to_sh_point(_CP), **solid)  # C' CP
plot_segment(to_sh_point(_CP), to_sh_point(_D), **solid)  # CP D
plot_segment(to_sh_point(_0), to_sh_point(_1), **solid)  # 0 1
plot_segment(to_sh_point(_0), to_sh_point(_1_s), **dashed)  # 0 1_s
plot_segment(to_sh_point(_1), to_sh_point(_SH), **solid)  # 1 SH

# ########## Расчет тепловой схемы ##########

eff = S(0.98)  # КПД теплообменника
a_1, a_2, a_w2, a_3, a_w3, a_4, a_w4, a_5, a_w5, a_D, a_6, a_7, a_8 = symbols(
    "a_1, a_2, a_w2, a_3, a_w3, a_4, a_w4, a_5, a_w5, a_D, a_6, a_7, a_8")
h_mp2, h_mp3, h_mp4 = symbols(
    "h_mp2, h_mp3, h_mp4")

h_1, h_2, h_3, h_4, h_5, h_D, h_6, h_7, h_8 = [S(s.h) for s in _s]
h_d1, h_d2, h_d3, h_d4, h_d5, h_dD, h_d6, h_d7, h_d8 = [S(d.h) for d in _d]
h_w1, h_w2, h_w3, h_w4, h_w5, h_wD, h_w6, h_w7, h_w8 = [S(w.h) for w in _h]
h_w0 = S(IAPWS97(x=0, t=_C_.t + delta_t_OE_PU).h)

# sl - steam leakage # s - seals # p - purge (только с барабаном)
a_sl, a_s, a_p = (0.02, 0.036, 0)
a_fw = (1 + a_s + a_sl + a_p)
K_r = 1.2  # 1.15..1.25
eff_em = 0.98
H = S((_0.h - _1.h) + (_SH.h - _C.h))
N_e = S(_.N * 1000)  # MW -> kW

_d.append(IAPWS97(x=0, t=t_sv + delta_tsp))
h_dT1 = S(_d[-1].h)
_d.append(IAPWS97(x=0, t=t_ps + delta_tsp))
h_dT2 = S(_d[-1].h)
h_os = S(IAPWS97(x=0, t=t_os).h)
_h.append(IAPWS97(x=0, t=t_sv))
h_sv = S(_h[-1].h)
_h.append(IAPWS97(x=0, t=t_ps))
h_ps = S(_h[-1].h)
_s.append(_s[1])
h_T1 = h_2
_s.append(_s[2])
h_T2 = h_3
y_T1 = (h_T1 - _C.h) / H  # (3.6) T - теплофикационный
y_T2 = (h_T2 - _C.h) / H  # (3.6)
D_T1 = Q * 1e3 / (2 * (h_T1 - h_dT1) * eff)  # (3.7) MW -> kW
D_T2 = Q * 1e3 / (2 * (h_T2 - h_dT2) * eff)  # (3.7) MW -> kW
D = K_r * (N_e / (H * eff_em) + y_T1 * D_T1 + y_T2 * D_T2)
a_os = S(Q * 1e3 / ((h_ps - h_os) * D))  # ? тепловая нагрузка MW -> kW
a_T1 = S(D_T1 / D)
a_T2 = S(D_T2 / D)

f_8 = a_8 * (h_8 - h_d8) * eff - (
    a_fw * (h_w8 - h_w7))
f_7 = a_7 * (h_7 - h_d7) * eff + a_8 * (h_d8 - h_d7) * eff - (
    a_fw * (h_w7 - h_w6))
f_6 = a_6 * (h_6 - h_d6) * eff + (a_8 + a_7) * (h_d7 - h_d6) * eff - (
    a_fw * (h_w6 - (h_wD + delta_h_FP)))  # h_wD с учетом подогрева в FP
fm_D = a_D + (a_8 + a_7 + a_6) + a_w5 - (
    a_fw)
f_D = a_D * h_D + (a_8 + a_7 + a_6) * h_d6 + a_w5 * h_w5 - (
    a_fw * h_wD)

s = solve([f_8, f_7, f_6, fm_D, f_D])
a_8 = s[a_8]
a_7 = s[a_7]
a_6 = s[a_6]
a_w5 = s[a_w5]
a_D = s[a_D]

f_5 = a_5 * (h_5 - h_d5) * eff - (
    a_w5 * (h_w5 - h_mp4))
fm_mp4 = a_w4 + (a_5 + a_4) - (
    a_w5)
f_mp4 = a_w4 * h_w4 + (a_5 + a_4) * h_d4 - (
    a_w5 * h_mp4)
f_4 = a_4 * (h_4 - h_d4) * eff + a_5 * (h_d5 - h_d4) * eff - (
    a_w4 * (h_w4 - h_mp3))
fm_mp3 = a_w3 + a_3 + a_T2 - (
    a_w4)
f_mp3 = a_w3 * h_w3 + a_3 * h_d3 + a_T2 * h_dT2 - (
    a_w4 * h_mp3)
f_3 = a_3 * (h_3 - h_d3) * eff - (
    a_w3 * (h_w3 - h_mp2))
f_T2 = a_T2 * (h_T2 - h_dT2) * eff - (
    a_os * (h_ps - h_sv))
fm_mp2 = a_w2 + a_2 + a_T1 - (
    a_w3)
f_mp2 = a_w2 * h_w2 + a_2 * h_d2 + a_T1 * h_dT1 - (
    a_w3 * h_mp2)
f_2 = a_2 * (h_2 - h_d2) * eff - (
    a_w2 * (h_w2 - h_w1))
f_T1 = a_T1 * (h_T1 - h_dT1) * eff - (
    a_os * (h_sv - h_os))
f_1 = a_1 * (h_1 - h_d1) * eff - (
    a_w2 * (h_w1 - h_w0))

# solve([f_5, fm_mp4, f_mp4, f_4, fm_mp3, f_mp3, f_3, fm_mp2, f_mp2, f_2, f_1])

def print_iapws97(l, label):
    print("\n********** %s **********" % label)
    print("%s\t%s\t%s\t%s\t" % ("t", "p", "h", "s"))
    for e in l:
        print("%.1f\t%.4g\t%.4g\t%.4g" % (e.t, e.p, e.h, e.s))


print_iapws97(_s, "Steam")
print_iapws97(_d, "Drainage")
print_iapws97(_h, "Water")


def iapws97_asarray(iapws_arrays, *args):
    rslt_arr = []
    for iapws_array in iapws_arrays:
        for a in args:
            rslt_arr.append([getattr(iapws, a) for iapws in iapws_array])
    return rslt_arr


np.savetxt(
    "out.csv",
    np.row_stack(iapws97_asarray([_s, _d, _h], "t", "p", "h", "s")),
    delimiter=",",
    fmt="%.4g")

system("libreoffice --calc out.csv 2>/dev/null &")
plt.show()  # (block=False)
