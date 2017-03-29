#!/bin/python3

f = open("result", "r")
lines = f.readlines()
f.close()
def get(symbol):
    for i in range(0, len(lines)):
        if lines[i] == ">>> " + symbol + "\n":
            return float(lines[i + 1])

D = get("D")
H = get("H")
_C_h = get("_C.h")
_SH_h = get("_SH.h")
eff_em = get("eff_em")
_1_h= get("_1.h")

_a = (

)

_y = [(_1_h - _h[0])]